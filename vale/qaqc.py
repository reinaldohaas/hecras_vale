# -*- coding: utf-8 -*-
"""
QA/QC do SIG-SC contra o Copernicus, e preenchimento de vazio com vies medido.

O QUE SABEMOS DOS DOIS, e que e o que torna a qualificacao possivel:

  SIG-SC   MDT, 1 m. Solo exposto: NAO tem arvore, nao tem telhado, nao tem a
           lamina d'agua. E o dado bom -- quando existe. Tem vazio (1,5 a 1,9%
           por tile) e os vazios NAO sao declarados como NoData: saem como
           0,00, que passa por cota valida.

  Copernicus  MDS, 30 m. Superficie: inclui copa de mata, construcao e a lamina
           d'agua gravada como um plano na cota do espelho. Cobre tudo, sem
           vazio.

Dai sai a relacao que interessa: sobre o mesmo ponto, Copernicus >= SIG-SC
quase sempre, e a diferenca E a altura do que esta em cima do solo. Medida no
Itajai do Sul, essa diferenca chegou a 25,1 m -- foi o que produziu um falso
"corcovo" de 9 m no leito, que por um dia se pensou ser a Barragem Sul.

PREENCHER VAZIO COM COPERNICUS CRU ESTA ERRADO, e e o erro que este modulo
existe para nao cometer: enxertar um MDS num buraco de MDT insere um degrau da
altura da vegetacao ali. O certo e medir o vies LOCAL -- a mediana de
(SIG-SC menos Copernicus) num anel em volta do vazio, onde os dois existem --
e desconta-lo do Copernicus antes de enxertar. O vies e local porque a
vegetacao e local: mata fechada no vale, campo na crista.

Uso:
    python -m vale.qaqc                      relatorio de qualificacao
    python -m vale.qaqc --preencher SAIDA    grava o MDT com vazio corrigido
"""
import argparse
import os

import numpy as np

JANELA_VIES = 15          # celulas do anel em volta do vazio
MIN_AMOSTRAS = 30         # pares validos para o vies local valer
VIES_MAX = 60.0           # m; acima disso nao e vegetacao, e erro de registro


def _abrir(caminho):
    import rasterio
    return rasterio.open(caminho)


def amostrar_par(mdt, mds, n=200_000, semente=20260818):
    """Pares (SIG-SC, Copernicus) sobre os mesmos pontos, para estatistica.

    Amostragem aleatoria e nao a grade inteira: a 10 m a bacia tem 105 milhoes
    de celulas, e a estatistica converge muito antes disso.
    """
    rng = np.random.default_rng(semente)
    b = mdt.bounds
    xs = rng.uniform(b.left, b.right, n)
    ys = rng.uniform(b.bottom, b.top, n)

    def ler(ds):
        inv = ~ds.transform
        col, lin = inv * (xs, ys)
        c = np.floor(col).astype(np.int64)
        l = np.floor(lin).astype(np.int64)
        ok = (c >= 0) & (c < ds.width) & (l >= 0) & (l < ds.height)
        out = np.full(n, np.nan)
        if ok.any():
            banda = ds.read(1)
            out[ok] = banda[l[ok], c[ok]]
            nod = ds.nodata
            if nod is not None:
                out[out == nod] = np.nan
            out[out < -1000.0] = np.nan
        return out

    return ler(mdt), ler(mds), xs, ys


def qualificar(caminho_mdt, caminho_mds, log=print):
    """Compara os dois modelos e devolve as medidas que decidem o uso."""
    mdt, mds = _abrir(caminho_mdt), _abrir(caminho_mds)
    try:
        a, b, _, _ = amostrar_par(mdt, mds)
    finally:
        mdt.close()
        mds.close()

    dentro = np.isfinite(b)                      # dentro da area do Copernicus
    tem_mdt = np.isfinite(a) & dentro
    vazio = dentro & ~np.isfinite(a)
    d = a[tem_mdt] - b[tem_mdt]                  # negativo = MDT abaixo do MDS

    r = {
        "amostras": int(dentro.sum()),
        "cobertura_mdt": float(tem_mdt.sum() / max(dentro.sum(), 1)),
        "vazio_mdt": float(vazio.sum() / max(dentro.sum(), 1)),
        "dif_mediana": float(np.median(d)) if len(d) else float("nan"),
        "dif_p5": float(np.percentile(d, 5)) if len(d) else float("nan"),
        "dif_p95": float(np.percentile(d, 95)) if len(d) else float("nan"),
        "fracao_mdt_abaixo": float(np.mean(d < 0.0)) if len(d) else float("nan"),
        "dif_max_negativa": float(d.min()) if len(d) else float("nan"),
    }

    log("=" * 68)
    log("QUALIFICACAO DO SIG-SC CONTRA O COPERNICUS")
    log("=" * 68)
    log(f"  amostras dentro da area           {r['amostras']:>12,}")
    log(f"  cobertura do MDT                  {100*r['cobertura_mdt']:>11.2f}%")
    log(f"  vazio do MDT (a preencher)        {100*r['vazio_mdt']:>11.2f}%")
    log("")
    log("  diferenca MDT menos MDS, em metros (negativo = MDT mais baixo):")
    log(f"     p5 {r['dif_p5']:>8.2f}    mediana {r['dif_mediana']:>8.2f}"
        f"    p95 {r['dif_p95']:>8.2f}")
    log(f"     MDT abaixo do MDS em {100*r['fracao_mdt_abaixo']:.1f}% dos pontos"
        f"   (maior diferenca: {r['dif_max_negativa']:.1f} m)")
    log("")

    # O teste que interessa: o MDS e um limite SUPERIOR do MDT?
    if r["fracao_mdt_abaixo"] > 0.80:
        log("  [OK] o Copernicus e limite SUPERIOR do SIG-SC na maior parte da")
        log("       area, como esperado de superficie contra terreno. A")
        log("       diferenca e a altura do que esta em cima do solo, e e ela")
        log("       que sera descontada ao preencher vazio.")
    else:
        log("  [ATENCAO] o Copernicus NAO e limite superior consistente. Isso")
        log("       nao se explica por vegetacao: pode ser diferenca de datum")
        log("       vertical ou erro de registro entre os dois. Investigue")
        log("       ANTES de usar um para preencher o outro.")
    if abs(r["dif_mediana"]) > 30.0:
        log(f"  [ATENCAO] mediana de {r['dif_mediana']:.1f} m e alta demais para")
        log("       ser so vegetacao -- suspeite de datum vertical diferente.")
    return r


def vies_local(a, b, mascara_vazio, janela=JANELA_VIES):
    """Vies (MDT menos MDS) medido em volta de cada vazio, nao global.

    Global nao serve: a vegetacao e local. Mata fechada no fundo do vale da
    20 m de diferenca; campo na crista da 1 m. Preencher os dois com a mesma
    correcao poe um degrau em um dos casos.
    """
    from scipy import ndimage
    d = np.where(np.isfinite(a) & np.isfinite(b), a - b, np.nan)
    valido = np.isfinite(d).astype(np.float32)
    soma = ndimage.uniform_filter(np.nan_to_num(d), size=janela) * janela ** 2
    cont = ndimage.uniform_filter(valido, size=janela) * janela ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        media = soma / cont
    media[cont < MIN_AMOSTRAS] = np.nan
    # onde nem no anel ha amostra suficiente, cai para o vies global
    global_ = float(np.nanmedian(d)) if np.isfinite(d).any() else 0.0
    media = np.where(np.isfinite(media), media, global_)
    return np.clip(media, -VIES_MAX, VIES_MAX)


def preencher(caminho_mdt, caminho_mds, destino, log=print):
    """Grava o MDT com os vazios preenchidos por Copernicus CORRIGIDO.

    Tres passos, e o do meio e o que faz sentido fisico:
      1. reamostra o Copernicus para a grade do MDT (interpolacao bilinear --
         a 30 m para 10 m, vizinho mais proximo poria degraus de 30 m);
      2. mede o vies local (MDT menos MDS) no anel em volta de cada vazio;
      3. preenche com (Copernicus mais vies), que e a melhor estimativa do
         SOLO naquele ponto -- e nao da copa das arvores.
    """
    import rasterio
    from rasterio.warp import Resampling, reproject

    with rasterio.open(caminho_mdt) as ds_a:
        a = ds_a.read(1).astype("float32")
        perfil = ds_a.profile.copy()
        nod_a = ds_a.nodata
        transform_a, crs_a = ds_a.transform, ds_a.crs
        forma = a.shape
    if nod_a is not None:
        a[a == nod_a] = np.nan
    a[a < -1000.0] = np.nan

    b = np.full(forma, np.nan, dtype="float32")
    with rasterio.open(caminho_mds) as ds_b:
        reproject(source=rasterio.band(ds_b, 1), destination=b,
                  src_transform=ds_b.transform, src_crs=ds_b.crs,
                  dst_transform=transform_a, dst_crs=crs_a,
                  resampling=Resampling.bilinear,
                  src_nodata=ds_b.nodata, dst_nodata=np.nan)
    b[b < -1000.0] = np.nan

    vazio = (~np.isfinite(a)) & np.isfinite(b)
    n_vazio = int(vazio.sum())
    if n_vazio == 0:
        log("   nenhum vazio a preencher")
    else:
        vies = vies_local(a, b, vazio)
        a_novo = a.copy()
        a_novo[vazio] = b[vazio] + vies[vazio]
        log(f"   vazios preenchidos: {n_vazio:,} celulas "
            f"({100*n_vazio/a.size:.2f}% da grade)")
        log(f"   vies aplicado: mediana {np.nanmedian(vies[vazio]):.2f} m, "
            f"faixa {np.nanmin(vies[vazio]):.2f} a {np.nanmax(vies[vazio]):.2f} m")
        a = a_novo

    perfil.update(dtype="float32", nodata=-9999.0, compress="deflate",
                  tiled=True, bigtiff="YES")
    saida = np.where(np.isfinite(a), a, -9999.0).astype("float32")
    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    with rasterio.open(destino, "w", **perfil) as ds:
        ds.write(saida, 1)
    log(f"   gravado: {destino}")
    return destino


def main(argv=None):
    p = argparse.ArgumentParser(
        description="QA/QC do SIG-SC contra o Copernicus")
    p.add_argument("--mdt", default="modelo/Terrain/vale_sigsc_10m.tif")
    p.add_argument("--mds", default="Terrain/Terreno_Copernicus.tif")
    p.add_argument("--preencher", default=None,
                   help="grava o MDT com vazio corrigido neste caminho")
    a = p.parse_args(argv)

    for c in (a.mdt, a.mds):
        if not os.path.exists(c):
            print(f"nao encontrei {c}")
            print("Gere o MDT primeiro: python -m vale 3 fonte=sigsc")
            return 2
    qualificar(a.mdt, a.mds)
    if a.preencher:
        print()
        print("PREENCHIMENTO DOS VAZIOS")
        preencher(a.mdt, a.mds, a.preencher)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
