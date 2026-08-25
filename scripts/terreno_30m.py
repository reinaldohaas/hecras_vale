# -*- coding: utf-8 -*-
"""Terreno do HEC-RAS a 30 m, do MDT SIG-SC, sobre o dominio de uma geometria.

    python scripts/terreno_30m.py modelo/so_mirim.g01 --nome mirim30

POR QUE 30 m, E POR QUE NAO E O TERRENO QUE JA EXISTE

  O terreno que o modelo carrega hoje TAMBEM e de 30 m, mas e Copernicus
  GLO-30: modelo de SUPERFICIE, com copa de mata e lamina d'agua dentro do
  dado. Serve de pano de fundo e nao serve para julgar leito.

  Este aqui e o SIG-SC -- terreno nu, levantado a 1 m -- reduzido a 30 m. E a
  mesma fonte contra a qual as secoes foram diagnosticadas, agora em forma que
  o RAS Mapper abre e desenha sobre os 141 km de uma vez.

O DOMINIO E O DO RIO, E NAO O DAS SECOES

  A envoltoria das cutlines NAO cobre o modelo: no Itajai-Mirim o eixo do rio
  segue 11,4 km ao sul e 4,6 km a oeste alem da secao mais extrema, porque ha
  trecho de centerline sem secao cortada. Um terreno feito so das cutlines
  deixa 9,9% dos vertices do eixo no vazio -- visivel no RAS Mapper como a
  ponta do rio pendurada fora do terreno. O dominio aqui e a uniao das
  cutlines COM os eixos (`Reach XY`), mais folga.

DE ONDE VEM, E O QUE ISSO CORRIGE

  Por padrao le as FOLHAS DE 1 m direto, por um VRT, e reduz a 30 m com media.

  O caminho alternativo (`--fonte <tif>`) parte do mosaico de 10 m que ja esta
  no disco. E muito mais rapido, mas aquele mosaico foi recortado num corredor
  em volta dos eixos: fora do corredor ele nao tem dado, e o terreno sai com um
  rombo no meio. Medido numa grade de 500 m sobre o dominio, 45,7% dos pontos
  estao vazios no mosaico de 10 m e 99,7% DESSES tem dado nas folhas de 1 m --
  ou seja, o rombo e recorte, e nao falta de levantamento.

ZERO E VAZIO -- MAS SO NA FONTE CERTA

  Nas folhas de 1 m o vazio e gravado como 0,00 e nao como NoData, e por isso
  o VRT nasce com `-srcnodata 0`. Ja o mosaico de 10 m JA declara -9999 e nao
  tem pixel algum em 0,00 (conferido na janela do Itajai-Mirim). Passar
  `-srcnodata 0` sobre ele seria erro, e nao redundancia: `-srcnodata`
  SUBSTITUI o nodata declarado, entao o -9999 viraria cota valida e a media de
  borda misturaria vazio com terreno -- mediram-se 9.458 celulas entre -9998 e
  -10 m assim, que no RAS Mapper viram fossos.
"""
import os
import sys

import numpy as np
import rasterio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qc_secoes import ler_secoes                            # noqa: E402
from qc_geometria import ler_eixos                          # noqa: E402
from mdt_sigsc import tiles_do_dominio, PASTA               # noqa: E402
from vale.terreno import _gdal, _rodar, piramides, tamanho  # noqa: E402
from vale.config import WKT                                 # noqa: E402

MARGEM = 2000.0          # m de folga alem do dominio do modelo
VAZIO = -9999.0


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def _e_valor(argv, x):
    """True se `x` e o VALOR de uma opcao `--algo`, e nao uma geometria."""
    i = argv.index(x)
    return i > 0 and argv[i - 1].startswith("--")


def dominio(geom, margem):
    """Envoltoria das cutlines E dos eixos, com folga.

    Aceita UMA geometria ou uma lista delas. Um terreno por rio nao serve
    quando os rios viram rede: o de 30 m que existia cobria o Mirim inteiro e
    deixava Norte, Sul e Oeste 100% de fora, com o Acu em 60% -- o RAS Mapper
    abria e nao desenhava nada. A uniao dos seis da 146 x 133 km.
    """
    geoms = [geom] if isinstance(geom, str) else list(geom)
    S, P, E = [], [], []
    for g in geoms:
        s = ler_secoes(g)
        S += s
        P += [np.asarray(d["cut"], float) for d in s]
        E += [np.asarray(e.coords, float) for e in ler_eixos(g).values()]
    T = np.vstack(P + E)
    C = np.vstack(P)
    if (T.min(0) < C.min(0)).any() or (T.max(0) > C.max(0)).any():
        d = np.abs(np.concatenate([C.min(0) - T.min(0), T.max(0) - C.max(0)]))
        print(f"   o eixo passa alem das cutlines em ate {d.max():.0f} m -- "
              "o dominio segue o eixo")
    return (T[:, 0].min() - margem, T[:, 1].min() - margem,
            T[:, 0].max() + margem, T[:, 1].max() + margem), S


def cobertura(tif, S, geom, log=print):
    """Mede a cobertura onde ela importa: nas secoes e ao longo do eixo."""
    with rasterio.open(tif) as s:
        a = s.read(1)
        tr = ~s.transform
        v = a > VAZIO + 1
        log(f"   {s.width} x {s.height} px   "
            f"cota {a[v].min():.2f} a {a[v].max():.2f} m")
        sujo = int((v & (a < -1)).sum())
        log(f"   celulas entre -9998 e -1 m (media contaminada): {sujo}")

        def tem(pts):
            c, l = tr * (pts[:, 0], pts[:, 1])
            d = (c >= 0) & (c < s.width) & (l >= 0) & (l < s.height)
            z = np.full(len(pts), VAZIO)
            ci = np.clip(c.astype(int), 0, s.width - 1)
            li = np.clip(l.astype(int), 0, s.height - 1)
            z[d] = a[li[d], ci[d]]
            return z > VAZIO + 1

        cs = []
        for d in S:
            A = np.asarray(d["cut"][0], float)
            B = np.asarray(d["cut"][-1], float)
            t = np.linspace(0, 1, 60)[:, None]
            cs.append(tem(A + (B - A) * t).mean())
        cs = np.array(cs)
        log(f"   cobertura nas secoes : mediana {100*np.median(cs):.0f}%   "
            f"min {100*cs.min():.0f}%   abaixo de 50%: {int((cs < 0.5).sum())}")
        E = np.vstack([np.asarray(e.coords, float)
                       for e in ler_eixos(geom).values()])
        fe = tem(E)
        log(f"   cobertura no eixo    : {100*fe.mean():.1f}%   "
            f"vertices no vazio: {int((~fe).sum())} de {len(E)}")


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    geom = [x for x in argv if not x.startswith("--")
            and not _e_valor(argv, x)]
    geom = geom if len(geom) > 1 else geom[0]
    nome = _arg(argv, "--nome", "mirim30")
    res = _arg(argv, "--res", 30.0, float)
    margem = _arg(argv, "--margem", MARGEM, float)
    fonte = _arg(argv, "--fonte")             # None = folhas de 1 m
    raiz = os.path.dirname(geom if isinstance(geom, str) else geom[0]) or "."
    if not isinstance(geom, str):
        raiz = os.path.dirname(raiz) or "."       # sobe da pasta do rio
    pasta = os.path.join(raiz, "Terrain")
    os.makedirs(pasta, exist_ok=True)

    (x0, y0, x1, y1), S = dominio(geom, margem)
    print(f"geometria : "
          + (geom if isinstance(geom, str)
             else ", ".join(os.path.basename(g) for g in geom))
          + f"   ({len(S)} secoes)")
    print(f"dominio   : {(x1-x0)/1000:.1f} x {(y1-y0)/1000:.1f} km  "
          f"(folga de {margem:g} m)")
    print(f"resolucao : {res:g} m")

    # ---- a entrada do warp
    if fonte:
        entrada = fonte
        print(f"fonte     : {fonte}  (mosaico pronto, nodata ja declarado)")
    else:
        tiles = tiles_do_dominio((x0, y0, x1, y1))
        if not tiles:
            raise SystemExit("nenhuma folha do SIG-SC sobre o dominio")
        b = sum(os.path.getsize(t) for t in tiles)
        print(f"fonte     : {len(tiles)} folhas do SIG-SC a 1 m "
              f"({tamanho(b)}) em {PASTA}")
        lista = os.path.join(pasta, f"{nome}_folhas.txt")
        with open(lista, "w", encoding="utf-8") as f:
            f.write(chr(10).join(tiles))
        entrada = os.path.join(pasta, f"{nome}_sigsc_1m.vrt")
        # -srcnodata 0 AQUI SIM: e a fonte em que zero e vazio
        _rodar([_gdal("gdalbuildvrt"),
                "-srcnodata", "0", "-vrtnodata", str(int(VAZIO)),
                "-input_file_list", lista, entrada], None, print)
        print(f"   VRT: {tamanho(os.path.getsize(entrada))} de indice")

    tif = os.path.join(pasta, f"MDT_SIGSC_{res:g}m.tif")
    # sem -srcnodata: o VRT ja declara -9999, e o mosaico de 10 m tambem.
    # Repetir aqui SUBSTITUIRIA essa declaracao -- ver o cabecalho.
    _rodar([_gdal("gdalwarp"),
            "-te", f"{x0:.0f}", f"{y0:.0f}", f"{x1:.0f}", f"{y1:.0f}",
            "-tr", str(res), str(res), "-tap",
            "-r", "average",             # media, e nao sorteio de um pixel
            "-dstnodata", str(int(VAZIO)),
            "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-wm", "2048",
            "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "ZLEVEL=6",
            "-co", "BIGTIFF=YES", "-co", "NUM_THREADS=ALL_CPUS",
            "-overwrite", entrada, tif], None, print)
    piramides(tif, log=print)
    print(f"\nraster: {tif}  ({tamanho(os.path.getsize(tif))})")
    cobertura(tif, S, geom if isinstance(geom, str) else geom[0])

    # ---- o .hdf, pelo RasProcess do proprio HEC-RAS
    from ras_commander import RasTerrain
    prj = os.path.join(pasta, f"{nome}.prj")
    with open(prj, "w", encoding="ascii") as f:
        f.write(WKT)
    destino = os.path.join(pasta, f"{nome}_Terreno.hdf")
    print(f"\nRasTerrain.create_terrain_hdf -> {os.path.basename(destino)}")
    RasTerrain.create_terrain_hdf(
        input_rasters=[tif], output_hdf=destino,
        projection_prj=prj, units="Meters", timeout_seconds=14400)

    saida = os.path.join(pasta,
                         f"{nome}_Terreno.{os.path.basename(tif)[:-4]}.tif")
    print("\nCONFERENCIA")
    for p in (destino, saida):
        ok = os.path.exists(p)
        print(f"   {'OK   ' if ok else 'FALTA'} {os.path.basename(p)}"
              + (f"   {tamanho(os.path.getsize(p))}" if ok else ""))
    return destino


if __name__ == "__main__":
    main(sys.argv[1:])
