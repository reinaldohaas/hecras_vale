# -*- coding: utf-8 -*-
"""QA/QC das folhas 1 m do SIG-SC com o Copernico de arbitro.

    python scripts/qaqc_sigsc.py --detectar   varre as folhas e classifica
    python scripts/qaqc_sigsc.py --mosaico    remonta o corredor 1 m (v2)
    python scripts/qaqc_sigsc.py --figuras    antes/depois nos piores miolos

A doenca do mosaico do corredor (descontinuidades que o professor viu no
QGIS) nasce nas folhas: colares de zero na borda, vazios internos e
folhas deslocadas verticalmente. Este programa:

  --detectar  le cada folha decimada (~10 m) e mede
                frac_zero      fracao de pixels <= ZERO_MAX (0 nao e mar)
                frac_nodata    fracao de nodata declarado
                colar          largura media (px) de moldura zero/nodata
                dif_med/p95    desvio contra o Copernico 30 m reamostrado
              e grava doc/qaqc_sigsc/relatorio.csv + mapa de saude.
  --corrigir  para folha doente: pixels invalidos recebem o Copernico
              reamostrado bilinear para 1 m, com esponja de ESPONJA m de
              mistura linear na fronteira valida/invalida. Escreve em
              PASTA_FIX (nunca sobrescreve a folha original).
  --figuras   gera antes/depois das N piores.

Aceite: folha corrigida sem zero/nodata; costura entre vizinhas com
degrau p95 < 0.5 m (medido em --detectar na rodada seguinte usando
PASTA_FIX por cima da original).
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)
os.chdir(RAIZ)

PASTA = r'C:\Users\haas\Downloads\sigsc'
PASTA_FIX = r'C:\Users\haas\Downloads\sigsc_corrigido'
COPERNICO = os.path.join('doc', 'qgis', 'relevo_copernicus_30m_utm.tif')
REL = os.path.join('doc', 'qaqc_sigsc', 'relatorio.csv')
ZERO_MAX = 0.05      # <= isto e zero de colar (mar de verdade e < 2 m
                     # so na foz; la o Copernico tambem e ~0, tanto faz)
DECIM = 10           # deteccao a ~10 m
ESPONJA = 20.0       # metros de mistura na fronteira do remendo
LIM_FRAC = 0.001     # folha doente se > 0.1% de pixels invalidos
LIM_DIF = 15.0       # ou se |mediana - Copernico| > 15 m (folha fora
                     # de nivel; vegetacao explica ate ~10 m)


def listar():
    return sorted(glob.glob(os.path.join(PASTA, 'MDT_*.tif')))


def invalido(a, nodata):
    m = ~np.isfinite(a) | (a <= ZERO_MAX)
    if nodata is not None and np.isfinite(nodata):
        m |= np.isclose(a, nodata)
    return m


def copernico_para(bounds, shape, vrt_src):
    """Copernico reamostrado bilinear para a janela/forma pedida."""
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_bounds
    dst = np.full(shape, np.nan, np.float32)
    tr = from_bounds(*bounds, shape[1], shape[0])
    reproject(source=vrt_src.read(1), destination=dst,
              src_transform=vrt_src.transform, src_crs=vrt_src.crs,
              dst_transform=tr, dst_crs=vrt_src.crs,
              src_nodata=vrt_src.nodata,
              resampling=Resampling.bilinear)
    return dst


def detectar():
    import rasterio
    cop = rasterio.open(COPERNICO)
    cop_arr = cop.read(1).astype(np.float32)
    if cop.nodata is not None:
        cop_arr[np.isclose(cop_arr, cop.nodata)] = np.nan
    os.makedirs(os.path.dirname(REL), exist_ok=True)
    linhas = []
    folhas = listar()
    for k, p in enumerate(folhas):
        with rasterio.open(p) as s:
            h = max(1, s.height // DECIM)
            w = max(1, s.width // DECIM)
            a = s.read(1, out_shape=(h, w)).astype(np.float32)
            nod, bounds = s.nodata, s.bounds
        inv = invalido(a, nod)
        frac = float(inv.mean())
        # colar: molduras concentricas totalmente invalidas -- NORMAL
        # (toda folha tem; o vizinho cobre no mosaico)
        colar = 0
        while (colar < min(h, w) // 2 and
               inv[colar, :].all() and inv[-1 - colar, :].all()
               and inv[:, colar].all() and inv[:, -1 - colar].all()):
            colar += 1
        # doenca de verdade: invalido NO MIOLO (alem do colar + 2 px)
        m = colar + 2
        miolo = inv[m:h - m, m:w - m]
        frac_miolo = float(miolo.mean()) if miolo.size else 1.0
        # desvio contra o Copernico na mesma grade decimada
        from rasterio.transform import from_bounds as fb
        from rasterio.warp import reproject, Resampling
        cop_dec = np.full((h, w), np.nan, np.float32)
        reproject(cop_arr, cop_dec, src_transform=cop.transform,
                  src_crs=cop.crs, dst_crs=cop.crs,
                  dst_transform=fb(*bounds, w, h),
                  resampling=Resampling.bilinear)
        ok = ~inv & np.isfinite(cop_dec)
        if ok.sum() > 50:
            d = (a - cop_dec)[ok]
            dmed = float(np.median(d))
            dp95 = float(np.percentile(np.abs(d - dmed), 95))
        else:
            dmed = dp95 = np.nan
        if frac > 0.99:
            classe = 'vazia'          # nada a aproveitar
        elif np.isfinite(dmed) and abs(dmed) > LIM_DIF:
            classe = 'fora_de_nivel'
        elif frac_miolo > LIM_FRAC:
            classe = 'furada'         # invalido no miolo
        else:
            classe = 'sa'             # so o colar, que e normal
        linhas.append({'folha': os.path.basename(p),
                       'frac_invalida': round(frac, 5),
                       'frac_miolo': round(frac_miolo, 5),
                       'colar_px': colar * DECIM,
                       'dif_mediana': round(dmed, 2),
                       'dif_p95': round(dp95, 2),
                       'classe': classe,
                       'doente': int(classe != 'sa')})
        if (k + 1) % 100 == 0:
            print(f'  {k + 1}/{len(folhas)}...', flush=True)
    with open(REL, 'w', newline='') as fh:
        wcsv = csv.DictWriter(fh, fieldnames=list(linhas[0]))
        wcsv.writeheader()
        wcsv.writerows(linhas)
    n = sum(r['doente'] for r in linhas)
    print(f'{len(linhas)} folhas, {n} doentes -> {REL}')
    mapa_saude(linhas)
    return linhas


def mapa_saude(linhas):
    import rasterio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(14, 11))
    CORES = {'sa': ('#2d6a4f', 0.25), 'furada': ('#d62828', 0.7),
             'vazia': ('#6c757d', 0.8), 'fora_de_nivel': ('#f77f00', 0.7)}
    for r in linhas:
        p = os.path.join(PASTA, r['folha'])
        with rasterio.open(p) as s:
            b = s.bounds
        cor, alfa = CORES[r['classe']]
        ax.add_patch(plt.Rectangle((b.left, b.bottom),
                     b.right - b.left, b.top - b.bottom,
                     fc=cor, alpha=alfa, ec='k', lw=0.2))
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_title('Saude das folhas SIG-SC: vermelho=furada no miolo, '
                 'cinza=vazia, laranja=fora de nivel, verde=sa')
    fig.tight_layout()
    fig.savefig('doc/figuras/qaqc_sigsc_saude.png', dpi=110)
    print('mapa: doc/figuras/qaqc_sigsc_saude.png')


def carregar_relatorio():
    with open(REL) as fh:
        return [dict(r) for r in csv.DictReader(fh)]


MOSAICO_VELHO = os.path.join('taha_ai_novo', 'Terrain',
                             'taha_ai_corredor_1m_completo.tif')
MOSAICO_NOVO = os.path.join('taha_ai_novo', 'Terrain',
                            'taha_ai_corredor_1m_v2.tif')
BLOCO = 4096
HALO = 64


def mosaico():
    """Remonta o corredor 1 m: zero e NODATA (o defeito raiz), folha
    valida ganha, furo recebe Copernico bilinear com esponja."""
    import rasterio
    from rasterio.windows import Window, from_bounds as win_de
    from rasterio.transform import from_bounds as tr_de
    from rasterio.warp import reproject, Resampling
    from scipy.ndimage import distance_transform_edt

    velho = rasterio.open(MOSAICO_VELHO)
    prof = velho.profile.copy()
    prof.update(nodata=-9999.0, compress='deflate', predictor=3,
                tiled=True, blockxsize=512, blockysize=512,
                bigtiff='YES')
    cop = rasterio.open(COPERNICO)
    cop_arr = cop.read(1).astype(np.float32)
    if cop.nodata is not None:
        cop_arr[np.isclose(cop_arr, cop.nodata)] = np.nan

    # indice espacial simples das folhas
    folhas = []
    for p in listar():
        with rasterio.open(p) as s:
            folhas.append((s.bounds, p))

    T = velho.transform
    W, H = velho.width, velho.height
    with rasterio.open(MOSAICO_NOVO, 'w', **prof) as out:
        ny = int(np.ceil(H / BLOCO))
        nx = int(np.ceil(W / BLOCO))
        for jy in range(ny):
            for jx in range(nx):
                r0, c0 = jy * BLOCO, jx * BLOCO
                r1 = min(r0 + BLOCO, H)
                c1 = min(c0 + BLOCO, W)
                # janela com halo (p/ esponja atravessar o bloco)
                hr0, hc0 = max(0, r0 - HALO), max(0, c0 - HALO)
                hr1, hc1 = min(H, r1 + HALO), min(W, c1 + HALO)
                hh, hw = hr1 - hr0, hc1 - hc0
                oeste, norte = T * (hc0, hr0)
                leste, sul = T * (hc1, hr1)
                # pegada do corredor = onde o mosaico velho tinha dado
                # (mesmo errado); fora dela nao se inventa terreno
                va = velho.read(1, window=Window(hc0, hr0, hw, hh)
                                ).astype(np.float32)
                pegada = np.isfinite(va)
                if velho.nodata is not None:
                    pegada &= ~np.isclose(va, velho.nodata)
                if not pegada.any():
                    out.write(np.full((r1 - r0, c1 - c0), -9999.0,
                                      np.float32), 1,
                              window=Window(c0, r0, c1 - c0, r1 - r0))
                    continue
                acc = np.full((hh, hw), np.nan, np.float32)
                for b, p in folhas:
                    if (b.right < oeste or b.left > leste
                            or b.top < sul or b.bottom > norte):
                        continue
                    with rasterio.open(p) as s:
                        try:
                            wj = win_de(oeste, sul, leste, norte,
                                        s.transform)
                            a = s.read(1, window=wj, boundless=True,
                                       fill_value=0.0,
                                       out_shape=(hh, hw)
                                       ).astype(np.float32)
                        except Exception:
                            continue
                    ok = np.isfinite(a) & (a > ZERO_MAX)
                    poe = ok & ~np.isfinite(acc)
                    acc[poe] = a[poe]
                acc[~pegada] = np.nan
                furo = ~np.isfinite(acc) & pegada
                if furo.any():
                    remendo = np.full((hh, hw), np.nan, np.float32)
                    reproject(cop_arr, remendo,
                              src_transform=cop.transform,
                              src_crs=cop.crs, dst_crs=cop.crs,
                              dst_transform=tr_de(oeste, sul, leste,
                                                  norte, hw, hh),
                              resampling=Resampling.bilinear)
                    dist = distance_transform_edt(furo)  # px = m (1 m)
                    peso = np.clip(dist / ESPONJA, 0, 1)
                    perto = suave_borda(acc, furo)
                    mix = peso * remendo + (1 - peso) * perto
                    mix[~np.isfinite(remendo)] = perto[
                        ~np.isfinite(remendo)]
                    acc[furo] = mix[furo]
                acc[~np.isfinite(acc)] = -9999.0
                out.write(acc[r0 - hr0:r0 - hr0 + (r1 - r0),
                              c0 - hc0:c0 - hc0 + (c1 - c0)],
                          1, window=Window(c0, r0, c1 - c0, r1 - r0))
            print(f'  faixa {jy + 1}/{ny}', flush=True)
    print(f'mosaico novo: {MOSAICO_NOVO} (o velho fica intacto)')


def suave_borda(a, inv):
    """Preenche invalidos com o valor valido mais proximo (nearest)."""
    from scipy.ndimage import distance_transform_edt
    idx = distance_transform_edt(inv, return_indices=True,
                                 return_distances=False)
    return a[tuple(idx)]


def figuras(n=6):
    """Antes/depois do mosaico nos miolos mais furados."""
    import rasterio
    from rasterio.windows import from_bounds as win_de
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    piores = sorted((r for r in carregar_relatorio()
                     if r['classe'] == 'furada'),
                    key=lambda r: -float(r['frac_miolo']))[:n]
    os.makedirs('doc/figuras/qaqc_sigsc', exist_ok=True)
    for r in piores:
        with rasterio.open(os.path.join(PASTA, r['folha'])) as s:
            b = s.bounds
        paineis = []
        for tif, tit in [(MOSAICO_VELHO, 'mosaico ANTES'),
                         (MOSAICO_NOVO, 'mosaico DEPOIS')]:
            if not os.path.exists(tif):
                continue
            with rasterio.open(tif) as s:
                wj = win_de(b.left, b.bottom, b.right, b.top,
                            s.transform)
                a = s.read(1, window=wj, boundless=True,
                           out_shape=(800, 1000)).astype(np.float32)
                if s.nodata is not None:
                    a[np.isclose(a, s.nodata)] = np.nan
            a[a <= ZERO_MAX] = np.nan
            paineis.append((a, tit))
        if not paineis:
            continue
        vmin = np.nanpercentile(paineis[-1][0], 2)
        vmax = np.nanpercentile(paineis[-1][0], 98)
        fig, axs = plt.subplots(1, len(paineis),
                                figsize=(7 * len(paineis), 6))
        axs = np.atleast_1d(axs)
        for ax, (dat, tit) in zip(axs, paineis):
            im = ax.imshow(dat, vmin=vmin, vmax=vmax, cmap='terrain')
            ax.set_title(f'{r["folha"][:-4]}\n{tit} (branco=invalido)')
            ax.axis('off')
        fig.colorbar(im, ax=list(axs), shrink=0.7, label='m')
        out = f'doc/figuras/qaqc_sigsc/{r["folha"][:-4]}.png'
        fig.savefig(out, dpi=100, bbox_inches='tight')
        plt.close(fig)
        print('figura:', out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--detectar', action='store_true')
    ap.add_argument('--mosaico', action='store_true')
    ap.add_argument('--figuras', action='store_true')
    args = ap.parse_args()
    if args.detectar:
        detectar()
    if args.mosaico:
        mosaico()
    if args.figuras:
        figuras()
    if not any(vars(args).values()):
        print(__doc__)


if __name__ == '__main__':
    main()
