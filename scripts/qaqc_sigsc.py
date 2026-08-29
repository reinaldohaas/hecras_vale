# -*- coding: utf-8 -*-
"""QA/QC das folhas 1 m do SIG-SC com o Copernico de arbitro.

    python scripts/qaqc_sigsc.py --detectar          varre as folhas
    python scripts/qaqc_sigsc.py --corrigir          conserta as doentes
    python scripts/qaqc_sigsc.py --figuras           antes/depois das piores

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
        # colar: molduras concentricas totalmente invalidas
        colar = 0
        while (colar < min(h, w) // 2 and
               inv[colar, :].all() and inv[-1 - colar, :].all()
               and inv[:, colar].all() and inv[:, -1 - colar].all()):
            colar += 1
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
        doente = (frac > LIM_FRAC
                  or (np.isfinite(dmed) and abs(dmed) > LIM_DIF))
        linhas.append({'folha': os.path.basename(p),
                       'frac_invalida': round(frac, 5),
                       'colar_px': colar * DECIM,
                       'dif_mediana': round(dmed, 2),
                       'dif_p95': round(dp95, 2),
                       'doente': int(doente)})
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
    for r in linhas:
        p = os.path.join(PASTA, r['folha'])
        with rasterio.open(p) as s:
            b = s.bounds
        if r['doente'] and r['frac_invalida'] > LIM_FRAC:
            cor, alfa = '#d62828', 0.7
        elif r['doente']:
            cor, alfa = '#f77f00', 0.7
        else:
            cor, alfa = '#2d6a4f', 0.25
        ax.add_patch(plt.Rectangle((b.left, b.bottom),
                     b.right - b.left, b.top - b.bottom,
                     fc=cor, alpha=alfa, ec='k', lw=0.2))
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_title('Saude das folhas SIG-SC: vermelho=furos/colar, '
                 'laranja=fora de nivel, verde=sa')
    fig.tight_layout()
    fig.savefig('doc/figuras/qaqc_sigsc_saude.png', dpi=110)
    print('mapa: doc/figuras/qaqc_sigsc_saude.png')


def carregar_relatorio():
    with open(REL) as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def corrigir():
    import rasterio
    from scipy.ndimage import distance_transform_edt
    cop = rasterio.open(COPERNICO)
    os.makedirs(PASTA_FIX, exist_ok=True)
    doentes = [r for r in carregar_relatorio() if r['doente'] == '1'
               and float(r['frac_invalida']) > LIM_FRAC]
    print(f'{len(doentes)} folhas a corrigir -> {PASTA_FIX}')
    for k, r in enumerate(doentes):
        src = os.path.join(PASTA, r['folha'])
        dst = os.path.join(PASTA_FIX, r['folha'])
        if os.path.exists(dst):
            continue
        with rasterio.open(src) as s:
            a = s.read(1).astype(np.float32)
            prof, nod = s.profile, s.nodata
            bounds = s.bounds
            res = s.transform.a
        inv = invalido(a, nod)
        if not inv.any():
            continue
        remendo = copernico_para(bounds, a.shape, cop)
        # esponja: peso do Copernico cai a zero a ESPONJA m da fronteira
        dist = distance_transform_edt(inv) * res
        w = np.clip(dist / ESPONJA, 0.0, 1.0).astype(np.float32)
        base = a.copy()
        base[inv] = 0.0
        mix = np.where(inv, w * remendo + (1 - w) * suave_borda(a, inv),
                       a)
        # onde o Copernico tambem e nan, herda o vizinho valido
        furo = inv & ~np.isfinite(remendo)
        if furo.any():
            mix[furo] = suave_borda(a, inv)[furo]
        prof.update(nodata=None, compress='deflate', predictor=3,
                    tiled=True)
        with rasterio.open(dst, 'w', **prof) as o:
            o.write(mix.astype(np.float32), 1)
        print(f'  [{k + 1}/{len(doentes)}] {r["folha"]}  '
              f'inv={float(r["frac_invalida"]):.3f}', flush=True)
    print('corrigidas. Aceite: rode --detectar de novo apontando '
          'PASTA_FIX (ou confira as figuras).')


def suave_borda(a, inv):
    """Preenche invalidos com o valor valido mais proximo (nearest)."""
    from scipy.ndimage import distance_transform_edt
    idx = distance_transform_edt(inv, return_indices=True,
                                 return_distances=False)
    return a[tuple(idx)]


def figuras(n=8):
    import rasterio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    doentes = sorted((r for r in carregar_relatorio()
                      if r['doente'] == '1'),
                     key=lambda r: -float(r['frac_invalida']))[:n]
    os.makedirs('doc/figuras/qaqc_sigsc', exist_ok=True)
    for r in doentes:
        fx = os.path.join(PASTA_FIX, r['folha'])
        if not os.path.exists(fx):
            continue
        with rasterio.open(os.path.join(PASTA, r['folha'])) as s:
            antes = s.read(1, out_shape=(s.height // 8, s.width // 8))
            nod = s.nodata
        with rasterio.open(fx) as s:
            depois = s.read(1, out_shape=(s.height // 8, s.width // 8))
        antes = antes.astype(np.float32)
        antes[invalido(antes, nod)] = np.nan
        vmin = np.nanpercentile(depois, 2)
        vmax = np.nanpercentile(depois, 98)
        fig, axs = plt.subplots(1, 2, figsize=(14, 7))
        for ax, dat, tit in [(axs[0], antes, 'antes (furos=branco)'),
                             (axs[1], depois, 'depois')]:
            im = ax.imshow(dat, vmin=vmin, vmax=vmax, cmap='terrain')
            ax.set_title(f'{r["folha"]}  {tit}')
            ax.axis('off')
        fig.colorbar(im, ax=axs, shrink=0.7, label='m')
        out = f'doc/figuras/qaqc_sigsc/{r["folha"][:-4]}.png'
        fig.savefig(out, dpi=100, bbox_inches='tight')
        plt.close(fig)
        print('figura:', out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--detectar', action='store_true')
    ap.add_argument('--corrigir', action='store_true')
    ap.add_argument('--figuras', action='store_true')
    args = ap.parse_args()
    if args.detectar:
        detectar()
    if args.corrigir:
        corrigir()
    if args.figuras:
        figuras()
    if not any(vars(args).values()):
        print(__doc__)


if __name__ == '__main__':
    main()
