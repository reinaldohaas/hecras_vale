# -*- coding: utf-8 -*-
"""Fecha o MDT ao longo dos rios do modelo: acha MUROS no talvegue do
mosaico v2 e classifica cada um (ponte OSM / represa SNISB-OSM / ?).

    python scripts/fechar_mdt.py --detectar     varre e classifica
    python scripts/fechar_mdt.py --figuras      perfil de cada muro

METODO

  O perfil do terreno ao longo da centerline, andando para a foz, so
  pode descer. O envelope condicionado q (minimo acumulado de montante
  para jusante) revela os muros: onde p - q > LIMIAR ha um obstaculo
  que a agua real nao ve -- tabuleiro de ponte que o MDT engoliu,
  aterro, ou uma represa DE VERDADE.

  Classificacao (fontes so para decidir duvida, ordem do professor):
    represa  ha barragem SNISB (oficial) ou dam/weir OSM a < 150 m
             -> NAO escavar; e candidata a estrutura no modelo
    ponte    ha ponte OSM a < 60 m -> escavar (brecha curta)
    ?        sem explicacao -> lista para inspecao no painel

Sai: doc/qaqc_sigsc/muros.csv + doc/painel/muros.geojson (p/ painel)
e figuras em doc/figuras/muros/.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)
os.chdir(RAIZ)

G01_HDF = 'taha_ai.g01.hdf'
MOSAICO = os.path.join('taha_ai_novo', 'Terrain',
                       'taha_ai_corredor_1m_v2.tif')
PASSO = 2.0          # m ao longo da centerline
LIMIAR = 1.0         # m acima do envelope = muro
MIN_LARG = 4.0       # m: muro mais curto que isto e ruido
CSV_MUROS = os.path.join('doc', 'qaqc_sigsc', 'muros.csv')
GJ_MUROS = os.path.join('doc', 'painel', 'muros.geojson')


def eixos_utm():
    """Centerlines em 31982, por rio/reach, de montante p/ jusante."""
    import h5py
    f = h5py.File(G01_HDF, 'r')
    rc = f['Geometry/River Centerlines']
    at = rc['Attributes'][:]
    pl_info = rc['Polyline Info'][:]
    pl_pts = rc['Polyline Points'][:]
    out = []
    for k in range(len(at)):
        rio = at['River Name'][k].decode().strip()
        reach = at['Reach Name'][k].decode().strip()
        j0, m = pl_info[k][0], pl_info[k][1]
        pts = pl_pts[j0:j0 + m][:, :2].astype(float)
        out.append((rio, reach, pts))
    return out


def amostrar(pts):
    """Reamostra a polyline em PASSO m; retorna (xy, s_acum)."""
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)])
    alvo = np.arange(0, s[-1], PASSO)
    x = np.interp(alvo, s, pts[:, 0])
    y = np.interp(alvo, s, pts[:, 1])
    return np.column_stack([x, y]), alvo


def pontos_de(caminho, tipos=None):
    """Centroides das feicoes de um geojson (4326 -> 31982)."""
    from pyproj import Transformer
    from shapely.geometry import shape
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    if not os.path.exists(caminho):
        return np.zeros((0, 2)), []
    gj = json.load(open(caminho, encoding='utf-8'))
    xs, meta = [], []
    for f in gj.get('features', []):
        if not f.get('geometry'):
            continue
        if tipos and f.get('properties', {}).get('tipo') not in tipos:
            pass
        c = shape(f['geometry']).centroid
        x, y = tr.transform(c.x, c.y)
        xs.append((x, y))
        meta.append(f.get('properties', {}))
    return np.array(xs), meta


def rs_max_por_rio():
    """Maior RS (m de estacionamento) com secao no modelo, por rio."""
    import h5py
    f = h5py.File(G01_HDF, 'r')
    at = f['Geometry/Cross Sections/Attributes'][:]
    out = {}
    for k in range(len(at)):
        rio = at['River'][k].decode().strip()
        try:
            rs = float(at['RS'][k].decode())
        except ValueError:
            continue
        out[rio] = max(out.get(rio, 0.0), rs)
    return out


def detectar_rede():
    """Deteccao NA REDE DA LAMINA (nao na centerline esquematica, que
    corta morros): grade 8 m, distancia geodesica desde a foz dentro
    da mascara d'agua, minimo de montante propagado pela rede;
    muro = z acima da agua que chega de montante."""
    import rasterio
    from rasterio.features import rasterize
    from rasterio.transform import from_bounds
    from scipy.spatial import cKDTree
    from scipy.ndimage import label as nlabel
    from shapely.geometry import shape
    from pyproj import Transformer
    import heapq

    RES = 8.0
    src = rasterio.open(MOSAICO)
    b = src.bounds
    W = int((b.right - b.left) / RES)
    H = int((b.top - b.bottom) / RES)
    t = from_bounds(b.left, b.bottom, b.right, b.top, W, H)
    print(f'grade {W}x{H} a {RES:.0f} m')

    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    lam = json.load(open(os.path.join('doc', 'painel',
                                      'lamina_sigsc.geojson')))
    formas = []
    from shapely.ops import transform as stransform
    for f in lam['features']:
        g = shape(f['geometry'])
        formas.append(stransform(lambda x, y: tr.transform(x, y), g))
    mask = rasterize(((g, 1) for g in formas), out_shape=(H, W),
                     transform=t, fill=0, dtype='uint8') == 1
    print(f'mascara d\'agua: {int(mask.sum())} celulas')

    z = src.read(1, out_shape=(H, W)).astype(np.float32)
    if src.nodata is not None:
        z[np.isclose(z, src.nodata)] = np.nan
    mask &= np.isfinite(z)

    # semente da foz: celula d'agua mais a leste na linha da barra
    ys, xs = np.where(mask)
    k0 = int(np.argmax(xs))
    print(f'foz na celula ({ys[k0]},{xs[k0]})')

    # Dijkstra geodesico (distancia em celulas) dentro da mascara
    INF = np.float32(np.inf)
    dist = np.full((H, W), INF, np.float32)
    dist[ys[k0], xs[k0]] = 0
    fila = [(0.0, ys[k0], xs[k0])]
    VIZ = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1),
           (1, -1), (1, 0), (1, 1)]
    while fila:
        d, i, j = heapq.heappop(fila)
        if d > dist[i, j]:
            continue
        for di, dj in VIZ:
            ni, nj = i + di, j + dj
            if 0 <= ni < H and 0 <= nj < W and mask[ni, nj]:
                nd = d + (1.414 if di and dj else 1.0)
                if nd < dist[ni, nj]:
                    dist[ni, nj] = nd
                    heapq.heappush(fila, (nd, ni, nj))
    lig = np.isfinite(dist) & mask
    print(f'ligadas a foz: {int(lig.sum())} de {int(mask.sum())}')

    # minimo de montante: processa em ordem DECRESCENTE de distancia
    ordem = np.argsort(dist[lig])[::-1]
    ii, jj = np.where(lig)
    ii, jj = ii[ordem], jj[ordem]
    q = np.where(lig, z, np.nan).astype(np.float32)
    for i, j in zip(ii, jj):
        m = q[i, j]
        for di, dj in VIZ:
            ni, nj = i + di, j + dj
            if (0 <= ni < H and 0 <= nj < W and lig[ni, nj]
                    and dist[ni, nj] > dist[i, j]
                    and q[ni, nj] < m):
                m = q[ni, nj]
        q[i, j] = min(q[i, j], m)
    exc = np.where(lig, z - q, 0.0)
    muro_mask = exc > LIMIAR
    print(f'celulas-muro: {int(muro_mask.sum())}')

    # agrupa muros e classifica
    rot, n = nlabel(muro_mask)
    p_pontes, _ = pontos_de('doc/osm/osm_pontes.geojson')
    p_rep1, m_rep1 = pontos_de('doc/osm/snisb_barragens.geojson')
    p_rep2, m_rep2 = pontos_de('doc/osm/osm_represas.geojson')
    p_rep = (np.vstack([p for p in (p_rep1, p_rep2) if len(p)])
             if (len(p_rep1) + len(p_rep2)) else np.zeros((0, 2)))
    m_rep = list(m_rep1) + list(m_rep2)
    kd_pontes = cKDTree(p_pontes) if len(p_pontes) else None
    kd_rep = cKDTree(p_rep) if len(p_rep) else None
    muros = []
    for lab in range(1, n + 1):
        sel = rot == lab
        npix = int(sel.sum())
        if npix * RES * RES < 200:      # < 200 m2 e ruido
            continue
        k = np.argmax(np.where(sel, exc, 0))
        i, j = np.unravel_index(k, exc.shape)
        x = b.left + (j + 0.5) * RES
        y = b.top - (i + 0.5) * RES
        classe, nome = '?', ''
        if kd_rep is not None:
            d, ki = kd_rep.query([x, y])
            if d < 150:
                classe = 'represa'
                nome = (m_rep[ki].get('BAR_NM_NOME')
                        or m_rep[ki].get('nome') or '')
        if classe == '?' and kd_pontes is not None:
            d, _ = kd_pontes.query([x, y])
            if d < 60:
                classe = 'ponte'
        muros.append({'x': round(float(x), 1),
                      'y': round(float(y), 1),
                      'altura_m': round(float(exc[i, j]), 2),
                      'area_m2': int(npix * RES * RES),
                      'classe': classe, 'nome': nome})
    muros.sort(key=lambda m: -m['altura_m'])
    os.makedirs(os.path.dirname(CSV_MUROS), exist_ok=True)
    with open(CSV_MUROS, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(muros[0]) if muros
                           else ['x'])
        w.writeheader()
        w.writerows(muros)
    inv = Transformer.from_crs(31982, 4326, always_xy=True)
    feats = []
    for m in muros:
        lon, lat = inv.transform(m['x'], m['y'])
        feats.append({'type': 'Feature', 'properties': m,
                      'geometry': {'type': 'Point',
                                   'coordinates': [round(lon, 6),
                                                   round(lat, 6)]}})
    json.dump({'type': 'FeatureCollection', 'features': feats},
              open(GJ_MUROS, 'w'), separators=(',', ':'))
    cont = {'ponte': 0, 'represa': 0, '?': 0}
    for m in muros:
        cont[m['classe']] += 1
    print(f'{len(muros)} muros: {cont["ponte"]} pontes, '
          f'{cont["represa"]} represas, {cont["?"]} sem explicacao')
    return muros


def detectar():
    import rasterio
    from scipy.spatial import cKDTree
    src = rasterio.open(MOSAICO)
    rs_max = rs_max_por_rio()
    p_pontes, _ = pontos_de('doc/osm/osm_pontes.geojson')
    p_rep1, m_rep1 = pontos_de('doc/osm/snisb_barragens.geojson')
    p_rep2, m_rep2 = pontos_de('doc/osm/osm_represas.geojson')
    p_rep = (np.vstack([p for p in (p_rep1, p_rep2) if len(p)])
             if (len(p_rep1) + len(p_rep2)) else np.zeros((0, 2)))
    m_rep = list(m_rep1) + list(m_rep2)
    kd_pontes = cKDTree(p_pontes) if len(p_pontes) else None
    kd_rep = cKDTree(p_rep) if len(p_rep) else None

    muros = []
    for rio, reach, pts in eixos_utm():
        xy, s = amostrar(pts)
        if len(xy) < 10:
            continue
        z = np.fromiter(
            (v[0] for v in src.sample([tuple(p) for p in xy])),
            dtype=np.float32, count=len(xy))
        z[np.isclose(z, src.nodata or -9999.0)] = np.nan
        ok = np.isfinite(z)
        if ok.sum() < 10:
            continue
        zi = np.interp(s, s[ok], z[ok])
        # garante ordem montante -> jusante (mediana dos 10% de cada
        # ponta decide; se a centerline veio da foz p/ cima, vira tudo)
        dez = max(3, len(zi) // 10)
        invertido = (np.median(zi[:dez]) < np.median(zi[-dez:]))
        if invertido:
            zi, xy, s = zi[::-1], xy[::-1], s[-1] - s[::-1]
        # a centerline guarda o rio ORIGINAL; o modelo esta amputado.
        # So vale o trecho de jusante coberto por secoes (RS maximo).
        alcance = rs_max.get(rio, 0.0) + 500.0
        if s[-1] - s[0] > alcance:
            corte = np.searchsorted(s, s[-1] - alcance)
            zi, xy, s = zi[corte:], xy[corte:], s[corte:]
        # envelope condicionado: minimo acumulado montante -> jusante
        q = np.minimum.accumulate(zi)
        exc = zi - q
        acima = exc > LIMIAR
        # segmentos contiguos
        i = 0
        while i < len(acima):
            if not acima[i]:
                i += 1
                continue
            j = i
            while j < len(acima) and acima[j]:
                j += 1
            larg = (j - i) * PASSO
            if larg >= MIN_LARG:
                k = i + int(np.argmax(exc[i:j]))
                x, y = xy[k]
                alt = float(exc[k])
                classe, nome = '?', ''
                if kd_rep is not None:
                    d, ki = kd_rep.query([x, y])
                    if d < 150:
                        classe = 'represa'
                        nome = (m_rep[ki].get('BAR_NM_NOME')
                                or m_rep[ki].get('nome') or '')
                if classe == '?' and kd_pontes is not None:
                    d, _ = kd_pontes.query([x, y])
                    if d < 60:
                        classe = 'ponte'
                muros.append({'rio': rio, 'reach': reach,
                              's_m': round(float(s[k]), 1),
                              'x': round(float(x), 1),
                              'y': round(float(y), 1),
                              'altura_m': round(alt, 2),
                              'largura_m': round(larg, 1),
                              'classe': classe, 'nome': nome})
            i = j
    os.makedirs(os.path.dirname(CSV_MUROS), exist_ok=True)
    with open(CSV_MUROS, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(muros[0]) if muros
                           else ['rio'])
        w.writeheader()
        w.writerows(muros)
    # geojson p/ o painel
    from pyproj import Transformer
    inv = Transformer.from_crs(31982, 4326, always_xy=True)
    feats = []
    for m in muros:
        lon, lat = inv.transform(m['x'], m['y'])
        feats.append({'type': 'Feature', 'properties': m,
                      'geometry': {'type': 'Point',
                                   'coordinates': [round(lon, 6),
                                                   round(lat, 6)]}})
    os.makedirs(os.path.dirname(GJ_MUROS), exist_ok=True)
    json.dump({'type': 'FeatureCollection', 'features': feats},
              open(GJ_MUROS, 'w'), separators=(',', ':'))
    n = {'ponte': 0, 'represa': 0, '?': 0}
    for m in muros:
        n[m['classe']] += 1
    print(f'{len(muros)} muros: {n["ponte"]} pontes, '
          f'{n["represa"]} represas, {n["?"]} sem explicacao')
    print(f'-> {CSV_MUROS} e {GJ_MUROS}')
    return muros


def figuras(n_max=12):
    import rasterio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    src = rasterio.open(MOSAICO)
    muros = [dict(r) for r in csv.DictReader(open(CSV_MUROS))]
    muros.sort(key=lambda m: -float(m['altura_m']))
    os.makedirs('doc/figuras/muros', exist_ok=True)
    eixos = {(r, rc): pts for r, rc, pts in eixos_utm()}
    for m in muros[:n_max]:
        pts = eixos.get((m['rio'], m['reach']))
        if pts is None:
            continue
        xy, s = amostrar(pts)
        s0 = float(m['s_m'])
        j = (s > s0 - 600) & (s < s0 + 600)
        z = np.fromiter(
            (v[0] for v in src.sample([tuple(p) for p in xy[j]])),
            dtype=np.float32, count=int(j.sum()))
        z[np.isclose(z, src.nodata or -9999.0)] = np.nan
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(s[j] / 1000, z, '-', lw=1.5, color='#1864ab')
        ax.axvline(s0 / 1000, color='#d62828', ls='--')
        ax.set_title(f"{m['rio']} {m['reach']} @ {s0:.0f} m -- "
                     f"{m['classe']} {m['nome']}  "
                     f"(+{m['altura_m']} m, {m['largura_m']} m)")
        ax.set_xlabel('s (km na centerline)')
        ax.set_ylabel('z (m)')
        ax.grid(alpha=0.3)
        out = (f"doc/figuras/muros/{m['rio']}_{m['reach']}_"
               f"{int(float(m['s_m']))}.png")
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print('figura:', out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--detectar', action='store_true')
    ap.add_argument('--detectar-rede', action='store_true')
    ap.add_argument('--figuras', action='store_true')
    args = ap.parse_args()
    if args.detectar_rede:
        detectar_rede()
    if args.detectar:
        detectar()
    if args.figuras:
        figuras()
    if not any(vars(args).values()):
        print(__doc__)


if __name__ == '__main__':
    main()
