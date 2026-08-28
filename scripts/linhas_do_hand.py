# -*- coding: utf-8 -*-
"""Gera as linhas INICIAIS de edge e bank a partir do terreno, por HAND.

    python scripts/linhas_do_hand.py taha_ai.g01.estado_mes_relevo1m \
        --hand-edge 15 --saida doc/linhas_hand

SO GERA LINHAS (GeoJSON por grupo + JSON do editor). Nada de g01.

COMO (ideia do Reinaldo, 27/08)

  EDGES  a pintura HAND vira o limite: altura acima do talvegue do rio
         mais proximo (Copernicus 45 m + talvegues do modelo). Cada
         pixel pertence ao rio mais proximo, entao onde duas planicies
         se tocam (Acu-Mirim) as edges dos dois rios se ENCONTRAM na
         linha media -- o divisor sai de graca. O contorno da mancha
         HAND <= --hand-edge de cada rio, separado em margem N e S pelo
         lado do eixo e ordenado pela quilometragem, e a edge inicial.

  BANKS  comecam no rio DESENHADO no SIG-SC/FBDS: o contorno dos
         poligonos RIOS_DUPLOS perto do eixo, tambem separado N/S e
         ordenado. E a margem d'agua real -- o editor ajusta dali.

  As linhas vao para o editor interativo (artifact) para ajuste a mouse
  e aprovacao; so depois viram secoes.
"""
import glob
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes                       # noqa: E402
from qc_geometria import ler_eixos                     # noqa: E402

BB = (580000, 6933000, 740000, 7070000)
RES = 45.0
CORES = {'Itajai_Acu': '#d62828', 'Itajai_Mirim': '#f77f00',
         'Itajai_Norte': '#5f0f40', 'Itajai_Oeste': '#9d4edd',
         'Itajai_Sul': '#0077b6', 'Rio_Benedito': '#2d6a4f',
         'Rio_dos_Cedros': '#74a57f', 'Rio_do_Testo': '#b5838d',
         'Rio_Iraputa': '#6c757d', 'Rio_Taio': '#c9184a',
         'Rio_Trombudo': '#7f5539', 'Rio_das_Pombas': '#3a5a40'}
ROTULOS = {'Itajai_Acu': 'Açu', 'Itajai_Mirim': 'Mirim',
           'Itajai_Norte': 'Norte', 'Itajai_Oeste': 'Oeste',
           'Itajai_Sul': 'Sul', 'Rio_Benedito': 'Benedito',
           'Rio_dos_Cedros': 'Cedros', 'Rio_do_Testo': 'Testo',
           'Rio_Iraputa': 'Iraputá', 'Rio_Taio': 'Taió',
           'Rio_Trombudo': 'Trombudo', 'Rio_das_Pombas': 'Pombas'}


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def eixos_e_talvegues(g01):
    from shapely.geometry import LineString
    S = ler_secoes(g01)
    E = ler_eixos(g01)
    eixos = {}
    for (rio, reach), lsg in sorted(E.items()):
        eixos.setdefault(rio, []).extend(list(lsg.coords))
    saida = {}
    for rio, coords in eixos.items():
        ls = LineString(coords)
        secs = sorted([(d['rs'], float(np.asarray(d['z'], float).min()))
                       for d in S if d['rio'] == rio and d['tipo'] == '1'])
        if not secs:
            continue
        rss = np.array([s[0] for s in secs])
        tals = np.array([s[1] for s in secs])
        pontos, cotas = [], []
        for s in np.arange(0, ls.length, 100):
            p = ls.interpolate(s)
            pontos.append((p.x, p.y))
            cotas.append(float(np.interp(ls.length - s, rss, tals)))
        saida[rio] = (ls, np.array(pontos), np.array(cotas))
    return saida


def lado_e_km(ls, P):
    """para cada ponto: lado (+1 N/esq, -1 S/dir) e quilometragem."""
    from shapely.geometry import Point
    lados, kms = [], []
    for x, y in P:
        s = ls.project(Point(x, y))
        a = ls.interpolate(max(s - 30, 0))
        b = ls.interpolate(min(s + 30, ls.length))
        cz = (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x)
        lados.append(1 if cz > 0 else -1)
        kms.append(ls.length - s)
    return np.array(lados), np.array(kms)


def em_partes(P, lados, kms, lado, simpl, salto=1500.0):
    """pontos de um lado -> TODAS as polilinhas continuas (ordenadas pela
    km, quebradas onde ha salto). Antes so o maior fragmento sobrevivia e
    os rios ficavam pelados; costurar fragmentos gerava cordas retas."""
    from shapely.geometry import LineString
    m = lados == lado
    if m.sum() < 4:
        return []
    o = np.argsort(-kms[m])
    pts = P[m][o]
    partes, atual = [], [pts[0]]
    for p in pts[1:]:
        if np.hypot(p[0] - atual[-1][0], p[1] - atual[-1][1]) > salto:
            partes.append(atual)
            atual = [p]
        else:
            atual.append(p)
    partes.append(atual)
    saida = []
    for pts in partes:
        if len(pts) < 4:
            continue
        ls = LineString(pts)
        if ls.length < 1500:      # fiapos nao valem edicao
            continue
        s = simpl
        ls = ls.simplify(s)
        while len(ls.coords) > 700:
            s *= 1.6
            ls = ls.simplify(s)
        saida.append([[round(x, 1), round(y, 1)] for x, y in ls.coords])
    return saida


def preencher_lacunas(ls, PB, lados, kms, mdt_src, passo=150.0,
                      km_min=0.0, km_max=1e12):
    """onde a FBDS nao tem margem num raio de 400 m de quilometragem,
    mede a LAMINA no MDT 1 m (transecto + medir de largura_do_sigsc) e
    acrescenta os dois pontos de borda d'agua."""
    if mdt_src is None:
        return PB, lados, kms
    import largura_do_sigsc as L
    T = mdt_src.transform
    tem = {1: np.sort(kms[lados == 1]), -1: np.sort(kms[lados == -1])}

    def coberto(lado, km):
        v = tem[lado]
        if len(v) == 0:
            return False
        i = np.searchsorted(v, km)
        for j in (i - 1, i):
            if 0 <= j < len(v) and abs(v[j] - km) < 400:
                return True
        return False

    novos_p, novos_l, novos_k = [], [], []
    n_novo = 0
    for s in np.arange(200, ls.length - 200, passo):
        km = ls.length - s
        if not (km_min <= km <= km_max):
            continue
        falta = [ld for ld in (1, -1) if not coberto(ld, km)]
        if not falta:
            continue
        P0 = np.asarray(ls.interpolate(s).coords[0])
        P1 = np.asarray(ls.interpolate(min(s + 30, ls.length)).coords[0])
        t = P1 - P0
        t = t / max(np.hypot(*t), 1e-9)
        nvec = np.array([-t[1], t[0]])
        z = L.transecto(mdt_src, T, P0, nvec)
        m = L.medir(z)
        if m is None or m == "solto":
            continue
        # bordas da lamina: reconstoi a faixa d'agua como em medir()
        c = len(z) // 2
        jan = z[c - 25:c + 26]
        if np.all(np.isnan(jan)):
            continue
        esp = np.nanmin(jan)
        i0 = c - 25 + int(np.nanargmin(jan))
        agua = np.abs(z - esp) <= L.TOL_LAMINA
        a = i0
        while a > 0 and agua[a - 1]:
            a -= 1
        b = i0
        while b < len(z) - 1 and agua[b + 1]:
            b += 1
        for off, lado in [(b - L.MEIA, 1), (a - L.MEIA, -1)]:
            if lado in falta:
                p = P0 + off * nvec
                novos_p.append(p)
                novos_l.append(lado)
                novos_k.append(km)
                n_novo += 1
    if novos_p:
        PB = np.vstack([PB, np.array(novos_p)])
        lados = np.concatenate([lados, np.array(novos_l)])
        kms = np.concatenate([kms, np.array(novos_k)])
        print(f'      lacunas da FBDS: +{n_novo} pontos da lamina SIG-SC')
    return PB, lados, kms


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    g01 = argv[0]
    hand_edge = _arg(argv, '--hand-edge', 15.0, float)
    pasta = _arg(argv, '--saida', 'doc/linhas_hand')
    os.makedirs(pasta, exist_ok=True)
    import rasterio
    from rasterio.vrt import WarpedVRT
    from rasterio.windows import from_bounds
    from scipy.spatial import cKDTree
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import LineString

    rios = eixos_e_talvegues(g01)
    S = ler_secoes(g01)
    todos_p, todos_c, todos_r = [], [], []
    for k, (rio, (ls, P, C)) in enumerate(sorted(rios.items())):
        todos_p.append(P)
        todos_c.append(C)
        todos_r += [rio] * len(P)
    arv = cKDTree(np.vstack(todos_p))
    cotas = np.concatenate(todos_c)
    rotulo = np.array(todos_r)

    src = rasterio.open('dem_bacia_itajai.tif')
    vrt = WarpedVRT(src, crs='EPSG:31982',
                    resampling=rasterio.enums.Resampling.bilinear)
    w = from_bounds(*BB, transform=vrt.transform)
    NY, NX = int((BB[3] - BB[1]) / RES), int((BB[2] - BB[0]) / RES)
    Z = vrt.read(1, window=w, out_shape=(NY, NX)).astype('float32')
    Z = np.where(Z < -10, np.nan, Z)
    gx = np.linspace(BB[0], BB[2], NX)
    gy = np.linspace(BB[3], BB[1], NY)
    GX, GY = np.meshgrid(gx, gy)
    dist, idx = arv.query(np.column_stack([GX.ravel(), GY.ravel()]),
                          workers=-1)
    hand = (Z.ravel() - cotas[idx]).reshape(Z.shape)
    perto = (dist < 12000).reshape(Z.shape)
    dono = rotulo[idx].reshape(Z.shape)
    print(f'HAND {NX}x{NY} pronto')

    # rio desenhado (FBDS) para os banks
    partes = []
    for shp in glob.glob('doc/fbds/*/SC_*_RIOS_DUPLOS.shp'):
        try:
            g = gpd.read_file(shp)
            if g.crs is None:
                g = g.set_crs(31982)
            partes.append(g.to_crs(31982)[['geometry']])
        except Exception:
            pass
    fbds = gpd.GeoDataFrame(pd.concat(partes, ignore_index=True),
                            crs=31982)
    import largura_do_sigsc as L
    try:
        mdt_src = rasterio.open(L.MDT)
    except Exception as e:
        print(f'MDT do SIG-SC indisponivel ({e}) -- lacunas ficam abertas')
        mdt_src = None

    linhas = []
    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name', 'properties': {'name': 'EPSG:31982'}},
           'features': []}
    for rio, (ls, P, C) in sorted(rios.items()):
        # ---- EDGES pelo HAND, por TRANSECTO sobre a mancha LIMPA por
        # MORFOLOGIA (ideia do Reinaldo: reconhecimento de imagem):
        # ABERTURA poda os dedos dendriticos dos vales laterais (mais
        # estreitos que ~500 m) e FECHAMENTO tapa buracos -- e o
        # contorno vira uma curva de verdade, nao uma teia
        from scipy import ndimage
        mancha = (hand <= hand_edge) & perto & (dono == rio)
        el = ndimage.generate_binary_structure(2, 1)
        mancha = ndimage.binary_closing(mancha, el, iterations=3)
        mancha = ndimage.binary_opening(mancha, el, iterations=6)
        mancha = ndimage.binary_closing(mancha, el, iterations=3)

        def hand_ok(x, y):
            j = int(round((x - BB[0]) / (BB[2] - BB[0]) * (NX - 1)))
            i = int(round((BB[3] - y) / (BB[3] - BB[1]) * (NY - 1)))
            if 0 <= i < NY and 0 <= j < NX:
                return bool(mancha[i, j])
            return False

        # so a QUILOMETRAGEM COM SECOES no modelo: alem da amputacao
        # (Mirim > km 40 etc.) o eixo continua mas nao ha modelo -- era
        # de la que vinham os rabiscos rio acima
        rss_rio = [d['rs'] for d in S
                   if d['rio'] == rio and d['tipo'] == '1']
        s_ini = max(100.0, ls.length - max(rss_rio))
        s_fim = min(ls.length - 100.0, ls.length - min(rss_rio))
        for lado, tag in [(1, 'N'), (-1, 'S')]:
            estacoes, afastas = [], []
            for s in np.arange(s_ini, s_fim, 250.0):
                P0 = np.asarray(ls.interpolate(s).coords[0])
                P1 = np.asarray(ls.interpolate(min(s + 30,
                                                   ls.length)).coords[0])
                t = P1 - P0
                t = t / max(np.hypot(*t), 1e-9)
                nv = np.array([-t[1], t[0]]) * lado
                d = RES
                while d <= 12000:
                    p = P0 + d * nv
                    if not hand_ok(p[0], p[1]):
                        break
                    d += RES
                estacoes.append((P0, nv))
                afastas.append(d - RES if d > 2 * RES else np.nan)
            afastas = np.array(afastas, float)
            # mediana movel (5 estacoes) tira o zigue-zague dos buracos
            lisas = np.array([np.nanmedian(afastas[max(0, i - 2):i + 3])
                              for i in range(len(afastas))])
            pontos_e = [P0 + a * nv for (P0, nv), a in
                        zip(estacoes, lisas) if np.isfinite(a)]
            if len(pontos_e) < 4:
                continue
            PE = np.array(pontos_e)
            pedacos = em_partes(PE, np.full(len(PE), lado),
                                np.arange(len(PE), 0, -1.0) * 250, lado,
                                80.0)
            for k, lin in enumerate(pedacos):
                suf = f' {k+1}/{len(pedacos)}' if len(pedacos) > 1 else ''
                nome = (f'{ROTULOS.get(rio, rio)} — edge {tag}{suf} '
                        f'(HAND {hand_edge:.0f} m)')
                linhas.append({'nome': nome, 'cor': CORES.get(rio, '#333'),
                               'grupo': ROTULOS.get(rio, rio),
                               'pontos': lin})
                geo['features'].append({'type': 'Feature',
                                        'properties': {'nome': nome},
                                        'geometry': {'type': 'LineString',
                                                     'coordinates': lin}})
        # ---- BANKS pelo rio desenhado
        try:
            faixa = ls.buffer(600)
            sel = fbds[fbds.intersects(faixa)]
            if len(sel):
                uniao = sel.union_all() if hasattr(sel, 'union_all') \
                    else sel.unary_union
                uniao = uniao.intersection(faixa)
                borda = []
                geoms = getattr(uniao, 'geoms', [uniao])
                for gpoly in geoms:
                    if gpoly.geom_type != 'Polygon':
                        continue
                    borda += list(gpoly.exterior.coords)
                if len(borda) > 10:
                    PB = np.array(borda)
                    lados, kms = lado_e_km(ls, PB)
                    # LACUNAS da FBDS preenchidas pela LAMINA do MDT 1 m
                    # (ordem do Reinaldo: "no FBDS e onde ele nao estiver
                    # disponivel usar o SIG-SC")
                    PB, lados, kms = preencher_lacunas(
                        ls, PB, lados, kms, mdt_src,
                        km_min=min(rss_rio), km_max=max(rss_rio))
                    # bank fora da quilometragem do modelo tambem sai
                    m_ok = (kms >= min(rss_rio) - 400) \
                        & (kms <= max(rss_rio) + 400)
                    PB, lados, kms = PB[m_ok], lados[m_ok], kms[m_ok]
                    for lado, tag in [(1, 'N'), (-1, 'S')]:
                        pedacos = em_partes(PB, lados, kms, lado, 25.0,
                                            salto=800.0)
                        for k, lin in enumerate(pedacos):
                            suf = (f' {k+1}/{len(pedacos)}'
                                   if len(pedacos) > 1 else '')
                            nome = (f'{ROTULOS.get(rio, rio)} — bank '
                                    f'{tag}{suf} (FBDS+SIG-SC)')
                            linhas.append({'nome': nome,
                                           'cor': CORES.get(rio, '#333'),
                                           'grupo': ROTULOS.get(rio, rio),
                                           'banco': True, 'pontos': lin})
                            geo['features'].append(
                                {'type': 'Feature',
                                 'properties': {'nome': nome},
                                 'geometry': {'type': 'LineString',
                                              'coordinates': lin}})
        except Exception as e:
            print(f'   {rio}: banks falharam ({e})')
        print(f'   {rio}: ok')

    # divisor Acu-Mirim: linha media entre os dois eixos
    if 'Itajai_Acu' in rios and 'Itajai_Mirim' in rios:
        acu_ls = rios['Itajai_Acu'][0]
        mir_ls = rios['Itajai_Mirim'][0]
        div = []
        for s in np.linspace(0, mir_ls.length, 500):
            pm = mir_ls.interpolate(mir_ls.length - s)
            pa = acu_ls.interpolate(acu_ls.project(pm))
            if np.hypot(pm.x - pa.x, pm.y - pa.y) > 45000:
                break
            div.append(((pm.x + pa.x) / 2, (pm.y + pa.y) / 2))
        div = [[round(x, 1), round(y, 1)] for x, y in
               LineString(div).simplify(60).coords]
        linhas.append({'nome': 'Divisor Açu–Mirim', 'cor': '#111111',
                       'grupo': 'Divisor', 'banco': True, 'pontos': div})

    dados = {'bbox': list(BB), 'linhas': linhas}
    json.dump(dados, open(os.path.join(pasta, 'linhas_editor.json'), 'w'))
    json.dump(geo, open(os.path.join(pasta, 'linhas_hand.geojson'), 'w'))
    n_v = sum(len(l['pontos']) for l in linhas)
    print(f'\n{len(linhas)} linhas, {n_v} vertices -> {pasta}/')


if __name__ == '__main__':
    main(sys.argv[1:])
