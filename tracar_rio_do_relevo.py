# -*- coding: utf-8 -*-
"""
Traca o eixo dos rios A PARTIR DO RELEVO, nao do desenho cartografico.

MOTIVO
------
Ate aqui o eixo vinha da base da ANA (BHO 2017), que e um tracado cartografico
a 1:100.000. Ele nao segue o talvegue do DEM: passa perto do rio, nao no fundo
do vale. Cortar secoes perpendiculares a um eixo que nao esta no talvegue faz
o corte atravessar a calha de esguelha, encontrar o leito antigo ou o meandro
vizinho, e as bank lines saem cruzadas em estrela -- foi o que apareceu no RAS
Mapper no Itajai-Mirim e no Acu, e nenhuma janela de busca conserta isso,
porque o defeito esta no eixo.

Um tracado derivado do proprio relevo resolve por construcao: ele passa onde a
agua acumula, entao o talvegue esta no eixo por definicao, e a secao
perpendicular cruza a calha em angulo reto.

COMO
----
  1. recorta o DEM na bacia e reamostra (o tracado e suavizado depois, entao
     90 m basta e deixa o processamento viavel em numpy puro);
  2. PREENCHE DEPRESSOES por priority-flood (Barnes et al. 2014) -- sem isso o
     caminhamento cai num pit e para. Impoe um gradiente minimo nas areas
     planas, que e o que permite atravessar a lamina achatada do Copernicus;
  3. direcao de fluxo D8 -- para cada celula, a vizinha de maior declive;
  4. caminha de JUSANTE PARA MONTANTE nao: caminha da cabeceira ate a foz
     seguindo D8, o que sempre chega ao mar;
  5. suaviza e devolve o eixo em EPSG:31982.

Uso:
    python tracar_rio_do_relevo.py                 traca todos os rios
    python tracar_rio_do_relevo.py mirim acu       so estes
Saida:
    eixos_do_relevo.geojson    um LineString por rio, em UTM 22S
"""
import heapq
import os
import sys

import numpy as np
import rasterio
from rasterio.enums import Resampling
import geopandas as gpd
from shapely.geometry import LineString, Point

DEM = os.path.join("Terrain", "Terreno_Copernicus.tif")
SAIDA = "eixos_do_relevo.geojson"
RES = 90.0          # m; o eixo e suavizado depois, entao isto basta
EPS = 0.001         # gradiente imposto nas areas planas (m por celula)
SUAVIZA = 5         # janela da media movel, em celulas

# vizinhanca D8: (di, dj, distancia relativa)
VIZ = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
       (-1, -1, 1.4142), (-1, 1, 1.4142), (1, -1, 1.4142), (1, 1, 1.4142)]


def ler_dem():
    with rasterio.open(DEM) as d:
        fator = RES / abs(d.transform.a)
        h = int(d.height / fator)
        w = int(d.width / fator)
        # Resampling.min nao vale para leitura no rasterio, so para warp. E o
        # minimo importa: com a media o talvegue some sob as encostas.
        # Entao le inteiro e reduz por blocos com np.nanmin.
        a = d.read(1).astype(np.float64)
        f = int(round(fator))
        hh, ww = (a.shape[0] // f) * f, (a.shape[1] // f) * f
        with np.errstate(all="ignore"):
            z = np.nanmin(a[:hh, :ww].reshape(hh // f, f, ww // f, f),
                          axis=(1, 3))
        h, w = z.shape
        tr = d.transform * d.transform.scale(f, f)
        nod = d.nodata
    if nod is not None:
        z[z == nod] = np.nan
    z[z < -500] = np.nan
    print(f"DEM {w}x{h} a {RES:.0f} m  ({w*h/1e6:.1f} M celulas)")
    return z, tr


def preencher(z):
    """Priority-flood com gradiente minimo (Barnes, Lehman & Mulla 2014).

    Cada celula sai da fila uma vez so, entao e O(n log n) e roda em um passe
    -- o preenchimento iterativo vetorizado precisaria de uma iteracao por
    celula de comprimento da area plana, o que no baixo Itajai seriam milhares.
    """
    h, w = z.shape
    cheio = np.full((h, w), np.inf)
    visto = np.zeros((h, w), dtype=bool)
    fila = []
    val = np.isfinite(z)
    # semente: a borda do dominio valido (inclui o mar e as bordas do recorte)
    borda = val & ~(
        np.pad(val[1:, :], ((0, 1), (0, 0))) &
        np.pad(val[:-1, :], ((1, 0), (0, 0))) &
        np.pad(val[:, 1:], ((0, 0), (0, 1))) &
        np.pad(val[:, :-1], ((0, 0), (1, 0))))
    for i, j in zip(*np.nonzero(borda)):
        heapq.heappush(fila, (float(z[i, j]), int(i), int(j)))
        cheio[i, j] = z[i, j]
        visto[i, j] = True
    n = 0
    while fila:
        zc, i, j = heapq.heappop(fila)
        n += 1
        for di, dj, _ in VIZ:
            a, b = i + di, j + dj
            if a < 0 or a >= h or b < 0 or b >= w or visto[a, b]:
                continue
            if not val[a, b]:
                continue
            cheio[a, b] = max(z[a, b], zc + EPS)
            visto[a, b] = True
            heapq.heappush(fila, (float(cheio[a, b]), a, b))
    print(f"  depressoes preenchidas: {n} celulas processadas, "
          f"{int((cheio > z + 1e-9).sum())} elevadas")
    cheio[~val] = np.nan
    return cheio


def direcao(zc):
    """D8: indice do vizinho de maior declive, -1 se nao houver."""
    h, w = zc.shape
    melhor = np.full((h, w), -1, dtype=np.int8)
    decl = np.zeros((h, w))
    for k, (di, dj, dist) in enumerate(VIZ):
        viz = np.full((h, w), np.nan)
        i0, i1 = max(0, -di), h - max(0, di)
        j0, j1 = max(0, -dj), w - max(0, dj)
        viz[i0:i1, j0:j1] = zc[i0 + di:i1 + di, j0 + dj:j1 + dj]
        with np.errstate(invalid="ignore"):
            s = (zc - viz) / (dist * RES)
            m = np.isfinite(s) & (s > decl)
        decl[m] = s[m]
        melhor[m] = k
    return melhor


def caminhar(dirs, i, j, limite=200000):
    """Segue D8 de (i,j) ate onde nao houver mais para onde ir."""
    h, w = dirs.shape
    cam = [(i, j)]
    for _ in range(limite):
        k = dirs[i, j]
        if k < 0:
            break
        di, dj, _ = VIZ[k]
        i, j = i + di, j + dj
        if i < 0 or i >= h or j < 0 or j >= w:
            break
        cam.append((i, j))
    return cam


def encaixar(zc, i, j, raio=8):
    """Move o ponto para a celula mais BAIXA numa vizinhanca -- o talvegue.

    A cabeceira da ANA cai na encosta com frequencia; comecando o caminhamento
    ali, o D8 desce a encosta ate o rio, o que poe um trecho de ladeira no
    inicio do eixo.
    """
    h, w = zc.shape
    i0, i1 = max(0, i - raio), min(h, i + raio + 1)
    j0, j1 = max(0, j - raio), min(w, j + raio + 1)
    bl = zc[i0:i1, j0:j1]
    if not np.isfinite(bl).any():
        return i, j
    a, b = np.unravel_index(np.nanargmin(bl), bl.shape)
    return i0 + a, j0 + b


def suavizar(xy, n=SUAVIZA):
    """Media movel: o caminho D8 e serrilhado (so anda em 8 direcoes)."""
    if len(xy) < n * 2:
        return xy
    a = np.asarray(xy, dtype=float)
    k = np.ones(n) / n
    s = np.column_stack([np.convolve(a[:, 0], k, mode="valid"),
                         np.convolve(a[:, 1], k, mode="valid")])
    return np.vstack([a[0], s, a[-1]])


def main():
    import gerar_rede_hecras as G
    alvos = [a for a in sys.argv[1:] if not a.startswith("-")]
    z, tr = ler_dem()
    zc = preencher(z)
    dirs = direcao(zc)
    inv = ~tr

    rede = G.montar_rede()
    feats = []
    for k, v in rede.items():
        if alvos and k not in alvos:
            continue
        p0 = Point(list(v["linha"].coords)[0])          # cabeceira da ANA
        col, lin = inv * (p0.x, p0.y)
        i, j = int(lin), int(col)
        if not (0 <= i < zc.shape[0] and 0 <= j < zc.shape[1]):
            print(f"  {k}: cabeceira fora do DEM")
            continue
        i, j = encaixar(zc, i, j)
        cam = caminhar(dirs, i, j)
        if len(cam) < 10:
            print(f"  {k}: caminho curto demais ({len(cam)} celulas)")
            continue
        xy = [tr * (c + 0.5, r + 0.5) for r, c in cam]
        # O caminhamento segue ate o MAR: um afluente passa pela confluencia e
        # desce a calha principal inteira (o Sul saia com 272 km em vez de 87).
        # Corta na celula mais proxima da foz que a ANA declara.
        foz = Point(list(v["linha"].coords)[-1])
        dfoz = [foz.distance(Point(p)) for p in xy]
        corte = int(np.argmin(dfoz))
        if corte >= 10:
            xy = xy[:corte + 1]
        ln = LineString(suavizar(xy))
        # quanto o tracado do relevo se afasta do da ANA
        d = np.array([ln.distance(Point(p)) for p in
                      list(v["linha"].coords)[::20]])
        print(f"  {k:<10} {ln.length/1000:7.1f} km  "
              f"(ANA: {v['linha'].length/1000:.1f} km)   "
              f"afastamento da ANA: mediana {np.median(d):5.0f} m, "
              f"max {d.max():6.0f} m")
        feats.append({"rio": k, "nome": v["nome"],
                      "km": round(ln.length / 1000, 2), "geometry": ln})

    if not feats:
        raise SystemExit("nenhum eixo tracado")
    gpd.GeoDataFrame(feats, crs="EPSG:31982").to_file(SAIDA, driver="GeoJSON")
    print(f"\n[OK] {SAIDA}  ({len(feats)} eixos)")


if __name__ == "__main__":
    main()
