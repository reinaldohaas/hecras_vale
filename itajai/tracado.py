# -*- coding: utf-8 -*-
"""
Eixo dos rios derivado do RELEVO.

O tracado da ANA (BHO 2017) e cartografico, a 1:100.000: passa PERTO do rio,
nao no fundo do vale. Seccao perpendicular a um eixo torto cruza a calha de
esguelha e encontra o leito antigo ou o meandro vizinho -- e o que embaralhava
as bank lines, e nao ha ajuste de janela que conserte, porque o defeito esta no
eixo.

O tracado do relevo passa onde a agua acumula: o talvegue fica no eixo por
construcao, e a secao perpendicular cruza a calha em angulo reto.

Da ANA continua vindo o que ela tem de bom -- a TOPOLOGIA (quem desagua em
quem) e a area de drenagem. Ver itajai/topologia.py.
"""
import heapq

import numpy as np
import rasterio
from shapely.geometry import LineString, Point

from . import terreno

RES = 90.0        # o eixo e suavizado depois; 90 m basta e cabe em numpy puro
EPS = 0.001       # gradiente imposto nas areas planas (m por celula)
SUAVIZA = 5       # janela da media movel, em celulas

VIZ = [(-1, 0), (1, 0), (0, -1), (0, 1),
       (-1, -1), (-1, 1), (1, -1), (1, 1)]
DIST = [1.0, 1.0, 1.0, 1.0, 1.4142, 1.4142, 1.4142, 1.4142]


def _grade():
    """DEM em UTM reamostrado por MINIMO em blocos.

    O minimo, nao a media: com a media o talvegue some sob as encostas e o
    caminhamento sobe para o divisor. (Resampling.min do rasterio so vale para
    warp, nao para leitura.)
    """
    caminho = terreno.preparar_utm()
    with rasterio.open(caminho) as d:
        f = int(round(RES / abs(d.transform.a)))
        a = d.read(1).astype(np.float64)
        h = (a.shape[0] // f) * f
        w = (a.shape[1] // f) * f
        with np.errstate(all="ignore"):
            z = np.nanmin(a[:h, :w].reshape(h // f, f, w // f, f), axis=(1, 3))
        tr = d.transform * d.transform.scale(f, f)
        nod = d.nodata
    if nod is not None:
        z[z == nod] = np.nan
    z[z < -500] = np.nan
    return z, tr


def preencher(z):
    """Priority-flood com gradiente minimo (Barnes, Lehman & Mulla 2014).

    Cada celula sai da fila uma unica vez -- O(n log n) em um passe. O
    preenchimento iterativo vetorizado precisaria de uma iteracao por celula de
    comprimento da area plana, e no baixo Itajai isso seriam milhares.

    O gradiente minimo e o que permite atravessar a lamina achatada do
    Copernicus: sem ele o caminhamento para no primeiro plano d'agua.
    """
    h, w = z.shape
    cheio = np.full((h, w), np.inf)
    visto = np.zeros((h, w), bool)
    val = np.isfinite(z)
    borda = val & ~(
        np.pad(val[1:, :], ((0, 1), (0, 0))) &
        np.pad(val[:-1, :], ((1, 0), (0, 0))) &
        np.pad(val[:, 1:], ((0, 0), (0, 1))) &
        np.pad(val[:, :-1], ((0, 0), (1, 0))))
    fila = []
    for i, j in zip(*np.nonzero(borda)):
        cheio[i, j] = z[i, j]
        visto[i, j] = True
        heapq.heappush(fila, (float(z[i, j]), int(i), int(j)))
    while fila:
        zc, i, j = heapq.heappop(fila)
        for di, dj in VIZ:
            a, b = i + di, j + dj
            if a < 0 or a >= h or b < 0 or b >= w or visto[a, b] or not val[a, b]:
                continue
            cheio[a, b] = max(z[a, b], zc + EPS)
            visto[a, b] = True
            heapq.heappush(fila, (float(cheio[a, b]), a, b))
    cheio[~val] = np.nan
    return cheio


def direcao(zc):
    """D8: indice do vizinho de maior declive; -1 se nao houver."""
    h, w = zc.shape
    melhor = np.full((h, w), -1, np.int8)
    decl = np.zeros((h, w))
    for k, ((di, dj), dist) in enumerate(zip(VIZ, DIST)):
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


def _encaixar(zc, i, j, raio=8):
    """Puxa o ponto para a celula mais baixa da vizinhanca -- o talvegue.

    A cabeceira declarada cai na encosta com frequencia; comecando ali, o D8
    desce a ladeira ate o rio e o eixo ganha um trecho que nao e rio.
    """
    h, w = zc.shape
    i0, i1 = max(0, i - raio), min(h, i + raio + 1)
    j0, j1 = max(0, j - raio), min(w, j + raio + 1)
    bl = zc[i0:i1, j0:j1]
    if not np.isfinite(bl).any():
        return i, j
    a, b = np.unravel_index(np.nanargmin(bl), bl.shape)
    return i0 + a, j0 + b


def _suavizar(xy, n=SUAVIZA):
    """Media movel: o caminho D8 e serrilhado, so anda em 8 direcoes."""
    a = np.asarray(xy, float)
    if len(a) < n * 2:
        return a
    k = np.ones(n) / n
    s = np.column_stack([np.convolve(a[:, 0], k, "valid"),
                         np.convolve(a[:, 1], k, "valid")])
    return np.vstack([a[0], s, a[-1]])


class Tracador:
    """Prepara a grade uma vez; traca quantos rios forem pedidos."""

    def __init__(self):
        z, self.tr = _grade()
        self.zc = preencher(z)
        self.dirs = direcao(self.zc)
        self.inv = ~self.tr

    def eixo(self, cabeceira, foz, limite=200000):
        """Caminha da cabeceira ate a foz seguindo o relevo.

        O caminhamento so para no mar, entao um afluente desceria a calha
        principal inteira depois da confluencia (o Itajai do Sul saia com
        272 km em vez de 87). Por isso o corte na celula mais proxima da foz.
        """
        col, lin = self.inv * (cabeceira.x, cabeceira.y)
        i, j = _encaixar(self.zc, int(lin), int(col))
        h, w = self.zc.shape
        cam = [(i, j)]
        for _ in range(limite):
            k = self.dirs[i, j]
            if k < 0:
                break
            di, dj = VIZ[k]
            i, j = i + di, j + dj
            if i < 0 or i >= h or j < 0 or j >= w:
                break
            cam.append((i, j))
        if len(cam) < 10:
            return None
        xy = [self.tr * (c + 0.5, r + 0.5) for r, c in cam]
        d = [foz.distance(Point(p)) for p in xy]
        corte = int(np.argmin(d))
        if corte >= 10:
            xy = xy[:corte + 1]
        return LineString(_suavizar(xy))
