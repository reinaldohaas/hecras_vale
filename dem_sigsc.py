# -*- coding: utf-8 -*-
"""
Amostrador do MDT do SIG-SC (1 m) para a geracao das secoes.

O levantamento aerofotogrametrico de Santa Catarina cobre a bacia do Itajai
inteira com Modelo Digital de TERRENO a 1 m, em EPSG:31982 -- o mesmo CRS do
modelo, sem reprojecao. Contra o Copernicus GLO-30 usado ate aqui:

    Copernicus GLO-30   30 m   modelo de SUPERFICIE   ~155 MB
    SIG-SC MDT           1 m   modelo de TERRENO      ~98 GB (829 tiles)

Sao 900x mais celulas por area, e terreno em vez de superficie -- some a
vegetacao, que no Copernicus levanta as margens artificialmente.

Como sao 98 GB, nao da para carregar tudo. Esta classe monta um indice das
extensoes pelos arquivos .tfw e carrega os tiles sob demanda, com cache LRU.
Como as secoes seguem o rio, tiles consecutivos se repetem e o cache resolve.

Nao usa leitura por janela nem ds.sample(): as duas derrubam o processo no
rasterio 1.5.x deste ambiente. Le o tile inteiro (0,2 s para 126 MB) e
indexa com numpy.
"""
import glob
import os
from collections import OrderedDict

import numpy as np
import rasterio

DIR_PADRAO = r"C:\Users\haas\Downloads\sigsc"
EPSG = 31982


class DemHibrido:
    """SIG-SC onde ha cobertura, DEM de reserva no resto.

    Os 829 tiles baixados nao formam um retangulo solido: o Itajai do Sul,
    por exemplo, esta apenas 5% coberto (falta a folha de Ituporanga),
    enquanto os outros cinco rios estao a 100%. Sem reserva, aquele rio fica
    sem secoes e some do modelo -- foi o que aconteceu: 5 secoes em 87,3 km.
    """

    def __init__(self, sigsc, reserva):
        self.sig = sigsc
        self.res = reserva
        self.n_sig = 0
        self.n_res = 0
        self.nome = f"{sigsc.nome} + reserva"

    def sample(self, xs, ys):
        z = self.sig.sample(xs, ys)
        falta = np.isnan(z)
        if falta.any():
            z2 = self.res.sample(np.asarray(xs)[falta], np.asarray(ys)[falta])
            z[falta] = z2
            self.n_res += int(falta.sum())
        self.n_sig += int(z.size - falta.sum())
        return z


class DemSIGSC:
    def __init__(self, pasta=DIR_PADRAO, cache=6):
        self.tiles = []
        for tfw in glob.glob(os.path.join(pasta, "*.tfw")):
            tif = tfw[:-4] + ".tif"
            if not os.path.exists(tif):
                continue
            v = [float(x) for x in open(tfw).read().split()]
            px, py, x0, y0 = v[0], v[3], v[4], v[5]
            # .tfw da o CENTRO do primeiro pixel; a borda fica meio pixel antes
            x0 -= px / 2.0
            y0 -= py / 2.0
            self.tiles.append({"tif": tif, "px": px, "py": py, "x0": x0, "y0": y0,
                               "w": None, "h": None})
        if not self.tiles:
            raise FileNotFoundError(f"nenhum tile em {pasta}")
        # As dimensoes VARIAM entre tiles (6463 x 4862, 6445 x ...). Assumir
        # iguais faz o ponto cair no tile errado e estourar o indice. Abrir o
        # .tif so le o cabecalho, entao indexar os 829 e barato.
        cache_idx = os.path.join(pasta, "_indice_tiles.npz")
        if os.path.exists(cache_idx):
            z = np.load(cache_idx, allow_pickle=True)
            dims = {str(k): v for k, v in zip(z["nomes"], z["dims"])}
        else:
            dims = {}
            for t in self.tiles:
                with rasterio.open(t["tif"]) as d:
                    dims[os.path.basename(t["tif"])] = (d.width, d.height)
            np.savez(cache_idx,
                     nomes=np.array(list(dims.keys())),
                     dims=np.array(list(dims.values())))
        for t in self.tiles:
            w, h = dims[os.path.basename(t["tif"])]
            t["w"], t["h"] = int(w), int(h)
            t["x1"] = t["x0"] + t["w"] * t["px"]
            t["y1"] = t["y0"] + t["h"] * t["py"]      # py e negativo
        self.cache = OrderedDict()
        self.cache_max = cache
        self.lidos = 0
        self.nome = f"SIG-SC MDT 1 m ({len(self.tiles)} tiles)"

    # ---------------------------------------------------------------- interno
    def _tile_de(self, x, y):
        for t in self.tiles:
            if t["x0"] <= x < t["x1"] and t["y1"] <= y < t["y0"]:
                return t
        return None

    def _array(self, t):
        k = t["tif"]
        if k in self.cache:
            self.cache.move_to_end(k)
            return self.cache[k]
        with rasterio.open(k) as d:
            a = d.read(1)
            nod = d.nodata
        a = a.astype(np.float32)
        if nod is not None:
            a[a == nod] = np.nan
        a[a < -100] = np.nan
        # ZERO EXATO E VAZIO, nao cota. O laser nao retorna da lamina d'agua,
        # e o MDT preenche esses pixels com 0,00. Tratar como cota faz o
        # talvegue despencar: a 145,9 km da foz o Copernicus da 114,00 m e o
        # SIG-SC da 0,00 -- 114 m de erro no meio do vale. Marcados como
        # ausentes, sao preenchidos por interpolacao na secao (e a calha e
        # escavada em seguida, que e o tratamento correto para o canal).
        a[a == 0.0] = np.nan
        self.cache[k] = a
        self.lidos += 1
        if len(self.cache) > self.cache_max:
            self.cache.popitem(last=False)
        return a

    # ----------------------------------------------------------------- publico
    def sample(self, xs, ys):
        """Cota nas coordenadas UTM 22S dadas. NaN fora da cobertura."""
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        out = np.full(xs.shape, np.nan, dtype=float)
        restante = np.ones(xs.shape, dtype=bool)
        while restante.any():
            i = int(np.flatnonzero(restante)[0])
            t = self._tile_de(xs[i], ys[i])
            if t is None:
                restante[i] = False
                continue
            # todos os pontos que caem NESTE tile, de uma vez
            m = (restante & (xs >= t["x0"]) & (xs < t["x1"])
                 & (ys >= t["y1"]) & (ys < t["y0"]))
            a = self._array(t)
            col = ((xs[m] - t["x0"]) / t["px"]).astype(int)
            row = ((ys[m] - t["y0"]) / t["py"]).astype(int)
            np.clip(col, 0, t["w"] - 1, out=col)
            np.clip(row, 0, t["h"] - 1, out=row)
            out[m] = a[row, col]
            restante[m] = False
        return out

    def cobre(self, x, y):
        return self._tile_de(float(x), float(y)) is not None


if __name__ == "__main__":
    d = DemSIGSC()
    print(d.nome)
    xs = np.linspace(660000, 660500, 20)
    ys = np.full(20, 6980000.0)
    z = d.sample(xs, ys)
    print("amostra:", np.round(z, 2))
    print(f"tiles lidos: {d.lidos}")
