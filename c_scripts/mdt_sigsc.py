# -*- coding: utf-8 -*-
"""Amostrador do MDT do SIG-SC a 1 m, direto sobre os tiles, sem mosaico.

O SIG-SC vem em folhas de carta (`MDT_SG-22-*.tif`), 1019 delas, 121 GB. Nao
da para juntar num raster so -- o dominio do Itajai-Mirim sozinho tem 64 x 44
km, que a 1 m sao 2,8 bilhoes de pixels. Aqui os tiles ficam fechados e cada
consulta e despachada para a folha que contem o ponto.

DUAS ARMADILHAS DO DADO, conferidas:

  NODATA NAO DECLARADO. 983 dos 1019 tiles trazem `nodata=None` e 36 trazem
  0.0. Medido nas folhas do dominio, 94 a 96% dos pixels sao exatamente 0.00:
  o vazio E o zero, e sem tratar isso o terreno "existe" em toda parte e vale
  zero metro. Por isso `0.0` e sempre tratado como vazio -- o que custa perder
  o nivel do mar exato, irrelevante num rio cujo leito esta acima de 4 m.

  CRS E RESOLUCAO SAO UNIFORMES. Os 1019 tiles estao em EPSG:31982 a 1 x 1 m,
  float32 -- o mesmo CRS do modelo. Nao ha reprojecao no caminho.
"""
import glob
import os

import numpy as np
import rasterio

PASTA = r"C:\Users\haas\Downloads\sigsc"


class MosaicoSigsc:
    def __init__(self, tiles=None, pasta=PASTA, zero_e_vazio=True):
        self.zero_e_vazio = zero_e_vazio
        self.caminhos = tiles or sorted(glob.glob(os.path.join(pasta, "*.tif")))
        self.bounds = []
        for p in self.caminhos:
            with rasterio.open(p) as s:
                self.bounds.append((s.bounds.left, s.bounds.bottom,
                                    s.bounds.right, s.bounds.top))
        self.bounds = np.array(self.bounds)
        self._abertos = {}

    def _src(self, i):
        if i not in self._abertos:
            self._abertos[i] = rasterio.open(self.caminhos[i])
        return self._abertos[i]

    def cobertura(self):
        b = self.bounds
        return (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())

    def cota(self, xs, ys):
        """Cota de cada ponto, NaN onde nao ha dado."""
        xs = np.atleast_1d(np.asarray(xs, float))
        ys = np.atleast_1d(np.asarray(ys, float))
        z = np.full(xs.shape, np.nan)
        b = self.bounds
        for i in range(len(self.caminhos)):
            m = ((xs >= b[i, 0]) & (xs < b[i, 2])
                 & (ys >= b[i, 1]) & (ys < b[i, 3]) & ~np.isfinite(z))
            if not m.any():
                continue
            s = self._src(i)
            v = np.array([q[0] for q in s.sample(
                list(zip(xs[m], ys[m])), indexes=1)], float)
            if s.nodata is not None:
                v[v == s.nodata] = np.nan
            if self.zero_e_vazio:
                v[v == 0.0] = np.nan
            v[v < -1000] = np.nan
            z[m] = v
        return z

    def fechar(self):
        for s in self._abertos.values():
            s.close()
        self._abertos.clear()


def tiles_do_dominio(bbox, pasta=PASTA):
    """Folhas que intersectam (xmin, ymin, xmax, ymax)."""
    fora = []
    for p in sorted(glob.glob(os.path.join(pasta, "*.tif"))):
        with rasterio.open(p) as s:
            b = s.bounds
            if not (b.right < bbox[0] or b.left > bbox[2]
                    or b.top < bbox[1] or b.bottom > bbox[3]):
                fora.append(p)
    return fora
