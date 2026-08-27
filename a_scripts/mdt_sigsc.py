# -*- coding: utf-8 -*-
"""
Módulo de Amostragem do MDT SIG-SC 1m (Direct GeoTIFF Sampler).
Parte dos scripts independentes do Antigravity (a_scripts).
"""
import glob
import math
import os
import rasterio
import numpy as np

PASTA_SIGSC = r"C:\Users\haas\Downloads\sigsc"

class MosaicoSigsc:
    def __init__(self, pasta=PASTA_SIGSC):
        self.pasta = pasta
        self._abertos = {}
        self._indexar_folhas()

    def _indexar_folhas(self):
        arquivos = glob.glob(os.path.join(self.pasta, "*.tif"))
        self.folhas = {}
        for arq in arquivos:
            nome = os.path.basename(arq)
            try:
                partes = nome.replace(".tif", "").split("_")
                easting = float(partes[1])
                northing = float(partes[2])
                self.folhas[(int(easting), int(northing))] = arq
            except (IndexError, ValueError):
                continue

    def _obter_dataset(self, arq):
        if arq not in self._abertos:
            self._abertos[arq] = rasterio.open(arq)
        return self._abertos[arq]

    def cota(self, xs, ys):
        """Amostra elevações para listas de coordenadas UTM (SIRGAS 2000 UTM 22S)."""
        xs = np.asarray(xs, dtype=np.float64)
        ys = np.asarray(ys, dtype=np.float64)
        cotas = np.full(xs.shape, np.nan, dtype=np.float64)

        if not self.folhas:
            return cotas

        blocos = {}
        for i, (x, y) in enumerate(zip(xs, ys)):
            k_x = int(math.floor(x / 5000.0) * 5000)
            k_y = int(math.floor(y / 5000.0) * 5000)
            blocos.setdefault((k_x, k_y), []).append(i)

        for (k_x, k_y), indices in blocos.items():
            arq = self.folhas.get((k_x, k_y))
            if not arq:
                # Tenta blocos vizinhos próximos se estiver na borda de 5km
                continue
            ds = self._obter_dataset(arq)
            pts = [(xs[idx], ys[idx]) for idx in indices]
            amostras = list(ds.sample(pts))
            for local_i, idx in enumerate(indices):
                val = amostras[local_i][0]
                if val != ds.nodata and np.isfinite(val) and val > -50.0:
                    cotas[idx] = val

        return cotas

    def fechar(self):
        for ds in self._abertos.values():
            try:
                ds.close()
            except Exception:
                pass
        self._abertos.clear()
