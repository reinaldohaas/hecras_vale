# -*- coding: utf-8 -*-
"""
Leitura e amostragem do terreno.

Duas decisoes que nao sao obvias:

BANDA INTEIRA EM MEMORIA. O caminho natural seria ds.sample() ou leitura por
janela, um ponto de cada vez. No rasterio 1.5 sobre os GeoTIFF grandes desta
bacia isso derruba o processo em codigo nativo, sem traceback em Python. Ler a
banda de uma vez e interpolar em numpy e mais rapido e nao quebra. Ha um limite
de tamanho: acima dele o programa avisa em vez de tentar e morrer.

AMOSTRAGEM NO CRS DO PROPRIO DEM. A geometria toda trabalha num CRS metrico,
mas os pontos sao convertidos para o CRS do raster na hora de amostrar. Assim
o DEM nunca precisa ser reprojetado -- reprojetar terreno para consultar cota
introduz uma reamostragem a mais e perde os NoData nas bordas.
"""
import numpy as np
import rasterio
from pyproj import CRS, Transformer

MAX_CELULAS = 400_000_000        # ~1,6 GB em float32


class DEM:
    """Um GeoTIFF de terreno, com amostragem bilinear e NoData preservado."""

    def __init__(self, caminho, banda=1):
        self.caminho = str(caminho)
        with rasterio.open(self.caminho) as ds:
            if ds.width * ds.height > MAX_CELULAS:
                raise MemoryError(
                    f"DEM com {ds.width}x{ds.height} celulas e grande demais "
                    f"para carregar inteiro ({MAX_CELULAS} celulas e o limite). "
                    f"Gere um recorte ou uma piramide antes.")
            self.crs = CRS.from_user_input(ds.crs) if ds.crs else None
            self.transform = ds.transform
            self.bounds = ds.bounds
            self.nodata = ds.nodata
            self.largura, self.altura = ds.width, ds.height
            self.res = (abs(ds.transform.a), abs(ds.transform.e))
            z = ds.read(banda, masked=True)
        self.z = np.ma.filled(z.astype("float32"), np.nan)
        if self.nodata is not None:
            self.z[self.z == self.nodata] = np.nan
        # o inverso da affine, para ir de coordenada a indice de celula
        self._inv = ~self.transform
        self._transformadores = {}

    # ------------------------------------------------------------------ CRS
    @property
    def crs_metrico(self):
        """CRS de trabalho, em metros.

        Se o DEM ja e projetado, e o dele. Se e geografico, estima o UTM do
        centro -- medir largura de secao em graus nao tem sentido.
        """
        if self.crs is None:
            return None
        if self.crs.is_projected:
            return self.crs
        x = 0.5 * (self.bounds.left + self.bounds.right)
        y = 0.5 * (self.bounds.bottom + self.bounds.top)
        zona = int((x + 180.0) // 6) + 1
        epsg = (32600 if y >= 0 else 32700) + zona
        return CRS.from_epsg(epsg)

    def _transformador(self, crs_origem):
        chave = crs_origem.to_string()
        if chave not in self._transformadores:
            self._transformadores[chave] = Transformer.from_crs(
                crs_origem, self.crs, always_xy=True)
        return self._transformadores[chave]

    # ------------------------------------------------------------ amostragem
    def cota(self, xs, ys, crs=None):
        """Cota nos pontos dados, por interpolacao bilinear.

        Devolve NaN onde ha NoData ou fora do raster -- nunca um valor
        inventado. Quem chama decide o que fazer com o buraco; o requisito e
        justamente nao mascarar terreno ausente.
        """
        xs = np.asarray(xs, dtype="float64").ravel()
        ys = np.asarray(ys, dtype="float64").ravel()
        if crs is not None and self.crs is not None and \
                CRS.from_user_input(crs) != self.crs:
            xs, ys = self._transformador(CRS.from_user_input(crs)).transform(xs, ys)

        col, lin = self._inv * (xs, ys)
        col = np.asarray(col) - 0.5           # centro da celula
        lin = np.asarray(lin) - 0.5

        c0 = np.floor(col).astype("int64")
        l0 = np.floor(lin).astype("int64")
        fc = col - c0
        fl = lin - l0

        out = np.full(xs.shape, np.nan, dtype="float64")
        ok = (c0 >= 0) & (l0 >= 0) & (c0 + 1 < self.largura) & (l0 + 1 < self.altura)
        if not np.any(ok):
            return out

        c0o, l0o, fco, flo = c0[ok], l0[ok], fc[ok], fl[ok]
        z00 = self.z[l0o, c0o]
        z01 = self.z[l0o, c0o + 1]
        z10 = self.z[l0o + 1, c0o]
        z11 = self.z[l0o + 1, c0o + 1]
        # com qualquer vizinho NoData o resultado e NoData: interpolar em cima
        # de um buraco produz uma rampa falsa exatamente na borda do dado
        v = (z00 * (1 - fco) * (1 - flo) + z01 * fco * (1 - flo) +
             z10 * (1 - fco) * flo + z11 * fco * flo)
        out[ok] = v
        return out

    def perfil_linha(self, linha, espacamento=2.0, crs=None):
        """Amostra o DEM ao longo de uma LineString.

        Devolve (estacas, cotas, xs, ys) com espacamento uniforme em metros.
        """
        L = float(linha.length)
        if L <= 0:
            return (np.zeros(0),) * 4
        n = max(int(round(L / float(espacamento))) + 1, 2)
        sta = np.linspace(0.0, L, n)
        pts = [linha.interpolate(float(s)) for s in sta]
        xs = np.array([p.x for p in pts])
        ys = np.array([p.y for p in pts])
        return sta, self.cota(xs, ys, crs=crs), xs, ys

    def resumo(self):
        val = np.isfinite(self.z)
        return {
            "caminho": self.caminho,
            "crs": self.crs.to_string() if self.crs else "(sem CRS)",
            "crs_metrico": self.crs_metrico.to_string() if self.crs else "",
            "dimensoes": f"{self.largura} x {self.altura}",
            "resolucao": f"{self.res[0]:.2f} x {self.res[1]:.2f}",
            "cotas": (f"{np.nanmin(self.z):.2f} a {np.nanmax(self.z):.2f}"
                      if val.any() else "(vazio)"),
            "nodata": f"{100.0 * (1 - val.mean()):.1f}% do raster",
        }
