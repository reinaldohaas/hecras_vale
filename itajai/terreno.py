# -*- coding: utf-8 -*-
"""
Relevo da bacia. UMA fonte: Copernicus GLO-30.

Misturar fontes desencontra a geometria -- o eixo sai do fundo de vale que uma
enxerga e a secao e amostrada noutra, cujo talvegue esta em lugar diferente. O
SIG-SC de 1 m e melhor terreno, mas volta depois, como troca isolada e
testavel contra esta referencia.

Duas representacoes, e a distincao importa:

    bruto()   EPSG:4326, o arquivo como baixado. So para amostrar cotas.
    utm()     EPSG:31982, 30 m. Para tracar rio e para o RAS Mapper.
"""
import os

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
from pyproj import Transformer

EPSG = 31982                   # SIRGAS 2000 / UTM 22S
BRUTO = "dem_itajai.tif"       # Copernicus GLO-30, como baixado (EPSG:4326)
PASTA = "Terrain"
UTM = os.path.join(PASTA, "Terreno_Copernicus.tif")
RES = 30.0


class Amostrador:
    """Le a banda inteira e indexa com numpy.

    ds.sample() e a leitura por janela derrubam o processo no rasterio 1.5.x
    deste ambiente -- crash nativo, sem traceback. Ler tudo custa 150 MB e
    resolve.
    """

    def __init__(self, caminho=None):
        self.caminho = caminho or BRUTO
        self.ds = rasterio.open(self.caminho)
        self.arr = self.ds.read(1)
        self.linhas, self.colunas = self.arr.shape
        self.nodata = self.ds.nodata
        self.inv = ~self.ds.transform
        self.tr = Transformer.from_crs(EPSG, self.ds.crs.to_epsg(),
                                       always_xy=True)
        self.nome = f"Copernicus GLO-30 ({os.path.basename(self.caminho)})"

    def cota(self, xs, ys):
        """Cota nas coordenadas UTM dadas. NaN fora da cobertura."""
        lon, lat = self.tr.transform(np.asarray(xs, float), np.asarray(ys, float))
        a, b, c, d, e, f = (self.inv.a, self.inv.b, self.inv.c,
                            self.inv.d, self.inv.e, self.inv.f)
        col = np.floor(a * lon + b * lat + c).astype(int)
        lin = np.floor(d * lon + e * lat + f).astype(int)
        ok = ((lin >= 0) & (lin < self.linhas)
              & (col >= 0) & (col < self.colunas))
        out = np.full(np.shape(lon), np.nan)
        out[ok] = self.arr[lin[ok], col[ok]]
        if self.nodata is not None:
            out[out == self.nodata] = np.nan
        out[out < -500] = np.nan
        return out

    def talvegue(self, xs, ys, raio=45.0, n=3):
        """Cota MINIMA numa janela em volta de cada ponto.

        O tracado nao passa exatamente no fundo; a cota pontual pega a margem.
        """
        xs = np.asarray(xs, float)
        ys = np.asarray(ys, float)
        d = np.linspace(-raio, raio, n)
        z = np.full((len(d) * len(d), xs.size), np.nan)
        k = 0
        for dx in d:
            for dy in d:
                z[k] = self.cota(xs + dx, ys + dy)
                k += 1
        with np.errstate(all="ignore"):
            return np.nanmin(z, axis=0)


def preparar_utm(forcar=False):
    """Reprojeta o Copernicus para UTM 22S a 30 m, com piramides.

    O RAS Mapper exige o terreno na projecao do projeto; em graus o import
    falha ou sai deslocado. As piramides evitam que ele redesenhe o raster
    inteiro a cada zoom.
    """
    os.makedirs(PASTA, exist_ok=True)
    if os.path.exists(UTM) and not forcar:
        return UTM
    with rasterio.open(BRUTO) as src:
        dst = f"EPSG:{EPSG}"
        tr, w, h = calculate_default_transform(
            src.crs, dst, src.width, src.height, *src.bounds,
            resolution=(RES, RES))
        perfil = src.profile.copy()
        perfil.update(crs=dst, transform=tr, width=w, height=h,
                      dtype="float32", count=1, nodata=-9999.0,
                      compress="deflate", predictor=3, tiled=True,
                      blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
        with rasterio.open(UTM, "w", **perfil) as out:
            reproject(rasterio.band(src, 1), rasterio.band(out, 1),
                      src_transform=src.transform, src_crs=src.crs,
                      dst_transform=tr, dst_crs=dst,
                      src_nodata=src.nodata, dst_nodata=-9999.0,
                      resampling=Resampling.bilinear)
    with rasterio.open(UTM, "r+") as d:
        d.build_overviews([2, 4, 8, 16, 32], Resampling.average)
    return UTM


def preparar_hdf(projection_wkt, nome="Terreno"):
    """Terreno no formato do RAS Mapper, via RasProcess.exe CreateTerrain.

    Nao precisa da interface: e o mesmo construtor que Project > New Terrain
    usa. Sem este .hdf o HEC-RAS nao calcula profundidade, porque profundidade
    e cota d'agua menos terreno.

    O RasProcess exige um .prj ESRI -- apontar um .projection roda ate 100% e
    falha com "Referencia de objeto nao definida", erro do .NET que nao diz
    nada sobre a causa.
    """
    from ras_commander import RasTerrain
    preparar_utm()
    prj = os.path.join(PASTA, "projecao.prj")
    with open(prj, "w", encoding="utf-8") as f:
        f.write(projection_wkt)
    saida = os.path.join(PASTA, f"{nome}.hdf")
    RasTerrain.create_terrain_hdf(
        input_rasters=[os.path.abspath(UTM)],
        output_hdf=os.path.abspath(saida),
        projection_prj=os.path.abspath(prj), units="Meters")
    return saida
