# -*- coding: utf-8 -*-
"""
Prepara o relevo Copernicus para o RAS Mapper.

O dem_itajai.tif vem do Copernicus em EPSG:4326 -- graus. O projeto e
EPSG:31982 (SIRGAS 2000 / UTM 22S), em metros. O RAS Mapper exige o terreno na
MESMA projecao do projeto: com o raster em graus, o import de terreno falha ou
sai deslocado, que e o motivo de nunca ter aparecido relevo nenhum.

Aqui o DEM e reprojetado para UTM 22S a 30 m (a resolucao nativa do GLO-30),
gravado com compressao e piramides internas -- sem elas o RAS Mapper redesenha
o raster inteiro a cada zoom.

Uso:  python preparar_terreno.py
Saida: Terrain/Terreno_Copernicus.tif
"""
import os

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject

ORIGEM = "dem_itajai.tif"
PASTA = "Terrain"
SAIDA = os.path.join(PASTA, "Terreno_Copernicus.tif")
EPSG = 31982
RES = 30.0


def main():
    os.makedirs(PASTA, exist_ok=True)
    with rasterio.open(ORIGEM) as src:
        dst_crs = f"EPSG:{EPSG}"
        transform, w, h = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds,
            resolution=(RES, RES))
        perfil = src.profile.copy()
        perfil.update(crs=dst_crs, transform=transform, width=w, height=h,
                      dtype="float32", count=1, nodata=-9999.0,
                      compress="deflate", predictor=3, tiled=True,
                      blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
        print(f"{ORIGEM}: {src.width}x{src.height} em {src.crs}")
        print(f"  -> {w}x{h} em {dst_crs} a {RES:.0f} m")
        with rasterio.open(SAIDA, "w", **perfil) as dst:
            reproject(source=rasterio.band(src, 1),
                      destination=rasterio.band(dst, 1),
                      src_transform=src.transform, src_crs=src.crs,
                      dst_transform=transform, dst_crs=dst_crs,
                      src_nodata=src.nodata, dst_nodata=-9999.0,
                      resampling=Resampling.bilinear)
    # piramides: o RAS Mapper redesenha o raster inteiro a cada zoom sem elas
    with rasterio.open(SAIDA, "r+") as d:
        d.build_overviews([2, 4, 8, 16, 32], Resampling.average)
        d.update_tags(ns="rio_overview", resampling="average")
    with rasterio.open(SAIDA) as d:
        a = d.read(1, out_shape=(1, d.height // 20, d.width // 20))
        v = a[a > -9998]
        print(f"\n[OK] {SAIDA}  ({os.path.getsize(SAIDA)/1e6:.0f} MB)")
        print(f"     extensao {[round(x) for x in d.bounds]}")
        print(f"     cotas de {v.min():.1f} a {v.max():.1f} m")


if __name__ == "__main__":
    main()
