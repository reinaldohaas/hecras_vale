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
    terreno_hdf()


def terreno_hdf():
    """Converte para o .hdf do RAS Mapper, sem passar pela interface.

    Eu vinha dizendo que este passo so existia em Project > New Terrain. Nao e
    o caso: o ras-commander chama o RasProcess.exe CreateTerrain, que e o mesmo
    construtor que a GUI usa -- gera as 7 piramides, o TIN e o armazenamento em
    tiles. Sem este .hdf o HEC-RAS nao calcula profundidade nenhuma, porque
    profundidade e cota d'agua menos terreno.

    Pega: o RasProcess exige um .prj ESRI. Apontar o nosso .projection falha em
    100% do progresso com "Referencia de objeto nao definida" -- erro do .NET
    que nao diz nada sobre a causa. Por isso a copia com a extensao certa.
    """
    try:
        from ras_commander import RasTerrain
    except ImportError:
        print("\n  (ras-commander ausente: pip install ras-commander)")
        return
    prj = os.path.join(PASTA, "projecao.prj")
    origem = None
    for c in ("Itajai_Rede_1983.projection", "Itajai_Rede.projection"):
        if os.path.exists(c):
            origem = c
            break
    if origem is None:
        print("\n  (nenhum .projection encontrado; rode o gerador antes)")
        return
    with open(origem, encoding="utf-8") as f:
        wkt = f.read()
    with open(prj, "w", encoding="utf-8") as f:
        f.write(wkt)
    saida = os.path.join(PASTA, "Terreno.hdf")
    try:
        RasTerrain.create_terrain_hdf(
            input_rasters=[os.path.abspath(SAIDA)],
            output_hdf=os.path.abspath(saida),
            projection_prj=os.path.abspath(prj), units="Meters")
    except Exception as e:
        print(f"\n  [ERRO] CreateTerrain: {e}")
        return
    print(f"[OK] {saida}  ({os.path.getsize(saida)/1e6:.2f} MB)")
    print("     o .rasmap passa a declarar este terreno automaticamente")


if __name__ == "__main__":
    main()
