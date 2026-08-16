#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT: baixar_dem_online.py
DESCRIÇÃO: Obtém dinamicamente o Modelo Digital de Elevação (DEM 30m) da bacia
           a partir da extensão geográfica da rede de rios extraída da ANA,
           adicionando uma folga de segurança para cobrir todo o vale.

REQUISITOS: requests, json, os, sys
SAÍDA: dem_<nome_bacia>.tif (GeoTIFF em WGS84 EPSG:4326)
================================================================================
"""

import os
import sys
import json
import urllib.parse
import urllib.request
import argparse
import unicodedata

# Chaves de API de fallback do OpenTopography
OPENTOPO_API_KEYS = [
    os.environ.get("OPENTOPO_API_KEY", "cbe30884d238441ae080d36328459286"),
    "demo_key",
    "guest"
]

DEM_TYPE = "COP30"  # Copernicus GLO-30 (~30m de resolução espacial)

def sanitizar_nome(nome):
    nfkd = unicodedata.normalize('NFKD', nome)
    return "".join([c for c in nfkd if not unicodedata.combining(c)]).lower().replace(" ", "_").replace("-", "_")

def calcular_bbox_geojson(caminho_geojson, folga_graus=0.15):
    """
    Lê o GeoJSON da rede de rios e calcula a Bounding Box [Sul, Oeste, Norte, Leste]
    acrescentando uma folga de segurança (folga_graus ~ 15 km).
    """
    with open(caminho_geojson, "r", encoding="utf-8") as f:
        data = json.load(f)

    lons, lats = [], []
    for feat in data.get("features", []):
        geom = feat.get("geometry")
        if not geom: continue
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])
        
        lines = [coords] if gtype == "LineString" else (coords if gtype == "MultiLineString" else [])
        for line in lines:
            for pt in line:
                if len(pt) >= 2:
                    lons.append(pt[0])
                    lats.append(pt[1])

    if not lons or not lats:
        raise ValueError(f"Nenhuma coordenada encontrada no GeoJSON: {caminho_geojson}")

    min_lon, max_lon = min(lons), max(lons)
    min_lat, max_lat = min(lats), max(lats)

    south = round(min_lat - folga_graus, 4)
    north = round(max_lat + folga_graus, 4)
    west = round(min_lon - folga_graus, 4)
    east = round(max_lon + folga_graus, 4)

    return south, west, north, east

def baixar_dem_opentopography(south, west, north, east, caminho_saida_tif):
    """Baixa o arquivo GeoTIFF de elevação do OpenTopography."""
    print(f"  i Bounding Box do Vale com folga: S={south}, W={west}, N={north}, E={east}")
    
    for api_key in OPENTOPO_API_KEYS:
        url = (
            "https://portal.opentopography.org/API/globaldem"
            f"?demtype={DEM_TYPE}&south={south}&north={north}"
            f"&west={west}&east={east}&outputFormat=GTiff&API_Key={api_key}"
        )
        print(f"  i Baixando modelo de elevação Copernicus DEM 30m ({DEM_TYPE})...")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Antigravity/2.0"})
            with urllib.request.urlopen(req, timeout=300) as response:
                if response.status == 200:
                    data = response.read()
                    with open(caminho_saida_tif, "wb") as f:
                        f.write(data)
                    tamanho_mb = len(data) / (1024 * 1024)
                    print(f"  ✓ DEM GeoTIFF baixado com sucesso: {os.path.basename(caminho_saida_tif)} ({tamanho_mb:.2f} MB)")
                    return True
        except Exception as e:
            print(f"  ⚠ Tentativa de download falhou (chave={api_key[:6]}...): {e}")

    print(f"  ⚠ Não foi possível baixar o DEM via OpenTopography. Verifique a conexão ou a chave API.")
    return False

def obter_dem_bacia(nome_bacia, pasta_trabalho, folga_graus=0.15):
    nome_limpo = sanitizar_nome(nome_bacia)
    caminho_geojson = os.path.join(pasta_trabalho, f"rios_{nome_limpo}.geojson")
    
    if not os.path.exists(caminho_geojson):
        caminho_geojson_orig = os.path.join(pasta_trabalho, f"rios_{nome_limpo}_original.geojson")
        if os.path.exists(caminho_geojson_orig):
            caminho_geojson = caminho_geojson_orig
        else:
            raise FileNotFoundError(f"Arquivo de rios não encontrado para calcular o DEM: {caminho_geojson}")

    caminho_tif = os.path.join(pasta_trabalho, f"dem_{nome_limpo}.tif")
    
    if os.path.exists(caminho_tif):
        print(f"  ✓ DEM já existe em disco: {os.path.basename(caminho_tif)}")
        return caminho_tif

    south, west, north, east = calcular_bbox_geojson(caminho_geojson, folga_graus=folga_graus)
    sucesso = baixar_dem_opentopography(south, west, north, east, caminho_tif)
    
    if sucesso:
        return caminho_tif
    else:
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Baixa o DEM 30m cobrindo a bacia com folga.")
    parser.add_argument("-b", "--bacia", type=str, default="Itajaí", help="Nome do rio/bacia (ex: Itajaí, Tijucas, Tubarão)")
    parser.add_argument("-o", "--saida", type=str, default=r"C:\Users\haas\github\hecras_vale", help="Pasta de trabalho")
    parser.add_argument("-f", "--folga", type=float, default=0.15, help="Folga de margem em graus (padrão: 0.15 ~ 15km)")
    
    args = parser.parse_args()
    obter_dem_bacia(args.bacia, args.saida, folga_graus=args.folga)
