#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, MultiLineString, Point

# ======================================================
# Área aproximada da bacia do Itajaí
# (S,W,N,E)
# ======================================================

south = -28.1
west  = -50.9
north = -26.3
east  = -48.4

# ======================================================
# Consulta Overpass
# ======================================================

query = f"""
[out:json][timeout:900];

(
  way["waterway"~"river|stream"]({south},{west},{north},{east});
);

out geom;
"""

print("Baixando dados...")

r = requests.post(
    "https://overpass.kumi.systems/api/interpreter",
    data=query
)

r.raise_for_status()

data = r.json()

print("Número de segmentos:", len(data["elements"]))

# ======================================================
# Converte para GeoDataFrame
# ======================================================

features = []

for e in data["elements"]:

    if "geometry" not in e:
        continue

    coords = [(p["lon"], p["lat"]) for p in e["geometry"]]

    if len(coords) < 2:
        continue

    features.append({
        "osm_id": e["id"],
        "name": e.get("tags", {}).get("name", ""),
        "geometry": LineString(coords)
    })

gdf = gpd.GeoDataFrame(features, crs="EPSG:4326")

print("Linhas:", len(gdf))

# ======================================================
# Monta o grafo
# ======================================================

G = nx.Graph()

for idx, row in gdf.iterrows():

    line = row.geometry

    a = tuple(line.coords[0])
    b = tuple(line.coords[-1])

    G.add_edge(a, b, index=idx)

print("Nós:", G.number_of_nodes())
print("Arestas:", G.number_of_edges())

# ======================================================
# Procura o Rio Itajaí-Açu
# ======================================================

principal = gdf[
    gdf["name"].str.contains("Itaja", case=False, na=False)
]

if len(principal) == 0:
    raise Exception("Rio Itajaí não encontrado.")

print(principal[["name"]])

# usa o primeiro segmento encontrado

linha = principal.iloc[0].geometry

inicio = tuple(linha.coords[0])
fim = tuple(linha.coords[-1])

# ======================================================
# Componente conectada
# ======================================================

comp1 = nx.node_connected_component(G, inicio)
comp2 = nx.node_connected_component(G, fim)

component = comp1 if len(comp1) > len(comp2) else comp2

print("Nós conectados:", len(component))

# ======================================================
# Seleciona linhas
# ======================================================

indices = []

for u, v, d in G.edges(data=True):

    if u in component and v in component:
        indices.append(d["index"])

rede = gdf.loc[sorted(set(indices))]

print("Segmentos finais:", len(rede))

# ======================================================
# Salva
# ======================================================

rede.to_file(
    "itajai_rede.geojson",
    driver="GeoJSON"
)

print("Arquivo salvo:")
print("itajai_rede.geojson")