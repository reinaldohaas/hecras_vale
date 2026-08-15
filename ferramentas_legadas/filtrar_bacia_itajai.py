"""
FILTRAGEM TOPOLÓGICA DA BACIA DO ITAJAÍ PELA EXUTÓRIA NO PORTO DE ITAJAÍ
========================================================================
Este script seleciona APENAS os rios e tributários pertencentes à Bacia
do Rio Itajaí, cuja foz/exutória deságua no Porto de Itajaí (lat -26.91, lon -48.65).

Elimina rios das bacias vizinhas (Itapocu ao norte, Tijucas ao sul, Canoas a oeste).
Gera a figura 'figuras/figura_rede_filtrada_itajai.png' para conferência visual.
"""

import json
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from shapely.geometry import LineString, MultiLineString, Point, shape, mapping
from shapely.strtree import STRtree

# Coordenada aproximada da Foz / Porto de Itajaí (Exutória)
FOZ_PORTO_ITAJAI = Point(-48.652, -26.910)

def main():
    print("=" * 60)
    print("FILTRANDO BACIA DO ITAJAÍ PELA EXUTÓRIA DO PORTO DE ITAJAÍ")
    print("=" * 60)

    input_file = "rios_itajai.geojson"
    if not os.path.exists(input_file):
        input_file = "vale_itajai_full_network.geojson"

    with open(input_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    print(f"Total de segmentos na área bruta (incluindo bacias vizinhas): {len(data['features'])}")

    # Prepara geometrias e propriedades
    records = []
    geoms = []
    for feat in data["features"]:
        g = shape(feat["geometry"])
        if g.is_empty:
            continue
        records.append({
            "name": feat["properties"].get("name", "Sem Nome"),
            "properties": feat["properties"],
            "geom": g
        })
        geoms.append(g)

    # Constrói árvore espacial para busca rápida de conectividade
    tree = STRtree(geoms)

    # 1. Encontra o segmento do Rio Itajaí-Açu na Foz / Porto de Itajaí
    best_dist = 1.0
    start_idx = None

    for idx, rec in enumerate(records):
        d = rec["geom"].distance(FOZ_PORTO_ITAJAI)
        name_low = rec["name"].lower()
        if "itajai" in name_low or "acu" in name_low:
            if d < best_dist:
                best_dist = d
                start_idx = idx

    if start_idx is None:
        # Fallback para o segmento mais próximo da foz
        dists = [g.distance(FOZ_PORTO_ITAJAI) for g in geoms]
        start_idx = int(np.argmin(dists))

    print(f"Segmento de Foz identificado: '{records[start_idx]['name']}' (distância da foz: {best_dist*111.0:.2f} km)")

    # 2. Busca em largura (BFS) para selecionar APENAS os rios conectados à foz
    visited_indices = set([start_idx])
    queue = [start_idx]

    # Tolerância de desconexão espacial (~250 metros em graus)
    SNAP_TOL_DEG = 0.0025

    while queue:
        curr_idx = queue.pop(0)
        curr_geom = geoms[curr_idx]
        
        # Encontra candidatos vizinhos via bounding box
        possible_matches = tree.query(curr_geom)
        
        for match_idx in possible_matches:
            if match_idx not in visited_indices:
                other_geom = geoms[match_idx]
                if curr_geom.distance(other_geom) <= SNAP_TOL_DEG:
                    visited_indices.add(match_idx)
                    queue.append(match_idx)

    filtered_features = []
    for idx in sorted(visited_indices):
        rec = records[idx]
        filtered_features.append({
            "type": "Feature",
            "properties": rec["properties"],
            "geometry": mapping(rec["geom"])
        })

    print(f"✓ Selecionados APENAS {len(filtered_features)} rios/segmentos pertencentes à Bacia do Itajaí!")

    # Salva o GeoJSON filtrado limpo
    output_geojson = "bacia_itajai_exutoria_porto.geojson"
    out_fc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
        "features": filtered_features
    }

    with open(output_geojson, "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)
    
    # Atualiza rios_itajai.geojson e vale_itajai_full_network.geojson com a bacia limpa
    with open("rios_itajai.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)
    with open("vale_itajai_full_network.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)

    print(f"✓ Arquivos 'rios_itajai.geojson' e 'vale_itajai_full_network.geojson' atualizados com a bacia limpa!")

    # 3. GERA A FIGURA DE VERIFICAÇÃO VISUAL DA BACIA FILTRADA
    fig, ax = plt.subplots(figsize=(11, 9), dpi=300)
    
    # Plota rios filtrados em azul hidrológico
    main_rivers_plotted = set()
    for feat in filtered_features:
        name = feat["properties"].get("name", "")
        name_low = name.lower()
        
        g = shape(feat["geometry"])
        
        lw = 1.0
        c = '#4a90e2' # Azul afluentes
        
        if "itajai-acu" in name_low or "itajai acu" in name_low or (name == "Rio Itajaí" and g.length > 0.5):
            c = '#003366' # Azul escuro principal
            lw = 2.5
        elif "mirim" in name_low:
            c = '#8e44ad' # Roxo Mirim
            lw = 1.8
        elif "sul" in name_low:
            c = '#27ae60' # Verde Sul
            lw = 1.8
        elif "oeste" in name_low:
            c = '#e67e22' # Laranja Oeste
            lw = 1.8
        elif "norte" in name_low or "hercilio" in name_low:
            c = '#c0392b' # Vermelho Norte
            lw = 1.8

        if g.geom_type == "LineString":
            xs, ys = g.xy
            ax.plot(xs, ys, color=c, linewidth=lw, alpha=0.85)
        elif g.geom_type == "MultiLineString":
            for ls in g.geoms:
                xs, ys = ls.xy
                ax.plot(xs, ys, color=c, linewidth=lw, alpha=0.85)

    # Marca a Exutória do Porto de Itajaí com um destaque vermelho especial
    ax.plot(FOZ_PORTO_ITAJAI.x, FOZ_PORTO_ITAJAI.y, '*', color='#e74c3c', markersize=18, markeredgecolor='black', zorder=10, label='Exutória: Porto de Itajaí')
    
    ax.annotate('FOZ DO RIO ITAJAÍ-AÇU\n(Porto de Itajaí / Oceano)', 
                (FOZ_PORTO_ITAJAI.x, FOZ_PORTO_ITAJAI.y), 
                xytext=(FOZ_PORTO_ITAJAI.x + 0.08, FOZ_PORTO_ITAJAI.y - 0.05),
                arrowprops=dict(facecolor='#e74c3c', shrink=0.08, width=2, headwidth=8),
                fontsize=10, fontweight='bold', color='#900c3f',
                bbox=dict(boxstyle="round,pad=0.3", fc="#fbeee6", ec="#e74c3c", lw=1.5))

    # Adiciona cidades de referência da bacia
    cities = [
        (-49.60, -27.28, 'Taió (Oeste)'),
        (-49.60, -27.18, 'Ituporanga (Sul)'),
        (-49.65, -26.92, 'Doutor Pedrinho / Norte'),
        (-49.64, -27.21, 'Rio do Sul (Confluência)'),
        (-49.07, -26.92, 'Blumenau'),
        (-48.90, -27.11, 'Brusque (Mirim)'),
    ]
    for cx, cy, cname in cities:
        ax.plot(cx, cy, 'o', color='#f39c12', markersize=6, markeredgecolor='black', zorder=8)
        ax.text(cx + 0.02, cy + 0.01, cname, fontsize=8, fontweight='bold', color='#2c3e50',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#bdc3c7", alpha=0.8))

    ax.set_title('Rede Hidrográfica Filtrada: Apenas Rios da Bacia do Itajaí (Exutória no Porto de Itajaí)', fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Longitude (Graus WGS84)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Latitude (Graus WGS84)', fontsize=10, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Salva a figura em figuras/figura_rede_filtrada_itajai.png
    output_dir = Path("figuras").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    target_fig = output_dir / "figura_rede_filtrada_itajai.png"
    fig.savefig(str(target_fig), dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"✓ FIGURA DE CONFERÊNCIA VISUAL SALVA EM: {target_fig}")

if __name__ == "__main__":
    main()
