"""
FILTRAGEM E VISUALIZAÇÃO DA BACIA DO RIO ITAJAÍ (EXUTÓRIA NO PORTO DE ITAJAÍ)
"""
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from shapely.geometry import shape, mapping, Point

FOZ_PORTO = Point(-48.652, -26.910)

def main():
    print("Iniciando filtragem da bacia do Itajaí...")
    with open("rios_itajai.geojson", "r", encoding="utf-8") as f:
        data = json.load(f)

    # Identifica os rios pertencentes à bacia hidrográfica do Itajaí pelas palavras-chave e conexão
    keywords = [
        "itajai", "acu", "mirim", "hercilio", "benedito", "cedros", "garuva",
        "salto", "testo", "blumenau", "brusque", "taio", "ituporanga", "boiteux",
        "pombas", "garuva", "canoinhas", "louro", "luiz", "indial", "perimbó"
    ]

    filtered_features = []
    excluded_count = 0

    for feat in data["features"]:
        name = feat["properties"].get("name", "")
        name_ascii = feat["properties"].get("name_ascii", "").lower()
        
        # Filtra rios conhecidos das bacias vizinhas (Itapocu, Tijucas, Canoas, etc.)
        if any(ex in name_ascii for ex in ["itapocu", "tijucas", "canoas", "cubatao", "biguacu"]):
            excluded_count += 1
            continue

        # Seleciona rios da rede do Itajaí
        g = shape(feat["geometry"])
        if g.is_empty:
            continue
            
        # Se contiver nome da bacia ou estiver dentro da caixa delimitadora principal do Itajaí
        bounds = g.bounds # (minx, miny, maxx, maxy)
        # BBOX da Bacia do Itajaí: lon [-49.9, -48.6], lat [-27.6, -26.7]
        if bounds[0] >= -50.1 and bounds[2] <= -48.5 and bounds[1] >= -27.65 and bounds[3] <= -26.55:
            filtered_features.append(feat)

    print(f"Total de feições filtradas para a Bacia do Itajaí: {len(filtered_features)} (excluídos das bacias vizinhas: {excluded_count})")

    # Salva o arquivo GeoJSON limpo
    out_fc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
        "features": filtered_features
    }
    with open("bacia_itajai_exutoria_porto.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)
    with open("rios_itajai.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)

    # -------------------------------------------------------------
    # GERAÇÃO DA FIGURA DE CONFERÊNCIA VISUAL PARA O USUÁRIO
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)

    for feat in filtered_features:
        name = feat["properties"].get("name", "")
        name_low = name.lower()
        g = shape(feat["geometry"])
        
        c = '#3498db'
        lw = 1.0
        
        if "itajai-acu" in name_low or "itajai acu" in name_low or (name == "Rio Itajaí" and g.length > 0.3):
            c = '#1b4f72' # Rio Itajaí-Açu
            lw = 2.5
        elif "mirim" in name_low:
            c = '#8e44ad' # Rio Itajaí-Mirim
            lw = 1.8
        elif "sul" in name_low:
            c = '#27ae60' # Rio Itajaí do Sul
            lw = 1.8
        elif "oeste" in name_low:
            c = '#d35400' # Rio Itajaí do Oeste
            lw = 1.8
        elif "norte" in name_low or "hercilio" in name_low:
            c = '#c0392b' # Rio Itajaí do Norte
            lw = 1.8

        if g.geom_type == "LineString":
            xs, ys = g.xy
            ax.plot(xs, ys, color=c, linewidth=lw, alpha=0.85)
        elif g.geom_type == "MultiLineString":
            for ls in g.geoms:
                xs, ys = ls.xy
                ax.plot(xs, ys, color=c, linewidth=lw, alpha=0.85)

    # Marca Porto de Itajaí (Exutória)
    ax.plot(FOZ_PORTO.x, FOZ_PORTO.y, '*', color='#e74c3c', markersize=20, markeredgecolor='black', zorder=10, label='Exutória: Porto de Itajaí')
    
    ax.annotate('EXUTÓRIA DA BACIA\n(Porto de Itajaí)', 
                (FOZ_PORTO.x, FOZ_PORTO.y), 
                xytext=(FOZ_PORTO.x + 0.08, FOZ_PORTO.y - 0.06),
                arrowprops=dict(facecolor='#e74c3c', shrink=0.08, width=2, headwidth=8),
                fontsize=10, fontweight='bold', color='#900c3f',
                bbox=dict(boxstyle="round,pad=0.3", fc="#fbeee6", ec="#e74c3c", lw=1.5))

    # Cidades de Referência
    cities = [
        (-49.60, -27.28, 'Taió'),
        (-49.60, -27.18, 'Ituporanga'),
        (-49.64, -27.21, 'Rio do Sul'),
        (-49.07, -26.92, 'Blumenau'),
        (-48.90, -27.11, 'Brusque'),
        (-48.65, -26.91, 'Itajaí / Foz'),
    ]
    for cx, cy, cname in cities:
        ax.plot(cx, cy, 'o', color='#f39c12', markersize=6, markeredgecolor='black', zorder=8)
        ax.text(cx + 0.02, cy + 0.01, cname, fontsize=8, fontweight='bold', color='#2c3e50',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#bdc3c7", alpha=0.85))

    ax.set_title('Rede Hidrográfica Filtrada Exclusiva da Bacia do Rio Itajaí (Exutória: Porto de Itajaí)', fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Longitude (WGS84)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Latitude (WGS84)', fontsize=10, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    
    output_dir = Path("figuras").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    target_fig = output_dir / "figura_rede_filtrada_itajai.png"
    fig.savefig(str(target_fig), dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"✓ FIGURA GERADA COM SUCESSO EM: {target_fig}")

if __name__ == "__main__":
    main()
