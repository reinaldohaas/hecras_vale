"""
FILTRAGEM E ILUSTRAÇÃO DA BACIA DO RIO ITAJAÍ (SEM DEPENDÊNCIA DE C-EXTENSIONS)
"""
import json
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FOZ_LON = -48.652
FOZ_LAT = -26.910

def main():
    print("Filtrando rios exclusivos da Bacia do Itajaí...")
    with open("rios_itajai.geojson", "r", encoding="utf-8") as f:
        data = json.load(f)

    filtered = []
    excluded = 0

    for feat in data["features"]:
        name = feat["properties"].get("name", "")
        name_ascii = feat["properties"].get("name_ascii", "").lower()
        
        # Filtra bacias vizinhas externas
        if any(b in name_ascii for b in ["itapocu", "tijucas", "canoas", "cubatao", "biguacu", "negro", "iguacu"]):
            excluded += 1
            continue
            
        geom = feat.get("geometry", {})
        coords = geom.get("coordinates", [])
        if not coords:
            continue
            
        if geom["type"] == "LineString":
            lines = [coords]
        elif geom["type"] == "MultiLineString":
            lines = coords
        else:
            continue
            
        lons = [pt[0] for line in lines for pt in line]
        lats = [pt[1] for line in lines for pt in line]
        
        min_lon, max_lon = min(lons), max(lons)
        min_lat, max_lat = min(lats), max(lats)
        
        # Limites geográficos estritos da Bacia do Itajaí
        if min_lat < -27.45 and max_lon > -49.1:  # Bacia do Tijucas
            excluded += 1
            continue
        if max_lat > -26.55:  # Bacia do Itapocu
            excluded += 1
            continue
        if min_lon < -50.15:  # Bacia do Canoas
            excluded += 1
            continue
            
        filtered.append(feat)

    print(f"✓ Selecionados {len(filtered)} rios pertencentes à Bacia do Itajaí (removidos {excluded} de outras bacias).")

    out_fc = {
        "type": "FeatureCollection",
        "features": filtered
    }
    with open("bacia_itajai_exutoria_porto.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)
    with open("rios_itajai.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)
    with open("vale_itajai_full_network.geojson", "w", encoding="utf-8") as f:
        json.dump(out_fc, f, ensure_ascii=False, indent=2)

    # -------------------------------------------------------------
    # GERAÇÃO DA FIGURA DE CONFERÊNCIA VISUAL
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)

    for feat in filtered:
        name = feat["properties"].get("name", "")
        name_low = name.lower()
        geom = feat.get("geometry", {})
        coords = geom.get("coordinates", [])
        
        if geom["type"] == "LineString":
            lines = [coords]
        else:
            lines = coords
            
        c = '#3498db'
        lw = 1.0
        
        if "itajai-acu" in name_low or "itajai acu" in name_low or name == "Rio Itajaí":
            c = '#1b4f72'; lw = 2.5
        elif "mirim" in name_low:
            c = '#8e44ad'; lw = 1.8
        elif "sul" in name_low:
            c = '#27ae60'; lw = 1.8
        elif "oeste" in name_low:
            c = '#d35400'; lw = 1.8
        elif "norte" in name_low or "hercilio" in name_low:
            c = '#c0392b'; lw = 1.8

        for line in lines:
            xs = [pt[0] for pt in line]
            ys = [pt[1] for pt in line]
            ax.plot(xs, ys, color=c, linewidth=lw, alpha=0.85)

    # Marca Exutória no Porto de Itajaí
    ax.plot(FOZ_LON, FOZ_LAT, '*', color='#e74c3c', markersize=20, markeredgecolor='black', zorder=10)
    ax.annotate('EXUTÓRIA DA BACIA\n(Porto de Itajaí / Oceano)', 
                (FOZ_LON, FOZ_LAT), 
                xytext=(FOZ_LON + 0.08, FOZ_LAT - 0.06),
                arrowprops=dict(facecolor='#e74c3c', shrink=0.08, width=2, headwidth=8),
                fontsize=10, fontweight='bold', color='#900c3f',
                bbox=dict(boxstyle="round,pad=0.3", fc="#fbeee6", ec="#e74c3c", lw=1.5))

    cities = [
        (-49.60, -27.28, 'Taió (Oeste)'),
        (-49.60, -27.18, 'Ituporanga (Sul)'),
        (-49.65, -26.92, 'José Boiteux (Norte)'),
        (-49.64, -27.21, 'Rio do Sul (Confluência)'),
        (-49.07, -26.92, 'Blumenau'),
        (-48.90, -27.11, 'Brusque (Mirim)'),
        (-48.65, -26.91, 'Itajaí / Foz'),
    ]
    for cx, cy, cname in cities:
        ax.plot(cx, cy, 'o', color='#f39c12', markersize=6, markeredgecolor='black', zorder=8)
        ax.text(cx + 0.015, cy + 0.01, cname, fontsize=8, fontweight='bold', color='#2c3e50',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#bdc3c7", alpha=0.85))

    ax.set_title('Rede Hidrográfica Filtrada Exclusiva da Bacia do Rio Itajaí (Exutória: Porto de Itajaí)', fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Longitude (Graus WGS84)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Latitude (Graus WGS84)', fontsize=10, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)

    os.makedirs("figuras", exist_ok=True)
    fig_path = os.path.abspath("figuras/figura_rede_filtrada_itajai.png")
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ FIGURA GERADA COM SUCESSO EM: {fig_path}")

if __name__ == "__main__":
    main()
