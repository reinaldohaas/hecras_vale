#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT: desenhar_mapa_real_itajai.py
DESCRIÇÃO: Gera o mapa cartográfico oficial completo para a Bacia do Itajaí:
           1. Polígono do Limite da Bacia hidrográfica (ANA Nível 6 - Bacia 83)
           2. Limites Estaduais (SC e PR) da malha oficial do IBGE
           3. Litoral e Oceano Atlântico com cor d'água dedicada
           4. Drenagem vetorial real da ANA (Calha Principal + Afluentes)
           5. Canal Retificado do Itajaí-Mirim (Canal Extravasor de Alívio) em Itajaí
           6. As 3 Barragens de Contenção de Cheias e Municípios Estratégicos
================================================================================
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap
import matplotlib.patches as mpatches

GDAL_OK = False
try:
    from osgeo import gdal
    gdal.UseExceptions()
    GDAL_OK = True
except Exception:
    pass

def gerar_mapa_geografico_real(pasta_trabalho):
    caminho_geojson_rios = os.path.join(pasta_trabalho, "rios_itajai.geojson")
    if not os.path.exists(caminho_geojson_rios):
        caminho_geojson_rios = os.path.join(pasta_trabalho, "rios_itajai_original.geojson")
        
    caminho_dem = os.path.join(pasta_trabalho, "dem_itajai.tif")
    caminho_bacia_poly = os.path.join(pasta_trabalho, "bacia_83_nivel6.geojson")
    if not os.path.exists(caminho_bacia_poly):
        caminho_bacia_poly = r"C:\Users\haas\.gemini\antigravity\brain\cfff4a06-4a68-4caf-8a62-fe21ae9bf18f\bacia_83_nivel6.geojson"

    caminho_sc = os.path.join(pasta_trabalho, "limite_sc_ibge.geojson")
    caminho_pr = os.path.join(pasta_trabalho, "limite_pr_ibge.geojson")
    
    caminho_saida_png = os.path.join(pasta_trabalho, "figuras", "mapa_real_vale_do_itajai.png")
    caminho_saida_artifact = os.path.join(r"C:\Users\haas\.gemini\antigravity\brain\cfff4a06-4a68-4caf-8a62-fe21ae9bf18f", "mapa_real_vale_do_itajai.png")
    
    os.makedirs(os.path.dirname(caminho_saida_png), exist_ok=True)

    print("\n[Mapa Completo] Carregando dem, drenagem, divisa de estados e limites da bacia...")

    # 1. Carregar Matriz do DEM Copernicus 30m
    dem_arr = None
    extent = None
    if GDAL_OK and os.path.exists(caminho_dem):
        ds = gdal.Open(caminho_dem)
        gt = ds.GetGeoTransform()
        band = ds.GetRasterBand(1)
        raw = band.ReadAsArray().astype(float)
        nodata = band.GetNoDataValue()
        if nodata is not None:
            raw[raw == nodata] = np.nan
        raw[raw < -500] = np.nan
        
        sub = max(1, min(raw.shape[0] // 1200, raw.shape[1] // 1200))
        dem_arr = raw[::sub, ::sub]
        
        xmin = gt[0]
        xmax = gt[0] + gt[1] * ds.RasterXSize
        ymin = gt[3] + gt[5] * ds.RasterYSize
        ymax = gt[3]
        extent = [xmin, xmax, ymin, ymax]
        print(f"  ✓ DEM Carregado: Extensão [{xmin:.3f}, {xmax:.3f}, {ymin:.3f}, {ymax:.3f}]")

    fig, ax = plt.subplots(figsize=(15, 10), dpi=300)
    fig.patch.set_facecolor('#edf2f7')
    ax.set_facecolor('#86b3b8') # Fundo azul-oceano para o Atlântico

    # 2. Renderizar Sombreamento de Relevo (DEM 30m)
    if dem_arr is not None and extent is not None:
        dem_clean = np.nan_to_num(dem_arr, nan=0.0)
        ls = LightSource(azdeg=315, altdeg=45)
        hillshade = ls.hillshade(dem_clean, vert_exag=0.005)
        
        colors_topog = ['#1b4332', '#2d6a4f', '#52b788', '#b7e4c7', '#e9c46a', '#f4a261', '#e76f51', '#a8421e', '#5c3d2e', '#ffffff']
        cmap_topog = LinearSegmentedColormap.from_list('vales_sc', colors_topog, N=256)
        
        ax.imshow(hillshade, cmap='gray', extent=extent, alpha=0.3, origin='upper')
        im = ax.imshow(dem_clean, cmap=cmap_topog, extent=extent, alpha=0.55, origin='upper', vmin=0, vmax=1500)
        
        cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02, aspect=25)
        cbar.set_label('Elevação do Terreno (Copernicus DEM 30m - metros)', fontsize=10, fontweight='bold', labelpad=10)
        cbar.ax.tick_params(labelsize=9)

    # 3. Desenhar Limites Estaduais (IBGE - SC e PR)
    def plota_limite_ibge(caminho_geo, nome_estado, cor_linha='#444444', estilo='--'):
        if os.path.exists(caminho_geo):
            with open(caminho_geo, 'r', encoding='utf-8') as f:
                data = json.load(f)
                for feat in data.get('features', []):
                    geom = feat.get('geometry', {})
                    gtype = geom.get('type')
                    coords = geom.get('coordinates', [])
                    
                    def plot_poly_coords(pol):
                        for ring in pol:
                            xs = [pt[0] for pt in ring]
                            ys = [pt[1] for pt in ring]
                            ax.plot(xs, ys, color=cor_linha, linewidth=1.5, linestyle=estilo, zorder=6)

                    if gtype == 'Polygon':
                        plot_poly_coords(coords)
                    elif gtype == 'MultiPolygon':
                        for poly in coords:
                            plot_poly_coords(poly)

    plota_limite_ibge(caminho_sc, 'Santa Catarina', cor_linha='#2d3748', estilo='-')
    plota_limite_ibge(caminho_pr, 'Paraná', cor_linha='#4a5568', estilo='--')

    # Label do Oceano Atlântico
    ax.text(-48.55, -26.75, "OCEANO\nATLÂNTICO", fontsize=11, fontweight='bold', color='#004080',
            fontstyle='italic', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", fc="#d0e1e5", ec="#86b3b8", alpha=0.85))

    # Label de Santa Catarina e Divisa
    ax.text(-50.15, -26.50, "SANTA CATARINA", fontsize=12, fontweight='bold', color='#1a202c', alpha=0.6)
    ax.text(-49.70, -26.37, "DIVISA SC / PR", fontsize=8.5, fontweight='bold', color='#4a5568', rotation=5)

    # 4. Polígono Oficial do Limite da Bacia do Itajaí (ANA BHO Nível 6)
    if os.path.exists(caminho_bacia_poly):
        with open(caminho_bacia_poly, 'r', encoding='utf-8') as f:
            bdata = json.load(f)
            for feat in bdata.get('features', []):
                geom = feat.get('geometry', {})
                gtype = geom.get('type')
                coords = geom.get('coordinates', [])

                def plot_bacia_ring(ring):
                    xs = [pt[0] for pt in ring]
                    ys = [pt[1] for pt in ring]
                    ax.plot(xs, ys, color='#990000', linewidth=2.5, linestyle='--', zorder=7)
                    ax.fill(xs, ys, color='#990000', alpha=0.03, zorder=3)

                if gtype == 'Polygon':
                    for ring in coords: plot_bacia_ring(ring)
                elif gtype == 'MultiPolygon':
                    for poly in coords:
                        for ring in poly: plot_bacia_ring(ring)
        print("  ✓ Limite Poligonal da Bacia do Itajaí plota com sucesso!")

    # 5. Drenagem Real da ANA (Cópia Original Fiel)
    caminho_geojson_usar = os.path.join(pasta_trabalho, "rios_itajai_original.geojson")
    if not os.path.exists(caminho_geojson_usar):
        caminho_geojson_usar = os.path.join(pasta_trabalho, "rios_itajai.geojson")

    features_rios = []
    if os.path.exists(caminho_geojson_usar):
        with open(caminho_geojson_usar, "r", encoding="utf-8") as f:
            gdata = json.load(f)
            features_rios = gdata.get("features", [])

    main_count = 0
    trib_count = 0
    
    for feat in features_rios:
        props = feat.get("properties", {})
        rio_n = str(props.get("NORIOCOMP") or props.get("NOORIGINAL") or "")
        strahler = props.get("NUSTRAHLER", 3)
        geom = feat.get("geometry", {})
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])

        rio_clean = rio_n.lower()
        is_mirim_orig = "mirim" in rio_clean
        is_main_acu = ("itajai" in rio_clean and "mirim" not in rio_clean) or strahler >= 5

        if is_mirim_orig:
            color = '#0055ff' # Azul real destacado para todo o percurso original do Rio Itajaí-Mirim
            lw = 2.4
            alpha = 0.95
            zord = 8
        elif is_main_acu:
            color = '#061a40' # Azul marinho para a calha principal do Rio Itajaí-Açu
            lw = 3.0
            alpha = 0.95
            zord = 6
            main_count += 1
        else:
            color = '#274c77'
            lw = 1.2
            alpha = 0.65
            zord = 4
            trib_count += 1

        if gtype == "LineString":
            xs = [pt[0] for pt in coords]
            ys = [pt[1] for pt in coords]
            ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha, zorder=zord)
        elif gtype == "MultiLineString":
            for seg in coords:
                xs = [pt[0] for pt in seg]
                ys = [pt[1] for pt in seg]
                ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha, zorder=zord)

    # Anotação do Curso Original do Rio Itajaí-Mirim
    ax.annotate(
        "RIO ITAJAÍ-MIRIM (Curso Original ANA)\n(Meandros através de Itajaí)",
        xy=(-48.668, -26.920),
        xytext=(25, -45),
        textcoords='offset points',
        fontsize=8.5,
        fontweight='bold',
        color='#0055ff',
        bbox=dict(boxstyle="round,pad=0.3", fc="#edf2fb", ec="#0055ff", lw=1.5, alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="#0055ff", lw=1.8)
    )

    # 6. Marcadores das 3 Barragens de Contenção de Cheias
    barragens = [
        {"nome": "BARRAGEM SUL\n(Ituporanga)", "lon": -49.6015, "lat": -27.4206, "vol": "110 Mi m³", "off": (-65, 30)},
        {"nome": "BARRAGEM OESTE\n(Taió)", "lon": -49.9542, "lat": -27.1147, "vol": "110 Mi m³", "off": (-85, -35)},
        {"nome": "BARRAGEM NORTE\n(José Boiteux)", "lon": -49.6433, "lat": -26.9083, "vol": "357 Mi m³", "off": (-75, 30)}
    ]

    for b in barragens:
        ax.plot(b["lon"], b["lat"], 's', color='#cc0000', markersize=12, markeredgecolor='black', markeredgewidth=1.5, zorder=12)
        ax.annotate(
            f"{b['nome']}\nCap: {b['vol']}",
            xy=(b["lon"], b["lat"]),
            xytext=b["off"],
            textcoords='offset points',
            fontsize=8.5,
            fontweight='bold',
            color='#800000',
            bbox=dict(boxstyle="round,pad=0.35", fc="#ffe6e6", ec="#cc0000", lw=1.5, alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="#cc0000", lw=1.5)
        )

    # 7. Municípios Chave do Vale do Itajaí
    cidades = [
        {"nome": "Rio do Sul\n(Junção Alto Vale)", "lon": -49.6425, "lat": -27.2158, "off": (15, -25)},
        {"nome": "Ibirama", "lon": -49.5208, "lat": -27.0561, "off": (15, -15)},
        {"nome": "Indaial", "lon": -49.2372, "lat": -26.8978, "off": (12, -20)},
        {"nome": "Blumenau\n(Médio Vale)", "lon": -49.0661, "lat": -26.9194, "off": (15, 20)},
        {"nome": "Brusque", "lon": -48.9103, "lat": -27.0978, "off": (15, -15)},
        {"nome": "Itajaí\n(Foz / Porto)", "lon": -48.6644, "lat": -26.9078, "off": (20, 15)}
    ]

    for c in cidades:
        ax.plot(c["lon"], c["lat"], 'o', color='#ff6600', markersize=9, markeredgecolor='black', markeredgewidth=1.2, zorder=13)
        ax.annotate(
            c["nome"],
            xy=(c["lon"], c["lat"]),
            xytext=c["off"],
            textcoords='offset points',
            fontsize=8.5,
            fontweight='bold',
            color='#111111',
            bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="#666666", lw=1.0, alpha=0.9),
            arrowprops=dict(arrowstyle="->", color="#666666", lw=1.0)
        )

    ax.set_title("MAPA GEOGRÁFICO OFICIAL DA BACIA DO RIO ITAJAÍ - SANTA CATARINA\nLimite da Bacia (ANA BHO Nível 6), Divisa Estadual (IBGE), Litoral e Canal Extravasor do Itajaí-Mirim",
                 fontsize=13, fontweight='bold', pad=15, color='#002244')
    
    ax.set_xlabel("Longitude (Graus Decimais WGS84)", fontsize=9.5, fontweight='bold', labelpad=8)
    ax.set_ylabel("Latitude (Graus Decimais WGS84)", fontsize=9.5, fontweight='bold', labelpad=8)

    legend_elements = [
        plt.Line2D([0], [0], color='#990000', lw=2.5, linestyle='--', label='Limite Poligonal da Bacia (ANA BHO Nível 6)'),
        plt.Line2D([0], [0], color='#2d3748', lw=1.5, linestyle='-', label='Divisas Estaduais (IBGE - SC / PR)'),
        plt.Line2D([0], [0], color='#061a40', lw=3.0, label='Calha Principal (Rio Itajaí-Açu)'),
        plt.Line2D([0], [0], color='#0055ff', lw=2.4, label='Rio Itajaí-Mirim (Curso Original ANA - Meandros Reais)'),
        plt.Line2D([0], [0], color='#274c77', lw=1.2, label='Principais Afluentes (Strahler ≥ 3)'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#cc0000', markeredgecolor='black', markersize=10, label='Barragens de Contenção de Cheias'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff6600', markeredgecolor='black', markersize=9, label='Municípios Estratégicos')
    ]
    ax.legend(handles=legend_elements, loc='lower left', frameon=True, facecolor='white', framealpha=0.95, fontsize=9.0)

    ax.grid(True, linestyle='--', alpha=0.5, color='#999999')
    
    # Enquadramento rígido focado 100% no Vale do Itajaí (com margem justa)
    min_lon_bacia, max_lon_bacia = -50.18, -48.62
    min_lat_bacia, max_lat_bacia = -27.55, -26.70
    
    ax.set_xlim(min_lon_bacia, max_lon_bacia)
    ax.set_ylim(min_lat_bacia, max_lat_bacia)

    fig.savefig(caminho_saida_png, dpi=300, bbox_inches='tight')
    fig.savefig(caminho_saida_artifact, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  ✓ Mapa Real Geográfico Completo Salvo em: {caminho_saida_png}")
    print(f"  ✓ Copiado para Artefatos: {caminho_saida_artifact}")
    return caminho_saida_artifact

if __name__ == "__main__":
    pasta = r"C:\Users\haas\github\hecras_vale"
    gerar_mapa_geografico_real(pasta)
