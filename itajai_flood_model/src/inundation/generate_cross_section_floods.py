"""
Gerador de Manchas 2D de Inundação por Seções Transversais (Cross-Section Delineation):
Calcula a largura de inundação W_left(x) e W_right(x) para cada uma das seções transversais
ao longo dos 10 rios e exporta polígonos GeoJSON contínuos de lâmina d'água.
"""

import sys
import json
from pathlib import Path
import numpy as np

repo_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(repo_root))

from itajai_flood_model.src.inundation.cross_sections import CrossSectionDelineator

def generate_cross_section_flood_layers():
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    dem_model_path = repo_root / "app" / "itajai_real_dem_model.json"
    out_geojson_path = repo_root / "app" / "manchas_inundacao_secoes_transversais.geojson"
    
    with open(dem_model_path, 'r', encoding='utf-8') as f:
        dem_data = json.load(f)
        
    profiles = dem_data.get('river_profiles', {})
    
    # Cotas de referência por evento em Blumenau e Rio do Sul
    events = {
        '1983': {'h_blu': 15.34, 'h_rs': 13.00, 'h_bq': 7.60},
        '2008': {'h_blu': 11.52, 'h_rs': 4.20, 'h_bq': 8.50},
        '2011': {'h_blu': 12.60, 'h_rs': 11.00, 'h_bq': 6.20},
        '2023': {'h_blu': 10.76, 'h_rs': 8.50, 'h_bq': 6.00}
    }
    
    features_out = []
    
    for ev_name, ev_h in events.items():
        for r_key, p_data in profiles.items():
            coords = p_data.get('coords', [])
            z_bed = p_data.get('z_dem') or p_data.get('elevations', [])
            if len(coords) < 2 or len(z_bed) < len(coords):
                continue
                
            # Construir seções transversais
            sections = CrossSectionDelineator.build_cross_sections_for_river(coords, z_bed, r_key)
            
            # Gerar perfil de linha d'água Z_water(x)
            z_water_prof = []
            for i, sec in enumerate(sections):
                frac = i / float(len(sections) - 1)
                
                # Cota de régua local
                if r_key == 'acu':
                    h_loc = (1.0 - frac) * ev_h['h_rs'] + frac * (ev_h['h_blu'] * 0.7)
                    z0_loc = (1.0 - frac) * 335.5 + frac * 4.88
                elif r_key == 'mirim':
                    h_loc = ev_h['h_bq'] * (1.0 - 0.3 * frac)
                    z0_loc = (1.0 - frac) * 18.5 + frac * 0.0
                elif r_key in ('oeste', 'sul', 'trombudo', 'perimbo', 'mirim_doce'):
                    h_loc = ev_h['h_rs'] * 0.85
                    z0_loc = sec.z_bed + 2.0
                else:
                    h_loc = ev_h['h_blu'] * 0.65
                    z0_loc = sec.z_bed + 2.0
                    
                z_w = z0_loc + h_loc
                z_water_prof.append(float(z_w))
                
            poly_feat = CrossSectionDelineator.delineate_flood_polygon(sections, z_water_prof)
            if poly_feat:
                poly_feat['properties']['event'] = ev_name
                poly_feat['properties']['river_name'] = p_data.get('name', r_key)
                poly_feat['properties']['fill_color'] = '#00f0ff' if ev_name == '2023' else ('#ff0055' if ev_name == '1983' else '#f59e0b')
                features_out.append(poly_feat)
                
    geojson_final = {
        'type': 'FeatureCollection',
        'features': features_out
    }
    
    with open(out_geojson_path, 'w', encoding='utf-8') as f:
        json.dump(geojson_final, f, indent=2)
        
    print(f"✅ Manchas 2D baseadas em seções transversais exportadas para {out_geojson_path} ({len(features_out)} polígonos)!")

if __name__ == '__main__':
    generate_cross_section_flood_layers()
