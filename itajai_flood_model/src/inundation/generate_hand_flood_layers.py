"""
Gerador de Camadas de Inundação 2D HAND em Alta Fidelidade Topográfica:
Executa o modelo HAND sobre o DEM real GeoTIFF e a hidrografia da ANA, gerando
os polígonos fiéis à topografia e curvas de nível para os eventos históricos.
"""

import sys
import json
from pathlib import Path
import numpy as np

repo_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(repo_root))

from itajai_flood_model.src.inundation.hand_engine import HANDModel

def build_all_hand_flood_layers():
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    
    # Preferir dem_bacia_itajai.tif se disponível, senão dem_blumenau_itajai.tif
    dem_path = repo_root / "dem_blumenau_itajai.tif"
    river_geojson = repo_root / "app" / "itajai_ana_rios_alta_resolucao.geojson"
    out_geojson = repo_root / "app" / "manchas_inundacao_hand_real.geojson"

    print(f"🌊 Inicializando Modelo HAND com DEM: {dem_path.name}...")
    hand_model = HANDModel(str(dem_path), str(river_geojson))
    
    events_config = {
        '1983': {
            'name': 'Grande Cheia de 1983 (15.34m)',
            'stages': {'blumenau': 15.34, 'rio_sul': 13.00, 'brusque': 8.50, 'baixo_vale': 4.80},
            'color': '#ff0055'
        },
        '2008': {
            'name': 'Cheia de 2008 (11.52m / Brusque 8.50m)',
            'stages': {'blumenau': 11.52, 'rio_sul': 4.20, 'brusque': 8.50, 'baixo_vale': 3.50},
            'color': '#f97316'
        },
        '2011': {
            'name': 'Cheia de 2011 (12.60m)',
            'stages': {'blumenau': 12.60, 'rio_sul': 11.00, 'brusque': 6.20, 'baixo_vale': 3.80},
            'color': '#eab308'
        },
        '2023': {
            'name': 'Cheia de 2023 (10.76m)',
            'stages': {'blumenau': 10.76, 'rio_sul': 8.50, 'brusque': 6.00, 'baixo_vale': 3.20},
            'color': '#00f0ff'
        }
    }
    
    all_features = []
    
    for ev_id, ev_info in events_config.items():
        print(f"📊 Processando evento {ev_id}: {ev_info['name']}...")
        res = hand_model.generate_flood_inundation(ev_info['stages'], max_dist_cells=180, min_depth_m=0.15)
        print(f"   -> Área Inundada Calculada: {res['area_km2']} km² | Volume: {res['volume_hm3']} hm³ | Profundidade Média: {res['mean_depth_m']}m")
        
        # Vetorizar faixas de severidade
        depth = res['depth_raster']
        classes = [
            {'name': 'Lâmina Baixa (0.15m - 1.0m)', 'mask': (depth >= 0.15) & (depth < 1.0), 'opacity': 0.50, 'color': '#38bdf8'},
            {'name': 'Lâmina Média (1.0m - 2.5m)', 'mask': (depth >= 1.0) & (depth < 2.5), 'opacity': 0.70, 'color': '#0284c7'},
            {'name': 'Lâmina Severa (> 2.5m)', 'mask': (depth >= 2.5), 'opacity': 0.85, 'color': '#0369a1' if ev_id != '1983' else '#ff0055'}
        ]
        
        from rasterio.features import shapes
        from shapely.geometry import shape, mapping
        
        for cls in classes:
            m = cls['mask'].astype(np.uint8)
            for geom_dict, val in shapes(m, mask=(m == 1), transform=hand_model.transform):
                if val == 1:
                    poly_shp = shape(geom_dict)
                    if poly_shp.area < 1e-6:
                        continue
                    poly_shp = poly_shp.simplify(0.0003, preserve_topology=True)
                    feat = {
                        'type': 'Feature',
                        'geometry': mapping(poly_shp),
                        'properties': {
                            'event': ev_id,
                            'event_name': ev_info['name'],
                            'class_name': cls['name'],
                            'fill_color': cls['color'],
                            'fill_opacity': cls['opacity'],
                            'area_total_km2': res['area_km2'],
                            'volume_total_hm3': res['volume_hm3']
                        }
                    }
                    all_features.append(feat)

    final_fc = {
        'type': 'FeatureCollection',
        'properties': {
            'modelo': 'HAND (Height Above Nearest Drainage)',
            'base_topografica': 'Copernicus DEM 30m Real',
            'hidrografia': 'ANA 1:5.000 BHO',
            'total_features': len(all_features)
        },
        'features': all_features
    }
    
    with open(out_geojson, 'w', encoding='utf-8') as f:
        json.dump(final_fc, f, indent=2)
        
    print(f"🎉 SUCESSO TOTAL! {len(all_features)} polígonos HAND gerados e salvos em {out_geojson}!")

if __name__ == '__main__':
    build_all_hand_flood_layers()
