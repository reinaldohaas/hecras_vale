"""
Gerador de Polígonos de Inundação HAND Espaço-Temporais Sincronizados com Marés e Batimetria:
Produz os polígonos GeoJSON das manchas de inundação hora a hora:
- Em nível normal (t=0h): vazão contida na calha profunda (Área = 0.0 km², sem manchas espúrias).
- Em cheia (t=12h a 48h): transbordo na planície sincronizado com a passagem da onda e o ciclo das marés na foz.
"""

import sys
import json
from pathlib import Path
import numpy as np
from rasterio.features import shapes
from shapely.geometry import shape, mapping

repo_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(repo_root))

from itajai_flood_model.src.inundation.unified_hand_engine import (
    TopographicHANDModel,
    SynchronizedFloodEngine
)

def generate_all_synchronized_flood_layers():
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    dem_path = repo_root / "dem_blumenau_itajai.tif"
    dem_model_json = repo_root / "app" / "itajai_real_dem_model.json"
    out_geojson = repo_root / "app" / "manchas_inundacao_hand_dinamico.geojson"

    print("🌊 1. Inicializando TopographicHANDModel (Batimetria Real Z_bed + Cotas de Margem Z_bank)...")
    topo_model = TopographicHANDModel(str(dem_path), str(dem_model_json))
    flood_engine = SynchronizedFloodEngine(topo_model, min_flood_depth_m=0.10)
    
    with open(dem_model_json, 'r', encoding='utf-8') as f:
        dem_data = json.load(f)
        
    profiles = dem_data.get('river_profiles', {})
    
    # Parâmetros de pico dos eventos históricos e maré meteorológica (ressaca)
    events = {
        '1983': {'h_blu_peak': 15.34, 'h_rs_peak': 13.00, 'h_bq_peak': 8.50, 't_peak_rs': 18.0, 't_peak_blu': 24.0, 't_peak_foz': 30.0, 'surge_m': 1.40},
        '2008': {'h_blu_peak': 11.52, 'h_rs_peak': 4.20, 'h_bq_peak': 8.50, 't_peak_rs': 14.0, 't_peak_blu': 22.0, 't_peak_foz': 28.0, 'surge_m': 1.20},
        '2011': {'h_blu_peak': 12.60, 'h_rs_peak': 11.00, 'h_bq_peak': 6.20, 't_peak_rs': 16.0, 't_peak_blu': 24.0, 't_peak_foz': 30.0, 'surge_m': 0.90},
        '2023': {'h_blu_peak': 10.76, 'h_rs_peak': 8.50, 'h_bq_peak': 6.00, 't_peak_rs': 16.0, 't_peak_blu': 24.0, 't_peak_foz': 30.0, 'surge_m': 0.80}
    }
    
    time_steps = [0, 12, 18, 24, 30, 36, 48]
    all_features = []
    
    for ev_id, ev_cfg in events.items():
        print(f"\n📊 Processando Evento {ev_id} com Batimetria e Marés...")
        
        for t_h in time_steps:
            # 1. Construir perfis de cota absoluta Z_water(s, t) para os 10 rios
            water_profiles = {}
            for r_key, prof in profiles.items():
                coords = prof.get('coords', [])
                n_sec = len(coords)
                if n_sec < 2:
                    continue
                z_dem_prof = np.asarray(prof.get('z_dem') or prof.get('elevations', [10.0]*n_sec), dtype=float)
                
                z_water_sec = np.zeros(n_sec, dtype=float)
                
                for i in range(n_sec):
                    frac = i / float(n_sec - 1)
                    
                    if r_key == 'acu':
                        t_peak_local = (1.0 - frac) * ev_cfg['t_peak_rs'] + frac * ev_cfg['t_peak_foz']
                        h_peak_local = (1.0 - frac) * ev_cfg['h_rs_peak'] + frac * ev_cfg['h_blu_peak']
                        
                        pulse = np.exp(-0.5 * ((t_h - t_peak_local) / 9.0) ** 2)
                        
                        # Em t=0 (normal), H = 1.5m (abaixo da cota de transbordo de 8.0m em Blumenau e 7.0m em RS)
                        h_stage = 1.50 + (h_peak_local - 1.50) * pulse
                        
                        # Cota de zero da régua
                        if frac < 0.68:
                            w = frac / 0.68
                            z_zero_regua = (1.0 - w) * 332.0 + w * 4.88
                        else:
                            w = (frac - 0.68) / 0.32
                            z_zero_regua = (1.0 - w) * 4.88 + w * 0.0
                            
                        z_water_sec[i] = z_zero_regua + h_stage

                    elif r_key == 'mirim':
                        pulse_bq = np.exp(-0.5 * ((t_h - (ev_cfg['t_peak_blu'] - 2.0)) / 8.5) ** 2)
                        h_stage = 1.20 + (ev_cfg['h_bq_peak'] - 1.20) * pulse_bq
                        z_zero = (1.0 - frac) * 180.0 + frac * 0.5
                        z_water_sec[i] = z_zero + h_stage

                    else:
                        pulse_gen = np.exp(-0.5 * ((t_h - ev_cfg['t_peak_blu']) / 9.0) ** 2)
                        h_stage = 1.00 + (ev_cfg['h_blu_peak'] * 0.50 - 1.00) * pulse_gen
                        z_ref = z_dem_prof[i]
                        z_water_sec[i] = (z_ref - 3.5) + h_stage
                        
                water_profiles[r_key] = {'z_water_m': z_water_sec}

            # 2. Executar acoplamento no SynchronizedFloodEngine (com marés e transbordo)
            res = flood_engine.compute_instantaneous_flood(
                water_profiles,
                t_hour=t_h,
                storm_surge_peak_m=ev_cfg['surge_m'] if t_h > 12 else 0.0,
                max_corridor_cells=150
            )

            depth = res.depth_connected_2d
            print(f"   t = {t_h:02d}h -> Área Inundada: {res.area_km2:6.2f} km² | Vol: {res.volume_hm3:7.2f} hm³ | Maré: {res.ocean_level_z:+.2f}m")

            if res.area_km2 < 0.05:
                # Nível normal dentro da calha -> sem polígonos na planície
                continue

            classes = [
                {'name': 'Lâmina Baixa (0.1 - 1.0m)', 'mask': (depth >= 0.10) & (depth < 1.0), 'color': '#38bdf8', 'opacity': 0.50},
                {'name': 'Lâmina Média (1.0 - 2.5m)', 'mask': (depth >= 1.0) & (depth < 2.5), 'color': '#0284c7', 'opacity': 0.70},
                {'name': 'Lâmina Severa (> 2.5m)', 'mask': (depth >= 2.5), 'color': '#ff0055' if ev_id == '1983' else '#0369a1', 'opacity': 0.85}
            ]

            for cls in classes:
                m = cls['mask'].astype(np.uint8)
                if not np.any(m):
                    continue
                for geom_dict, val in shapes(m, mask=(m == 1), transform=topo_model.transform):
                    if val == 1:
                        poly_shp = shape(geom_dict)
                        if poly_shp.area < 5e-6:
                            continue
                        poly_shp = poly_shp.simplify(0.0006, preserve_topology=True)
                        feat = {
                            'type': 'Feature',
                            'geometry': mapping(poly_shp),
                            'properties': {
                                'event': ev_id,
                                'time_hour': t_h,
                                'class_name': cls['name'],
                                'fill_color': cls['color'],
                                'fill_opacity': cls['opacity'],
                                'area_km2': res.area_km2,
                                'volume_hm3': res.volume_hm3,
                                'ocean_tide_m': res.ocean_level_z
                            }
                        }
                        all_features.append(feat)

    final_fc = {
        'type': 'FeatureCollection',
        'properties': {
            'modelo': 'Synchronized Flood Engine + Topographic HAND + Ocean Tides',
            'calha_profunda': 'Batimetria real com Z_bed < 0m na foz e estuario',
            'mares': 'Astronômica Semidiurna + Ressaca Meteorológica',
            'transbordo': 'Bankfull Stage (Área = 0 km² em nivel normal)',
            'total_features': len(all_features)
        },
        'features': all_features
    }

    with open(out_geojson, 'w', encoding='utf-8') as f:
        json.dump(final_fc, f, indent=2)

    print(f"\n🎉 SUCESSO! {len(all_features)} polígonos sincronizados exportados para {out_geojson}!")

if __name__ == '__main__':
    generate_all_synchronized_flood_layers()
