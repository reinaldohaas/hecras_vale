"""
Script de Geração e Exportação das 8 Camadas de Debug do Modelo HAND Sincronizado:
Gera e salva no diretório `app/debug_rasters/` os seguintes arquivos:
1. 01_dem.tif
2. 02_drainage_id.tif
3. 03_drainage_elevation.tif
4. 04_hand.tif
5. 05_water_surface_t24.tif
6. 06_relative_water_level_t24.tif
7. 07_depth_t24.tif
8. 08_connected_flood_t24.tif
E exporta também a grade JSON correspondente para visualização imediata no navegador (app/debug_hand_layers.html).
"""

import sys
import json
from pathlib import Path
import numpy as np

repo_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(repo_root))

from itajai_flood_model.src.inundation.unified_hand_engine import (
    TopographicHANDModel,
    SynchronizedFloodEngine
)

def export_all_debug_layers():
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    dem_path = repo_root / "dem_blumenau_itajai.tif"
    net_json = repo_root / "app" / "itajai_real_dem_model.json"
    out_dir = repo_root / "app" / "debug_rasters"
    out_json = repo_root / "app" / "debug_hand_matrices.json"

    print("🌊 1. Inicializando Modelo Topográfico HAND Estático...")
    topo_model = TopographicHANDModel(dem_path, net_json)
    
    print("📁 2. Exportando Rasters Topográficos (01 a 04)...")
    topo_model.export_debug_rasters(out_dir)

    print("⚡ 3. Executando Acoplamento Hidráulico Sincronizado para t=24h (Pico de 1983)...")
    flood_engine = SynchronizedFloodEngine(topo_model, min_flood_depth_m=0.05)

    with open(net_json, 'r', encoding='utf-8') as f:
        net_data = json.load(f)
    profiles = net_data.get('river_profiles', {})

    # Perfil longitudinal da grande cheia de 1983 no pico (t=24h)
    # Blumenau = 15.34m (Z=20.22m), Rio do Sul = 13.0m (Z=348.5m), Foz = 4.5m (Z=1.0m)
    water_profiles = {}
    for r_key, prof in profiles.items():
        coords = prof.get('coords', [])
        z_bed = np.asarray(prof.get('z_dem') or prof.get('elevations', [10.0]*len(coords)), dtype=float)
        n_pts = len(coords)
        if n_pts < 2:
            continue

        if r_key == 'acu':
            # Rio do Sul (km 0) -> Blumenau (km 105) -> Foz (km 153)
            s_vec = np.linspace(0, 1, n_pts)
            z_w = np.zeros(n_pts)
            for i in range(n_pts):
                s = s_vec[i]
                if s < 0.68: # Rio do Sul até Blumenau
                    f = s / 0.68
                    z_w[i] = (1.0 - f) * (335.5 + 13.0) + f * (4.88 + 15.34)
                else: # Blumenau até Foz
                    f = (s - 0.68) / 0.32
                    z_w[i] = (1.0 - f) * (4.88 + 15.34) + f * (0.0 + 3.80)
            water_profiles['acu'] = {'z_water_m': z_w}
        elif r_key == 'mirim':
            z_w = z_bed + 8.50 # Brusque em 8.50m
            water_profiles['mirim'] = {'z_water_m': z_w}
        else:
            z_w = z_bed + 6.00
            water_profiles[r_key] = {'z_water_m': z_w}

    flood_raster_t24 = flood_engine.compute_instantaneous_flood(water_profiles, t_hour=24.0, max_corridor_cells=160)

    print("📁 4. Exportando Rasters Hidráulicos (05 a 08)...")
    flood_engine.export_flood_debug_rasters(flood_raster_t24, out_dir)

    print("📊 5. Exportando Dados JSON Amostrados para o Visualizador Web...")
    # Amostrar grade para exibição leve no navegador (ex: step de 6 células -> ~120x300 pontos)
    step = 6
    sample_dem = topo_model.dem[::step, ::step]
    sample_drain_id = topo_model.drainage_id[::step, ::step]
    sample_drain_z = topo_model.drainage_elevation[::step, ::step]
    sample_hand = topo_model.hand[::step, ::step]
    sample_zw = flood_raster_t24.z_water_2d[::step, ::step]
    sample_eta = flood_raster_t24.eta_2d[::step, ::step]
    sample_depth = flood_raster_t24.depth_raw_2d[::step, ::step]
    sample_conn = flood_raster_t24.depth_connected_2d[::step, ::step]

    # Bounds e coordenadas
    nrows_s, ncols_s = sample_dem.shape
    lons = np.linspace(topo_model.bounds.left, topo_model.bounds.right, ncols_s)
    lats = np.linspace(topo_model.bounds.top, topo_model.bounds.bottom, nrows_s)

    debug_json_data = {
        'bounds': [topo_model.bounds.left, topo_model.bounds.bottom, topo_model.bounds.right, topo_model.bounds.top],
        'shape': [nrows_s, ncols_s],
        'lons': np.round(lons, 4).tolist(),
        'lats': np.round(lats, 4).tolist(),
        'metrics': {
            'area_connected_km2': flood_raster_t24.area_km2,
            'volume_connected_hm3': flood_raster_t24.volume_hm3,
            'max_depth_m': flood_raster_t24.max_depth_m,
            'mean_depth_m': flood_raster_t24.mean_depth_m
        },
        'layers': {
            'dem': np.round(sample_dem, 1).tolist(),
            'drainage_id': sample_drain_id.tolist(),
            'drainage_elevation': np.round(sample_drain_z, 1).tolist(),
            'hand': np.round(sample_hand, 2).tolist(),
            'water_surface': np.round(sample_zw, 2).tolist(),
            'relative_water_level_eta': np.round(sample_eta, 2).tolist(),
            'depth_raw': np.round(sample_depth, 2).tolist(),
            'depth_connected': np.round(sample_conn, 2).tolist()
        }
    }

    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(debug_json_data, f)

    print(f"🎉 SUCESSO! 8 Rasters GeoTIFF exportados em {out_dir} e matrizes JSON salvas em {out_json}!")

if __name__ == '__main__':
    export_all_debug_layers()
