import json
import numpy as np
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge
from pyproj import Transformer
import os
import gerar_geometria_hecras
from ras_commander import init_ras_project, RasCmdr, RasPreprocess

dem_path = r'C:\Users\haas\github\hecras_vale\dem_itajai.tif'
geojson_path = r'C:\Users\haas\github\hecras_vale\rios_itajai.geojson'

dem_sampler = gerar_geometria_hecras.DemSampler(dem_path)
to_utm = Transformer.from_crs(4326, 31982, always_xy=True)

with open(geojson_path, 'r', encoding='utf-8') as f:
    data = json.load(f)

lines_utm = []
for feat in data['features']:
    geom = feat.get('geometry')
    if not geom: continue
    gtype = geom.get('type')
    coords = geom.get('coordinates', [])
    if gtype == 'LineString':
        lines_utm.append(LineString([to_utm.transform(pt[0], pt[1]) for pt in coords]))
    elif gtype == 'MultiLineString':
        for seg in coords:
            lines_utm.append(LineString([to_utm.transform(pt[0], pt[1]) for pt in seg]))

merged = linemerge(lines_utm)
if isinstance(merged, MultiLineString):
    main_line = max(merged.geoms, key=lambda l: l.length)
else:
    main_line = merged

print('Comprimento do canal principal:', main_line.length / 1000.0, 'km')

# Extrair coordenadas do canal principal para o Reach Header
reach_coords = list(main_line.coords)
# Simplificar pontos para evitar mais de 500 pontos no Reach Header
if len(reach_coords) > 200:
    idx = np.linspace(0, len(reach_coords) - 1, 200).round().astype(int)
    reach_coords = [reach_coords[i] for i in idx]

n_sec = 60
dists = np.linspace(0, main_line.length, n_sec)
xs_list = []

for d in dists:
    eps = 5.0
    p = main_line.interpolate(d)
    pa = main_line.interpolate(max(d - eps, 0))
    pb = main_line.interpolate(min(d + eps, main_line.length))
    tx, ty = pb.x - pa.x, pb.y - pa.y
    tl = np.hypot(tx, ty) or 1.0
    rx, ry = ty / tl, -tx / tl
    
    p_left = (p.x - 600.0 * rx, p.y - 600.0 * ry)
    p_right = (p.x + 600.0 * rx, p.y + 600.0 * ry)
    
    res = gerar_geometria_hecras.cortar_secao_dem(main_line, d, dem_sampler)
    if res:
        sta, z = res
        lb, rb = gerar_geometria_hecras.estacoes_margem(sta, z)
        xs_list.append({
            'station_m': round(main_line.length - d, 1),
            'sta': sta,
            'z': z,
            'lb': lb,
            'rb': rb,
            'cutline': [p_left[0], p_left[1], p_right[0], p_right[1]]
        })

prj_name = 'Itajai_Calha_Georef'
g01_file = r'C:\Users\haas\github\hecras_vale\Itajai_Calha_Georef.g01'
prj_file = r'C:\Users\haas\github\hecras_vale\Itajai_Calha_Georef.prj'
u01_file = r'C:\Users\haas\github\hecras_vale\Itajai_Calha_Georef.u01'
p01_file = r'C:\Users\haas\github\hecras_vale\Itajai_Calha_Georef.p01'

with open(g01_file, 'w', encoding='ascii', errors='replace') as f:
    f.write(f'Geom Title={prj_name}\n')
    f.write('Program Version=6.10\n\n')
    f.write('River Title=ITAJAI\n')
    f.write('Reach Title=PRINCIPAL\n')
    
    # Write Reach Header centerline coordinates
    f.write(f'Reach Header=ITAJAI          ,PRINCIPAL       ,{len(reach_coords)}\n')
    rc_flat = []
    for cx, cy in reach_coords:
        rc_flat.extend([cx, cy])
    f.write(gerar_geometria_hecras.fixed_series(rc_flat) + '\n')

    for xs in xs_list:
        f.write('Type Header= 1 \n')
        f.write(f'River Reach Header=ITAJAI          ,PRINCIPAL       ,{xs["station_m"]:8.1f}\n')
        f.write('BEGIN DESCRIPTION:\nSecao amostrada DEM 30m\nEND DESCRIPTION:\n')
        f.write(f'NODE ID={xs["station_m"]:8.1f}\n')
        f.write('Node Last Edited Time=Aug/05/2026 00:00:00\n')
        f.write(f'#XS STA & ELEV= {len(xs["sta"])}\n')
        pts = []
        for s_val, z_val in zip(xs["sta"], xs["z"]):
            pts.extend([s_val, z_val])
        f.write(gerar_geometria_hecras.fixed_series(pts) + '\n')
        f.write("Manning's n Values= 3\n")
        f.write('    0.00    0.060    0.00    0.035    0.00    0.060\n')
        f.write(f'Bank Sta={xs["lb"]:8.2f},{xs["rb"]:8.2f}\n')
        f.write('Length Values= 1000.00,1000.00,1000.00\n')
        f.write('GIS Cut Line= 2\n')
        f.write(f'  {xs["cutline"][0]:12.2f}{xs["cutline"][1]:12.2f}{xs["cutline"][2]:12.2f}{xs["cutline"][3]:12.2f}\n')
        f.write('Type Header= 0 \n\n')

with open(prj_file, 'w', encoding='ascii', errors='replace') as f:
    f.write(f'Proj Title={prj_name}\n')
    f.write('Current Plan=p01\n')
    f.write('Default Directory=\n')
    f.write('Geom File=g01\n')
    f.write('Unsteady File=u01\n')
    f.write('Plan File=p01\n')

with open(u01_file, 'w', encoding='ascii', errors='replace') as f:
    f.write('Flow Title=Cheia_Sintetica\n')
    f.write('Program Version=6.10\n\n')
    hidro = gerar_geometria_hecras.gerar_hidrograma(1800.0)
    f.write(f'Boundary Location=ITAJAI          ,PRINCIPAL       ,{xs_list[0]["station_m"]:8.1f},, , , , \n')
    f.write('Interval=1HOUR\n')
    f.write(f'Flow Hydrograph= {len(hidro)}\n')
    f.write(gerar_geometria_hecras.fixed_series(hidro) + '\n')

with open(p01_file, 'w', encoding='ascii', errors='replace') as f:
    f.write('Plan Title=Plano_Simulacao\n')
    f.write('Program Version=6.10\n')
    f.write('Short Identifier=01\n')
    f.write('Simulation Date=01AUG2026,0000,03AUG2026,0100\n')
    f.write('Geom File=g01\n')
    f.write('Flow File=u01\n')
    f.write('Unsteady File=u01\n')

print('✓ Geometria georreferenciada HEC-RAS gerada com sucesso!')

ver = r'C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe'
print('Testando ras-commander no projeto georreferenciado...')
ras = init_ras_project(prj_file, ver)
res = RasPreprocess.preprocess_plan('01', ras_object=ras)
print('Preprocess Result:', res)
