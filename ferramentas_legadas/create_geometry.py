import os
import math
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge

def process_river_centerline(geojson_path):
    print(f"Lendo o arquivo GeoJSON: {geojson_path}")
    gdf = gpd.read_file(geojson_path)
    
    gdf_utm = gdf.to_crs(epsg=31982)
    lines = gdf_utm.geometry.tolist()
    merged_line = linemerge(lines)
    
    if isinstance(merged_line, MultiLineString):
        longest_line = max(merged_line.geoms, key=lambda line: line.length)
        merged_line = longest_line
        
    coords = list(merged_line.coords)
    
    start_point = coords[0]
    end_point = coords[-1]
    
    if start_point[0] > end_point[0]:
        coords.reverse()
        merged_line = LineString(coords)
        print("Linha invertida para fluir de Montante para Jusante.")
        
    return coords, merged_line

def format_16(val):
    """Formata para o padrão HEC-RAS de 16 caracteres para dados GIS"""
    s = f"{val:.4f}"
    return f"{s:>16}"

def format_8(val):
    """Formata para o padrão HEC-RAS de 8 caracteres para dados de seção"""
    s = f"{val:.2f}"
    return f"{s:>8}"

def create_hecras_files():
    project_name = "Itajai_Blumenau"
    prj_file = f"{project_name}.prj"
    geom_file = f"{project_name}.g01"
    flow_file = f"{project_name}.u01"
    plan_file = f"{project_name}.p01"
    geojson_path = r"C:\Users\haas\Downloads\export.geojson"
    
    try:
        river_coords, merged_line = process_river_centerline(geojson_path)
        total_length = merged_line.length
        print(f"Comprimento total do rio processado: {total_length/1000:.2f} km")
    except Exception as e:
        print(f"Erro ao processar GeoJSON: {e}")
        return

    # 1. Project File
    with open(prj_file, "w") as f:
        f.write(f"Proj Title={project_name}\n")
        f.write("Program Version=7.01\n")
        f.write(f"Geom File=g01\n")
        f.write(f"Unsteady File=u01\n")
        f.write(f"Plan File=p01\n")
        f.write("Y Axis Title=Elevation\n")
        f.write("X Axis Title(PR)=Distance\n")
        f.write("X Axis Title(CS)=Station\n")
        f.write("SI Units\n")

    # 2. Geometry File
    spacing = 1000
    stations = np.arange(int(total_length), -spacing, -spacing)
    if stations[-1] != 0:
        stations = np.append(stations, 0)
        
    with open(geom_file, "w") as f:
        f.write("Geom Title=Geometria Real 1D\n")
        f.write("Program Version=7.01\n")
        
        f.write("River Reach= Itajai_Acu,BaixoVale\n")
        
        # Reach XY usa formato de 16 caracteres, 2 pares por linha
        f.write(f"Reach XY= {len(river_coords)} \n")
        for i in range(0, len(river_coords), 2):
            line = ""
            for j in range(2):
                if i+j < len(river_coords):
                    x, y = river_coords[i+j]
                    line += format_16(x) + format_16(y)
            f.write(line + "\n")
            
        for st in stations:
            dist_from_upstream = total_length - st
            p_center = merged_line.interpolate(dist_from_upstream)
            
            p1 = merged_line.interpolate(max(0, dist_from_upstream - 1.0))
            p2 = merged_line.interpolate(min(total_length, dist_from_upstream + 1.0))
            dx = p2.x - p1.x
            dy = p2.y - p1.y
            length = math.hypot(dx, dy)
            if length == 0:
                nx, ny = 0, 1
            else:
                nx = -dy / length
                ny = dx / length
                
            left_x = p_center.x - nx * (-100)
            left_y = p_center.y - ny * (-100)
            right_x = p_center.x - nx * 100
            right_y = p_center.y - ny * 100
            
            z_bottom = -15.0 + (st / total_length) * 15.0
            
            st_str = f"{st:.2f}"
            
            # Reach length (distância até a próxima seção de jusante)
            # Se for a última seção (st == 0), a distância é 0.
            if st == 0:
                rl = 0.0
            else:
                # Distância real até a próxima seção na array
                current_idx = np.where(stations == st)[0][0]
                next_st = stations[current_idx + 1]
                rl = st - next_st
                
            f.write(f"Type RM Length L Ch R = 1 , {st_str} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
            f.write("Node Last Edited Time= Jan/01/2026 00:00:00\n")
            
            f.write("XS GIS Cut Line= 2 \n")
            f.write(format_16(left_x) + format_16(left_y) + format_16(right_x) + format_16(right_y) + "\n")
            
            # AQUI ESTAVA O ERRO! O HEC-RAS exige a tag #Sta/Elev= e formato de 8 colunas
            f.write("#Sta/Elev= 5 \n")
            pts = [
                (-100, z_bottom + 10),
                (-75, z_bottom),
                (0, z_bottom),
                (75, z_bottom),
                (100, z_bottom + 10)
            ]
            line = ""
            for px, py in pts:
                line += format_8(px) + format_8(py)
            f.write(line + "\n")
            
            f.write("Bank Sta=-75, 75\n")
            # A sintaxe do Manning também deve usar as tags corretas. Para simplificar, usamos n uniforme no canal.
            # O formato correto básico para Manning pode ser dado com a tag abaixo
            f.write("#Mann= 3 , -1 , 0 \n")
            # Station, n_value, 0
            f.write(format_8(-100) + format_8(0.05) + format_8(0) + format_8(-75) + format_8(0.035) + format_8(0) + format_8(75) + format_8(0.05) + format_8(0) + "\n")
            f.write("Bank Sta= -75, 75 \n")
            
    # 3. Unsteady Flow File
    with open(flow_file, "w") as f:
        f.write("Flow Title=Dummy_Data\n")
        f.write("Program Version=7.01\n")
        
        upstream_station = f"{stations[0]:.2f}"
        f.write(f"Boundary Location= Itajai_Acu,BaixoVale, {upstream_station} \n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 24 \n")
        flows = ["1000"] * 6 + ["3000"] * 6 + ["5000"] * 6 + ["2000"] * 6
        f.write(" ".join(flows) + "\n")
        
        f.write("Boundary Location= Itajai_Acu,BaixoVale, 0.00 \n")
        f.write("Interval= 1HOUR\n")
        f.write("Stage Hydrograph= 24 \n")
        stages = ["0"] * 12 + ["1"] * 12
        f.write(" ".join(stages) + "\n")
        
        f.write("Initial Stage= 2 \n")
        f.write("Initial Flow= 1000 \n")

    # 4. Plan File
    with open(plan_file, "w") as f:
        f.write("Plan Title=Simulacao_Previsao\n")
        f.write("Program Version=7.01\n")
        f.write("Geom File=g01\n")
        f.write("Unsteady File=u01\n")
        f.write("Simulation Date= 01Jan2026, 0000, 02Jan2026, 0000\n")
        f.write("Computation Interval= 1MIN\n")
        f.write("Output Interval= 1HOUR\n")
        f.write("Unsteady Routing= 1\n")

    print(f"Arquivos HEC-RAS gerados com sintaxe nativa estrita corrigida!")

if __name__ == "__main__":
    create_hecras_files()
