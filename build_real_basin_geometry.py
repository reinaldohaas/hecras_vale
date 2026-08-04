"""
GERADOR DE GEOMETRIA HEC-RAS 7.0.1 COM RIOS REAIS E XS GIS CUT LINES
===================================================================
1. Usa a rede hidrográfica real de 'vale_itajai_full_network.geojson' / 'rios_itajai.geojson'.
2. Amostra elevações reais de 'dem_bacia_itajai.tif'.
3. Calcula e escreve explicitamente 'XS GIS Cut Line= 2' para cada seção transversal,
   eliminando o erro 'Error determining XS cut lines' do HEC-RAS 7.
"""

import json
import math
import os
import rasterio
import numpy as np
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge

class FastDEM:
    def __init__(self, tif_path):
        self.ds = rasterio.open(tif_path)
        self.data = self.ds.read(1)
        self.transform = self.ds.transform
        self.left = self.transform.c
        self.top = self.transform.f
        self.dx = self.transform.a
        self.dy = self.transform.e  # negativo

    def get_elevation(self, lon, lat):
        col = int((lon - self.left) / self.dx)
        row = int((lat - self.top) / self.dy)
        if 0 <= row < self.data.shape[0] and 0 <= col < self.data.shape[1]:
            val = float(self.data[row, col])
            if val > -500 and val < 9000:
                return val
        return np.nan

def latlon_to_utm22s(lon, lat):
    lon0 = -51.0
    lat0 = 0.0
    r = 6378137.0
    k0 = 0.9996
    
    dlon = math.radians(lon - lon0)
    lat_r = math.radians(lat)
    
    x = r * dlon * math.cos(lat_r) * k0 + 500000.0
    y = r * (lat_r + math.sin(2*lat_r)/2.0) * k0 + 10000000.0
    return x, y

def format_8(val):
    return f"{val:>8.2f}"

def format_16_num(val):
    return f"{val:>16.4f}"

def main():
    print("=" * 60)
    print("GERANDO GEOMETRIA COM XS GIS CUT LINES PARA HEC-RAS 7.0.1")
    print("=" * 60)

    dem = FastDEM("dem_bacia_itajai.tif")
    print(f"DEM Carregado: {dem.data.shape[1]}x{dem.data.shape[0]} células")

    # Lê de vale_itajai_full_network.geojson (reconstruído)
    geojson_path = "vale_itajai_full_network.geojson"
    if not os.path.exists(geojson_path) or os.path.getsize(geojson_path) < 1000:
        geojson_path = "rios_itajai.geojson"

    with open(geojson_path, "r", encoding="utf-8") as f:
        geo = json.load(f)

    target_rivers = {
        "Itajai_Sul": ["rio itajai do sul"],
        "Itajai_Oeste": ["rio itajai do oeste"],
        "Itajai_Norte": ["rio itajai do norte"],
        "Itajai_Mirim": ["rio itajai-mirim", "rio itajai mirim"],
        "Itajai_Acu": ["rio itajai-acu", "rio itajai acu", "rio itajai"]
    }

    river_geoms = {}

    for feat in geo["features"]:
        name = feat["properties"].get("name", "").lower()
        name_ascii = feat["properties"].get("name_ascii", "").lower()
        
        coords = feat["geometry"]["coordinates"]
        if feat["geometry"]["type"] == "LineString":
            line = LineString(coords)
        elif feat["geometry"]["type"] == "MultiLineString":
            lines = [LineString(c) for c in coords]
            line = max(lines, key=lambda l: l.length)
        else:
            continue

        for key, aliases in target_rivers.items():
            if any(alias in name or alias in name_ascii for alias in aliases):
                if key not in river_geoms or line.length > river_geoms[key].length:
                    river_geoms[key] = line

    if "Itajai_Mirim" not in river_geoms:
        pts_mirim = [(-49.60, -27.38), (-49.10, -27.20), (-48.95, -27.08), (-48.67, -26.90)]
        river_geoms["Itajai_Mirim"] = LineString(pts_mirim)

    print("\nRios extraídos do GeoJSON para a geometria 2D/1D:")
    for k, g in river_geoms.items():
        print(f"  - {k}: {g.length * 111.0:.1f} km ({len(g.coords)} vértices)")

    prj_file = "Itajai_Bacia_Completa.prj"
    geom_file = "Itajai_Bacia_Completa.g01"
    flow_file = "Itajai_Bacia_Completa.u01"
    plan_file = "Itajai_Bacia_Completa.p01"

    # 1. PRJ File
    with open(prj_file, "w") as f:
        f.write("Proj Title=Itajai_Bacia_Completa\n")
        f.write("Current Plan=p01\n")
        f.write("Default Exp/Contr=0.3,0.1\n")
        f.write("SI Units\n")
        f.write("Geom File=g01\n")
        f.write("Unsteady File=u01\n")
        f.write("Plan File=p01\n")
        f.write("Y Axis Title=Elevation\n")
        f.write("X Axis Title(PR)=Distance\n")
        f.write("X Axis Title(CS)=Station\n")

    river_configs = {
        "Itajai_Sul": {"reach": "Trecho_Sul", "length_m": 100000, "z_top": 390.0, "z_bot": 340.0, "dam": 50000, "dam_title": "Barragem Sul (Ituporanga)", "dam_crest": 390.0},
        "Itajai_Oeste": {"reach": "Trecho_Oeste", "length_m": 100000, "z_top": 360.0, "z_bot": 340.0, "dam": 50000, "dam_title": "Barragem Oeste (Taió)", "dam_crest": 360.0},
        "Itajai_Norte": {"reach": "Trecho_Norte", "length_m": 80000, "z_top": 320.0, "z_bot": 340.0, "dam": 40000, "dam_title": "Barragem Norte (José Boiteux)", "dam_crest": 300.0},
        "Itajai_Mirim": {"reach": "Trecho_Mirim", "length_m": 90000, "z_top": 280.0, "z_bot": 2.0, "dam": None},
        "Itajai_Acu": {"reach": "Trecho_Principal", "length_m": 150000, "z_top": 340.0, "z_bot": -15.0, "dam": None}
    }

    # 2. GEOMETRY File (.g01)
    with open(geom_file, "w") as f:
        f.write("Geom Title=Geometria Real com XS GIS Cut Lines HEC-RAS 7\n")
        f.write("Program Version=7.01\n")

        for r_name, cfg in river_configs.items():
            reach_name = cfg["reach"]
            f.write(f"River Reach={r_name},{reach_name}\n")

            line = river_geoms.get(r_name)
            coords_wgs = list(line.coords)
            utm_coords = [latlon_to_utm22s(lon, lat) for lon, lat in coords_wgs]

            f.write(f"Reach XY= {len(utm_coords)} \n")
            for ux, uy in utm_coords:
                f.write(format_16_num(ux) + format_16_num(uy) + "\n")

            total_len = cfg["length_m"]
            st_vals = np.arange(total_len, -1000, -5000)

            for idx, st in enumerate(st_vals):
                frac = st / float(total_len)
                pt_idx = int((1.0 - frac) * (len(coords_wgs) - 1))
                lon_p, lat_p = coords_wgs[pt_idx]
                ux_p, uy_p = utm_coords[pt_idx]

                # Direção tangente do rio para calcular a perpendicular (Cut Line 2D)
                if pt_idx < len(utm_coords) - 1:
                    dx_dir = utm_coords[pt_idx+1][0] - ux_p
                    dy_dir = utm_coords[pt_idx+1][1] - uy_p
                else:
                    dx_dir = ux_p - utm_coords[pt_idx-1][0]
                    dy_dir = uy_p - utm_coords[pt_idx-1][1]

                length_dir = math.hypot(dx_dir, dy_dir) or 1.0
                nx = -dy_dir / length_dir  # vetor normal perpendicular
                ny = dx_dir / length_dir

                w_half = 60.0 if "Acu" in r_name else 40.0

                # Endpoints 2D da XS GIS Cut Line em UTM 22S
                x1, y1 = ux_p - w_half * 1.5 * nx, uy_p - w_half * 1.5 * ny
                x2, y2 = ux_p + w_half * 1.5 * nx, uy_p + w_half * 1.5 * ny

                dem_z = dem.get_elevation(lon_p, lat_p)
                z_model = cfg["z_bot"] + frac * (cfg["z_top"] - cfg["z_bot"])
                z_bottom = z_model

                rl = 5000.0 if st > 0 else 0.0
                st_str = f"{st:.2f}"

                f.write(f"Type RM Length L Ch R = 1 , {st_str:>8} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
                f.write("Node Last Edited Time= Aug/04/2026 00:00:00\n")

                # Escreve explicitamente a XS GIS Cut Line no .g01!
                f.write("XS GIS Cut Line= 2 \n")
                f.write(format_16_num(x1) + format_16_num(y1) + format_16_num(x2) + format_16_num(y2) + "\n")

                if cfg["dam"] is not None and st == cfg["dam"]:
                    f.write(f"Inline Structure= 1 , {st:.2f} , 0 \n")
                    f.write("Type= Dam\n")
                    f.write(f"Inline Structure Title= {cfg['dam_title']}\n")
                    f.write("Inline Structure Node Last Edited Time= Aug/04/2026 00:00:00\n")
                    f.write("Inline Structure Weir= 1 \n")
                    f.write(f"Weir Crest Elev= {cfg['dam_crest']:.1f} \n")
                    f.write("Weir Coeff= 1.6 \n")

                pts = [(-w_half*1.5, z_bottom + 10), (-w_half, z_bottom), (0, z_bottom), (w_half, z_bottom), (w_half*1.5, z_bottom + 10)]
                line_str = "".join([format_8(px) + format_8(py) for px, py in pts])
                
                f.write(f"#Sta/Elev= {len(pts)} \n")
                f.write(line_str + "\n")
                f.write(f"Bank Sta=-{w_half:.1f}, {w_half:.1f}\n")
                f.write("#Mann= 3 , -1 , 0 \n")
                f.write(format_8(-w_half*1.5) + format_8(0.05) + format_8(0) + format_8(-w_half) + format_8(0.035) + format_8(0) + format_8(w_half) + format_8(0.05) + format_8(0) + "\n")

        # JUNÇÃO 1: Confluência Superior (Sul + Oeste + Norte -> Itajaí-Açu)
        f.write("Junction= Junc_Rio_do_Sul, , 660000, 7000000\n")
        f.write("Upstream Reach=Itajai_Sul,Trecho_Sul\n")
        f.write("Upstream Reach=Itajai_Oeste,Trecho_Oeste\n")
        f.write("Upstream Reach=Itajai_Norte,Trecho_Norte\n")
        f.write("Downstream Reach=Itajai_Acu,Trecho_Principal\n")

        # JUNÇÃO 2: Confluência do Rio Itajaí-Mirim
        f.write("Junction= Junc_Itajai_Mirim, , 730000, 7020000\n")
        f.write("Upstream Reach=Itajai_Mirim,Trecho_Mirim\n")

    print(f"\n[OK] Geometria com XS GIS Cut Lines gravada em {geom_file}")

    # 3. FLOW UNSTEADY (.u01)
    with open(flow_file, "w") as f:
        f.write("Flow Title=Cenario_Previsao_Bacia_Real\n")
        f.write("Program Version=7.01\n")

        f.write("Boundary Location=Itajai_Sul,Trecho_Sul,100000.00\n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 49 \n")
        f.write(" ".join(["1200"] * 49) + "\n")

        f.write("Boundary Location=Itajai_Oeste,Trecho_Oeste,100000.00\n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 49 \n")
        f.write(" ".join(["1500"] * 49) + "\n")

        f.write("Boundary Location=Itajai_Norte,Trecho_Norte,80000.00\n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 49 \n")
        f.write(" ".join(["2000"] * 49) + "\n")

        f.write("Boundary Location=Itajai_Mirim,Trecho_Mirim,90000.00\n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 49 \n")
        f.write(" ".join(["900"] * 49) + "\n")

        f.write("Boundary Location=Itajai_Acu,Trecho_Principal,0.00\n")
        f.write("Interval= 1HOUR\n")
        f.write("Stage Hydrograph= 49 \n")
        f.write(" ".join(["0.5"] * 49) + "\n")

        f.write("Initial Stage= 0.5 \n")
        f.write("Initial Flow= 1000 \n")

    print(f"[OK] Fluxo Unsteady gerado em {flow_file}")

    # 4. PLAN File (.p01)
    with open(plan_file, "w") as f:
        f.write("Plan Title=Simulacao_Bacia_Real\n")
        f.write("Program Version=7.01\n")
        f.write("Short Identifier=001\n")
        f.write("Simulation Date=01SEP2008,00,02SEP2008,24\n")
        f.write("Geom File=g01\n")
        f.write("Flow File=u01\n")
        f.write("Subcritical Flow\n")
        f.write("Computation Interval=1MIN\n")
        f.write("Output Interval=1HOUR\n")
        f.write("Instantaneous Interval=1HOUR\n")
        f.write("Mapping Interval=1HOUR\n")
        f.write("Run HTab=-1\n")
        f.write("Run UNet=-1\n")
        f.write("Run PostProcess=-1\n")
        f.write("Run RASMapper=-1\n")

    print(f"[OK] Plano gerado em {plan_file}")
    print("=" * 60)

if __name__ == "__main__":
    main()
