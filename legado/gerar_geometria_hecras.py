#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT: gerar_geometria_hecras.py
DESCRIÇÃO: Amostra elevações do DEM ao longo das seções transversais da rede de
           rios extraída da ANA e gera os arquivos válidos do projeto HEC-RAS:
           - .prj (Arquivo de Projeto HEC-RAS - ASCII/ANSI)
           - .g01 (Geometria 1D Georreferenciada com GIS Cut Lines)
           - .u01 (Condições de Contorno de Vazão Não-Permanente / Hidrograma)
           - .p01 (Plano de Simulação Hidrodinâmica)

REQUISITOS: osgeo.gdal (ou rasterio), shapely, pyproj, numpy
================================================================================
"""

import os
import sys
import json
import numpy as np
import argparse
import unicodedata

# Suporte GDAL (nativo QGIS/Miniforge) e rasterio
GDAL_OK = False
RASTERIO_OK = False

try:
    from osgeo import gdal
    gdal.UseExceptions()
    GDAL_OK = True
except Exception:
    pass

if not GDAL_OK:
    try:
        import rasterio
        RASTERIO_OK = True
    except Exception:
        pass

from shapely.geometry import LineString, MultiLineString, Point, mapping
from shapely.ops import linemerge
from pyproj import Transformer

UTM_EPSG = 31982  # SIRGAS 2000 / UTM 22S (metros)
SPACING = 1000.0   # Espaçamento entre seções transversais (m)
HALFWIDTH = 600.0  # Meia-largura da seção transversal (m)
STEP = 30.0        # Passo de amostragem no DEM (m)
MAXPTS = 120       # Máximo de pontos por seção (limite HEC-RAS = 450)
MIN_SLOPE = 1e-4   # Declividade mínima para monotonicidade do talvegue
N_HORAS = 49       # ordinatas horárias (h=0..48) — deve cobrir TODA a simulação


def fmt_rs(station_m):
    """RS em campo de 8 caracteres alinhado à ESQUERDA (formato oficial HEC-RAS)."""
    return f"{station_m:.1f}"[:8].ljust(8)


def bc_line(river16, reach16, station_m):
    """Linha 'Boundary Location' no formato oficial de 6 campos:
       River(16), Reach(16), RS(8), + 3 campos vazios (8, 16, 16).
       Validado contra os projetos-exemplo oficiais do HEC-RAS; escrever
       menos campos (ou sem padding) faz o RAS associar a BC a uma storage
       area inexistente -> 'needs a upstream boundary condition'."""
    return (f"Boundary Location={river16[:16]:<16},{reach16[:16]:<16},"
            f"{fmt_rs(station_m)},        ,                ,                ")

def sanitizar_ascii(texto):
    """Remove totalmente acentos e caracteres não-ASCII para total compatibilidade HEC-RAS."""
    if not texto: return "RIO"
    nfkd = unicodedata.normalize('NFKD', str(texto))
    ascii_str = "".join([c for c in nfkd if not unicodedata.combining(c)])
    clean = "".join([c if (c.isalnum() or c in "_- ") else "_" for c in ascii_str])
    return clean.strip().replace(" ", "_")

def fixed_series_8(vals):
    """Retorna valores formatados em colunas fixas de 8 caracteres (10 por linha) para estações/elevações e hidrogramas."""
    return "\n".join("".join(f"{v:8.2f}" for v in vals[i:i+10]) for i in range(0, len(vals), 10))

def fixed_series_16(vals):
    """Retorna coordenadas UTM formatadas em colunas fixas de 16 caracteres (5 por linha) para Reach Header e Cut Lines."""
    return "\n".join("".join(f"{v:16.2f}" for v in vals[i:i+5]) for i in range(0, len(vals), 5))

def gerar_hidrograma(vazao_pico, base=None, n_horas=49, tp=18):
    base = base if base is not None else max(vazao_pico * 0.15, 50)
    horas = np.arange(n_horas)
    # Curva suave de hidrograma gaussiana/gama horária (48h)
    sigma = tp / 2.5
    q_raw = np.exp(-((horas - tp) ** 2) / (2 * sigma ** 2))
    q = base + (vazao_pico - base) * q_raw
    return [round(float(v), 2) for v in q]

class DemSampler:
    def __init__(self, dem_path):
        self.dem_path = dem_path
        if GDAL_OK:
            self.mode = "gdal"
            self.ds = gdal.Open(dem_path)
            self.gt = self.ds.GetGeoTransform()
            self.band = self.ds.GetRasterBand(1)
            self.arr = self.band.ReadAsArray()
            self.nodata = self.band.GetNoDataValue()
            self.rows, self.cols = self.arr.shape
            self.to_dem = Transformer.from_crs(UTM_EPSG, 4326, always_xy=True)
        elif RASTERIO_OK:
            # NOTA: NAO usar ds.sample() — em rasterio 1.5.x essa chamada
            # derruba o processo (crash nativo, sem traceback). Carregamos a
            # banda inteira em memoria (~155 MB p/ o DEM da bacia) e indexamos
            # direto pela transform afim: mais robusto e muito mais rapido.
            self.mode = "rasterio"
            self.ds = rasterio.open(dem_path)
            self.nodata = self.ds.nodata
            self.arr = self.ds.read(1)
            self.rows, self.cols = self.arr.shape
            self.inv = ~self.ds.transform
            self.to_dem = Transformer.from_crs(UTM_EPSG, self.ds.crs.to_epsg(), always_xy=True)
        else:
            raise RuntimeError("Nem GDAL nem Rasterio estão disponíveis para amostragem do DEM.")

    def sample(self, xs_utm, ys_utm):
        lons, lats = self.to_dem.transform(xs_utm, ys_utm)
        if self.mode == "gdal":
            vals = []
            gt = self.gt
            for lon, lat in zip(lons, lats):
                px = int((lon - gt[0]) / gt[1])
                py = int((lat - gt[3]) / gt[5])
                if 0 <= px < self.cols and 0 <= py < self.rows:
                    val = float(self.arr[py, px])
                    if self.nodata is not None and val == self.nodata:
                        val = np.nan
                    vals.append(val)
                else:
                    vals.append(np.nan)
            return np.array(vals)
        else:
            lons = np.asarray(lons, dtype=float)
            lats = np.asarray(lats, dtype=float)
            a, b, c, d, e, f = (self.inv.a, self.inv.b, self.inv.c,
                                self.inv.d, self.inv.e, self.inv.f)
            cols = np.floor(a * lons + b * lats + c).astype(int)
            rows = np.floor(d * lons + e * lats + f).astype(int)
            dentro = ((rows >= 0) & (rows < self.rows) &
                      (cols >= 0) & (cols < self.cols))
            vals = np.full(lons.shape, np.nan, dtype=float)
            vals[dentro] = self.arr[rows[dentro], cols[dentro]]
            if self.nodata is not None:
                vals[vals == self.nodata] = np.nan
            vals[vals < -1000] = np.nan          # sentinelas de oceano/void
            return vals

def cortar_secao_dem(line, dist_m, dem_sampler):
    pt = line.interpolate(dist_m)
    d_delta = 1.0
    pt_prev = line.interpolate(max(0, dist_m - d_delta))
    pt_next = line.interpolate(min(line.length, dist_m + d_delta))
    
    dx = pt_next.x - pt_prev.x
    dy = pt_next.y - pt_prev.y
    length = np.hypot(dx, dy)
    if length == 0: return None
    
    nx = -dy / length
    ny = dx / length
    
    offsets = np.linspace(-HALFWIDTH, HALFWIDTH, MAXPTS)
    xs_utm = pt.x + nx * offsets
    ys_utm = pt.y + ny * offsets
    
    z_raw = dem_sampler.sample(xs_utm, ys_utm)
    if np.isnan(z_raw).all(): return None
    
    valid_idx = np.where(~np.isnan(z_raw))[0]
    if len(valid_idx) < 5: return None
    
    sta = offsets[valid_idx] - offsets[valid_idx[0]]
    z = z_raw[valid_idx]
    z = np.where(np.isnan(z), np.nan_to_num(z, nan=np.nanmean(z)), z)
    
    # Garantir talvegue bem definido
    mid = len(z) // 2
    z[mid-2:mid+3] = np.minimum(z[mid-2:mid+3], np.min(z) - 1.5)
    
    cutline = (xs_utm[valid_idx[0]], ys_utm[valid_idx[0]], xs_utm[valid_idx[-1]], ys_utm[valid_idx[-1]])
    return sta, z, cutline

def estacoes_margem(sta, z, altura_margem=3.0):
    """Margens TOPOGRAFICAS: parte do talvegue (ponto mais baixo) e sobe ate
    'altura_margem' metros de cada lado.

    IMPORTANTE: o valor devolvido tem de ser EXATAMENTE um dos valores da
    tabela #Sta/Elev, com a MESMA precisao com que ela e escrita (2 casas),
    senao o HEC-RAS recusa: "Left bank station not in station elevation data".
    """
    z = np.asarray(z, dtype=float)
    i0 = int(np.nanargmin(z))
    limiar = z[i0] + altura_margem
    li = i0
    while li > 0 and z[li] < limiar:
        li -= 1
    ri = i0
    while ri < len(z) - 1 and z[ri] < limiar:
        ri += 1
    # margens estritamente interiores e distintas
    li = min(max(li, 1), len(sta) - 3)
    ri = max(min(ri, len(sta) - 2), li + 1)
    # arredonda para a MESMA precisao usada em #Sta/Elev (fixed_series_8 -> .2f)
    return round(float(sta[li]), 2), round(float(sta[ri]), 2)

def processar_geometria_bacia(nome_bacia, pasta_trabalho):
    bacia_ascii = sanitizar_ascii(nome_bacia)
    prj_name = f"{bacia_ascii}_Oficial"
    nome_limpo_arq = sanitizar_ascii(nome_bacia).lower()
    
    caminho_geojson = os.path.join(pasta_trabalho, f"rios_{nome_limpo_arq}.geojson")
    caminho_dem = os.path.join(pasta_trabalho, f"dem_{nome_limpo_arq}.tif")
    
    if not os.path.exists(caminho_geojson):
        caminho_geojson = os.path.join(pasta_trabalho, f"rios_{nome_limpo_arq}_original.geojson")
        
    if not os.path.exists(caminho_geojson) or not os.path.exists(caminho_dem):
        print(f"  ⚠ Não foi possível localizar o GeoJSON ({caminho_geojson}) ou o DEM ({caminho_dem}).")
        return False

    print(f"\n[1/3] Amostrando seções transversais do DEM 30m para o {bacia_ascii}...")
    dem_sampler = DemSampler(caminho_dem)

    with open(caminho_geojson, "r", encoding="utf-8") as f:
        geojson_data = json.load(f)

    to_utm = Transformer.from_crs(4326, UTM_EPSG, always_xy=True)

    lines_utm = []
    for feat in geojson_data.get("features", []):
        geom = feat.get("geometry")
        if not geom: continue
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])
        if gtype == "LineString":
            lines_utm.append(LineString([to_utm.transform(pt[0], pt[1]) for pt in coords]))
        elif gtype == "MultiLineString":
            for seg in coords:
                lines_utm.append(LineString([to_utm.transform(pt[0], pt[1]) for pt in seg]))

    merged = linemerge(lines_utm)
    if isinstance(merged, MultiLineString):
        main_line = max(merged.geoms, key=lambda l: l.length)
    else:
        main_line = merged

    reach_coords = list(main_line.coords)
    if len(reach_coords) > 200:
        idx = np.linspace(0, len(reach_coords) - 1, 200).round().astype(int)
        reach_coords = [reach_coords[i] for i in idx]

    n_sec = max(10, int(main_line.length // SPACING))
    dists = np.linspace(0, main_line.length, n_sec)
    
    xs_list = []
    for d in dists:
        res = cortar_secao_dem(main_line, d, dem_sampler)
        if res:
            sta, z, cutline = res
            lb, rb = estacoes_margem(sta, z)
            xs_list.append({
                "station_m": round(main_line.length - d, 1),
                "sta": sta,
                "z": z,
                "lb": lb,
                "rb": rb,
                "cutline": cutline
            })

    print(f"  ✓ Total de {len(xs_list)} seções transversais georreferenciadas no canal principal.")

    # WKT oficial para SIRGAS 2000 / UTM Zone 22S (EPSG:31982)
    wkt_crs = 'PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1]]'

    # Calcular min/max UTM reais das seções e da calha para o Viewing Rectangle do HEC-RAS
    all_x = [pt[0] for pt in reach_coords] + [xs['cutline'][0] for xs in xs_list] + [xs['cutline'][2] for xs in xs_list]
    all_y = [pt[1] for pt in reach_coords] + [xs['cutline'][1] for xs in xs_list] + [xs['cutline'][3] for xs in xs_list]
    min_x, max_x = min(all_x) - 1000.0, max(all_x) + 1000.0
    min_y, max_y = min(all_y) - 1000.0, max(all_y) + 1000.0

    # 2. Gerar Arquivo de Geometria HEC-RAS (.g01)
    g01_file = os.path.join(pasta_trabalho, f"{prj_name}.g01")
    print(f"\n[2/3] Gerando arquivo de geometria HEC-RAS: {os.path.basename(g01_file)}...")
    
    with open(g01_file, "w", encoding="ascii", errors="replace") as f:
        f.write(f"Geom Title={prj_name}\n")
        f.write("Program Version=7.01\n")
        f.write("Viewing Rectangle= 0 , 1 , 1 , 0 \n")
        f.write(f"Spatial Reference System={wkt_crs}\n\n")
        
        rv, rc = f"{bacia_ascii.upper()[:16]:<16}", "CALHA_PRINCIPAL "
        f.write(f"River Reach={rv},{rc}\n")
        f.write(f"Reach XY= {len(reach_coords)} \n")
        for i in range(0, len(reach_coords), 2):
            pair1 = reach_coords[i]
            if i + 1 < len(reach_coords):
                pair2 = reach_coords[i+1]
                f.write(f"{pair1[0]:16.4f}{pair1[1]:16.4f}{pair2[0]:16.4f}{pair2[1]:16.4f}\n")
            else:
                f.write(f"{pair1[0]:16.4f}{pair1[1]:16.4f}\n")
        f.write("\n")

        for i_xs, xs in enumerate(xs_list):
            # Comprimento ATE a proxima secao de jusante; a ultima secao deve ser 0,0,0
            if i_xs < len(xs_list) - 1:
                dx = round(abs(xs['station_m'] - xs_list[i_xs + 1]['station_m']), 1)
            else:
                dx = 0.0
            f.write(f"Type RM Length L Ch R = 1 ,{xs['station_m']:.1f},{dx:.1f},{dx:.1f},{dx:.1f}\n")
            f.write("XS GIS Cut Line=2\n")
            cl = xs['cutline']
            f.write(f"     {cl[0]:11.4f}     {cl[1]:11.4f}     {cl[2]:11.4f}     {cl[3]:11.4f}\n")
            f.write("Node Last Edited Time=Aug/06/2026 00:00:00\n")
            sta, z = xs["sta"], xs["z"]
            f.write(f"#Sta/Elev= {len(sta)} \n")
            pts_pair = []
            for s_val, z_val in zip(sta, z):
                pts_pair.extend([s_val, z_val])
            f.write(fixed_series_8(pts_pair) + "\n")
            f.write("#Mann= 3 ,-1,0\n")
            # 9 campos de EXATAMENTE 8 caracteres (sta, n, 0) x3
            f.write("".join(f"{v:>8}" for v in
                            [f"{float(sta[0]):.2f}", ".05", "0",
                             f"{xs['lb']:.2f}", ".035", "0",
                             f"{xs['rb']:.2f}", ".05", "0"]) + "\n")
            # MESMA precisao de #Sta/Elev (.2f) — obrigatorio p/ o RAS casar
            f.write(f"Bank Sta={xs['lb']:.2f},{xs['rb']:.2f}\n")
            f.write("XS Rating Curve= 0 ,0\n")
            f.write("XS HTab Horizontal Distribution= 5 , 5 , 5 \n")
            f.write("Exp/Cntr=0.3,0.1\n\n")

    # 3. Gerar Arquivo de Projeto (.prj), Vazão (.u01), Plano (.p01) e RASMap (.rasmap) em ASCII puro
    prj_file = os.path.join(pasta_trabalho, f"{prj_name}.prj")
    u01_file = os.path.join(pasta_trabalho, f"{prj_name}.u01")
    p01_file = os.path.join(pasta_trabalho, f"{prj_name}.p01")
    rasmap_file = os.path.join(pasta_trabalho, f"{prj_name}.rasmap")
    
    with open(prj_file, "w", encoding="ascii", errors="replace") as f:
        f.write(f"Proj Title={prj_name}\n")
        f.write("Current Plan=p01\n")
        f.write("Default Directory=\n")
        f.write("Default Exp/Contr=0.3,0.1\n")
        f.write("SI Units\n")          # CRITICO: sem isto o HEC-RAS assume pes
        f.write(f"Geom File=g01\n")
        f.write(f"Unsteady File=u01\n")
        f.write(f"Plan File=p01\n")
        f.write("Y Axis Title=Elevation\n")
        f.write("X Axis Title(PR)=Distance\n")
        f.write("X Axis Title(CS)=Station\n")
        f.write(f"RASMap Filename={prj_name}.rasmap\n")

    # Criar arquivo de projeção WKT para RASMapper (.projection sidecar)
    proj_sidecar = os.path.join(pasta_trabalho, f"{prj_name}.projection")
    with open(proj_sidecar, "w", encoding="utf-8") as f:
        f.write(wkt_crs)

    # Gerar arquivo .rasmap válido para o RASMapper do HEC-RAS com a projeção espacial
    with open(rasmap_file, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n')
        f.write('<RASMapper>\n')
        f.write('  <Version>2.00</Version>\n')
        f.write(f'  <RASProjectionFilename Filename="{prj_name}.projection" />\n')
        f.write('  <RASData>\n')
        f.write(f'    <Layer Name="Terrain" Filename="dem_{nome_limpo_arq}.tif" />\n')
        f.write('  </RASData>\n')
        f.write('</RASMapper>\n')
        
    with open(u01_file, "w", encoding="ascii", errors="replace") as f:
        f.write(f"Flow Title=Cheia_Sintetica_{bacia_ascii}\n")
        f.write("Program Version=7.01\n")
        f.write("Use Restart= 0 \n")
        if xs_list:
            vazao_pico = 1500.0
            hidro = gerar_hidrograma(vazao_pico, n_horas=N_HORAS)
            base_ini = hidro[0]

            # Condicao inicial (obrigatoria no unsteady) no RS de montante
            f.write(f"Initial Flow Loc={rv},{rc},"
                    f"{fmt_rs(xs_list[0]['station_m'])},{base_ini:.0f}\n")

            # Contorno de MONTANTE: hidrograma de vazao
            f.write(bc_line(rv, rc, xs_list[0]['station_m']) + "\n")
            f.write("Interval=1HOUR\n")
            f.write(f"Flow Hydrograph= {len(hidro)} \n")
            f.write(fixed_series_8(hidro) + "\n")

            # Contorno de JUSANTE: profundidade normal (Friction Slope)
            f.write(bc_line(rv, rc, xs_list[-1]['station_m']) + "\n")
            f.write("Friction Slope=0.001\n")

    with open(p01_file, "w", encoding="ascii", errors="replace") as f:
        f.write(f"Plan Title=Plano_Simulacao_Inundacao\n")
        f.write("Program Version=7.01\n")
        f.write(f"Short Identifier=01\n")
        f.write("Geom File=g01\n")
        f.write("Flow File=u01\n")
        f.write("Unsteady File=u01\n")
        # A serie tem N_HORAS ordinatas (h=0..N_HORAS-1); a simulacao deve
        # terminar em h=N_HORAS-1, senao: "Time series data ends before the
        # end of the simulation."
        f.write("Simulation Date=01AUG2026,0000,03AUG2026,0000\n")
        f.write("Subcritical Flow\n")
        f.write("Computation Interval=1MIN\n")
        f.write("Output Interval=1HOUR\n")
        f.write("Instantaneous Interval=1HOUR\n")
        f.write("Mapping Interval=1HOUR\n")
        f.write("Run HTab=-1\n")
        f.write("Run UNet=-1\n")
        f.write("Run Sediment= 0 \n")
        f.write("Run PostProcess=-1\n")
        f.write("Run RASMapper=0\n")

    print(f"\n[3/3] Projeto HEC-RAS criado com sucesso: {prj_name}.prj")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gera geometria HEC-RAS 1D a partir de rios e DEM.")
    parser.add_argument("-b", "--bacia", type=str, default="Itajaí", help="Nome do rio/bacia")
    parser.add_argument("-o", "--saida", type=str, default=r"C:\Users\haas\github\hecras_vale", help="Pasta de trabalho")
    
    args = parser.parse_args()
    processar_geometria_bacia(args.bacia, args.saida)
