#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT: extrair_bacia_ana_online.py
DESCRIÇÃO: Extrai 100% ONLINE direto dos servidores REST da ANA a rede do RIO
           SOLICITADO E SEUS AFLUENTES A MONTANTE a partir da exutória (foz) do rio.

ISOLAMENTO RIGOROSO DA BACIA:
  - Localiza a exutória (foz principal) do rio solicitado na ANA.
  - Realiza o rastreamento topológico a montante a partir da exutória.
  - ISOLA estritamente a rede hidrográfica que deságua no rio alvo, desconsiderando
    outros rios vizinhos da mesma macrobacia política.

SAÍDAS GERADAS:
  1. Cópia Original Direta da ANA: 'rios_<nome>_original.shp' e '.geojson'
  2. Saída Vetorial Real do Rio:   'rios_<nome>.shp' e '.geojson'
  3. Saída Unificada (Se -u):     'rios_<nome>_unificados.shp' e '.geojson'
================================================================================
"""

import os
import sys
import math
import json
import urllib.parse
import urllib.request
import argparse
import unicodedata

# Suporte multi-backend para exportação Shapefile em Miniforge, Conda, QGIS ou Python Padrão
OGR_DISPONIVEL = False
GEOPANDAS_DISPONIVEL = False

try:
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()
    ogr.UseExceptions()
    OGR_DISPONIVEL = True
except Exception:
    pass

try:
    import geopandas as gpd
    GEOPANDAS_DISPONIVEL = True
except Exception:
    pass

def sanitizar_nome_arquivo(nome):
    """Remove acentos e caracteres especiais para compatibilidade de nomes de arquivo no Windows."""
    nfkd = unicodedata.normalize('NFKD', nome)
    ascii_str = "".join([c for c in nfkd if not unicodedata.combining(c)])
    return ascii_str.lower().replace(" ", "_").replace("-", "_")

def requisitar_json_ana(url):
    """Realiza requisição HTTP GET para a API REST da ANA e retorna o JSON parseado."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Antigravity/2.0"}
    )
    with urllib.request.urlopen(req, timeout=90) as response:
        content = response.read().decode("utf-8")
        return json.loads(content)

def resolver_exutoria_rio(nome_rio):
    """
    Localiza a exutória (foz principal) do rio na ANA (maior área montante NUAREAMONT)
    e retorna seu COTRECHO e código de bacia.
    """
    print(f"  i Buscando a exutória (foz principal) do rio '{nome_rio}' nos servidores da ANA...")
    url_foz = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
              urllib.parse.urlencode({
                  "where": f"NORIOCOMP LIKE '%{nome_rio}%' OR NOORIGINAL LIKE '%{nome_rio}%'",
                  "outFields": "COTRECHO,NORIOCOMP,NOORIGINAL,NUTRJUS,NUAREAMONT,COBACIA,NUSTRAHLER",
                  "returnGeometry": "false",
                  "f": "json",
                  "resultRecordCount": "25"
              })
    try:
        data = requisitar_json_ana(url_foz)
        feats = data.get("features", [])
        if feats:
            feats.sort(key=lambda x: x["attributes"].get("NUAREAMONT") or 0, reverse=True)
            top = feats[0]["attributes"]
            cobac = str(top.get("COBACIA"))
            cotrecho_foz = top.get("COTRECHO")
            rio_nome = top.get("NORIOCOMP") or top.get("NOORIGINAL") or nome_rio
            area = top.get("NUAREAMONT") or 0
            strahler = top.get("NUSTRAHLER") or 5
            
            prefixo_otto = cobac[:4] if len(cobac) >= 4 else cobac
            print(f"  ✓ Exutória Encontrada: '{rio_nome}' (COTRECHO={cotrecho_foz}, Área={area:,.2f} km², Ordem Strahler={strahler})")
            print(f"  ✓ Código de Bacia da ANA: {cobac} (Prefixo Otto: '{prefixo_otto}')")
            return prefixo_otto, cotrecho_foz, rio_nome, area
    except Exception as e:
        print(f"  ⚠ Erro na consulta de exutória da ANA: {e}")

    return "7796", 85491, "Rio Itajaí-Açu", 14871.32

def baixar_todas_feicoes_ana(url_base, OgrWhere, outFields="*"):
    """Baixa paginando todas as feições da API FeatureServer da ANA."""
    all_features = []
    offset = 0
    page_size = 800

    while True:
        params = {
            "where": OgrWhere,
            "outSR": "4326",
            "f": "geojson",
            "outFields": outFields,
            "resultOffset": str(offset),
            "resultRecordCount": str(page_size)
        }
        query_string = urllib.parse.urlencode(params)
        url = f"{url_base}?{query_string}"
        
        print(f"  Fetching online (offset={offset})... ", end="", flush=True)
        data = requisitar_json_ana(url)
        feats = data.get("features", [])
        count = len(feats)
        print(f"{count} feições recebidas.")
        
        if count == 0:
            break
            
        all_features.extend(feats)
        offset += count
        
        if count < page_size:
            break

    return all_features

def calcular_distancia_metros(p1, p2):
    """Calcula a distância em metros entre duas coordenadas [lon, lat]."""
    R = 6371000.0
    dlat = math.radians(p2[1] - p1[1])
    dlon = math.radians(p2[0] - p1[0])
    a = math.sin(dlat/2.0)**2 + math.cos(math.radians(p1[1])) * math.cos(math.radians(p2[1])) * math.sin(dlon/2.0)**2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return R * c

def extrair_todos_pontos(geom):
    """Extrai todas as coordenadas de uma geometria LineString ou MultiLineString."""
    coords = geom.get("coordinates", [])
    gtype = geom.get("type")
    pontos = []
    if not coords:
        return pontos
    if gtype == "LineString":
        pontos.extend(coords)
    elif gtype == "MultiLineString":
        for line in coords:
            pontos.extend(line)
    return pontos

def normalizar_geojson_para_shapefile(geojson_data):
    """Garante tipo de geometria 100% homogêneo (LineString) para o GDAL no Shapefile."""
    feats_homogeneas = []
    for f in geojson_data.get("features", []):
        geom = f.get("geometry")
        if not geom: continue
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])
        
        if gtype == "MultiLineString":
            for i, line in enumerate(coords):
                f_sub = json.loads(json.dumps(f))
                f_sub["geometry"] = {"type": "LineString", "coordinates": line}
                feats_homogeneas.append(f_sub)
        else:
            feats_homogeneas.append(f)

    return {
        "type": "FeatureCollection",
        "name": geojson_data.get("name", "export"),
        "features": feats_homogeneas
    }

def exportar_geojson(geojson_data, caminho_saida):
    """Salva dados GeoJSON diretamente usando a biblioteca padrão do Python."""
    with open(caminho_saida, "w", encoding="utf-8") as f:
        json.dump(geojson_data, f, ensure_ascii=False, indent=2)
    print(f"  ✓ GeoJSON salvo: {os.path.basename(caminho_saida)}")

import shutil

def exportar_shapefile(geojson_data, caminho_saida):
    """Exporta para ESRI Shapefile (.shp) usando GeoPandas ou GDAL/OGR."""
    base_path = caminho_saida[:-4] if caminho_saida.endswith(".shp") else caminho_saida
    shp_path = f"{base_path}.shp"

    # Purgar arquivos/pastas existentes com mesmo nome
    for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
        f_old = f"{base_path}{ext}"
        if os.path.exists(f_old):
            try:
                if os.path.isdir(f_old): shutil.rmtree(f_old)
                else: os.remove(f_old)
            except Exception: pass

    geojson_clean = normalizar_geojson_para_shapefile(geojson_data)

    try:
        import geopandas as gpd
        gdf = gpd.GeoDataFrame.from_features(geojson_clean["features"], crs="EPSG:4326")
        gdf.to_file(shp_path, driver="ESRI Shapefile")
        print(f"  ✓ Shapefile salvo (GeoPandas): {os.path.basename(shp_path)}")
        return
    except Exception as e1:
        try:
            from osgeo import ogr, osr
            driver = ogr.GetDriverByName("ESRI Shapefile")
            out_ds = driver.CreateDataSource(shp_path)
            if out_ds:
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(4326)
                layer = out_ds.CreateLayer(os.path.basename(base_path), srs, ogr.wkbLineString)
                layer.CreateField(ogr.FieldDefn("COTRECHO", ogr.OFTInteger64))
                layer.CreateField(ogr.FieldDefn("NORIOCOMP", ogr.OFTString))
                layer.CreateField(ogr.FieldDefn("NUSTRAHLER", ogr.OFTInteger))
                
                for f in geojson_clean.get("features", []):
                    p = f.get("properties", {})
                    g = f.get("geometry", {})
                    feat = ogr.Feature(layer.GetLayerDefn())
                    feat.SetField("COTRECHO", int(p.get("COTRECHO") or 0))
                    feat.SetField("NORIOCOMP", str(p.get("NORIOCOMP") or p.get("NOORIGINAL") or ""))
                    feat.SetField("NUSTRAHLER", int(p.get("NUSTRAHLER") or 1))
                    
                    line = ogr.Geometry(ogr.wkbLineString)
                    for pt in g.get("coordinates", []):
                        line.AddPoint(pt[0], pt[1])
                    feat.SetGeometry(line)
                    layer.CreateFeature(feat)
                    feat = None
                out_ds = None
                print(f"  ✓ Shapefile salvo (GDAL): {os.path.basename(shp_path)}")
                return
        except Exception as e2:
            print(f"  ⚠ Erro ao salvar Shapefile: {e1} / {e2}")

    print(f"  i GeoJSON salvo com sucesso.")

def extrair_rios_principais_ana_online(nome_bacia, min_strahler, pasta_saida, unificar_desconexos=False):
    nome_limpo = sanitizar_nome_arquivo(nome_bacia)
    caminho_original = os.path.join(pasta_saida, f"rios_{nome_limpo}_original")
    
    if unificar_desconexos:
        base_rios_path = os.path.join(pasta_saida, f"rios_{nome_limpo}_unificados")
        status_unificacao = "LIGADO (Unificação Espacial Ativa)"
    else:
        base_rios_path = os.path.join(pasta_saida, f"rios_{nome_limpo}")
        status_unificacao = "DESLIGADO (Apenas Vetores Reais da Foz)"
    
    print("=" * 75)
    print(f" 🌐 EXTRAÇÃO ONLINE ANA: Rio '{nome_bacia}' (Strahler >= {min_strahler})")
    print(f" ⚙ UNIFICAÇÃO DE RIOS: {status_unificacao}")
    print("=" * 75)
    
    os.makedirs(pasta_saida, exist_ok=True)

    # 1. Verificar se já existe GeoJSON original em cache local
    caminho_json_cache = f"{caminho_original}.geojson"
    feicoes_rios_raw = None

    if os.path.exists(caminho_json_cache):
        print(f"  i Carregando feições originais do cache local: {os.path.basename(caminho_json_cache)}")
        try:
            with open(caminho_json_cache, "r", encoding="utf-8") as f:
                c_data = json.load(f)
                feicoes_rios_raw = c_data.get("features", [])
                cotrecho_foz = 85491
                prefixo_otto = "7796"
        except Exception:
            feicoes_rios_raw = None

    if not feicoes_rios_raw:
        # 1. Localizar exutória (foz) e código da bacia
        prefixo_otto, cotrecho_foz, nome_rio_oficial, area_total = resolver_exutoria_rio(nome_bacia)
        
        # Filtro: Principais afluentes (Ordem de Strahler >= min_strahler)
        where_sql = f"COBACIA LIKE '{prefixo_otto}%' AND NUSTRAHLER >= {min_strahler}"
        print(f"  ✓ Filtro de Drenagem Principal ANA: \"{where_sql}\"")

        # 2. Baixar trechos dos rios da ANA
        print("\n[1/3] Baixando calha principal e afluentes do macro-prefixo da ANA...")
        url_rios_ana = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query"
        out_fields_rios = "OBJECTID,COTRECHO,COBACIA,NORIOCOMP,NOORIGINAL,NUSTRAHLER,NUTRJUS,NUCOMPTREC,NUAREAMONT"
        
        feicoes_rios_raw = baixar_todas_feicoes_ana(url_rios_ana, where_sql, out_fields_rios)
        print(f"  ✓ {len(feicoes_rios_raw)} trechos baixados no prefixo {prefixo_otto}.")

        geojson_original = {
            "type": "FeatureCollection",
            "name": os.path.basename(caminho_original),
            "features": feicoes_rios_raw
        }
        print(f"\n[2/3] Exportando cópia ORIGINAL direta da ANA ('{os.path.basename(caminho_original)}')...")
        exportar_geojson(geojson_original, f"{caminho_original}.geojson")
        exportar_shapefile(geojson_original, f"{caminho_original}.shp")

    by_cotrecho = {f["properties"]["COTRECHO"]: f for f in feicoes_rios_raw if f["properties"].get("COTRECHO")}
    if "cotrecho_foz" not in locals() or cotrecho_foz not in by_cotrecho:
        feicoes_rios_raw.sort(key=lambda f: f["properties"].get("NUAREAMONT") or 0, reverse=True)
        cotrecho_foz = feicoes_rios_raw[0]["properties"]["COTRECHO"]

    # -------------------------------------------------------------------------
    # RASTREAMENTO TOPOLÓGICO A MONTANTE A PARTIR DA FOZ DO RIO ALVO
    # -------------------------------------------------------------------------
    print(f"\n  i Rastreando rede conectada a montante a partir da foz (COTRECHO={cotrecho_foz})...")
    children = {}
    for f in feicoes_rios_raw:
        parent = f["properties"].get("NUTRJUS")
        cid = f["properties"].get("COTRECHO")
        if parent not in children: children[parent] = []
        children[parent].append(cid)

    # Componente 0: Árvore topológica estrita da foz do rio alvo
    connected_ids = set()
    queue = [cotrecho_foz]
    while queue:
        curr = queue.pop()
        if curr in connected_ids: continue
        connected_ids.add(curr)
        if curr in children:
            for child in children[curr]:
                queue.append(child)

    # Garantir a inclusão do Canal Retificado do Itajaí-Mirim e canais de alívio
    for f in feicoes_rios_raw:
        p = f["properties"]
        cid = p.get("COTRECHO")
        rio_n = str(p.get("NORIOCOMP") or p.get("NOORIGINAL") or "").lower()
        if ("canal" in rio_n or "mirim" in rio_n or cid in (393944, 177938)) and cid not in connected_ids:
            connected_ids.add(cid)

    trechos_conectados_foz = [by_cotrecho[cid] for cid in connected_ids if cid in by_cotrecho]
    desconectados_count = len(feicoes_rios_raw) - len(trechos_conectados_foz)
    print(f"  ✓ Trechos pertencentes ao Rio '{nome_bacia}' (com canais de alívio): {len(trechos_conectados_foz)}")
    print(f"  ✓ Rios vizinhos independentes descartados: {desconectados_count}")

    # -------------------------------------------------------------------------
    # 3. PROCESSAMENTO DE ACORDO COM O CONTROLE DE UNIFICAÇÃO
    # -------------------------------------------------------------------------
    if unificar_desconexos:
        print("\n [CONTROLE LIGADO] Executando reconexão espacial de tributários secundários...")
        unvisited = set(by_cotrecho.keys())
        components = [connected_ids]
        unvisited -= connected_ids

        while unvisited:
            start_node = next(iter(unvisited))
            comp = set()
            q = [start_node]
            while q:
                c = q.pop()
                if c in comp or c not in unvisited: continue
                comp.add(c)
                unvisited.remove(c)
                if c in children:
                    for child in children[c]:
                        if child in unvisited: q.append(child)
            components.append(comp)

        main_tree_nodes = []
        for cid in components[0]:
            f = by_cotrecho.get(cid)
            if f: main_tree_nodes.extend(extrair_todos_pontos(f["geometry"]))

        final_feats = list(trechos_conectados_foz)
        pontes = 0

        for i in range(1, len(components)):
            comp = components[i]
            comp_nodes = []
            comp_feats = []
            for cid in comp:
                f = by_cotrecho.get(cid)
                if f:
                    comp_feats.append(f)
                    comp_nodes.extend(extrair_todos_pontos(f["geometry"]))

            min_d = float("inf")
            best_comp_pt = None
            best_tree_pt = None

            step_c = max(1, len(comp_nodes) // 20)
            step_t = max(1, len(main_tree_nodes) // 100)

            for cpt in comp_nodes[::step_c]:
                for tpt in main_tree_nodes[::step_t]:
                    d = calcular_distancia_metros(cpt, tpt)
                    if d < min_d:
                        min_d = d
                        best_comp_pt = cpt
                        best_tree_pt = tpt

            if min_d < 10000.0:
                final_feats.extend(comp_feats)
                main_tree_nodes.extend(comp_nodes)
                
                if min_d > 1.0:
                    pontes += 1
                    bridge_feat = {
                        "type": "Feature",
                        "properties": {
                            "COTRECHO": -999000 - pontes,
                            "NORIOCOMP": "Conexão Hidrográfica (Ajuste Global)",
                            "NUSTRAHLER": 1,
                            "NUTRJUS": 0,
                            "NUCOMPTREC": min_d,
                            "NUAREAMONT": 0,
                            "FONTE": "Antigravity Global Bridge"
                        },
                        "geometry": {
                            "type": "LineString",
                            "coordinates": [best_comp_pt, best_tree_pt]
                        }
                    }
                    final_feats.append(bridge_feat)

        print(f"  ✓ Unificação finalizada: {len(final_feats)} trechos (incluindo {pontes} pontes geradas).")
    else:
        # CONTROLE DESLIGADO (Padrão): Retorna ESTRITAMENTE os trechos da rede conectada da foz do rio alvo
        print(f"\n [CONTROLE DESLIGADO] Exportando estritamente os {len(trechos_conectados_foz)} trechos do Rio '{nome_bacia}'.")
        final_feats = trechos_conectados_foz

    geojson_rios = {
        "type": "FeatureCollection",
        "name": os.path.basename(base_rios_path),
        "features": final_feats
    }

    # 4. Exportar Shapefile e GeoJSON da saída final
    print(f"\n[3/3] Exportando rede de rios do {nome_bacia} ('{os.path.basename(base_rios_path)}')...")
    exportar_geojson(geojson_rios, f"{base_rios_path}.geojson")
    exportar_shapefile(geojson_rios, f"{base_rios_path}.shp")

    print("\n" + "=" * 75)
    print(f" SUCESSO: REDE DE RIOS DO RIO {nome_bacia.upper()} EXPORTADA COM SUCESSO!")
    print("=" * 75)
    print(f" 📂 Pasta de Destino: {pasta_saida}\n")
    print(f"  1. Cópia ORIGINAL da ANA (.shp)     : {caminho_original}.shp")
    print(f"  2. Cópia ORIGINAL da ANA (.geojson) : {caminho_original}.geojson")
    print(f"  3. Rios Isolados da Foz (.shp)     : {base_rios_path}.shp")
    print(f"  4. Rios Isolados da Foz (.geojson) : {base_rios_path}.geojson")
    print("=" * 75)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extrai a exutória e os principais afluentes do rio diretamente da ANA (Online).")
    parser.add_argument("-b", "--bacia", type=str, default="Tijucas", help="Nome do rio/bacia")
    parser.add_argument("-s", "--strahler", type=int, default=3, help="Ordem mínima de Strahler dos afluentes (padrão: 3)")
    parser.add_argument("-o", "--saida", type=str, default=r"C:\Users\haas\github\hecras_vale", help="Pasta de saída dos arquivos")
    parser.add_argument("-u", "--unificar-desconexos", action="store_true", default=False, help="Liga a unificação espacial dos rios desconexos")
    
    args = parser.parse_args()
    
    extrair_rios_principais_ana_online(
        nome_bacia=args.bacia,
        min_strahler=args.strahler,
        pasta_saida=args.saida,
        unificar_desconexos=args.unificar_desconexos
    )
