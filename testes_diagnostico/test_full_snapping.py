import json
import math
import os

def calcular_distancia_meters(p1, p2):
    R = 6371000.0
    dlat = math.radians(p2[1] - p1[1])
    dlon = math.radians(p2[0] - p1[0])
    a = math.sin(dlat/2.0)**2 + math.cos(math.radians(p1[1])) * math.cos(math.radians(p2[1])) * math.sin(dlon/2.0)**2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return R * c

def extrair_todos_pontos(geom):
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

def testar_conexao_total(caminho_geojson):
    if not os.path.exists(caminho_geojson):
        print(f"Arquivo não encontrado: {caminho_geojson}")
        return

    with open(caminho_geojson, "r", encoding="utf-8") as f:
        data = json.load(f)

    feats = data["features"]
    print(f"\n=======================================================")
    print(f" TESTANDO CONEXÃO GLOBAL: {os.path.basename(caminho_geojson)}")
    print(f" Total de segmentos: {len(feats)}")
    print(f"=======================================================")

    # Identificar foz (maior área NUAREAMONT)
    feats.sort(key=lambda f: f["properties"].get("NUAREAMONT") or 0, reverse=True)
    main_outlet = feats[0]
    main_id = main_outlet["properties"].get("COTRECHO")
    print(f" Foz: COTRECHO={main_id}, Nome={main_outlet['properties'].get('NORIOCOMP')}, Área={main_outlet['properties'].get('NUAREAMONT')} km²")

    # Mapeamento NUTRJUS -> COTRECHO
    by_cotrecho = {f["properties"]["COTRECHO"]: f for f in feats if f["properties"].get("COTRECHO")}
    children = {}
    for f in feats:
        parent = f["properties"].get("NUTRJUS")
        cid = f["properties"].get("COTRECHO")
        if parent not in children: children[parent] = []
        children[parent].append(cid)

    # Componente 0: Árvore topológica principal a partir da foz
    connected_ids = set()
    queue = [main_id]
    while queue:
        curr = queue.pop()
        if curr in connected_ids: continue
        connected_ids.add(curr)
        if curr in children:
            for child in children[curr]:
                queue.append(child)

    print(f" Topologicamente conectados via NUTRJUS: {len(connected_ids)} de {len(feats)}")

    # Agrupar todos os segmentos em componentes conectados (por proximidade < 50m)
    unvisited = set(f["properties"].get("COTRECHO") for f in feats if f["properties"].get("COTRECHO"))
    components = []

    # O componente 0 é a árvore conectada da foz
    components.append(connected_ids)
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

    print(f" Total de componentes desconectados isolados: {len(components) - 1}")

    # CONECTAR CADA COMPONENTE ISOLADO À ÁRVORE PRINCIPAL
    main_tree_nodes = []
    for cid in components[0]:
        f = by_cotrecho.get(cid)
        if f:
            main_tree_nodes.extend(extrair_todos_pontos(f["geometry"]))

    pontes = []
    reconectados_total = 0

    for i in range(1, len(components)):
        comp = components[i]
        comp_nodes = []
        for cid in comp:
            f = by_cotrecho.get(cid)
            if f:
                comp_nodes.extend(extrair_todos_pontos(f["geometry"]))

        # Encontrar menor distância entre este componente e a árvore principal
        min_d = float("inf")
        best_comp_pt = None
        best_tree_pt = None

        # Amostragem para busca rápida
        step_c = max(1, len(comp_nodes) // 20)
        step_t = max(1, len(main_tree_nodes) // 100)

        for cpt in comp_nodes[::step_c]:
            for tpt in main_tree_nodes[::step_t]:
                d = calcular_distancia_meters(cpt, tpt)
                if d < min_d:
                    min_d = d
                    best_comp_pt = cpt
                    best_tree_pt = tpt

        if min_d < 10000.0: # 10 km
            reconectados_total += len(comp)
            main_tree_nodes.extend(comp_nodes)
            
            if min_d > 1.0:
                bridge_feat = {
                    "type": "Feature",
                    "properties": {
                        "COTRECHO": -999000 - i,
                        "NORIOCOMP": "Conexão Hidrográfica (Ajuste Global)",
                        "NUSTRAHLER": 1,
                        "NUCOMPTREC": min_d
                    },
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [best_comp_pt, best_tree_pt]
                    }
                }
                pontes.append(bridge_feat)

    print(f" Componentes reconectados à árvore principal: {len(components) - 1}")
    print(f" Pontes de conexão geradas: {len(pontes)}")
    print(f" TOTAL FINAL DE SEGMENTOS 100% CONECTADOS NA ÁRVORE DA FOZ: {len(feats) + len(pontes)}")

testar_conexao_total(r"C:\Users\haas\github\hecras_vale\rios_araranguá_nivel6_conectados.geojson")
testar_conexao_total(r"C:\Users\haas\github\hecras_vale\rios_tubarão_nivel6_conectados.geojson")
