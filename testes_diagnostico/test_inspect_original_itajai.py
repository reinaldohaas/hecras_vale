import json

caminho_orig = r"C:\Users\haas\github\hecras_vale\rios_itajai_original.geojson"
with open(caminho_orig, "r", encoding="utf-8") as f:
    data = json.load(f)

feats = data.get("features", [])
print(f"Total de feições em rios_itajai_original.geojson: {len(feats)}")

# Filtrar feições na caixa delimitadora de Itajaí
itajai_area_feats = []
for f in feats:
    geom = f.get("geometry", {})
    coords = geom.get("coordinates", [])
    gtype = geom.get("type")
    
    # Extrair todos os pontos para verificar BBox
    pts = []
    if gtype == "LineString": pts = coords
    elif gtype == "MultiLineString":
        for seg in coords: pts.extend(seg)
        
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    
    if lons and lats:
        min_x, max_x = min(lons), max(lons)
        min_y, max_y = min(lats), max(lats)
        
        # BBox de Itajaí (-48.75 a -48.60, -26.96 a -26.85)
        if max_x >= -48.75 and min_x <= -48.60 and max_y >= -26.96 and min_y <= -26.85:
            p = f["properties"]
            itajai_area_feats.append((p.get("COTRECHO"), p.get("NORIOCOMP") or p.get("NOORIGINAL"), p.get("NUSTRAHLER"), len(pts)))

print(f"\nFeições originais na região de Itajaí ({len(itajai_area_feats)} trechos):")
for cot, nm, st, npts in itajai_area_feats:
    print(f"  - COTRECHO={cot}, Nome='{nm}', Strahler={st}, n_pts={npts}")
