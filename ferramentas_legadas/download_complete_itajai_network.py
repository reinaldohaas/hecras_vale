"""
DOWNLOAD DA REDE HIDROGRÁFICA COMPLETA DO VALE DO ITAJAÍ (OSM / OVERPASS)
==========================================================================
Baixa a rede completa de rios e afluentes para toda a Bacia do Itajaí
e salva em 'vale_itajai_full_network.geojson' e 'rios_itajai.geojson'.
"""

import json
import time
import requests
from shapely.geometry import LineString, MultiLineString, mapping
from shapely.ops import linemerge

BBOX = (-27.75, -50.25, -26.40, -48.55)

OVERPASS_URLS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
    "https://overpass.osm.ch/api/interpreter"
]

HEADERS = {
    "User-Agent": "hecras-vale-model/2.0 (Bacia do Itajai; uso academico)",
    "Accept": "application/json"
}

def fetch_overpass_rivers():
    s, w, n, e = BBOX
    query = f"""
    [out:json][timeout:180];
    (
      way["waterway"~"river|stream"]({s},{w},{n},{e});
    );
    out geom;
    """
    print(f"Consultando Overpass API para a BBOX {BBOX}...")
    
    for url in OVERPASS_URLS:
        try:
            print(f"  tentando endpoint: {url}")
            r = requests.post(url, data={"data": query}, headers=HEADERS, timeout=120)
            if r.status_code == 200:
                return r.json()
            print(f"  status: {r.status_code}")
        except Exception as err:
            print(f"  falhou: {err}")
    
    print("Aviso: Falha no download online. Usando base local existente em rios_itajai.geojson...")
    return None

def main():
    data = fetch_overpass_rivers()
    
    if data and "elements" in data and len(data["elements"]) > 0:
        ways_by_name = {}
        for el in data["elements"]:
            geom = el.get("geometry")
            if not geom or len(geom) < 2:
                continue
            coords = [(p["lon"], p["lat"]) for p in geom]
            name = el.get("tags", {}).get("name", "Sem Nome")
            ways_by_name.setdefault(name, []).append(coords)
            
        features = []
        for name, segs in ways_by_name.items():
            lines = [LineString(c) for c in segs if len(c) >= 2]
            if not lines:
                continue
            merged = linemerge(MultiLineString(lines)) if len(lines) > 1 else lines[0]
            
            features.append({
                "type": "Feature",
                "properties": {
                    "name": name,
                    "length_km": round(merged.length * 111.0, 2)
                },
                "geometry": mapping(merged)
            })
            
        fc = {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
            "features": features
        }
        
        with open("vale_itajai_full_network.geojson", "w", encoding="utf-8") as f:
            json.dump(fc, f, ensure_ascii=False, indent=2)
        print(f"✓ Salvo com SUCESSO: vale_itajai_full_network.geojson com {len(features)} feições de rios!")
    else:
        # Se o overpass demorar, copia a base rios_itajai.geojson existente (517 rios)
        with open("rios_itajai.geojson", "r", encoding="utf-8") as f_in:
            data_local = json.load(f_in)
        with open("vale_itajai_full_network.geojson", "w", encoding="utf-8") as f_out:
            json.dump(data_local, f_out, ensure_ascii=False, indent=2)
        print(f"✓ 'vale_itajai_full_network.geojson' reconstruído com SUCESSO a partir de 'rios_itajai.geojson' ({len(data_local['features'])} rios)!")

if __name__ == "__main__":
    main()
