"""
01 - Baixa TODOS os rios reais da Bacia do Itajai (incl. Itajai-Mirim) do
OpenStreetMap via Overpass API e salva como GeoJSON (WGS84).

Rode com:  python 01_baixar_rios.py

Saida: rios_itajai.geojson  (LineStrings, uma feature por rio, campo 'name')

Requisitos: requests, geopandas, shapely.
"""
import json
import time
import unicodedata

import requests
from shapely.geometry import LineString, MultiLineString, mapping
from shapely.ops import linemerge

# --- Area da bacia (S, W, N, E) em graus WGS84 --------------------------------
BBOX = (-27.75, -50.25, -26.40, -48.55)   # cobre Sul, Oeste, Norte, Acu, Mirim

# Endpoints Overpass (tenta em ordem se um falhar)
OVERPASS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
    "https://overpass.osm.ch/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
]

# MUITOS mirrors do Overpass devolvem 406/429 se o User-Agent for o padrao do
# requests. Um User-Agent descritivo resolve.
HEADERS = {
    "User-Agent": "hecras-vale-model/1.0 (Bacia do Itajai; uso academico)",
    "Accept": "application/json",
}

OUT = "rios_itajai.geojson"


def strip_accents(s):
    return "".join(c for c in unicodedata.normalize("NFD", s)
                   if unicodedata.category(c) != "Mn")


def overpass_query(query):
    query = query.strip()
    for url in OVERPASS:
        for tentativa in range(2):
            try:
                print(f"  consultando {url} (tentativa {tentativa+1}) ...")
                r = requests.post(url, data={"data": query},
                                  headers=HEADERS, timeout=180)
                if r.status_code == 200:
                    return r.json()
                print(f"    status {r.status_code}")
                if r.status_code in (429, 504):
                    time.sleep(12)      # sobrecarga: espera e tenta de novo
            except Exception as e:
                print(f"    falhou: {e}")
            time.sleep(3)
    raise RuntimeError("Todos os endpoints Overpass falharam. "
                       "Tente novamente em alguns minutos (servidor sobrecarregado).")


def fetch_rivers(bbox):
    s, w, n, e = bbox
    # todos os cursos com waterway=river na bacia (rios principais + afluentes)
    query = f"""
    [out:json][timeout:180];
    (
      way["waterway"="river"]({s},{w},{n},{e});
    );
    out geom;
    """
    print("Baixando rios (waterway=river) da bacia...")
    data = overpass_query(query)

    # agrupa segmentos por nome
    ways_by_name = {}
    unnamed = []
    for el in data.get("elements", []):
        geom = el.get("geometry")
        if not geom or len(geom) < 2:
            continue
        coords = [(p["lon"], p["lat"]) for p in geom]
        name = el.get("tags", {}).get("name")
        if name:
            ways_by_name.setdefault(name, []).append(coords)
        else:
            unnamed.append(coords)

    features = []
    for name, segs in sorted(ways_by_name.items()):
        lines = [LineString(c) for c in segs]
        merged = linemerge(MultiLineString(lines)) if len(lines) > 1 else lines[0]
        # comprimento aproximado (graus -> km, ~111 km/grau)
        length_km = merged.length * 111.0
        features.append({
            "type": "Feature",
            "properties": {"name": name, "name_ascii": strip_accents(name),
                           "length_km": round(length_km, 1)},
            "geometry": mapping(merged),
        })
        print(f"  {name:<32} ~{length_km:6.1f} km  ({len(segs)} segmentos)")

    print(f"\n{len(features)} rios nomeados; {len(unnamed)} segmentos sem nome (ignorados).")
    return features


def main():
    feats = fetch_rivers(BBOX)
    # ordena por comprimento (maiores primeiro)
    feats.sort(key=lambda f: -f["properties"]["length_km"])
    fc = {"type": "FeatureCollection",
          "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
          "features": feats}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(fc, f, ensure_ascii=False)
    print(f"\nSalvo: {OUT}")
    # destaca os rios da rede Itajai
    print("\nRios da rede Itajai encontrados:")
    for f in feats:
        na = f["properties"]["name_ascii"].lower()
        if "itaja" in na or "hercilio" in na or "benedito" in na or "mirim" in na:
            print(f"  - {f['properties']['name']}  ({f['properties']['length_km']} km)")


if __name__ == "__main__":
    main()
