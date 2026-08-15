import requests
import json
import os
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge

def fetch_river(river_name, bbox):
    """Baixa as vias de um rio específico via Overpass API."""
    south, west, north, east = bbox
    overpass_url = "http://overpass-api.de/api/interpreter"
    
    query = f"""
    [out:json][timeout:25];
    (
      way["waterway"="river"]["name"~"{river_name}",i]({south},{west},{north},{east});
    );
    out geometry;
    """
    print(f"Buscando no OpenStreetMap: {river_name}...")
    try:
        response = requests.post(overpass_url, data={'data': query}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            features = []
            for elem in data.get('elements', []):
                if 'geometry' in elem:
                    coords = [(pt['lon'], pt['lat']) for pt in elem['geometry']]
                    if len(coords) >= 2:
                        features.append({
                            'type': 'Feature',
                            'properties': {'name': river_name},
                            'geometry': {
                                'type': 'LineString',
                                'coordinates': coords
                            }
                        })
            return features
    except Exception as e:
        print(f"Erro ao buscar {river_name}: {e}")
    return []

def main():
    bbox = (-27.4, -49.8, -26.6, -48.5) # Bounding Box do Vale do Itajaí
    
    rivers_to_fetch = [
        ("Rio Itajaí do Sul", "Itajai_do_Sul"),
        ("Rio Itajaí do Oeste", "Itajai_do_Oeste"),
        ("Rio Itajaí do Norte", "Itajai_do_Norte"),
        ("Rio Itajaí-Açu", "Itajai_Acu"),
        ("Rio Benedito", "Rio_Benedito"),
        ("Rio Itajaí-Mirim", "Itajai_Mirim")
    ]
    
    all_features = []
    
    for osm_name, code_name in rivers_to_fetch:
        feats = fetch_river(osm_name, bbox)
        if feats:
            print(f"-> {len(feats)} trechos encontrados para {osm_name}")
            all_features.extend(feats)
        else:
            print(f"-> Nenhum trecho retornado diretamente para {osm_name}.")
            
    geojson_data = {
        "type": "FeatureCollection",
        "features": all_features
    }
    
    output_path = "vale_itajai_full_network.geojson"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(geojson_data, f, ensure_ascii=False, indent=2)
        
    print(f"\nRede hidrográfica salva em: {output_path}")

if __name__ == "__main__":
    main()
