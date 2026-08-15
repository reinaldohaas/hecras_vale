import json
import urllib.parse
import urllib.request

url_ana = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "geometry": "-48.80,-27.00,-48.55,-26.80",
              "geometryType": "esriGeometryEnvelope",
              "inSR": "4326",
              "outSR": "4326",
              "spatialRel": "esriSpatialRelIntersects",
              "outFields": "COTRECHO,COBACIA,NORIOCOMP,NOORIGINAL,NUSTRAHLER,NUTRJUS,NUAREAMONT",
              "returnGeometry": "true",
              "f": "json"
          })

req = urllib.request.Request(url_ana, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    data = json.loads(resp.read().decode("utf-8"))
    feats = data.get("features", [])

print(f"Total de trechos em Itajaí (inSR=4326): {len(feats)}")
for f in feats:
    p = f["attributes"]
    print(f"  - COTRECHO={p.get('COTRECHO')}, Nome='{p.get('NORIOCOMP') or p.get('NOORIGINAL')}', Strahler={p.get('NUSTRAHLER')}, Area={p.get('NUAREAMONT')}")
