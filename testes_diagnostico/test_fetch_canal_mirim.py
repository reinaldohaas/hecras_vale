import json
import urllib.parse
import urllib.request

url_ana = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "where": "NORIOCOMP LIKE '%Itajaí-mirim%' OR NORIOCOMP LIKE '%Itajai-mirim%' OR NORIOCOMP LIKE '%Canal Itajaí%' OR NORIOCOMP LIKE '%Canal Itajai%'",
              "outFields": "COTRECHO,COBACIA,NORIOCOMP,NOORIGINAL,NUSTRAHLER,NUTRJUS,NUAREAMONT",
              "returnGeometry": "true",
              "f": "json"
          })

req = urllib.request.Request(url_ana, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    data = json.loads(resp.read().decode("utf-8"))
    feats = data.get("features", [])

print(f"Total de trechos do Rio/Canal Itajaí-Mirim: {len(feats)}")
for f in feats:
    p = f["attributes"]
    print(f"  - COTRECHO={p.get('COTRECHO')}, Nome='{p.get('NORIOCOMP') or p.get('NOORIGINAL')}', Strahler={p.get('NUSTRAHLER')}, TRJUS={p.get('NUTRJUS')}")
