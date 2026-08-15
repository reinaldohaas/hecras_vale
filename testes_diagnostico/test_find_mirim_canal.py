import json
import urllib.parse
import urllib.request
import os

print("Buscando trechos da ANA BHO no prefixo 7796...")
url_ana = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "where": "COBACIA LIKE '7796%'",
              "outFields": "COTRECHO,COBACIA,NORIOCOMP,NOORIGINAL,NUSTRAHLER,NUTRJUS,NUAREAMONT",
              "returnGeometry": "true",
              "f": "json",
              "resultRecordCount": "1000"
          })

req = urllib.request.Request(url_ana, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=60) as resp:
    data = json.loads(resp.read().decode("utf-8"))
    feats = data.get("features", [])

print(f"Total de trechos recebidos: {len(feats)}")

mirim_feats = []
for f in feats:
    attrs = f.get("attributes", {})
    rio = str(attrs.get("NORIOCOMP") or attrs.get("NOORIGINAL") or "").lower()
    if "mirim" in rio or "canal" in rio or "extravasor" in rio or attrs.get("NUSTRAHLER") == 1:
        mirim_feats.append(attrs)

print(f"Encontrados {len(mirim_feats)} trechos relacionados ao Itajaí-Mirim / canais / ordem Strahler < 3.")
for m in mirim_feats[:20]:
    print(f"  - COTRECHO={m.get('COTRECHO')}, Nome='{m.get('NORIOCOMP') or m.get('NOORIGINAL')}', Strahler={m.get('NUSTRAHLER')}, TRJUS={m.get('NUTRJUS')}")
