import json
import urllib.parse
import urllib.request

def requisitar(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Antigravity/2.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode('utf-8'))

bacias = ["Tubarão", "Araranguá", "Itajaí"]

for nome in bacias:
    print(f"\n--- RASTREAMENTO DO RIO '{nome}' E SEUS TRIBUTÁRIOS ---")
    url_foz = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
              urllib.parse.urlencode({
                  "where": f"NORIOCOMP LIKE '%{nome}%' OR NOORIGINAL LIKE '%{nome}%'",
                  "outFields": "COTRECHO,NORIOCOMP,NUTRJUS,NUAREAMONT,COBACIA",
                  "returnGeometry": "false",
                  "f": "json",
                  "resultRecordCount": "20"
              })
    data_foz = requisitar(url_foz)
    feats = data_foz.get("features", [])
    if not feats:
        print("  Nenhum segmento encontrado.")
        continue
    feats.sort(key=lambda x: x["attributes"].get("NUAREAMONT") or 0, reverse=True)
    foz = feats[0]["attributes"]
    print(f"  Foz encontrada: COTRECHO={foz['COTRECHO']}, Nome='{foz['NORIOCOMP']}', Área={foz['NUAREAMONT']} km², COBACIA={foz['COBACIA']}")
