import json
import urllib.parse
import urllib.request

def requisitar_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Antigravity/2.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode('utf-8'))

nome = "Tijucas"
url_foz = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "where": f"NORIOCOMP LIKE '%{nome}%' OR NOORIGINAL LIKE '%{nome}%'",
              "outFields": "COTRECHO,NORIOCOMP,NOORIGINAL,NUTRJUS,NUAREAMONT,COBACIA,NUSTRAHLER",
              "returnGeometry": "false",
              "f": "json"
          })
data = requisitar_json(url_foz)
feats = data.get("features", [])
feats.sort(key=lambda f: f["attributes"].get("NUAREAMONT") or 0, reverse=True)
foz = feats[0]["attributes"]
print("Foz Encontrada:", foz)
cobac = str(foz["COBACIA"])
print("COBACIA da Foz:", cobac)
print("Prefixo Otto (4 dígitos):", cobac[:4])
print("Prefixo Otto (5 dígitos):", cobac[:5])
print("Prefixo Otto (6 dígitos):", cobac[:6])

# Baixar todos os trechos com COBACIA LIKE cobac[:4]
url_all = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "where": f"COBACIA LIKE '{cobac[:4]}%' AND NUSTRAHLER >= 3",
              "outFields": "COTRECHO,NORIOCOMP,NOORIGINAL,NUTRJUS,NUAREAMONT,COBACIA,NUSTRAHLER",
              "returnGeometry": "false",
              "f": "json",
              "resultRecordCount": "1000"
          })
all_data = requisitar_json(url_all)
all_feats = all_data.get("features", [])
print(f"\nTotal trechos no prefixo {cobac[:4]} com Strahler >= 3: {len(all_feats)}")

# Tracing a montante a partir da foz
by_cotrecho = {f["attributes"]["COTRECHO"]: f["attributes"] for f in all_feats if f["attributes"].get("COTRECHO")}
children = {}
for f in all_feats:
    attr = f["attributes"]
    parent = attr.get("NUTRJUS")
    cid = attr.get("COTRECHO")
    if parent not in children: children[parent] = []
    children[parent].append(cid)

cotrecho_foz = foz["COTRECHO"]
connected = set()
queue = [cotrecho_foz]
while queue:
    curr = queue.pop()
    if curr in connected: continue
    connected.add(curr)
    if curr in children:
        for child in children[curr]:
            queue.append(child)

print(f"Trechos CONECTADOS À FOZ DO TIJUCAS: {len(connected)}")

# Mostrar alguns trechos desconectados que foram ignorados
disconnected = [f["attributes"] for f in all_feats if f["attributes"]["COTRECHO"] not in connected]
print(f"Trechos DESCONECTADOS (de outros rios como Camboriú, Biguaçu, Cubatão): {len(disconnected)}")
rios_desconectados = set(d.get("NORIOCOMP") or d.get("NOORIGINAL") for d in disconnected if d.get("NORIOCOMP") or d.get("NOORIGINAL"))
print("Rios desconectados ignorados na filtragem estrita da foz:", sorted(list(rios_desconectados)))
