import json
import urllib.parse
import urllib.request

def requisitar_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Antigravity/2.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode('utf-8'))

bacias = ["Itajaí", "Tubarão", "Araranguá"]

for nome in bacias:
    print(f"\n=======================================================")
    print(f" PROCESSANDO EXUTÓRIA E PRINCIPAIS AFLUENTES: {nome}")
    print(f"=======================================================")
    
    # 1. Localizar foz na ANA
    url_foz = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
              urllib.parse.urlencode({
                  "where": f"NORIOCOMP LIKE '%{nome}%' OR NOORIGINAL LIKE '%{nome}%'",
                  "outFields": "COTRECHO,NORIOCOMP,NOORIGINAL,NUTRJUS,NUAREAMONT,COBACIA,NUSTRAHLER",
                  "returnGeometry": "false",
                  "f": "json",
                  "resultRecordCount": "25"
              })
    data = requisitar_json(url_foz)
    feats = data.get("features", [])
    if not feats:
        print("  Nenhum resultado para foz.")
        continue
    feats.sort(key=lambda f: f["attributes"].get("NUAREAMONT") or 0, reverse=True)
    foz = feats[0]["attributes"]
    cobac_prefix = str(foz["COBACIA"])[:4]
    
    print(f" Exutória (Foz): '{foz.get('NORIOCOMP')}' (COTRECHO={foz.get('COTRECHO')}, Área={foz.get('NUAREAMONT')} km², Strahler={foz.get('NUSTRAHLER')})")
    print(f" Prefixo Otto ANA: '{cobac_prefix}'")

    # 2. Baixar trechos dos rios da bacia com Strahler >= 3 (Principais Afluentes)
    url_rios = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
               urllib.parse.urlencode({
                   "where": f"COBACIA LIKE '{cobac_prefix}%' AND NUSTRAHLER >= 3",
                   "outFields": "COTRECHO,NORIOCOMP,NUTRJUS,NUSTRAHLER,NUAREAMONT",
                   "returnGeometry": "false",
                   "f": "json",
                   "returnCountOnly": "true"
               })
    count_data = requisitar_json(url_rios)
    print(f" Total de Principais Afluentes (Strahler >= 3): {count_data.get('count')} trechos")
