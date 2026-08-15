import json
import urllib.parse
import urllib.request

bacias_teste = ["Itajaí", "Tubarão", "Araranguá", "Tietê", "Doce", "Paraíba do Sul", "São Francisco"]

for nome in bacias_teste:
    url = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_TRECHODRENAGEM/FeatureServer/0/query?" + \
          urllib.parse.urlencode({
              "where": f"NORIOCOMP LIKE '%{nome}%' OR NOORIGINAL LIKE '%{nome}%'",
              "outFields": "COBACIA,NORIOCOMP,NUAREAMONT",
              "returnGeometry": "false",
              "f": "json",
              "resultRecordCount": "5"
          })
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Antigravity/2.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            feats = data.get("features", [])
            if feats:
                feats.sort(key=lambda x: x["attributes"].get("NUAREAMONT") or 0, reverse=True)
                top = feats[0]["attributes"]
                cobac = str(top.get("COBACIA"))
                print(f"Bacia: {nome:15s} -> Rio: {top.get('NORIOCOMP'):30s} -> COBACIA: {cobac:10s} (Prefixo Otto: {cobac[:4]}) -> Área: {top.get('NUAREAMONT')} km²")
            else:
                print(f"Bacia: {nome:15s} -> Nenhum trecho encontrado com LIKE '%{nome}%'")
    except Exception as e:
        print(f"Bacia: {nome:15s} -> Erro: {e}")
