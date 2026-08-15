import json
import urllib.request

url_ana_otto = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_OTTOBACIA/FeatureServer/0/query?" + \
               "where=COBACIA+LIKE+%2783%25%27&outFields=COBACIA,NUAREABAC&outSR=4326&f=geojson&resultRecordCount=2000"

req = urllib.request.Request(url_ana_otto, headers={"User-Agent": "Mozilla/5.0"})
try:
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
        feats = data.get("features", [])
        print("  ✓ Ottobacias BHO 83 recebidas com sucesso! Feições:", len(feats))
        with open(r"C:\Users\haas\github\hecras_vale\bacia_itajai_poligono.geojson", "w", encoding="utf-8") as f:
            json.dump(data, f)
except Exception as e:
    print("  ⚠ Erro ANA Ottobacias 83:", e)
