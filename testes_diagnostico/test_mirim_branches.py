import json

caminho_original = r"C:\Users\haas\github\hecras_vale\rios_itajai_original.geojson"
with open(caminho_original, "r", encoding="utf-8") as f:
    data = json.load(f)

feats = data.get("features", [])
print(f"Total no cache original: {len(feats)}")

mirim_all = []
for f in feats:
    p = f["properties"]
    rio = str(p.get("NORIOCOMP") or p.get("NOORIGINAL") or "").lower()
    cotrecho = p.get("COTRECHO")
    # Buscar todos os trechos perto da foz do Itajaí-Mirim
    if "mirim" in rio or "canal" in rio or cotrecho in (393944, 177938, 344713, 24977, 21231, 177939, 344720, 177942, 344726):
        mirim_all.append((cotrecho, rio, p.get("NUSTRAHLER"), p.get("NUTRJUS")))

print(f"Encontrados {len(mirim_all)} trechos do sistema Itajaí-Mirim:")
for cot, nm, st, trj in mirim_all:
    print(f"  - COTRECHO={cot}, Nome='{nm}', Strahler={st}, TRJUS={trj}")
