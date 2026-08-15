import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

caminho_orig = r"C:\Users\haas\github\hecras_vale\rios_itajai_original.geojson"
with open(caminho_orig, "r", encoding="utf-8") as f:
    data = json.load(f)

feats = data.get("features", [])

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

for f in feats:
    p = f["properties"]
    cot = p.get("COTRECHO")
    rio = str(p.get("NORIOCOMP") or p.get("NOORIGINAL") or "")
    geom = f.get("geometry", {})
    gtype = geom.get("type")
    coords = geom.get("coordinates", [])

    is_mirim = "mirim" in rio.lower()
    is_acu = "itajai" in rio.lower() and not is_mirim

    color = 'blue' if is_acu else ('red' if is_mirim else 'gray')
    lw = 2.5 if (is_acu or is_mirim) else 0.8
    alpha = 0.9 if (is_acu or is_mirim) else 0.4

    if gtype == "LineString":
        xs = [pt[0] for pt in coords]
        ys = [pt[1] for pt in coords]
        ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha)
    elif gtype == "MultiLineString":
        for seg in coords:
            xs = [pt[0] for pt in seg]
            ys = [pt[1] for pt in seg]
            ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha)

ax.set_xlim(-48.75, -48.60)
ax.set_ylim(-26.96, -26.85)
ax.set_title("Trechos Originais em Itajaí (Vermelho = Rio Itajaí-Mirim Original)")

fig.savefig(r"C:\Users\haas\github\hecras_vale\figuras\teste_mirim_original.png", dpi=300)
print("Salvo em figuras/teste_mirim_original.png")
