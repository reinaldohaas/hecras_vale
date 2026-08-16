# -*- coding: utf-8 -*-
"""
Exporta o modelo HEC-RAS para KMZ (Google Earth).

Gera um KMZ com quatro camadas, organizadas em pastas:

  1. Lamina d'agua 3D  -- o perfil longitudinal como uma CORTINA vertical:
     a linha e desenhada na cota d'agua absoluta e extrudada ate o terreno,
     entao a altura da parede e a profundidade da cheia. E o perfil
     longitudinal visto em perspectiva, sobre o relevo real.
  2. Talvegue 3D       -- o leito escavado, na cota absoluta.
  3. Secoes transversais -- as cutlines, coloridas por extravasamento.
  4. Mancha de inundacao -- os poligonos por classe de profundidade.

Uso:  python exportar_kml.py [PROJETO]
      python exportar_kml.py Itajai_Rede_1983

Saida: <PROJETO>.kmz  (abre com duplo clique no Google Earth)
"""
import json
import os
import sys
import zipfile

import numpy as np
import h5py
from pyproj import Transformer

UTM_EPSG = 31982
MANCHA = os.path.join("app", "manchas_inundacao_hecras.geojson")

# ABGR (o KML inverte: alpha, blue, green, red)
COR = {
    "agua":     "a0f0a000",   # azul-claro translucido
    "talvegue": "ff503010",   # marrom escuro
    "dentro":   "ff50c878",   # verde
    "extravasa":"ff2020e0",   # vermelho
    "baixa":    "60f0d038",
    "media":    "70e08420",
    "severa":   "808b3a1e",
}


def esc(t):
    return (str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def ler_geometria(projeto):
    """Eixo de cada trecho e as cutlines, do .g01."""
    txt = open(f"{projeto}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    eixos, secoes = {}, []
    rio = rea = rs = None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, rea = p[0].strip(), p[1].strip()
            n = int(txt[i + 1].split("=")[1]); v = []; j = i + 2
            while len(v) < 2 * n and j < len(txt) and not txt[j].startswith("Type RM"):
                v += [float(x) for x in txt[j].split()]; j += 1
            eixos[(rio, rea)] = list(zip(v[0::2], v[1::2]))
            i = j
            continue
        if l.startswith("Type RM"):
            try: rs = float(l.split(",")[1])
            except ValueError: rs = None
        elif l.startswith("XS GIS Cut Line"):
            v = [float(x) for x in txt[i + 1].split()]
            secoes.append({"rio": rio, "reach": rea, "rs": rs, "cut": v[:4]})
            i += 1
        i += 1
    return eixos, secoes


def kml_linha(nome, coords, cor, largura=3, extrude=False, modo="absolute",
              desc=""):
    c = " ".join(f"{lo:.6f},{la:.6f},{al:.1f}" for lo, la, al in coords)
    return f"""<Placemark><name>{esc(nome)}</name>
<description><![CDATA[{desc}]]></description>
<Style><LineStyle><color>{cor}</color><width>{largura}</width></LineStyle>
<PolyStyle><color>{COR['agua']}</color></PolyStyle></Style>
<LineString><extrude>{1 if extrude else 0}</extrude><tessellate>1</tessellate>
<altitudeMode>{modo}</altitudeMode><coordinates>{c}</coordinates></LineString></Placemark>"""


def kml_poligono(coords, cor, nome="", desc=""):
    c = " ".join(f"{lo:.6f},{la:.6f},0" for lo, la in coords)
    return f"""<Placemark><name>{esc(nome)}</name>
<description><![CDATA[{desc}]]></description>
<Style><PolyStyle><color>{cor}</color></PolyStyle>
<LineStyle><width>0</width></LineStyle></Style>
<Polygon><altitudeMode>clampToGround</altitudeMode><outerBoundaryIs><LinearRing>
<coordinates>{c}</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>"""


def main():
    projeto = sys.argv[1] if len(sys.argv) > 1 else "Itajai_Rede_1983"
    eixos, secoes = ler_geometria(projeto)
    tr = Transformer.from_crs(UTM_EPSG, 4326, always_xy=True)

    with h5py.File(f"{projeto}.p01.hdf", "r") as f:
        g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
              "Unsteady Time Series/Cross Sections"]
        ws, q = g["Water Surface"][:], g["Flow"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = [x["River"].decode().strip() for x in at]
        rch = [x["Reach"].decode().strip() for x in at]
        rsh = np.array([float(x["RS"].decode()) for x in at])
    print(f"{projeto}: {ws.shape[0]} instantes, {ws.shape[1]} secoes")

    partes = []

    # --- 1 e 2. lamina d'agua e talvegue em 3D, por trecho -----------------
    agua, leito = [], []
    for (rio, rea), pts in eixos.items():
        idx = [j for j in range(len(riv)) if riv[j] == rio and rch[j] == rea]
        if len(idx) < 2:
            continue
        idx.sort(key=lambda j: -rsh[j])          # montante -> jusante
        rs_o = rsh[idx]
        z_ag = ws[:, idx].max(axis=0)            # cota de PICO
        z_le = ws[0, idx]                        # 1o instante ~ leito
        L = len(pts)
        # cada ponto do eixo recebe a cota interpolada pela estaca
        prog = np.linspace(rs_o[0], rs_o[-1], L)
        ca, cl = [], []
        for k, (x, y) in enumerate(pts):
            lo, la = tr.transform(x, y)
            ca.append((lo, la, float(np.interp(prog[k], rs_o[::-1], z_ag[::-1]))))
            cl.append((lo, la, float(np.interp(prog[k], rs_o[::-1], z_le[::-1]))))
        d = (f"<b>{rio} / {rea}</b><br>Cota de pico: "
             f"{z_ag.max():.2f} m<br>Vazao de pico: {q[:, idx].max():.0f} m&sup3;/s"
             f"<br>Extensao: {abs(rs_o[0]-rs_o[-1])/1000:.1f} km")
        agua.append(kml_linha(f"{rio} / {rea}", ca, COR["agua"], 2,
                              extrude=True, desc=d))
        leito.append(kml_linha(f"{rio} / {rea}", cl, COR["talvegue"], 2, desc=d))
    partes.append("<Folder><name>Lamina d'agua no pico (3D)</name>"
                  "<description><![CDATA[A altura da parede e a profundidade "
                  "da cheia sobre o terreno.]]></description>"
                  + "".join(agua) + "</Folder>")
    partes.append("<Folder><name>Talvegue (leito escavado)</name><visibility>0</visibility>"
                  + "".join(leito) + "</Folder>")

    # --- 3. secoes transversais -------------------------------------------
    sec = []
    for s in secoes:
        if s["rs"] is None:
            continue
        cand = [j for j in range(len(riv))
                if riv[j] == s["rio"] and rch[j] == s["reach"]]
        if not cand:
            continue
        j = min(cand, key=lambda x: abs(rsh[x] - s["rs"]))
        x1, y1, x2, y2 = s["cut"]
        lo1, la1 = tr.transform(x1, y1); lo2, la2 = tr.transform(x2, y2)
        hmax = float(ws[:, j].max()); hmin = float(ws[0, j])
        prof = hmax - hmin
        cor = COR["extravasa"] if prof > 3.0 else COR["dentro"]
        d = (f"<b>{s['rio']} RS {s['rs']/1000:.2f} km</b><br>"
             f"Cota de pico: {hmax:.2f} m<br>Subida: {prof:.2f} m<br>"
             f"Vazao de pico: {q[:, j].max():.0f} m&sup3;/s")
        sec.append(kml_linha(f"RS {s['rs']/1000:.2f} km",
                             [(lo1, la1, 0), (lo2, la2, 0)], cor, 2,
                             modo="clampToGround", desc=d))
    partes.append(f"<Folder><name>Secoes transversais ({len(sec)})</name>"
                  "<visibility>0</visibility>" + "".join(sec) + "</Folder>")

    # --- 4. mancha de inundacao -------------------------------------------
    if os.path.exists(MANCHA):
        d = json.load(open(MANCHA, encoding="utf-8"))
        hs = sorted(set(f["properties"]["time_hour"] for f in d["features"]))
        pico = hs[-1] if not hs else max(
            hs, key=lambda h: sum(f["properties"]["area_km2"]
                                  for f in d["features"]
                                  if f["properties"]["time_hour"] == h))
        cls = {"Baixa": COR["baixa"], "dia": COR["media"], "Severa": COR["severa"]}
        pol = []
        for f in d["features"]:
            p = f["properties"]
            if p["time_hour"] != pico:
                continue
            cor = next((v for k, v in cls.items() if k in p["class_name"]),
                       COR["media"])
            gs = ([f["geometry"]["coordinates"]]
                  if f["geometry"]["type"] == "Polygon"
                  else f["geometry"]["coordinates"])
            for g_ in gs:
                anel = g_[0] if isinstance(g_[0][0], (list, tuple)) else g_
                if len(anel) < 4:
                    continue
                pol.append(kml_poligono(anel, cor, p["class_name"],
                                        f"{p['class_name']}<br>"
                                        f"{p['area_km2']} km&sup2; na classe"))
        partes.append(f"<Folder><name>Mancha de inundacao (t={pico}h)</name>"
                      + "".join(pol) + "</Folder>")
        print(f"  mancha: {len(pol)} poligonos no instante de pico (t={pico}h)")

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"><Document>
<name>{esc(projeto)} - Bacia do Itajai</name>
<description><![CDATA[Modelo hidrodinamico 1D no HEC-RAS 7.0.1.<br>
Rede da ANA (BHO 2017) &middot; relevo Copernicus GLO-30 &middot;
calha escavada por geometria hidraulica.<br>
Gerado por exportar_kml.py]]></description>
{''.join(partes)}
</Document></kml>"""

    kmz = f"{projeto}.kmz"
    with zipfile.ZipFile(kmz, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("doc.kml", kml)
    with open(f"{projeto}.kml", "w", encoding="utf-8") as f:
        f.write(kml)
    print(f"\n[OK] {kmz}  ({os.path.getsize(kmz)/1e6:.1f} MB)")
    print(f"[OK] {projeto}.kml  ({len(kml)/1e6:.1f} MB)")
    print("\nAbra o .kmz com duplo clique. Incline a vista (Shift + arrastar)")
    print("para ver a lamina d'agua em perspectiva sobre o relevo.")


if __name__ == "__main__":
    main()
