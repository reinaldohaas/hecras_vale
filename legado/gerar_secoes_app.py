# -*- coding: utf-8 -*-
"""
Exporta as secoes transversais do modelo HEC-RAS para o app.

Substitui app/hecras_secoes_transversais.geojson, que era sintetico e tinha
tres inconsistencias detectaveis no proprio grafico da interface:

  1. z_bed ficava FORA do perfil desenhado em 95% das secoes (mediana 1,06 m
     abaixo do terreno, pior caso 5,48 m). A interface desenhava agua
     preenchendo um canal que nao existia na geometria mostrada ao lado.
  2. top_width_m nao era largura de superficie: 1.191 das 1.262 secoes
     declaravam exatamente 700 m, que e a largura da propria cutline.
  3. bank_h tinha 4 valores distintos em 1.262 secoes -- z_bank era
     z_bed + constante por rio, nao medido do terreno.

Aqui tudo vem da mesma fonte: a geometria do .g01 e a cota do solver no
.p01.hdf. Em particular:
  - z_bed  = minimo REAL do perfil (com a calha ja escavada);
  - z_bank = cota das margens topograficas gravadas em Bank Sta;
  - top_width = largura MOLHADA calculada na cota d'agua do instante;
  - water_stages = serie do solver, nao curva arbitrada.

Uso:  python gerar_secoes_app.py [PROJETO] [N_PASSOS]
      python gerar_secoes_app.py Itajai_Rede_1983 17
"""
import json
import os
import sys

import numpy as np
import h5py
from pyproj import Transformer

UTM_EPSG = 31982
SAIDA = os.path.join("app", "hecras_secoes_transversais.geojson")

# Cidades de referencia, para rotular a secao pela mais proxima.
CIDADES = [
    ("Rio do Sul", -27.2150, -49.6430), ("Ibirama", -27.0550, -49.5200),
    ("Apiúna", -27.0380, -49.3890), ("Indaial", -26.8980, -49.2320),
    ("Timbó", -26.8230, -49.2720), ("Blumenau", -26.9180, -49.0660),
    ("Gaspar", -26.9320, -48.9590), ("Ilhota", -26.9020, -48.8280),
    ("Itajaí", -26.9060, -48.6650), ("Navegantes", -26.8990, -48.6540),
    ("Brusque", -27.0980, -48.9120), ("Botuverá", -27.2010, -49.0700),
    ("Vidal Ramos", -27.3890, -49.3600),
]


def ler_geometria(projeto):
    """Secoes do .g01: estacas, cotas, margens e cutline georreferenciada."""
    txt = open(f"{projeto}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    secoes, rio, rea, rs, cut = [], None, None, None, None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, rea = p[0].strip(), p[1].strip()
        elif l.startswith("Type RM"):
            try:
                rs = float(l.split(",")[1])
            except ValueError:
                rs = None
        elif l.startswith("XS GIS Cut Line"):
            v = [float(x) for x in txt[i + 1].split()]
            cut = v[:4]
            i += 1
        elif l.startswith("#Sta/Elev="):
            n = int(l.split("=")[1])
            v = []
            i += 1
            while i < len(txt) and len(v) < 2 * n:
                s = txt[i]
                v += [float(s[j:j + 8]) for j in range(0, len(s.rstrip()), 8)
                      if s[j:j + 8].strip()]
                i += 1
            sta, z = np.array(v[0::2]), np.array(v[1::2])
            lb = rb = None
            for k in range(i, min(i + 6, len(txt))):
                if txt[k].startswith("Bank Sta="):
                    lb, rb = [float(x) for x in txt[k].split("=")[1].split(",")]
                    break
            secoes.append({"rio": rio, "reach": rea, "rs": rs, "sta": sta,
                           "z": z, "lb": lb, "rb": rb, "cut": cut})
            continue
        i += 1
    return secoes


def largura_molhada(sta, z, nivel):
    """Largura da superficie livre na cota dada -- a largura de verdade."""
    molhado = z <= nivel
    if not molhado.any():
        return 0.0
    dsta = np.gradient(sta)
    return float(np.sum(dsta[molhado]))


def main():
    projeto = sys.argv[1] if len(sys.argv) > 1 else "Itajai_Rede_1983"
    n_passos = int(sys.argv[2]) if len(sys.argv) > 2 else 17
    evento = projeto.split("_")[-1] if projeto.count("_") >= 2 else "sim"

    secoes = ler_geometria(projeto)
    print(f"{len(secoes)} secoes lidas de {projeto}.g01")

    with h5py.File(f"{projeto}.p01.hdf", "r") as f:
        g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
              "Unsteady Time Series/Cross Sections"]
        ws = g["Water Surface"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = [x["River"].decode().strip() for x in at]
        rch = [x["Reach"].decode().strip() for x in at]
        rs_h = np.array([float(x["RS"].decode()) for x in at])
    print(f"{ws.shape[0]} instantes no HDF")
    idx = np.linspace(0, ws.shape[0] - 1, n_passos).round().astype(int)

    tr = Transformer.from_crs(UTM_EPSG, 4326, always_xy=True)
    feats = []
    for k, s in enumerate(secoes):
        if s["rs"] is None or s["lb"] is None or s["cut"] is None:
            continue
        # casa a secao com a coluna correspondente no HDF
        cand = [j for j in range(len(riv))
                if riv[j] == s["rio"] and rch[j] == s["reach"]]
        if not cand:
            continue
        j = min(cand, key=lambda x: abs(rs_h[x] - s["rs"]))
        serie = [round(float(ws[t, j]), 2) for t in idx]

        sta, z = s["sta"], s["z"]
        z_bed = float(z.min())                     # minimo REAL do perfil
        i_lb = int(np.argmin(np.abs(sta - s["lb"])))
        i_rb = int(np.argmin(np.abs(sta - s["rb"])))
        z_bank = float(min(z[i_lb], z[i_rb]))      # margem topografica
        z_pico = max(serie)

        x1, y1, x2, y2 = s["cut"]
        lo1, la1 = tr.transform(x1, y1)
        lo2, la2 = tr.transform(x2, y2)
        lo_c, la_c = tr.transform((x1 + x2) / 2, (y1 + y2) / 2)
        cidade = min(CIDADES, key=lambda c: (c[1] - la_c) ** 2 + (c[2] - lo_c) ** 2)[0]

        # Decima o perfil para a interface, mas FORCANDO a inclusao do
        # talvegue e das duas margens. Sem isso a calha escavada (estreita)
        # cai fora da amostra e o grafico volta a mostrar agua num canal que
        # nao aparece -- que era exatamente o defeito do arquivo antigo.
        i_min = int(np.argmin(z))
        sel = set(np.linspace(0, len(sta) - 1, 34).round().astype(int).tolist())
        sel |= {i_min, max(i_min - 1, 0), min(i_min + 1, len(sta) - 1),
                i_lb, i_rb}
        sel = sorted(sel)
        feats.append({
            "type": "Feature",
            "properties": {
                "xs_id": f"{s['rio'].lower()}_XS_{k:04d}",
                "river": s["rio"], "reach": s["reach"],
                "rs_km": round(s["rs"] / 1000, 2),
                "city": cidade,
                "z_bed": round(z_bed, 2),
                "z_bank": round(z_bank, 2),
                "bank_h": round(z_bank - z_bed, 2),
                "top_width_m": round(largura_molhada(sta, z, z_pico), 1),
                "z_water_pico": round(z_pico, 2),
                "is_overtop": bool(z_pico > z_bank),
                # Aliases dos nomes que index.html e mapa_perfis_hecras.html
                # ja leem, para as duas paginas seguirem funcionando sem
                # edicao. O sufixo _t24 vinha de "t = 24 h" do arquivo antigo;
                # aqui o valor e a cota de PICO, que e o que interessa.
                "z_water_t24": round(z_pico, 2),
                "is_overtop_t24": bool(z_pico > z_bank),
                "river_key": s["rio"].split("_")[-1].lower(),
                "dist_km": None,          # preenchido abaixo (km da cabeceira)
                "offsets_m": [round(float(sta[i] - sta[len(sta)//2]), 1) for i in sel],
                "sta_elevs": [round(float(z[i]), 2) for i in sel],
                f"water_stages_{evento}": serie,
                "fonte": f"HEC-RAS {projeto}",
            },
            "geometry": {"type": "LineString",
                         "coordinates": [[round(lo1, 6), round(la1, 6)],
                                         [round(lo2, 6), round(la2, 6)]]},
        })

    # dist_km = distancia desde a CABECEIRA (rs_km e medido da foz). O app
    # antigo mostrava os dois; manter evita confusao com o painel lateral.
    for rio in set(x["properties"]["river"] for x in feats):
        rs_max = max(x["properties"]["rs_km"] for x in feats
                     if x["properties"]["river"] == rio)
        for x in feats:
            if x["properties"]["river"] == rio:
                x["properties"]["dist_km"] = round(rs_max - x["properties"]["rs_km"], 2)

    with open(SAIDA, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)

    # --- conferencia das tres inconsistencias do arquivo antigo ---
    P = [x["properties"] for x in feats]
    fora = sum(1 for p in P
               if abs(p["z_bed"] - min(p["sta_elevs"])) > 0.5)
    larg = len(set(p["top_width_m"] for p in P))
    bh = len(set(p["bank_h"] for p in P))
    over = sum(1 for p in P if p["is_overtop"])
    print(f"\n[OK] {SAIDA}  ({len(feats)} secoes)")
    print(f"  z_bed fora do perfil : {fora}/{len(P)}   (antes: 1.199/1.262)")
    print(f"  larguras distintas   : {larg}          (antes: 2)")
    print(f"  bank_h distintos     : {bh}          (antes: 4)")
    print(f"  secoes que extravasam: {over}/{len(P)} ({100*over/len(P):.0f}%)")


if __name__ == "__main__":
    main()
