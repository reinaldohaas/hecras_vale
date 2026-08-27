# -*- coding: utf-8 -*-
"""Auditoria da geometria HEC-RAS: espacial, perfil, e comparacao com o MDT.

    python scripts/qc_geometria.py modelo/so_mirim.prj

NAO ALTERA NADA. Le o .prj, resolve o `Geom File=`, le as cross sections e a
polilinha do eixo do proprio .gNN, amostra o MDT do SIG-SC a 1 m ao longo de
cada cutline, e grava:

    <projeto>_qc.csv       uma linha por River Station
    <projeto>_qc.geojson   as cutlines, com todos os campos e o status
    <projeto>_qc.html      mapa clicavel: perfil HEC-RAS x MDT por secao

PREMISSAS, EXPLICITAS

  O PERFIL DO HEC-RAS E O DADO. O MDT serve para COMPARAR, nunca para
  substituir. A batimetria nao esta no MDT -- ele ve a superficie da agua --
  entao canal abaixo do MDT e o esperado, e nao entra em nenhum criterio de
  reprovacao.

  A TANGENTE DO RIO E LOCAL E ADAPTATIVA. Nada de janela fixa de +-250 m: a
  janela sai da largura do canal daquela secao (2x, entre 20 e 150 m), porque
  num meandro fechado uma janela longa mede a CORDA e nao a tangente.

  NENHUM LIMIAR DEFINE FISICA. Largura, fracao do talvegue e angulo entram
  como INDICADORES; o status diz o que foi medido e por que, e a decisao de
  corrigir e de quem le.

CRITERIOS (todos com o numero ao lado, no relatorio)

  station_length_error  |L_cutline - (station[-1]-station[0])|
                        > 0,50 m  CRITICAL     > 0,05 m  WARNING
  river_intersections   cruzamentos com o eixo do PROPRIO reach
                        0  CRITICAL   >=2  WARNING (com o motivo)
  angle                 angulo entre a cutline e a tangente local do rio
                        desvio de 90 graus:  > 50  CRITICAL   > 30  WARNING
  neighbor_overlap      cutline cruzando a da vizinha imediata
                        qualquer  WARNING
  dog_leg               maior deflexao interna da polilinha da cutline
                        > 30 graus  WARNING
  spikes                ponto > 3 m fora da mediana movel, FORA do canal
                        qualquer  WARNING
  dem_difference        mediana (HEC-RAS - MDT) FORA do canal
                        |.| > 3,0 m  CRITICAL   > 1,0 m  WARNING
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import (_bloco, _mediana_movel, geometria_do_projeto,  # noqa: E402
                       ler_secoes)
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio  # noqa: E402

TOL_LEN_W, TOL_LEN_C = 0.05, 0.50
ANG_W, ANG_C = 30.0, 50.0
DEM_W, DEM_C = 1.0, 3.0
SPIKE = 3.0
DOGLEG = 30.0


def ler_eixos(g01):
    """Polilinhas dos reaches, do proprio .gNN ('Reach XY=')."""
    from shapely.geometry import LineString
    linhas = open(g01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    eixos, chave = {}, None
    for i, l in enumerate(linhas):
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            chave = (p[0].strip(), p[1].strip() if len(p) > 1 else "")
        elif l.startswith("Reach XY=") and chave:
            v = _bloco(linhas, i, 16)
            eixos[chave] = LineString(np.array(v).reshape(-1, 2))
    return eixos


def _az(v):
    return np.degrees(np.arctan2(v[1], v[0]))


def _ang_entre(a, b):
    d = abs((a - b + 180.0) % 360.0 - 180.0)
    return min(d, 180.0 - d)


def tangente_local(eixo, s, janela):
    """Tangente do eixo em s, com janela ADAPTATIVA (nao fixa)."""
    a = np.array(eixo.interpolate(max(0.0, s - janela)).coords[0])
    b = np.array(eixo.interpolate(min(eixo.length, s + janela)).coords[0])
    v = b - a
    if np.hypot(*v) < 1e-9:
        v = np.array([1.0, 0.0])
    return v


def auditar(d, eixo, ant, prox, mdt):
    from shapely.geometry import LineString, Point
    st, z = d["sta"], d["z"]
    lb, rb = d.get("lb"), d.get("rb")
    cut = LineString(d["cut"])
    L_cut = float(cut.length)
    L_sta = float(st[-1] - st[0])
    larg_canal = (rb - lb) if (lb is not None and rb is not None) else np.nan

    r = {"River": d["rio"].strip(), "Reach": d["reach"].strip(),
         "RiverStation": d["rs"], "n_pontos": len(st),
         "length": L_cut, "station_length": L_sta,
         "station_length_error": abs(L_cut - L_sta),
         "channel_width": larg_canal,
         "bank_left": lb, "bank_right": rb}

    # ---------------- geometria espacial
    A = np.array(d["cut"][0], float); B = np.array(d["cut"][-1], float)
    az_cut = _az(B - A)
    r["angle_cutline"] = az_cut

    inter = cut.intersection(eixo)
    if inter.is_empty:
        ncr = 0; pontos = []
    elif inter.geom_type == "Point":
        ncr = 1; pontos = [inter]
    elif inter.geom_type == "MultiPoint":
        pontos = list(inter.geoms); ncr = len(pontos)
    else:
        pontos = []; ncr = 1
    r["river_intersections"] = ncr

    if pontos:
        s_eixo = float(np.median([eixo.project(p) for p in pontos]))
        jan = float(np.clip(2.0 * (larg_canal if np.isfinite(larg_canal) else 50.0),
                            20.0, 150.0))
        r["tangent_window"] = jan
        tg = tangente_local(eixo, s_eixo, jan)
        r["angle_river"] = _az(tg)
        r["angle"] = _ang_entre(az_cut, r["angle_river"])
    else:
        r["tangent_window"] = np.nan
        r["angle_river"] = np.nan
        r["angle"] = np.nan
    r["angle_desvio_90"] = abs(90.0 - r["angle"]) if np.isfinite(r["angle"]) else np.nan

    ov = 0
    for o in (ant, prox):
        if o is not None and "cut" in o:
            if cut.intersects(LineString(o["cut"])):
                ov += 1
    r["neighbor_overlap"] = ov
    r["dist_ant"] = (abs(d["rs"] - ant["rs"]) if ant is not None else np.nan)
    r["dist_prox"] = (abs(d["rs"] - prox["rs"]) if prox is not None else np.nan)

    P = np.asarray(d["cut"], float)
    if len(P) > 2:
        v = np.diff(P, axis=0)
        a = np.array([_az(x) for x in v])
        r["dog_leg"] = float(max(_ang_entre(a[i], a[i + 1])
                                 for i in range(len(a) - 1)))
    else:
        r["dog_leg"] = 0.0

    # ---------------- perfil HEC-RAS
    k = int(np.argmin(z))
    r["talweg_station"] = float(st[k])
    r["talweg_elev"] = float(z[k])
    r["talweg_fraction"] = float((st[k] - st[0]) / max(L_sta, 1e-9))
    r["bank_elev_left"] = float(np.interp(lb, st, z)) if lb is not None else np.nan
    r["bank_elev_right"] = float(np.interp(rb, st, z)) if rb is not None else np.nan
    r["pontos_duplicados"] = int(np.sum(np.diff(st) == 0))
    r["estacas_fora_de_ordem"] = int(np.sum(np.diff(st) < 0))
    ds = np.diff(st); dz = np.diff(z)
    r["salto_max"] = float(np.max(np.abs(dz))) if dz.size else 0.0
    r["incl_max"] = float(np.max(np.abs(dz) / np.maximum(ds, 1e-6))) if dz.size else 0.0

    # ---------------- MDT
    f = np.clip((st - st[0]) / max(L_sta, 1e-9), 0.0, 1.0)
    Q = [cut.interpolate(float(x), normalized=True) for x in f]
    zm = mdt.cota([p.x for p in Q], [p.y for p in Q])
    dif = z - zm
    canal = (((st >= lb) & (st <= rb)) if (lb is not None and rb is not None
                                           and rb > lb) else np.zeros(len(st), bool))
    fora = ~canal
    ok = np.isfinite(dif)
    r["mdt_cobertura"] = float(np.mean(np.isfinite(zm)))
    r["dem_difference"] = (float(np.median(dif[fora & ok]))
                           if (fora & ok).any() else np.nan)
    r["dem_difference_canal"] = (float(np.median(dif[canal & ok]))
                                 if (canal & ok).any() else np.nan)
    r["dem_difference_max"] = (float(dif[fora & ok][np.argmax(np.abs(dif[fora & ok]))])
                               if (fora & ok).any() else np.nan)
    if np.isfinite(zm).any():
        km = int(np.nanargmin(zm))
        r["talweg_station_mdt"] = float(st[km])
        r["talweg_elev_mdt"] = float(zm[km])
        r["talweg_desloc"] = float(st[km] - st[k])
    else:
        r["talweg_station_mdt"] = r["talweg_elev_mdt"] = r["talweg_desloc"] = np.nan

    mm = _mediana_movel(z)
    r["spikes"] = int(np.sum(fora & (np.abs(z - mm) > SPIKE)))
    r["spike_max"] = float(np.max(np.abs(z - mm))) if len(z) else 0.0

    # ---------------- status
    motivos, sev = [], 0

    def marca(cond, nivel, texto):
        nonlocal sev
        if cond:
            motivos.append(texto)
            sev = max(sev, nivel)

    marca(r["station_length_error"] > TOL_LEN_C, 2,
          f"cutline x stations difere {r['station_length_error']:.2f} m")
    marca(TOL_LEN_W < r["station_length_error"] <= TOL_LEN_C, 1,
          f"cutline x stations difere {r['station_length_error']:.3f} m")
    marca(ncr == 0, 2, "cutline nao cruza o eixo do reach")
    marca(ncr >= 2, 1, f"cutline cruza o eixo do MESMO reach {ncr} vezes")
    marca(np.isfinite(r["angle_desvio_90"]) and r["angle_desvio_90"] > ANG_C, 2,
          f"cutline a {r['angle']:.0f} graus da tangente (quase paralela ao fluxo)")
    marca(np.isfinite(r["angle_desvio_90"])
          and ANG_W < r["angle_desvio_90"] <= ANG_C, 1,
          f"cutline a {r['angle']:.0f} graus da tangente")
    marca(ov > 0, 1, f"cruza a cutline de {ov} vizinha(s)")
    marca(r["dog_leg"] > DOGLEG, 1, f"dog-leg de {r['dog_leg']:.0f} graus")
    marca(r["estacas_fora_de_ordem"] > 0, 2, "estacas fora de ordem")
    marca(r["pontos_duplicados"] > 0, 1, "estacas repetidas")
    marca(lb is not None and rb is not None and rb <= lb, 2, "margens invertidas")
    marca(np.isfinite(r["dem_difference"]) and abs(r["dem_difference"]) > DEM_C, 2,
          f"overbank {r['dem_difference']:+.1f} m fora do MDT")
    marca(np.isfinite(r["dem_difference"])
          and DEM_W < abs(r["dem_difference"]) <= DEM_C, 1,
          f"overbank {r['dem_difference']:+.1f} m fora do MDT")
    marca(r["spikes"] > 0, 1, f"{r['spikes']} spike(s) fora do canal")
    marca(r["mdt_cobertura"] < 0.9, 1,
          f"MDT cobre so {100*r['mdt_cobertura']:.0f}% da secao")

    r["status"] = ("CRITICAL" if sev == 2 else "WARNING" if sev == 1 else "OK")
    r["reason"] = "; ".join(motivos)
    r["_z"] = z; r["_zm"] = zm; r["_st"] = st
    return r


CAMPOS = ["River", "Reach", "RiverStation", "n_pontos", "length",
          "station_length", "station_length_error", "angle", "angle_cutline",
          "angle_river", "angle_desvio_90", "tangent_window",
          "river_intersections", "neighbor_overlap", "dist_ant", "dist_prox",
          "dog_leg", "talweg_station", "talweg_fraction", "talweg_elev",
          "talweg_station_mdt", "talweg_elev_mdt", "talweg_desloc",
          "channel_width", "bank_left", "bank_right",
          "bank_elev_left", "bank_elev_right",
          "pontos_duplicados", "estacas_fora_de_ordem", "salto_max", "incl_max",
          "spikes", "spike_max", "mdt_cobertura",
          "dem_difference", "dem_difference_canal", "dem_difference_max",
          "status", "reason"]


def main(argv):
    prj = argv[0]
    g01 = geometria_do_projeto(prj)
    raiz = os.path.dirname(prj) or "."
    nome = os.path.splitext(os.path.basename(prj))[0]
    print(f"projeto  : {prj}")
    print(f"geometria: {g01}")
    secoes = ler_secoes(g01)
    eixos = ler_eixos(g01)
    print(f"secoes   : {len(secoes)}   reaches: {len(eixos)}")

    P = np.vstack([d["cut"] for d in secoes if "cut" in d])
    bbox = (P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max())
    lista = os.path.join(raiz, f"sigsc_tiles_{nome}.txt")
    if os.path.exists(lista):
        tiles = open(lista).read().split("\n")
    else:
        tiles = tiles_do_dominio(bbox)
        open(lista, "w").write("\n".join(tiles))
    print(f"MDT      : SIG-SC 1 m, {len(tiles)} folhas sobre o dominio")
    mdt = MosaicoSigsc(tiles=tiles)

    res = []
    for i, d in enumerate(secoes):
        ant = secoes[i - 1] if i > 0 and secoes[i - 1]["reach"] == d["reach"] else None
        prox = (secoes[i + 1] if i + 1 < len(secoes)
                and secoes[i + 1]["reach"] == d["reach"] else None)
        eixo = eixos.get((d["rio"].strip(), d["reach"].strip())) or \
            eixos.get((d["rio"], d["reach"])) or list(eixos.values())[0]
        res.append(auditar(d, eixo, ant, prox, mdt))
        if (i + 1) % 200 == 0:
            print(f"   {i+1}/{len(secoes)}")

    import csv
    csvp = os.path.join(raiz, f"{nome}_qc.csv")
    with open(csvp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(CAMPOS)
        for r in res:
            w.writerow([r.get(k) for k in CAMPOS])

    feats = []
    for d, r in zip(secoes, res):
        feats.append({"type": "Feature",
                      "geometry": {"type": "LineString",
                                   "coordinates": [[float(x), float(y)]
                                                   for x, y in d["cut"]]},
                      "properties": {k: (None if isinstance(r.get(k), float)
                                         and not np.isfinite(r.get(k))
                                         else r.get(k)) for k in CAMPOS}})
    gj = os.path.join(raiz, f"{nome}_qc.geojson")
    json.dump({"type": "FeatureCollection",
               "crs": {"type": "name",
                       "properties": {"name": "urn:ogc:def:crs:EPSG::31982"}},
               "features": feats}, open(gj, "w"), ensure_ascii=False)

    from collections import Counter
    c = Counter(r["status"] for r in res)
    print()
    print("=" * 74)
    print(f"STATUS   OK {c['OK']:5d}   WARNING {c['WARNING']:5d}   "
          f"CRITICAL {c['CRITICAL']:5d}   de {len(res)}")
    print("=" * 74)
    mot = Counter()
    for r in res:
        for m in r["reason"].split("; "):
            if m:
                mot[m.split(" ")[0] + " " + (m.split(" ")[1] if len(m.split(" ")) > 1 else "")] += 1
    print(f"\ntabela : {csvp}")
    print(f"geojson: {gj}")
    return res, secoes, eixos, mdt


if __name__ == "__main__":
    main(sys.argv[1:])
