# -*- coding: utf-8 -*-
"""
Topologia da rede, da base da ANA (BHO 2017).

Da ANA vem o que ela tem de confiavel -- QUEM desagua em QUEM e a AREA DE
DRENAGEM de cada rio. A geometria (por onde o rio passa) vem do relevo; ver
itajai/tracado.py.

A separacao nao e capricho. Deixar a topologia sair da geometria quebra a rede:
Sul, Oeste e Acu se encontram no MESMO ponto em Rio do Sul, e pelo criterio de
"rio maior mais proximo da foz" o Sul acaba pendurado no Oeste em vez do Acu.
"""
import geopandas as gpd
from shapely.geometry import LineString, Point

from .config import EPSG

BASE = "rios_itajai.geojson"
AREA_MINIMA = 200.0      # km2; abaixo disso e torrente de cabeceira

# chave -> (padrao no nome da ANA, nome no HEC-RAS)
RIOS = {
    "acu":      ("Itajaí-açu",               "Itajai_Acu"),
    "sul":      ("Itajaí do Sul",            "Itajai_Sul"),
    "oeste":    ("Itajaí do Oeste",          "Itajai_Oeste"),
    "norte":    ("Itajaí do Norte|Hercílio", "Itajai_Norte"),
    "benedito": ("Benedito",                 "Rio_Benedito"),
    "mirim":    ("Itajaí-mirim",             "Itajai_Mirim"),
    # afluentes de 2a ordem. Nenhum destes desagua no Acu: Trombudo, das
    # Pombas e Taio entram no Oeste; Iraputa no Norte (o Hercilio); dos Cedros
    # no Benedito. So o do Testo entra na calha principal.
    "trombudo": ("Trombudo",                 "Rio_Trombudo"),
    "pombas":   ("das Pombas",               "Rio_das_Pombas"),
    "taio":     ("Taió",                     "Rio_Taio"),
    "iraputa":  ("Iraputã",                  "Rio_Iraputa"),
    "cedros":   ("dos Cedros",               "Rio_dos_Cedros"),
    "testo":    ("do Testo",                 "Rio_do_Testo"),
}
PRINCIPAL = "acu"


def carregar(escopo=None):
    """Cadeia principal de cada rio: area, cabeceira e foz.

    Devolve dict chave -> {nome, area, cabeceira: Point, foz: Point,
                           linha_ana: LineString}
    """
    g = gpd.read_file(BASE).to_crs(EPSG)
    g["NORIOCOMP"] = g["NORIOCOMP"].astype(str)
    por_cod = {int(r.COTRECHO): r for r in g.itertuples()}
    montante = {}
    for r in g.itertuples():
        montante.setdefault(int(r.NUTRJUS), []).append(int(r.COTRECHO))

    def cadeia(padrao):
        """Da foz para montante, sempre pelo ramo de maior area."""
        sub = g[g["NORIOCOMP"].str.contains(padrao, case=False, na=False)]
        ch = [int(sub.loc[sub["NUAREAMONT"].idxmax(), "COTRECHO"])]
        alvos = [a.lower() for a in padrao.split("|")]
        while True:
            # str(): nem todo trecho da BHO tem nome, e o valor vem NaN (float)
            ups = [c for c in montante.get(ch[-1], [])
                   if c in por_cod
                   and any(a in str(por_cod[c].NORIOCOMP).lower() for a in alvos)]
            if not ups:
                break
            m = max(ups, key=lambda c: float(por_cod[c].NUAREAMONT or 0))
            if float(por_cod[m].NUAREAMONT or 0) < AREA_MINIMA:
                break
            ch.append(m)
        return ch[::-1]                       # cabeceira -> foz

    def eixo(ch):
        pts = []
        for c in ch:
            geo = por_cod[c].geometry
            for l in ([geo] if geo.geom_type == "LineString" else list(geo.geoms)):
                cc = list(l.coords)
                if pts and (Point(pts[-1]).distance(Point(cc[0])) >
                            Point(pts[-1]).distance(Point(cc[-1]))):
                    cc = cc[::-1]
                pts += cc if not pts else cc[1:]
        return LineString(pts)

    rede = {}
    for k, (padrao, nome) in RIOS.items():
        if escopo and k != PRINCIPAL and k not in escopo:
            continue
        ch = cadeia(padrao)
        ln = eixo(ch)
        rede[k] = {"nome": nome, "area": float(por_cod[ch[-1]].NUAREAMONT),
                   "cabeceira": Point(ln.coords[0]), "foz": Point(ln.coords[-1]),
                   "linha_ana": ln}
    return rede


def arvore(rede):
    """Quem desagua em quem, pelas linhas da ANA.

    O receptor e sempre um rio de area MAIOR: assim um afluente nunca e
    pendurado noutro afluente menor que por acaso passe perto da foz dele.
    """
    receptor, filhos = {}, {k: [] for k in rede}
    for k, v in rede.items():
        if k == PRINCIPAL:
            continue
        cand = [m for m in rede if m != k and rede[m]["area"] > v["area"]]
        if not cand:
            continue
        alvo = min(cand, key=lambda m: rede[m]["linha_ana"].distance(v["foz"]))
        if rede[alvo]["linha_ana"].distance(v["foz"]) > 500.0:
            continue
        receptor[k] = alvo
        filhos[alvo].append(k)
    return receptor, filhos
