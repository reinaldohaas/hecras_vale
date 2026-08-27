# -*- coding: utf-8 -*-
"""Gera do zero a geometria de qualquer rio da bacia, com o metodo do mirim_novo.

    python scripts/gerar_rio_do_zero.py --rio Itajai_Norte --saida modelo/norte

Generaliza `gerar_mirim_do_zero.py`, que estava amarrado ao Itajai-Mirim: o
eixo vem de `eixos_do_relevo.geojson` pelo nome, e as larguras deixam de ser
constantes escritas no codigo.

POR QUE SECAO MODERADA, E NAO SECAO QUE CONTEM A CHEIA

  Os dois extremos estao medidos nesta bacia, com o mesmo validador:

      legado/Itajai_Rede_1983   secoes de 1.474 a 4.358 m   1.963 erros (1.492 Fatal)
      modelo/mirim_novo         secoes de     205 m             2 erros (1 Fatal)

  Secao larga demais em rio meandrico atravessa o proprio eixo e embaraca as
  edge lines -- e a rede legada mostra o resultado: 730 "XS intersects < 2
  banklines" e 653 avisos de Bank Station descolada da bankline. Aqui a
  meia-largura fica presa entre `MEIA_MIN` e `MEIA_MAX`.

  O preco esta declarado: secao moderada NAO contem a cheia de projeto na
  varzea plana, e onde a planicie e larga isso pede area de armazenamento ou
  2D, e nao secao maior. O relatorio no fim diz quantas secoes ficam abaixo da
  lamina esperada.

DE ONDE VEM A LARGURA DE CALHA

  Da rede legada, que e a unica fonte medida para estes rios -- mediana por
  reach, em `CALHA`. Ela varia ao longo do rio (o rio alarga para jusante),
  entao entra como rampa entre a cota de montante e a de jusante.

  Nao ha batimetria para nenhum deles. O perfil sai do MDT SIG-SC 1 m, que ve
  a LAMINA e nao o fundo -- a mesma ressalva que vale para o canal do Mirim.
"""
import argparse
import json
import os
import sys

import numpy as np
from shapely.geometry import LineString

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))

from ras_io import escrever                              # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio      # noqa: E402
from gerar_mirim_do_zero import amostrar_terreno_e_perfil  # noqa: E402

EIXOS = "eixos_do_relevo.geojson"
DX = 150.0
JANELA_TANGENTE = 60.0
MEIA_MIN, MEIA_MAX = 100.0, 300.0
FATOR_MEIA = 2.5          # meia-largura alvo, em larguras de calha

# largura de calha (m) medida na rede legada: (montante, jusante).
# A mediana do reach e tomada como valor de JUSANTE; a montante entra metade,
# porque o rio estreita para a cabeceira e a rede legada nao resolve isso.
CALHA = {
    "Itajai_Acu":   (60.0, 450.0),
    "Itajai_Norte": (45.0,  90.0),
    "Itajai_Oeste": (45.0,  92.0),
    "Itajai_Sul":   (50.0, 100.0),
    "Rio_Benedito": (36.0,  72.0),
    "Itajai_Mirim": (37.0,  74.0),
}
WKT = ('PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",'
       'DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],'
       'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
       'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
       'PARAMETER["False_Northing",10000000.0],'
       'PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],'
       'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')


def eixo_do_rio(nome, caminho=EIXOS):
    d = json.load(open(caminho, encoding="utf-8"))
    for f in d["features"]:
        if f["properties"].get("nome") == nome:
            return LineString(np.asarray(f["geometry"]["coordinates"], float))
    disp = [f["properties"].get("nome") for f in d["features"]]
    raise SystemExit(f"rio '{nome}' nao esta em {caminho}. Ha: {disp}")


def cutlines(eixo, calha0, calha1, dx=DX):
    """Cutlines ortogonais, com calha e meia-largura em rampa."""
    L = eixo.length
    est = np.arange(0.0, L, dx)
    if L - est[-1] > 20.0:
        est = np.append(est, L)
    secoes = []
    for s in est:
        p = np.array(eixo.interpolate(s).coords[0])
        a = np.array(eixo.interpolate(max(s - JANELA_TANGENTE, 0.0)).coords[0])
        b = np.array(eixo.interpolate(min(s + JANELA_TANGENTE, L)).coords[0])
        t = b - a
        nt = float(np.hypot(*t))
        if nt < 1e-6:
            continue
        t /= nt
        n = np.array([-t[1], t[0]])
        f = s / L                       # 0 na cabeceira, 1 na foz
        calha = calha0 + (calha1 - calha0) * f
        meia = float(np.clip(FATOR_MEIA * calha, MEIA_MIN, MEIA_MAX))
        secoes.append({
            "s_eixo": float(s), "rs": round(float(L - s), 2),
            "cut": LineString([p + n * meia, p - n * meia]),
            "meia_largura": meia, "largura_canal": float(calha),
            "no_canal": False})
    return secoes


def escrever_g01(caminho, secoes, eixo, rio, reach="R1", titulo=None):
    l = [f"Geom Title={titulo or rio}", "Program Version=7.01"]
    P = np.vstack([np.asarray(s["cut"].coords) for s in secoes])
    l.append("Viewing Rectangle= %.2f , %.2f , %.2f , %.2f "
             % (P[:, 0].min(), P[:, 0].max(), P[:, 1].max(), P[:, 1].min()))
    l.append("Spatial Reference System=" + WKT)
    l.append("")
    l.append(f"River Reach={rio:<16.16},{reach:<16.16}")
    c = list(eixo.coords)
    l.append(f"Reach XY= {len(c)} ")
    s_ = [f"{x:16.4f}{y:16.4f}" for x, y in c]
    for k in range(0, len(s_), 2):
        l.append("".join(s_[k:k + 2]))
    l.append("Rch Text X Y=0,0,0,0")
    l.append("")
    for i, s in enumerate(secoes):
        d = (round(float(abs(secoes[i + 1]["s_eixo"] - s["s_eixo"])), 2)
             if i + 1 < len(secoes) else 0.0)
        l.append(f"Type RM Length L Ch R = 1 ,{s['rs']:.2f},"
                 f"{d:8.2f},{d:8.2f},{d:8.2f}")
        l.append(f"Bank Sta={s['lb']:.2f},{s['rb']:.2f}")
        cc = list(s["cut"].coords)
        l.append(f"XS GIS Cut Line= {len(cc)}")
        l.append("".join(f"{p[0]:16.2f}{p[1]:16.2f}" for p in cc))
        l.append(f"#Sta/Elev= {len(s['sta'])} ")
        pf = [f"{a:8.2f}{b:8.2f}" for a, b in zip(s["sta"], s["z"])]
        for k in range(0, len(pf), 5):
            l.append("".join(pf[k:k + 5]))
        l.append("#Mann= 3 , 0 , 0 ")
        l.append(f"{s['sta'][0]:8.2f}{0.055:8.3f}{0:8d}"
                 f"{s['lb']:8.2f}{0.032:8.3f}{0:8d}"
                 f"{s['rb']:8.2f}{0.055:8.3f}{0:8d}")
        l.append(f"XS HTab Starting El and Incr={s['z_min']+0.02:.2f},"
                 "0.100, 500 ")
        l.append("XS HTab Horizontal Distribution=-1,-1,-1")
        l.append("XS Rating Curve= 0 ,0")
        l.append("Exp/Cntr=0.3,0.1")
        l.append("")
    escrever(caminho, "\n".join(l))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rio", required=True)
    ap.add_argument("--saida", required=True)
    ap.add_argument("--reach", default="R1")
    ap.add_argument("--dx", type=float, default=DX)
    ap.add_argument("--calha-montante", type=float, default=None)
    ap.add_argument("--calha-jusante", type=float, default=None)
    a = ap.parse_args()

    if a.rio not in CALHA and (a.calha_montante is None):
        raise SystemExit(f"sem largura de calha para '{a.rio}'. "
                         "Informe --calha-montante e --calha-jusante.")
    c0, c1 = CALHA.get(a.rio, (None, None))
    c0 = a.calha_montante if a.calha_montante is not None else c0
    c1 = a.calha_jusante if a.calha_jusante is not None else c1

    eixo = eixo_do_rio(a.rio)
    print(f"rio    : {a.rio}   eixo {eixo.length/1000:.2f} km   "
          f"{len(eixo.coords)} vertices")
    print(f"calha  : {c0:.0f} m na cabeceira -> {c1:.0f} m na foz")
    print(f"secao  : {2*np.clip(FATOR_MEIA*c0, MEIA_MIN, MEIA_MAX):.0f} m -> "
          f"{2*np.clip(FATOR_MEIA*c1, MEIA_MIN, MEIA_MAX):.0f} m   "
          f"(teto de {2*MEIA_MAX:.0f} m)")
    S = cutlines(eixo, c0, c1, a.dx)
    print(f"secoes : {len(S)}   a cada {a.dx:g} m")

    P = np.vstack([np.asarray(s["cut"].coords) for s in S])
    bb = (P[:, 0].min() - 60, P[:, 1].min() - 60,
          P[:, 0].max() + 60, P[:, 1].max() + 60)
    tiles = tiles_do_dominio(bb)
    print(f"MDT    : {len(tiles)} folhas do SIG-SC")
    S = amostrar_terreno_e_perfil(S, MosaicoSigsc(tiles=tiles))

    os.makedirs(a.saida, exist_ok=True)
    nome = os.path.basename(a.saida.rstrip("/\\"))
    g = os.path.join(a.saida, f"{nome}.g01")
    escrever_g01(g, S, eixo, a.rio, a.reach, nome)

    z = np.array([float(s["z_min"]) for s in S])
    lc = np.array([float(s["rb"] - s["lb"]) for s in S])
    ls = np.array([float(s["sta"][-1]) for s in S])
    d = np.diff(z)
    print(f"\ngeometria: {g}")
    print(f"   leito       : {z.min():.2f} a {z.max():.2f} m   "
          f"declividade media {(z.max()-z.min())/eixo.length:.5f}")
    print(f"   sobem p/ jusante: {int((d > 1e-9).sum())}   "
          f"pares de cota igual: {int((np.abs(d) < 0.005).sum())}")
    print(f"   secao       : mediana {np.median(ls):.0f} m   "
          f"max {ls.max():.0f} m")
    print(f"   calha       : mediana {np.median(lc):.0f} m")
    print(f"   pontos/secao: mediana {np.median([len(s['sta']) for s in S]):.0f}"
          f"   max {max(len(s['sta']) for s in S)}   (limite 500)")
    return g


if __name__ == "__main__":
    main()
