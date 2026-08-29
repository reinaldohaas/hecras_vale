# -*- coding: utf-8 -*-
"""Analise de conflito geometrico dos pares de cutlines que se cruzam.

    python scripts/analisar_overlaps.py modelo/so_mirim.g07

SO LE. Nao altera nada.

Para cada PAR de secoes cuja cutline cruza a da vizinha, mede o que permite
distinguir erro de geometria de consequencia inevitavel da curvatura:

    RS das duas · distancia longitudinal · angulo de cada uma em relacao a
    normal do rio · ponto exato da intersecao · lado do rio · largura de cada
    uma · raio de curvatura local do eixo · se o cruzamento cai dentro do canal
    ou na planicie · se cada secao esta centrada no eixo

CLASSIFICACAO, com o criterio explicito:

  A  erro inequivoco de geometria
     ao menos uma das duas com desvio da normal > 30 graus, OU centro da
     secao a mais de 25% da largura de distancia do eixo. A secao esta torta
     ou fora do lugar -- o cruzamento e sintoma disso.

  B  consequencia inevitavel da curvatura
     as duas praticamente normais (desvio <= 15 graus), as duas centradas, e
     a meia-largura maior que o raio de curvatura local. Duas secoes normais
     a um arco de raio R se cruzam a distancia R do centro: e geometria, nao
     defeito. Confirma-se olhando se a intersecao cai do lado INTERNO da curva.

  C  secao excessivamente larga
     nao e B (meia-largura ainda cabe no raio) mas a largura de uma delas
     passa de 3x a mediana do modelo. O cruzamento vem do tamanho, nao da
     curva.

  D  problema de espacamento
     nao e A, B nem C, e a distancia longitudinal entre as duas e menor que a
     largura media do canal delas. Estao perto demais para a escala do rio.

  E  ambiguo
     o que sobra.

A ordem de teste importa e e essa: A antes de B porque secao torta explica o
cruzamento sozinha; B antes de C porque curvatura e razao fisica e largura e
escolha; D por ultimo entre as objetivas.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                     # noqa: E402
from qc_geometria import ler_eixos, tangente_local   # noqa: E402

DESV_A = 30.0      # graus de desvio da normal que ja acusam secao torta
DESV_B = 15.0      # graus abaixo dos quais a secao e "praticamente normal"
FORA_EIXO = 0.25   # fracao da largura: centro longe do eixo
LARGA_C = 3.0      # multiplo da largura mediana


def _az(v):
    return np.degrees(np.arctan2(v[1], v[0]))


def _ang(a, b):
    d = abs((a - b + 180.0) % 360.0 - 180.0)
    return min(d, 180.0 - d)


def raio_local(eixo, s, jan):
    """R = ds / dtheta, com o angulo medido em duas cordas de meia janela."""
    h = max(jan, 5.0)
    a = np.array(eixo.interpolate(max(0.0, s - h)).coords[0])
    m = np.array(eixo.interpolate(s).coords[0])
    b = np.array(eixo.interpolate(min(eixo.length, s + h)).coords[0])
    t1, t2 = m - a, b - m
    dth = np.radians(_ang(_az(t1), _az(t2)))
    ds = float(np.hypot(*t1) + np.hypot(*t2)) / 2.0
    if dth < 1e-6:
        return np.inf, 0.0
    return ds / dth, np.degrees(dth)


def analisar(g01, eixo=None):
    from shapely.geometry import LineString, Point
    S = ler_secoes(g01)
    S.sort(key=lambda d: -d["rs"])
    if eixo is None:
        eixo = list(ler_eixos(g01).values())[0]
    larg = np.array([float(d["sta"][-1] - d["sta"][0]) for d in S])
    larg_med = float(np.median(larg))
    cut = [LineString(d["cut"]) for d in S]

    info = []
    for i, d in enumerate(S):
        g = cut[i].intersection(eixo)
        if g.is_empty:
            info.append(None); continue
        p = g if g.geom_type == "Point" else list(g.geoms)[0]
        s_ = float(eixo.project(p))
        lc = float(d["rb"] - d["lb"])
        jan = float(np.clip(2.0 * lc, 20.0, 150.0))
        t = tangente_local(eixo, s_, jan)
        A = np.array(d["cut"][0], float); B = np.array(d["cut"][-1], float)
        desv = abs(90.0 - _ang(_az(B - A), _az(t)))
        C = 0.5 * (A + B)
        R, dth = raio_local(eixo, s_, jan)
        info.append({"s": s_, "p": np.asarray(p.coords[0]), "t": t,
                     "desv": desv, "centro": C,
                     "fora_eixo": float(np.hypot(*(C - np.asarray(p.coords[0])))),
                     "R": R, "dtheta": dth, "larg": larg[i], "lc": lc})

    pares = []
    for i in range(len(S) - 1):
        for j in (i + 1,):
            if info[i] is None or info[j] is None:
                continue
            if not cut[i].intersects(cut[j]):
                continue
            x = cut[i].intersection(cut[j])
            if x.is_empty:
                continue
            q = x if x.geom_type == "Point" else list(x.geoms)[0]
            Q = np.asarray(q.coords[0])
            di = float(S[i]["rs"] - S[j]["rs"])
            # lado do rio em que a intersecao cai, e se e o lado interno
            t = info[i]["t"]; t = t / max(float(np.hypot(*t)), 1e-9)
            v = Q - info[i]["p"]
            lado = float(t[0] * v[1] - t[1] * v[0])     # >0 esquerda do fluxo
            # o lado INTERNO da curva e o lado para onde o eixo vira
            ta = info[i]["t"]; tb = info[j]["t"]
            vira = float(ta[0] * tb[1] - ta[1] * tb[0])
            interno = (np.sign(lado) == np.sign(vira)) and abs(vira) > 1e-9
            # a intersecao cai dentro do canal de alguma das duas?
            def no_canal(k):
                A = np.array(S[k]["cut"][0], float); B = np.array(S[k]["cut"][-1], float)
                L = float(np.hypot(*(B - A)))
                if L < 1e-9:
                    return False
                u = (B - A) / L
                st = float(np.dot(Q - A, u))
                return S[k]["lb"] <= st <= S[k]["rb"]
            dist_eixo = float(np.hypot(*(Q - info[i]["p"])))
            pares.append({
                "rs_i": S[i]["rs"], "rs_j": S[j]["rs"], "dx": di,
                "desv_i": info[i]["desv"], "desv_j": info[j]["desv"],
                "x": float(Q[0]), "y": float(Q[1]),
                "lado": "esquerda" if lado > 0 else "direita",
                "interno": bool(interno),
                "larg_i": info[i]["larg"], "larg_j": info[j]["larg"],
                "lc_i": info[i]["lc"], "lc_j": info[j]["lc"],
                "R": info[i]["R"], "dtheta": info[i]["dtheta"],
                "fora_i": info[i]["fora_eixo"], "fora_j": info[j]["fora_eixo"],
                "dist_eixo": dist_eixo,
                "canal": no_canal(i) or no_canal(j),
                "meia_larg": 0.5 * max(info[i]["larg"], info[j]["larg"]),
            })

    # ------------------------------------------------------------ classificar
    for p in pares:
        torta = max(p["desv_i"], p["desv_j"]) > DESV_A
        deslocada = (p["fora_i"] > FORA_EIXO * p["larg_i"]
                     or p["fora_j"] > FORA_EIXO * p["larg_j"])
        normais = max(p["desv_i"], p["desv_j"]) <= DESV_B
        centradas = not deslocada
        curva_aperta = p["meia_larg"] > p["R"]
        larga = max(p["larg_i"], p["larg_j"]) > LARGA_C * larg_med
        perto = p["dx"] < 0.5 * (p["lc_i"] + p["lc_j"])
        if torta or deslocada:
            c = "A"
        elif normais and centradas and curva_aperta:
            c = "B"
        elif larga:
            c = "C"
        elif perto:
            c = "D"
        else:
            c = "E"
        p["classe"] = c
    return pares, S, larg_med


def main(argv):
    import csv
    g = argv[0] if argv else "modelo/so_mirim.g07"
    pares, S, lm = analisar(g)
    from collections import Counter
    c = Counter(p["classe"] for p in pares)
    print(f"geometria: {g}")
    print(f"largura mediana do modelo: {lm:.1f} m")
    print(f"pares de cutlines que se cruzam: {len(pares)}")
    print(f"secoes envolvidas: {len({p['rs_i'] for p in pares} | {p['rs_j'] for p in pares})}")
    print()
    rot = {"A": "erro inequivoco de geometria",
           "B": "consequencia inevitavel da curvatura",
           "C": "secao excessivamente larga",
           "D": "problema de espacamento",
           "E": "ambiguo"}
    for k in "ABCDE":
        print("   %s  %-40s %4d  (%4.1f%%)"
              % (k, rot[k], c.get(k, 0), 100 * c.get(k, 0) / max(len(pares), 1)))
    saida = os.path.join(os.path.dirname(g) or ".", "overlaps_g07.csv")
    campos = ["classe", "rs_i", "rs_j", "dx", "desv_i", "desv_j", "larg_i",
              "larg_j", "lc_i", "lc_j", "R", "dtheta", "meia_larg", "fora_i",
              "fora_j", "x", "y", "lado", "interno", "canal", "dist_eixo"]
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(campos)
        for p in sorted(pares, key=lambda p: (p["classe"], -p["rs_i"])):
            w.writerow([p[k] for k in campos])
    print(f"\ntabela: {saida}")
    return pares


if __name__ == "__main__":
    main(sys.argv[1:])
