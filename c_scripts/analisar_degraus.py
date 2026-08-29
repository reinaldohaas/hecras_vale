# -*- coding: utf-8 -*-
"""Classifica os degraus de fundo entre secoes vizinhas: real, artefato ou erro.

    python scripts/analisar_degraus.py modelo/so_mirim.g07 --limiar 1.0

SO LE. Nao altera nada, nao alisa nada.

O TESTE QUE DECIDE

  Um degrau entre duas secoes vizinhas pode ser tres coisas. O que as separa e
  comparar o degrau REGISTRADO NAS SECOES com a queda do terreno MEDIDA AO
  LONGO DO EIXO entre elas:

      d_secoes = |min(z_i) - min(z_j)|              o que o modelo diz
      d_eixo   = |MDT(eixo em s_i) - MDT(eixo em s_j)|   o que o terreno diz

  Se o terreno explica o degrau, ele e real. Se nao explica, o degrau nasceu de
  as duas secoes amostrarem minimos em lugares diferentes -- artefato.

  Validado na RS 127448.69 (doc/QC_so_mirim_validation.md): degrau de 2,77 m
  contra 1,33 m de queda real ao longo do eixo, com degrau interno maximo de
  0,11 m no terreno. Metade real, metade artefato.

CLASSIFICACAO

  C  erro de geometria      o minimo de alguma das duas cai FORA das margens,
                            ou a secao esta a mais de 25% da largura do eixo.
                            Antes de discutir o degrau, a secao esta errada.
  A  provavelmente real     d_eixo explica >= 70% de d_secoes
  B  provavelmente artefato d_eixo explica <= 30% de d_secoes
  D  inconclusivo           entre 30% e 70%

  A ordem importa: C antes de tudo, porque secao mal posta invalida a
  comparacao. Depois A/B pelo quanto o terreno explica.

  Registra-se tambem o maior degrau INTERNO do MDT ao longo do eixo entre as
  duas: um degrau real de verdade costuma aparecer ali como salto localizado,
  e nao como rampa suave.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                     # noqa: E402
from qc_geometria import ler_eixos, tangente_local   # noqa: E402
from mdt_sigsc import MosaicoSigsc                   # noqa: E402

FORA_EIXO = 0.25
EXPLICA_A = 0.70
EXPLICA_B = 0.30


def main(argv):
    from shapely.geometry import LineString
    import csv
    g = argv[0] if argv else "modelo/so_mirim.g07"
    lim = float(argv[argv.index("--limiar") + 1]) if "--limiar" in argv else 1.0
    raiz = os.path.dirname(g) or "."
    base = os.path.basename(g).split(".")[0]
    S = ler_secoes(g)
    S.sort(key=lambda d: -d["rs"])
    eixo = list(ler_eixos(g).values())[0]
    mdt = MosaicoSigsc(tiles=open(os.path.join(
        raiz, f"sigsc_tiles_{base}.txt")).read().split("\n"))

    # posicao de cada secao sobre o eixo, e o minimo
    s_ = np.full(len(S), np.nan)
    fora = np.zeros(len(S), bool)
    zmin = np.array([float(np.min(d["z"])) for d in S])
    st_min = np.array([float(d["sta"][int(np.argmin(d["z"]))]) for d in S])
    dentro = np.array([bool(d["lb"] <= d["sta"][int(np.argmin(d["z"]))] <= d["rb"])
                       for d in S])
    for i, d in enumerate(S):
        ln = LineString(d["cut"])
        x = ln.intersection(eixo)
        if x.is_empty:
            fora[i] = True; continue
        p = x if x.geom_type == "Point" else list(x.geoms)[0]
        s_[i] = float(eixo.project(p))
        C = 0.5 * (np.array(d["cut"][0], float) + np.array(d["cut"][-1], float))
        L = float(d["sta"][-1] - d["sta"][0])
        if float(np.hypot(*(C - np.asarray(p.coords[0])))) > FORA_EIXO * L:
            fora[i] = True

    d_sec = np.abs(np.diff(zmin))
    alvo = np.flatnonzero(d_sec > lim)
    print(f"geometria: {g}")
    print(f"degraus acima de {lim:.1f} m: {len(alvo)} de {len(S)-1} pares")

    linhas = []
    for k in alvo:
        i, j = k, k + 1
        if not (np.isfinite(s_[i]) and np.isfinite(s_[j])):
            cls = "C"; d_eixo = np.nan; explic = np.nan; salto = np.nan
        else:
            a, b = sorted((s_[i], s_[j]))
            ss = np.linspace(a, b, max(20, int((b - a) / 3) + 2))
            P = [eixo.interpolate(float(x)) for x in ss]
            zz = mdt.cota([p.x for p in P], [p.y for p in P])
            ok = np.isfinite(zz)
            if ok.sum() < 3:
                d_eixo = np.nan; explic = np.nan; salto = np.nan
            else:
                d_eixo = abs(float(zz[ok][0] - zz[ok][-1]))
                explic = d_eixo / max(d_sec[k], 1e-9)
                salto = float(np.abs(np.diff(zz[ok])).max())
            if fora[i] or fora[j] or not dentro[i] or not dentro[j]:
                cls = "C"
            elif not np.isfinite(explic):
                cls = "D"
            elif explic >= EXPLICA_A:
                cls = "A"
            elif explic <= EXPLICA_B:
                cls = "B"
            else:
                cls = "D"
        linhas.append({
            "classe": cls,
            "rs_ant": S[i - 1]["rs"] if i > 0 else np.nan,
            "rs_i": S[i]["rs"], "rs_j": S[j]["rs"],
            "rs_pos": S[j + 1]["rs"] if j + 1 < len(S) else np.nan,
            "dx": S[i]["rs"] - S[j]["rs"],
            "z_i": zmin[i], "z_j": zmin[j], "degrau": zmin[j] - zmin[i],
            "d_eixo": d_eixo, "explica_pct": 100 * explic if np.isfinite(explic) else np.nan,
            "salto_mdt_eixo": salto,
            "min_i_frac": st_min[i] / float(S[i]["sta"][-1] - S[i]["sta"][0]),
            "min_j_frac": st_min[j] / float(S[j]["sta"][-1] - S[j]["sta"][0]),
            "min_i_no_canal": bool(dentro[i]), "min_j_no_canal": bool(dentro[j]),
            "lc_i": float(S[i]["rb"] - S[i]["lb"]),
            "lc_j": float(S[j]["rb"] - S[j]["lb"]),
            "larg_i": float(S[i]["sta"][-1] - S[i]["sta"][0]),
            "larg_j": float(S[j]["sta"][-1] - S[j]["sta"][0]),
            "secao_fora_do_eixo": bool(fora[i] or fora[j]),
        })

    from collections import Counter
    c = Counter(l["classe"] for l in linhas)
    rot = {"A": "provavelmente real", "B": "provavelmente artefato de amostragem",
           "C": "provavelmente erro de geometria", "D": "inconclusivo"}
    print()
    for k in "ABCD":
        print("   %s  %-42s %4d  (%4.1f%%)"
              % (k, rot[k], c.get(k, 0), 100 * c.get(k, 0) / max(len(linhas), 1)))
    saida = os.path.join(raiz, f"degraus_{base}.csv")
    campos = list(linhas[0].keys()) if linhas else []
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(campos)
        for l in sorted(linhas, key=lambda l: (l["classe"], -abs(l["degrau"]))):
            w.writerow([l[k] for k in campos])
    print(f"\ntabela: {saida}")
    return linhas


if __name__ == "__main__":
    main(sys.argv[1:])
