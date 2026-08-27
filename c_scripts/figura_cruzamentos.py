# -*- coding: utf-8 -*-
"""Desenha os cruzamentos de cutline de uma geometria do HEC-RAS.

    python scripts/figura_cruzamentos.py modelo/so_mirim.g01

Existe porque o contador do "Validate Geometry" da um NUMERO e nao um lugar, e
numeros diferentes na mesma geometria. O que da para apontar no mapa e o
cruzamento entre linhas de corte -- duas secoes que se tocam representam a
mesma agua duas vezes, e isso e verificavel aqui, sem depender do contador.

Sai um PNG: o rio inteiro com as secoes envolvidas em vermelho, e ampliacoes
dos focos onde elas se concentram.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from shapely.strtree import STRtree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402

N_ZOOM = 3
RAIO = 900.0      # m de meia-janela nas ampliacoes


def cruzamentos(S):
    L = [LineString(np.asarray(d["cut"], float)) for d in S]
    tree = STRtree(L)
    pares = set()
    for i, g in enumerate(L):
        for j in tree.query(g):
            j = int(j)
            if j > i and g.intersects(L[j]):
                pares.add((i, j))
    return L, sorted(pares)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    geom = argv[0]
    saida = argv[argv.index("--saida") + 1] if "--saida" in argv else \
        os.path.join("doc", os.path.basename(geom).replace(".", "_")
                     + "_cruzamentos.png")
    os.makedirs(os.path.dirname(saida) or ".", exist_ok=True)

    S = ler_secoes(geom)
    S.sort(key=lambda d: -d["rs"])
    L, pares = cruzamentos(S)
    env = sorted({i for p in pares for i in p})
    print(f"{geom}: {len(S)} secoes   {len(pares)} pares cruzados   "
          f"{len(env)} secoes envolvidas ({100*len(env)/len(S):.0f}%)")

    # focos: agrupa por proximidade ao longo do rio (indice ordenado por RS)
    focos, atual = [], [env[0]] if env else []
    for a, b in zip(env, env[1:]):
        if b - a <= 12:
            atual.append(b)
        else:
            focos.append(atual); atual = [b]
    if atual:
        focos.append(atual)
    focos.sort(key=len, reverse=True)
    print(f"focos (aglomerados ao longo do rio): {len(focos)}   "
          f"maiores: {[len(f) for f in focos[:N_ZOOM]]}")

    eixos = ler_eixos(geom)
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, N_ZOOM, height_ratios=[2.0, 1.25], hspace=0.22,
                          wspace=0.18)

    ax = fig.add_subplot(gs[0, :])
    for e in eixos.values():
        x, y = np.asarray(e.coords).T
        ax.plot(x, y, "-", color="#3b7dd8", lw=1.4, zorder=1,
                label="eixo do rio")
    for i, g in enumerate(L):
        x, y = np.asarray(g.coords).T
        mau = i in set(env)
        ax.plot(x, y, "-", lw=1.5 if mau else 0.4,
                color="#d62728" if mau else "#999999",
                alpha=1.0 if mau else 0.45, zorder=3 if mau else 2)
    for k, f in enumerate(focos[:N_ZOOM]):
        c = np.vstack([np.asarray(L[i].coords) for i in f]).mean(0)
        ax.plot(*c, "o", ms=16, mfc="none", mec="#111111", mew=1.8, zorder=5)
        ax.annotate(str(k + 1), c, ha="center", va="center",
                    fontsize=11, fontweight="bold", zorder=6)
    ax.set_aspect("equal")
    ax.set_title(f"{os.path.basename(geom)} -- {len(pares)} pares de linhas de "
                 f"corte se cruzam, envolvendo {len(env)} das {len(S)} secoes "
                 f"({100*len(env)/len(S):.0f}%)", fontsize=12)
    h = [plt.Line2D([], [], color="#999999", lw=1.0),
         plt.Line2D([], [], color="#d62728", lw=1.8),
         plt.Line2D([], [], color="#3b7dd8", lw=1.4)]
    ax.legend(h, ["secao sem conflito", "secao que cruza outra", "eixo do rio"],
              loc="upper left", fontsize=9)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=7)

    for k, f in enumerate(focos[:N_ZOOM]):
        a = fig.add_subplot(gs[1, k])
        c = np.vstack([np.asarray(L[i].coords) for i in f]).mean(0)
        for e in eixos.values():
            x, y = np.asarray(e.coords).T
            a.plot(x, y, "-", color="#3b7dd8", lw=1.8, zorder=1)
        alvo = set(f)
        for i, g in enumerate(L):
            x, y = np.asarray(g.coords).T
            if max(abs(x.mean() - c[0]), abs(y.mean() - c[1])) > 2 * RAIO:
                continue
            mau = i in alvo
            a.plot(x, y, "-", lw=1.8 if mau else 0.8,
                   color="#d62728" if mau else "#999999", zorder=3 if mau else 2)
            if mau:
                a.annotate(f"{S[i]['rs']:.0f}", (x.mean(), y.mean()),
                           fontsize=6, color="#7a1216", ha="center")
        a.set_xlim(c[0] - RAIO, c[0] + RAIO)
        a.set_ylim(c[1] - RAIO, c[1] + RAIO)
        a.set_aspect("equal")
        a.set_title(f"foco {k+1}: {len(f)} secoes   "
                    f"RS {S[max(f)]['rs']:.0f} a {S[min(f)]['rs']:.0f}",
                    fontsize=9)
        a.ticklabel_format(style="plain", useOffset=False)
        a.tick_params(labelsize=6)

    fig.savefig(saida, dpi=135, bbox_inches="tight")
    print(f"figura: {saida}")
    return saida


if __name__ == "__main__":
    main(sys.argv[1:])
