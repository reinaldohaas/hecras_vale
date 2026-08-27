# -*- coding: utf-8 -*-
"""Ajusta curvas-chave Q(h) com os pares cota x vazao diarios da ANA.

    python scripts/ajustar_curva_chave.py --series doc/ana_8085 \
        --saida doc/curvas_chave

Para cada estacao com cota E vazao no periodo, casa os dias e ajusta

    Q = a * (h - h0)^b        (h em m; cota da ANA vem em cm)

com h0 varrido em grade e (a, b) por minimos quadrados em log. Sai:

    parametros.csv            a, b, h0, faixa valida, erro relativo
    doc/figuras/curvas_chave.png  dispersao + curva por estacao, com o
                              pico de julho/1983 marcado (estrela)

E a regua para conferir o modelo nas 11 estacoes: com o nivel simulado
tira-se a vazao que a ANA teria publicado, e vice-versa.
"""
import csv
import glob
import os
import re
import sys

import numpy as np


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def ler(arq):
    out = {}
    for r in csv.reader(open(arq, encoding="utf-8"), delimiter=";"):
        if r[0] == "data" or len(r) < 2 or not r[1]:
            continue
        out[r[0]] = float(r[1])
    return out


def ajustar(h, q):
    """(a, b, h0) minimizando erro quadratico em log Q."""
    melhor = None
    for h0 in np.arange(h.min() - 3.0, h.min() + 0.01, 0.05):
        x = h - h0
        if (x <= 0).any():
            continue
        A = np.vstack([np.log(x), np.ones_like(x)]).T
        b, la = np.linalg.lstsq(A, np.log(q), rcond=None)[0]
        e = float(((A @ [b, la] - np.log(q)) ** 2).mean())
        if melhor is None or e < melhor[0]:
            melhor = (e, float(np.exp(la)), float(b), float(h0))
    return melhor


def main(argv):
    pasta = _arg(argv, "--series", "doc/ana_8085")
    saida = _arg(argv, "--saida", "doc/curvas_chave")
    pasta83 = _arg(argv, "--evento", "doc/ana_1983")
    os.makedirs(saida, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    estacoes = sorted({re.match(r"(\d+_.+)_vazao\.csv",
                                os.path.basename(p)).group(1)
                       for p in glob.glob(os.path.join(pasta,
                                                       "*_vazao.csv"))})
    linhas = []
    n = len(estacoes)
    ncol = 4
    nlin = -(-n // ncol)
    fig, eixos = plt.subplots(nlin, ncol, figsize=(16, 3.6 * nlin))
    eixos = np.atleast_2d(eixos)
    k = 0
    for est in estacoes:
        arq_q = os.path.join(pasta, f"{est}_vazao.csv")
        arq_h = os.path.join(pasta, f"{est}_cota.csv")
        if not os.path.exists(arq_h):
            continue
        Q = ler(arq_q)
        H = ler(arq_h)
        datas = sorted(set(Q) & set(H))
        if len(datas) < 60:
            continue
        q = np.array([Q[d] for d in datas])
        h = np.array([H[d] for d in datas]) / 100.0
        ok = (q > 0) & np.isfinite(h)
        q, h = q[ok], h[ok]
        r = ajustar(h, q)
        if r is None:
            continue
        e, a, b, h0 = r
        pred = a * np.maximum(h - h0, 1e-9) ** b
        erro = float(np.median(np.abs(pred - q) / q)) * 100

        ax = eixos[k // ncol][k % ncol]
        ax.plot(h, q, ".", ms=2.5, color="royalblue", alpha=0.5,
                label="diario 1980-85")
        hh = np.linspace(h.min(), h.max() * 1.15, 200)
        ax.plot(hh, a * np.maximum(hh - h0, 1e-9) ** b, "-",
                color="crimson", lw=1.5,
                label=f"Q={a:.1f}(h-{h0:.2f})^{b:.2f}")
        # pico de julho/1983
        try:
            Q83 = ler(os.path.join(pasta83, f"{est}_vazao.csv"))
            H83 = ler(os.path.join(pasta83, f"{est}_cota.csv"))
            dpico = max(Q83, key=Q83.get)
            if dpico in H83:
                ax.plot(H83[dpico] / 100.0, Q83[dpico], "*", ms=14,
                        color="darkorange", label="pico jul/1983")
        except FileNotFoundError:
            pass
        ax.set_title(est.replace("_", " "), fontsize=9)
        ax.set_xlabel("cota (m)", fontsize=8)
        ax.set_ylabel("Q (m3/s)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.5)
        k += 1
        linhas.append([est, f"{a:.3f}", f"{b:.3f}", f"{h0:.2f}",
                       f"{h.min():.2f}", f"{h.max():.2f}", len(q),
                       f"{erro:.1f}"])
        print(f"   {est:28s} Q={a:7.1f}(h-{h0:5.2f})^{b:4.2f}  "
              f"[{h.min():.2f}..{h.max():.2f} m]  n={len(q)}  "
              f"erro med {erro:.0f}%")
    for j in range(k, nlin * ncol):
        eixos[j // ncol][j % ncol].axis("off")
    fig.suptitle("Curvas-chave ajustadas dos pares diarios da ANA "
                 "(1980-1985), pico de jul/1983 em estrela", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig_arq = "doc/figuras/curvas_chave.png"
    fig.savefig(fig_arq, dpi=130)

    with open(os.path.join(saida, "parametros.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["estacao", "a", "b", "h0_m", "h_min_m", "h_max_m",
                    "n_pares", "erro_mediano_pct"])
        w.writerows(linhas)
    print(f"\nparametros: {saida}/parametros.csv")
    print(f"figura    : {fig_arq}")


if __name__ == "__main__":
    main(sys.argv[1:])
