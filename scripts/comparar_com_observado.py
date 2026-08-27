# -*- coding: utf-8 -*-
"""Compara a rodada de 1983 com o OBSERVADO da ANA, estacao por estacao.

    python scripts/comparar_com_observado.py taha_ai.p01.hdf \
        --series doc/ana_1983 --figura doc/figuras/calibracao_1983.png

SO LE E DESENHA. Vazao simulada (serie horaria do HDF) contra a vazao
diaria observada em cada estacao que cai dentro do modelo (projetada no
eixo por levantar_estacoes_ana + pyproj; RS abaixo). Imprime pico a pico
o erro -- a regua da calibracao.
"""
import csv
import datetime
import os
import sys

import numpy as np

# estacao -> (rio do modelo, RS aproximado em m, nome)
ESTACOES = {
    "83800002": ("Itajai_Acu",   64800,  "Blumenau"),
    "83300200": ("Itajai_Acu",  172900,  "Rio do Sul"),
    "83440000": ("Itajai_Norte",  3800,  "Ibirama"),
    "83345000": ("Itajai_Norte", 73600,  "Barra do Prata"),
    "83250000": ("Itajai_Sul",   24200,  "Ituporanga"),
    "83105000": ("Itajai_Sul",   75500,  "Saltinho"),
    "83050000": ("Itajai_Oeste", 56800,  "Taio"),
    "83900000": ("Itajai_Mirim", 36300,  "Brusque"),
    "83660000": ("Rio_Benedito", 22400,  "Benedito Novo"),
    "83675000": ("Rio_dos_Cedros", 11000, "Arrozeira"),
}
INICIO = datetime.datetime(1983, 7, 1)


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def obs(pasta, cod):
    for arq in os.listdir(pasta):
        if arq.startswith(cod) and arq.endswith("_vazao.csv"):
            datas, vals = [], []
            for r in csv.reader(open(os.path.join(pasta, arq),
                                     encoding="utf-8"), delimiter=";"):
                if r[0] == "data" or not r[1]:
                    continue
                d = datetime.date.fromisoformat(r[0])
                if d.year == 1983 and d.month in (6, 7, 8):
                    datas.append(datetime.datetime(d.year, d.month, d.day))
                    vals.append(float(r[1]))
            return datas, vals
    return [], []


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import h5py
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    hdf = argv[0]
    pasta = _arg(argv, "--series", "doc/ana_1983")
    fig_arq = _arg(argv, "--figura", "doc/figuras/calibracao_1983.png")

    f = h5py.File(hdf, "r")
    base = ("Results/Unsteady/Output/Output Blocks/Base Output/"
            "Unsteady Time Series")
    atrib = f["Geometry/Cross Sections/Attributes"][()]
    rios = [x.decode().strip() for x in atrib["River"]]
    rss = [float(x.decode().strip()) for x in atrib["RS"]]
    Q = f[base + "/Cross Sections/Flow"][()]
    tempo_min = f[base + "/Time"][()]            # dias desde o inicio
    datas_sim = [INICIO + datetime.timedelta(days=float(t))
                 for t in tempo_min]

    n = len(ESTACOES)
    ncol = 2
    nlin = -(-n // ncol)
    fig, eixos = plt.subplots(nlin, ncol, figsize=(14, 2.8 * nlin),
                              sharex=True)
    eixos = np.atleast_2d(eixos)
    print(f"{'estacao':16s} {'pico obs':>9s} {'pico sim':>9s} {'erro':>7s}")
    k = 0
    for cod, (rio, rs_alvo, nome) in ESTACOES.items():
        cand = [(abs(r - rs_alvo), j) for j, (ri, r) in
                enumerate(zip(rios, rss)) if ri == rio]
        if not cand:
            continue
        _, j = min(cand)
        q_sim = Q[:, j]
        dt_obs, q_obs = obs(pasta, cod)
        ax = eixos[k // ncol][k % ncol]
        ax.plot(datas_sim, q_sim, "-", color="crimson", lw=1.2,
                label="simulado")
        if q_obs:
            ax.plot(dt_obs, q_obs, ".-", color="royalblue", lw=0.8, ms=3,
                    label="observado (ANA)")
        ax.set_title(f"{nome}  ({rio} RS {rss[j]:.0f})", fontsize=9)
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=7)
        if k == 0:
            ax.legend(fontsize=8)
        if q_obs:
            po = max(q_obs)
            ps = float(np.nanmax(q_sim))
            print(f"{nome:16s} {po:8.0f} {ps:8.0f} {100*(ps-po)/po:6.0f}%")
        k += 1
    for j in range(k, nlin * ncol):
        eixos[j // ncol][j % ncol].axis("off")
    fig.suptitle("Julho/1983: vazao simulada vs observada (ANA)",
                 fontsize=12)
    fig.autofmt_xdate()
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(fig_arq) or ".", exist_ok=True)
    fig.savefig(fig_arq, dpi=130)
    print(f"\nfigura: {fig_arq}")


if __name__ == "__main__":
    main(sys.argv[1:])
