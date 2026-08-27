# -*- coding: utf-8 -*-
"""Audita as larguras do SIG-SC: DESENHA os transectos mais largos.

    python scripts/auditar_larguras.py taha_ai.g01 \
        --medidas doc/larguras_sigsc --por-rio 3 \
        --figura doc/figuras/auditoria_larguras.png

SO DESENHA. Para cada rio, refaz os `--por-rio` transectos de MAIOR
lamina e plota o perfil do terreno com a faixa que o algoritmo chamou
de agua (azul) e a calha plena (verde). Se a "agua" for arrozeira,
banhado ou reservatorio, aparece: fundo chato quilometrico em vez de
calha com barrancos.
"""
import csv
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_geometria import ler_eixos                     # noqa: E402
import largura_do_sigsc as L                           # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import rasterio
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g01 = argv[0]
    pasta = argv[argv.index("--medidas") + 1] if "--medidas" in argv \
        else "doc/larguras_sigsc"
    por_rio = int(argv[argv.index("--por-rio") + 1]) if "--por-rio" in argv \
        else 3
    fig_arq = argv[argv.index("--figura") + 1] if "--figura" in argv \
        else "doc/figuras/auditoria_larguras.png"

    src = rasterio.open(L.MDT)
    T = src.transform

    E = ler_eixos(g01)
    por = {}
    for (rio, reach), ls in E.items():
        por.setdefault(rio, []).append((reach, ls))

    alvos = []          # (rio, dist_km, lamina_csv)
    for rio in sorted(por):
        arq = os.path.join(pasta, f"{rio}.csv")
        if not os.path.exists(arq):
            continue
        med = []
        for r in csv.reader(open(arq, encoding="utf-8"), delimiter=";"):
            if r[0] == "dist_foz_km" or not r[1]:
                continue
            med.append((float(r[1]), float(r[0])))
        med.sort(reverse=True)
        for la, dk in med[:por_rio]:
            alvos.append((rio, dk, la))

    n = len(alvos)
    ncol = 3
    nlin = -(-n // ncol)
    fig, eixos = plt.subplots(nlin, ncol, figsize=(15, 2.6 * nlin))
    eixos = np.atleast_2d(eixos)

    from shapely.geometry import LineString
    for k, (rio, dk, la_csv) in enumerate(alvos):
        partes = sorted(por[rio])
        coords = []
        for _, ls in partes:
            coords += list(ls.coords)
        eixo = LineString(coords)
        s = eixo.length - dk * 1000.0
        P0 = np.asarray(eixo.interpolate(s).coords[0])
        P1 = np.asarray(eixo.interpolate(min(s + 30, eixo.length)).coords[0])
        t = P1 - P0
        t = t / max(np.hypot(*t), 1e-9)
        nvec = np.array([-t[1], t[0]])
        z = L.transecto(src, T, P0, nvec)
        x = np.arange(-L.MEIA, L.MEIA + 1)
        ax = eixos[k // ncol][k % ncol]
        ax.plot(x, z, "-", color="0.4", lw=0.8)
        m = L.medir(z)
        if m == "solto":
            ax.set_title(f"{rio}  km {dk:.1f}  sem margem (csv {la_csv:.0f})",
                         fontsize=9)
            m = None
        if m:
            la, ca, esp = m
            agua = np.abs(z - esp) <= L.TOL_LAMINA
            ax.fill_between(x, z, esp + L.TOL_LAMINA,
                            where=agua, color="royalblue", alpha=0.5)
            ax.axhline(esp, color="royalblue", lw=0.6, ls="--")
            ax.set_title(f"{rio}  km {dk:.1f}  lamina {la} m "
                         f"(csv {la_csv:.0f})", fontsize=9)
        else:
            ax.set_title(f"{rio}  km {dk:.1f}  sem agua", fontsize=9)
        ax.set_xlim(-L.MEIA, L.MEIA)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for k in range(n, nlin * ncol):
        eixos[k // ncol][k % ncol].axis("off")
    fig.suptitle("Auditoria: os transectos de MAIOR lamina de cada rio "
                 "(o que o algoritmo chamou de agua, em azul)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(fig_arq) or ".", exist_ok=True)
    fig.savefig(fig_arq, dpi=130)
    print(f"figura: {fig_arq}")


if __name__ == "__main__":
    main(sys.argv[1:])
