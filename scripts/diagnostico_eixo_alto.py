# -*- coding: utf-8 -*-
"""Onde o eixo do relevo corre ALTO, e a batimetria do legado nao serve.

    python scripts/diagnostico_eixo_alto.py
    python scripts/diagnostico_eixo_alto.py --rios Rio_Benedito Itajai_Acu

Serve a qualquer rio. Compara, ao longo do rio, o talvegue que o `rio_do_relevo`
leu do MDT (a LAMINA d'agua) com o fundo levantado em 1983 (legado). Onde os
dois andam juntos, a batimetria e um rebaixamento sadio de alguns metros. Onde
o eixo esquematico do legado -- de onde saiu `eixos_do_relevo.geojson` -- corre
pela ENCOSTA de um vale encaixado em vez do talvegue, a secao cortada do MDT
pega o fundo dezenas a centenas de metros ACIMA do canal real, e o
rebaixamento dispara.

POR QUE ISTO IMPORTA

  Ancorar a batimetria nesse trecho manda o `aplicar` DESCER o leito dezenas de
  metros dentro de uma calha que nao tem essa profundidade: cava um canion,
  contra as regras. Medido, foi o que:

    - deixou o Rio Benedito sem g02 (rebaixamento mediana 104 m a montante);
    - desestabilizou o solver do Itajai-Acu, cujo trecho de cabeceira (R1,
      RS > ~143000) pede rebaixar 55 m na mediana e ate 118 m -- o Acu sozinho
      da "Solution Solver Failed", e a rede vai instavel por causa dele.

  O conserto NAO e filtrar ponto a ponto (a interpolacao entre o que sobra so
  piora o degrau) nem inventar cota: e REFAZER O EIXO desse trecho seguindo o
  talvegue do MDT (acumulacao de fluxo), recortar as secoes nele, e so entao
  ancorar no levantamento -- que a jusante ja bate.

A figura marca a faixa onde o rebaixamento passa de `--limiar` metros: e ali
que o eixo precisa ser refeito.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from qc_secoes import ler_secoes                     # noqa: E402
from batimetria_do_legado import secoes_levantadas   # noqa: E402

LEGADO = "legado/Itajai_Rede_1983.g01"
# rio (nome no legado) -> geometria MDT crua para ler o talvegue
GEOM = {
    "Rio_Benedito": "modelo/rio_benedito/rio_benedito.g01",
    "Itajai_Acu":   "modelo/itajai_acu/itajai_acu.g01",
    "Itajai_Mirim": "modelo/itajai_mirim/itajai_mirim.g01",
    "Itajai_Norte": "modelo/itajai_norte/itajai_norte.g01",
    "Itajai_Sul":   "modelo/itajai_sul/itajai_sul.g01",
    "Itajai_Oeste": "modelo/itajai_oeste/itajai_oeste.g01",
}


def perfil(g):
    S = ler_secoes(g)
    S.sort(key=lambda d: -d["rs"])
    rs = np.array([d["rs"] for d in S])
    z = np.array([float(np.asarray(d["z"], float).min()) for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    x = np.r_[0.0, np.cumsum(ch[:-1])]
    return rs, x, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rios", nargs="*",
                    default=["Rio_Benedito", "Itajai_Acu"])
    ap.add_argument("--limiar", type=float, default=25.0,
                    help="m; rebaixamento acima disto marca eixo alto")
    ap.add_argument("--saida", default="doc/figuras/eixo_alto.png")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.saida), exist_ok=True)

    fig, axes = plt.subplots(len(a.rios), 1, figsize=(11, 4.3 * len(a.rios)),
                             squeeze=False)
    for ax, rio in zip(axes[:, 0], a.rios):
        g = GEOM.get(rio)
        if not g or not os.path.exists(g):
            ax.set_title(f"{rio}: geometria ausente")
            continue
        rs, x, z = perfil(g)
        L = secoes_levantadas(LEGADO, rio)
        Lrs, Linv = L[:, 0], L[:, 3]
        o = np.argsort(-Lrs)
        Lrs, Linv = Lrs[o], Linv[o]
        fundo = np.interp(x, np.interp(-Lrs, -rs, x), Linv)
        reb = z - fundo
        alto = reb > a.limiar
        ax.plot(x / 1000, z, color="tab:red", lw=1.5,
                label="talvegue lido do MDT (lamina)")
        ax.plot(x / 1000, fundo, color="tab:blue", lw=1.5,
                label="fundo levantado 1983 (legado)")
        if alto.any():
            ax.fill_between(x / 1000, z, fundo, where=alto, color="tab:red",
                            alpha=.25,
                            label=f"eixo alto (rebaixamento > {a.limiar:.0f} m)")
            faixa = x[alto]
            ax.set_title(f"{rio} -- eixo alto de {faixa.min()/1000:.0f} a "
                         f"{faixa.max()/1000:.0f} km "
                         f"(rebaixamento ate {reb[alto].max():.0f} m)",
                         fontsize=10)
        else:
            ax.set_title(f"{rio} -- batimetria sadia em todo o rio "
                         f"(rebaixamento max {reb.max():.0f} m)", fontsize=10)
        ax.set_xlabel("distancia ao longo do rio (km, de montante)")
        ax.set_ylabel("cota (m)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig(a.saida, dpi=110)
    print(f"figura -> {a.saida}")


if __name__ == "__main__":
    main()
