# -*- coding: utf-8 -*-
"""Por que a batimetria do legado nao serve ao Benedito a MONTANTE.

    python scripts/diagnostico_benedito.py

Ao aplicar a batimetria do legado nos seis rios, cinco pedem rebaixar o leito
alguns metros (a mediana casa com a calha levantada, ~5 a 10 m). O Benedito
pede rebaixar 104 m na mediana e ate 265 m -- fisicamente impossivel para um
rio. Este script mede de onde vem isso e deixa a figura.

O QUE ACONTECE

  O eixo do Benedito em `eixos_do_relevo.geojson` E a linha esquematica do
  reach `Rio_Benedito` do legado (a distancia entre eles e ZERO -- o eixo veio
  de la). Sao 471 vertices em 44 km, ~93 m entre pontos.

  A JUSANTE (RS 0 a ~18 km) o vale e largo e essa linha cai sobre o canal: a
  lamina lida no MDT fica a poucos metros do fundo levantado -- batimetria
  saudavel, como nos outros rios.

  A MONTANTE (RS ~22 a 44 km) o vale e estreito e encaixado. A linha
  esquematica, grossa, corre pela ENCOSTA, nao pelo talvegue. A secao que o
  `rio_do_relevo` corta do MDT ali pega o ponto mais baixo da encosta, dezenas
  a centenas de metros ACIMA do canal real. A diferenca lamina-fundo cresce de
  ~10 m para 265 m.

  Aplicar a batimetria nesse trecho mandaria o `aplicar` DESCER o leito 240 m
  dentro de uma calha que nao tem essa profundidade: escavar um canion, que e
  exatamente o que as regras proibem. Por isso o Benedito fica SEM g02 ate o
  eixo de montante ser refeito seguindo o talvegue do MDT (e nao a linha
  esquematica do legado).

  A saida NAO e inventar cota: e derivar o eixo de montante do proprio MDT
  (acumulacao de fluxo), recortar as secoes nesse eixo, e so entao ancorar no
  levantamento -- que a jusante ja bate.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from qc_secoes import ler_secoes              # noqa: E402
from batimetria_do_legado import secoes_levantadas   # noqa: E402

LEGADO = "legado/Itajai_Rede_1983.g01"
SAIDA = "doc/figuras/benedito_eixo_montante.png"


def perfil(g):
    S = ler_secoes(g)
    S.sort(key=lambda d: -d["rs"])
    rs = np.array([d["rs"] for d in S])
    z = np.array([float(np.asarray(d["z"], float).min()) for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    x = np.r_[0.0, np.cumsum(ch[:-1])]
    return rs, x, z


def main():
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    casos = [("modelo/rio_benedito/rio_benedito.g01", "Rio_Benedito",
              "Rio Benedito -- a montante o eixo esquematico sobe a encosta"),
             ("modelo/itajai_mirim/itajai_mirim.g01", "Itajai_Mirim",
              "Itajai-Mirim (controle) -- lamina acompanha o fundo levantado")]
    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    for ax, (g, rio, titulo) in zip(axes, casos):
        if not os.path.exists(g):
            ax.set_title(f"{rio}: geometria ausente ({g})")
            continue
        rs, x, z = perfil(g)
        L = secoes_levantadas(LEGADO, rio)
        Lrs, Linv = L[:, 0], L[:, 3]
        o = np.argsort(-Lrs)
        Lrs, Linv = Lrs[o], Linv[o]
        xl = np.interp(-Lrs, -rs, x)
        fundo = np.interp(x, xl, Linv)
        ax.plot(x / 1000, z, color="tab:red", lw=1.6,
                label="talvegue lido do MDT (lamina d'agua)")
        ax.plot(xl / 1000, Linv, color="tab:blue", lw=1.6,
                label="fundo levantado 1983 (legado)")
        ax.fill_between(x / 1000, z, fundo, color="0.85", zorder=0)
        dif = z - fundo
        ax.set_title(titulo, fontsize=10)
        ax.set_xlabel("distancia ao longo do rio (km, de montante)")
        ax.set_ylabel("cota (m)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)
        ax.text(0.98, 0.06, f"diferenca maxima {np.nanmax(dif):.0f} m",
                transform=ax.transAxes, ha="right", fontsize=9,
                color="tab:red")
    plt.tight_layout()
    plt.savefig(SAIDA, dpi=110)
    print(f"figura -> {SAIDA}")


if __name__ == "__main__":
    main()
