# -*- coding: utf-8 -*-
"""Le TODOS os perfis de uma geometria e mede o que nao pode estar errado.

    python scripts/qc_perfis.py modelo/itajai_acu/itajai_acu.g01
    python scripts/qc_perfis.py modelo/*/*.g02

Existe porque os erros simples vinham sendo descobertos UM POR UM, por
inspecao visual do usuario no RAS Mapper -- a secao que termina dentro
d'agua, a que nao alcanca o outro lado do rio, o vao de 1500 m -- quando uma
leitura de todos os perfis pegaria tudo de uma vez, antes do solver. Roda em
segundos, sem HEC-RAS; entra no `construir_rio.py` como porteiro.

O QUE MEDE, POR SECAO

  barranco      a folga de cada PONTA acima do talvegue. Ponta a menos de
                1 m e secao terminando dentro d'agua ("sem barranco para
                segurar", como o usuario viu): agua vaza na primeira cheia.
                Entre 1 e `--folga` m e aviso -- contem pouca cheia.
  calha contida a secao tem de conter as proprias margens (lb/rb dentro do
                perfil, rb > lb) e a calha nao pode ser a secao inteira.
  vaos de RS    espacamento entre secoes vizinhas acima de `--vao` m.
  pontos        mais de 500 pontos o HEC-RAS recusa.
  rs repetido   duas secoes na mesma estaca.

Saida: uma linha por rio com os totais, os piores casos com RS, e codigo de
saida 1 se houver GRAVE (ponta n'agua, calha fora, rs repetido, > 500 pontos).
"""
import argparse
import glob
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from qc_secoes import ler_secoes    # noqa: E402


def qc(caminho, folga=4.0, vao=800.0):
    S = ler_secoes(caminho)
    S.sort(key=lambda d: -d["rs"])
    graves, avisos = [], []
    rs_vistos = {}
    for d in S:
        rs = d["rs"]
        z = np.asarray(d["z"], float)
        st = np.asarray(d["sta"], float)
        zt = float(z.min())
        f_esq = float(z[0] - zt)
        f_dir = float(z[-1] - zt)
        pior = min(f_esq, f_dir)
        if pior <= 1.0:
            graves.append((rs, f"ponta n'agua ({pior:.1f} m acima do "
                               f"talvegue)"))
        elif pior < folga:
            avisos.append((rs, f"barranco baixo ({pior:.1f} m)"))
        if not (st[0] - 1e-6 <= d["lb"] < d["rb"] <= st[-1] + 1e-6):
            graves.append((rs, f"margens fora do perfil (lb {d['lb']:.0f} "
                               f"rb {d['rb']:.0f} perfil {st[0]:.0f}-"
                               f"{st[-1]:.0f})"))
        if len(st) > 500:
            graves.append((rs, f"{len(st)} pontos (limite 500)"))
        r2 = round(rs, 2)
        if r2 in rs_vistos:
            graves.append((rs, "RS repetida"))
        rs_vistos[r2] = True
    rs_arr = np.array([d["rs"] for d in S])
    dif = -np.diff(rs_arr)
    for i in np.flatnonzero(dif > vao):
        avisos.append((rs_arr[i], f"vao de {dif[i]:.0f} m ate a proxima"))
    return S, graves, avisos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geoms", nargs="+")
    ap.add_argument("--folga", type=float, default=4.0,
                    help="m de barranco abaixo disso e aviso")
    ap.add_argument("--vao", type=float, default=800.0,
                    help="m entre secoes acima disso e aviso")
    ap.add_argument("--detalhe", type=int, default=6,
                    help="quantos piores casos listar por rio")
    ap.add_argument("--rigido", action="store_true",
                    help="codigo de saida 1 se houver GRAVE (para CI)")
    a = ap.parse_args()
    caminhos = []
    for g in a.geoms:
        caminhos += glob.glob(g) or [g]

    total_graves = 0
    for c in caminhos:
        if not os.path.exists(c):
            print(f"{c}: NAO EXISTE")
            total_graves += 1
            continue
        S, graves, avisos = qc(c, a.folga, a.vao)
        total_graves += len(graves)
        print(f"{c}: {len(S)} secoes   GRAVES {len(graves)}   "
              f"avisos {len(avisos)}")
        for rs, m in sorted(graves)[:a.detalhe]:
            print(f"   GRAVE  RS {rs:10.1f}  {m}")
        if len(graves) > a.detalhe:
            print(f"   ... e mais {len(graves)-a.detalhe} grave(s)")
        for rs, m in sorted(avisos)[:max(a.detalhe - len(graves), 2)]:
            print(f"   aviso  RS {rs:10.1f}  {m}")
    return 1 if (total_graves and a.rigido) else 0


if __name__ == "__main__":
    sys.exit(main())
