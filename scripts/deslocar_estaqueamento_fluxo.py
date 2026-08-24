# -*- coding: utf-8 -*-
"""Acompanha no .uNN o deslocamento de estaqueamento feito na geometria.

    python scripts/deslocar_estaqueamento_fluxo.py \
        modelo/mirim_novo/mirim_novo.u01 --desloc 999.83 --saida u02

A ENTRADA NAO E TOCADA. Sai um .uNN novo.

POR QUE

  Estender o reach ate a foz renumerou todo o estaqueamento (`RS = distancia
  ate a foz`, com a foz em zero). As condicoes de contorno sao endereçadas por
  RS, entao todas precisam andar a mesma distancia -- menos uma.

  O CONTORNO DE JUSANTE NAO ANDA. Ele e o hidrograma de mare do Itajai-Acu, e
  estava na secao que era a ultima; agora existe secao abaixo dela. O lugar
  dele e a foz, que passou a ser RS 0 -- ou seja, ele fica no mesmo NUMERO e
  muda de SECAO. Somar o deslocamento nele o deixaria 1 km rio acima, que e
  exatamente o defeito que a extensao veio corrigir.

  Isto e explicito e nao automatico: quem chama diz qual contorno e o de
  jusante, por `--jusante`, e a troca sai impressa.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever    # noqa: E402


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    desl = _arg(argv, "--desloc", None, float)
    jus = _arg(argv, "--jusante", "Stage Hydrograph")
    ext = _arg(argv, "--saida", "u02")
    if desl is None:
        raise SystemExit("informe --desloc com o deslocamento em metros")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    t = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"deslocamento: +{desl:.2f} m   contorno de jusante: '{jus}'\n")

    # quais blocos sao o de jusante
    e_jus = [False] * len(t)
    ini = [i for i, l in enumerate(t) if l.startswith("Boundary Location=")]
    for a, b in zip(ini, ini[1:] + [len(t)]):
        if any(jus in x for x in t[a:b]):
            e_jus[a] = True

    saida, n_bc, n_ic = [], 0, 0
    for i, l in enumerate(t):
        if l.startswith("Initial RS="):
            p = l.split("=", 1)[1].split(",")
            v = float(p[2])
            p[2] = f"{v + desl:.2f}"
            saida.append("Initial RS=" + ",".join(p))
            n_ic += 1
            print(f"   Initial RS  {v:10.2f} -> {v+desl:10.2f}")
            continue
        if l.startswith("Boundary Location="):
            p = l.split("=", 1)[1].split(",")
            a0 = p[2].strip()
            a1 = p[3].strip()
            if e_jus[i]:
                novo_a0 = 0.0
                print(f"   {jus:<26} {float(a0):10.2f} -> "
                      f"{novo_a0:10.2f}   <- NAO desloca: vai para a foz")
            else:
                novo_a0 = float(a0) + desl
                print(f"   contorno    {float(a0):10.2f} -> {novo_a0:10.2f}"
                      + (f"   ate {float(a1)+desl:.2f}" if a1 else ""))
            p[2] = f"{novo_a0:<8.2f}"
            if a1:
                p[3] = f"{float(a1)+desl:<8.2f}"
            saida.append("Boundary Location=" + ",".join(p))
            n_bc += 1
            continue
        saida.append(l)
    escrever(novo, "\n".join(saida))

    print("\nCONFERENCIA")
    t2 = open(novo, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    print(f"   contornos : {n_bc}   condicoes iniciais: {n_ic}")
    for chave in ("Flow Hydrograph=", "Uniform Lateral Inflow Hydrograph=",
                  "Stage Hydrograph="):
        a = open(entrada, encoding="latin-1",
                 errors="replace").read().count(chave)
        print(f"   {chave:<36} {a} -> {t2.count(chave)}   "
              f"{'ok' if a == t2.count(chave) else 'DIVERGIU'}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
