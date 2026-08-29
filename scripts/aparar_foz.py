# -*- coding: utf-8 -*-
"""Apara a secao de foz que atravessa o eixo do rio RECEPTOR.

    python scripts/aparar_foz.py taha_ai_novo/taha_ai.g02 --saida g03

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  As ultimas secoes de um tributario (RS 75-200, a passos da juncao) sao
  largas o bastante para cruzar o EIXO do rio que as recebe -- e o Validate
  Geometry acusa "XS must intersect exactly one Reach", fatal. Medido no g02
  do taha_ai_novo: 6 secoes, todas de foz (Oeste R4 75, Sul 75, Mirim 172.88
  e 75, Taio 75, Testo 75).

POR QUE O corrigir_cutlines NAO RESOLVE

  La o canal (lb..rb) e INVIOLAVEL -- correto para meandro, onde a travessia
  indesejada cai na planicie. Na foz nao: o cruzamento com o eixo vizinho cai
  DENTRO da faixa marcada como canal, porque a "calha" dessas secoes foi
  alargada ate cobrir a varzea do rio receptor. Aquele chao pertence ao
  OUTRO rio; manter a secao ali conta a mesma agua duas vezes.

O QUE SE FAZ

  So nas secoes cuja cutline cruza eixo de outro reach: corta-se a ponta que
  contem esse cruzamento, com RECUO m de folga antes dele -- escolhendo a
  ponta que NAO contem o cruzamento com o proprio eixo, que e obrigatorio.
  Se a margem (Bank Sta) ficar fora da faixa restante, ela e trazida para o
  novo extremo; as quebras do #Mann acompanham, como no corrigir_cutlines.

  Nenhuma cota entre as estacas que sobrevivem muda; o talvegue nao muda,
  porque o cruzamento com o proprio eixo -- onde o canal de verdade esta --
  fica sempre do lado preservado.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from qc_geometria import ler_eixos                     # noqa: E402
from corrigir_cutlines import (mapa_reaches, travessias, _col, _fmt,
                               _arg, TOL)              # noqa: E402
from ras_io import escrever                            # noqa: E402

FOLGA = 10.0    # m alem da meia-largura do canal receptor


def faixa_foz(d, eixo, outros):
    """(s0, s1) preservados, ou None se a secao nao cruza eixo alheio.

    `outros` e uma lista de (eixo, meia_largura): o recuo antes do eixo
    alheio e a MEIA-LARGURA DO CANAL DELE mais FOLGA -- cortar a 3 m do
    eixo receptor deixaria a ponta DENTRO do canal do outro rio, terminando
    n'agua (foi medido: Mirim RS 75 e Sul RS 75 viraram ponta n'agua com
    recuo fixo de 3 m).
    """
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    L = float(np.hypot(*(B - A)))
    u = (B - A) / max(L, 1e-9)
    t_out = sorted((x, meia) for e, meia in outros
                   for x in travessias(A, u, 0.0, L, e))
    if not t_out:
        return None
    t_prop = travessias(A, u, 0.0, L, eixo)
    ancora = t_prop[0] if t_prop else 0.5 * (d["lb"] + d["rb"])
    s0, s1 = 0.0, L
    for x, meia in t_out:
        rec = meia + FOLGA
        if x < ancora:
            s0 = max(s0, x + rec)
        else:
            s1 = min(s1, x - rec)
    return s0, s1


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g03")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    # meia-largura de canal por reach (mediana), para o recuo na foz
    largs = {}
    for d, ch in zip(S, mapa):
        largs.setdefault(ch, []).append(float(d["rb"]) - float(d["lb"]))
    meia = {ch: 0.5 * float(np.median(v)) for ch, v in largs.items()}

    novos = {}
    for i, d in enumerate(S):
        ch = mapa[i]
        outros = [(e, meia.get(k, 30.0)) for k, e in eixos.items()
                  if k != ch]
        r = faixa_foz(d, eixos[ch], outros)
        if r is None:
            continue
        s0, s1 = r
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        L = float(np.hypot(*(B - A)))
        u = (B - A) / max(L, 1e-9)
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= s0 - TOL) & (st <= s1 + TOL)
        ns, nz = list(st[m]), list(z[m])
        if not ns or ns[0] > s0 + TOL:
            ns.insert(0, s0)
            nz.insert(0, float(np.interp(s0, st, z)))
        if ns[-1] < s1 - TOL:
            ns.append(s1)
            nz.append(float(np.interp(s1, st, z)))
        ns = np.array(ns) - s0
        nz = np.array(nz)
        # a comparacao e no VALOR JA ARREDONDADO a 2 casas (o formato do
        # arquivo): duas estacas a 0,006 m viram duplicata depois do %8.2f
        keep = [0]
        for k in range(1, len(ns)):
            if round(ns[k], 2) > round(ns[keep[-1]], 2):
                keep.append(k)
        ns, nz = ns[keep], nz[keep]
        # margens: entram na faixa se ficaram fora
        lb = float(d["lb"]) - s0
        rb = float(d["rb"]) - s0
        lb = min(max(lb, float(ns[0])), float(ns[-1]) - 0.01)
        rb = min(max(rb, lb + 0.01), float(ns[-1]))
        novos[i] = {"sta": ns, "z": nz, "lb": lb, "rb": rb, "desl": s0,
                    "cut": (A + s0 * u, A + s1 * u)}
        print(f"   {ch[0]:14s} {ch[1]:3s} RS {d['rs']:9.2f}  "
              f"largura {L:6.0f} -> {s1-s0:6.0f} m   "
              f"canal {d['rb']-d['lb']:5.0f} -> {rb-lb:5.0f} m")

    if not novos:
        print("nenhuma secao cruza eixo alheio -- nada a fazer")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
            if l.startswith("XS GIS Cut Line"):
                saida.append("XS GIS Cut Line= 2")
                saida.append("".join("%16.2f" % x for x in
                                     (nv["cut"][0][0], nv["cut"][0][1],
                                      nv["cut"][1][0], nv["cut"][1][1])))
                j += 1
                while j < len(linhas) and linhas[j].strip() and \
                        linhas[j][:1] in " -0123456789":
                    j += 1
                continue
            if l.startswith("#Sta/Elev"):
                v = []
                for a, b in zip(nv["sta"], nv["z"]):
                    v += [a, b]
                saida.append("#Sta/Elev= %d " % len(nv["sta"]))
                saida += _col(v, 8, 2)
                cnt = int(l.split("=")[1])
                j += 1
                lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(x), 8)
                                  if x[c:c + 8].strip()])
                    j += 1
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%s,%s"
                             % (_fmt(nv["lb"]), _fmt(nv["rb"])))
                j += 1
                continue
            if l.startswith("#Mann="):
                cnt = int(l.split("=")[1].split(",")[0])
                bruto, k2 = [], j + 1
                while k2 < len(linhas) and len(bruto) < 3 * cnt:
                    x = linhas[k2]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                              if x[c:c + 8].strip()]
                    k2 += 1
                val = [float(x) for x in bruto[:3 * cnt]]
                topo = float(nv["sta"][-1])
                for t in range(0, 3 * cnt, 3):
                    val[t] = min(max(val[t] - nv["desl"], 0.0), topo)
                saida.append(l)
                lin, corpo = "", []
                for t, x in enumerate(val):
                    lin += ("%8.2f" % x if t % 3 == 0 else
                            "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                    if (t + 1) % 9 == 0:
                        corpo.append(lin)
                        lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                j = k2
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    # -------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B2 = ler_secoes(novo)
    eixos2 = ler_eixos(novo)
    mapa2 = mapa_reaches(novo)
    resto = 0
    for i, d in enumerate(B2):
        A = np.asarray(d["cut"][0], float)
        Bp = np.asarray(d["cut"][-1], float)
        L = float(np.hypot(*(Bp - A)))
        u = (Bp - A) / max(L, 1e-9)
        n_out = sum(len(travessias(A, u, 0.0, L, e))
                    for k, e in eixos2.items() if k != mapa2[i])
        if n_out:
            resto += 1
            print(f"   AINDA cruza eixo alheio: {mapa2[i]} RS {d['rs']:.2f}")
    print(f"   secoes cruzando eixo alheio: {resto}   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
