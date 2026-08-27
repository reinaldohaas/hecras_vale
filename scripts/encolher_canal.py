# -*- coding: utf-8 -*-
"""Encolhe o canal ESCAVADO para a largura MEDIDA no SIG-SC, secao a secao.

    python scripts/encolher_canal.py taha_ai.g01 --saida g95 \
        --medidas doc/larguras_sigsc --fator 2.0

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  O gerador escavou canais por formula de regime (w = 5*A^0.4): mediana de
  115 m no Benedito onde a lamina medida no MDT 1 m e 16 m. Vazao baixa
  espalhada num prisma 7x largo vira lamina de centimetros -- exatamente
  onde o solver 1D oscila.

O QUE SE FAZ

  Largura-alvo por secao = lamina do SIG-SC interpolada na quilometragem
  da secao (mediana movel de `--janela-km`), vezes `--fator` (lamina e o
  PISO; 2x aproxima margens plenas), nunca abaixo de `--minimo`.

  So age se o canal atual for mais de `--somente-acima` vezes o alvo.
  O canal novo e centrado no talvegue da secao:

    - dentro do canal novo   : cotas intactas (talvegue e profundidade ficam)
    - resto do prisma antigo : ATERRADO de volta a linha do terreno
                               (reta entre as cotas das margens antigas;
                               so se levanta, nunca se cava)
    - pontos novos nas duas margens novas, `Bank Sta` e #Mann acompanham

  A secao total (planicie) nao muda de extensao. VETO respeitado: a foz do
  Itajai Mirim (`--poupar Itajai_Mirim:10`, km da foz) nao e tocada.

  CONFERENCIA relendo o gravado: mediana de largura por rio antes/depois,
  talvegues identicos, nenhum ponto rebaixado.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def serie_medida(pasta, rio, janela_km, cada_km=0.5):
    """(dist_km ordenado, lamina suavizada) ou None."""
    arq = os.path.join(pasta, f"{rio}.csv")
    if not os.path.exists(arq):
        return None
    pares = []
    for r in csv.reader(open(arq, encoding="utf-8"), delimiter=";"):
        if r[0] == "dist_foz_km" or not r[1]:
            continue
        pares.append((float(r[0]), float(r[1])))
    if len(pares) < 3:
        return None
    pares.sort()
    d = np.array([p[0] for p in pares])
    w = np.array([p[1] for p in pares])
    meia = max(1, int(round(janela_km / cada_km / 2)))
    suave = np.array([np.median(w[max(0, i - meia):i + meia + 1])
                      for i in range(len(w))])
    return d, suave


def encolher(sta, z, lb, rb, alvo):
    """Nova (sta, z, nb_l, nb_r) com canal de largura `alvo` no talvegue."""
    dentro = (sta >= lb - 1e-6) & (sta <= rb + 1e-6)
    i_tal = np.flatnonzero(dentro)[int(np.argmin(z[dentro]))]
    c = sta[i_tal]
    nb_l = max(lb, c - alvo / 2.0)
    nb_r = min(rb, c + alvo / 2.0)
    # se bateu numa borda, completa a largura pela outra
    falta = alvo - (nb_r - nb_l)
    if falta > 0:
        if nb_l <= lb + 1e-6:
            nb_r = min(rb, nb_r + falta)
        else:
            nb_l = max(lb, nb_l - falta)

    # Bank Sta TEM de coincidir com uma estacao gravada (exigencia do
    # RAS): ancora na estacao vizinha se ela nao deformar o canal em
    # mais de 20%; senao insere ponto novo ja arredondado como a
    # gravacao (2 casas) para o casamento ser exato
    def ancorar(xb):
        k = int(np.argmin(np.abs(sta - xb)))
        if abs(float(sta[k]) - xb) <= 0.1 * alvo:
            return float(sta[k]), False
        return round(xb, 2), True

    nb_l, ins_l = ancorar(nb_l)
    nb_r, ins_r = ancorar(nb_r)

    z_lb = float(np.interp(lb, sta, z))
    z_rb = float(np.interp(rb, sta, z))

    def chao(x):
        return z_lb + (z_rb - z_lb) * (x - lb) / max(rb - lb, 1e-9)

    st2, z2 = list(sta), list(z)
    for xb, ins in ((nb_l, ins_l), (nb_r, ins_r)):
        if ins:
            k = int(np.searchsorted(st2, xb))
            st2.insert(k, xb)
            z2.insert(k, float(np.interp(xb, sta, z)))
    st2 = np.asarray(st2)
    z2 = np.asarray(z2)
    aterro = (st2 >= lb - 1e-6) & (st2 <= rb + 1e-6) \
        & ((st2 < nb_l - 1e-6) | (st2 > nb_r + 1e-6))
    z2[aterro] = np.maximum(z2[aterro], [chao(x) for x in st2[aterro]])
    return st2, z2, nb_l, nb_r


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g95")
    pasta = _arg(argv, "--medidas", "doc/larguras_sigsc")
    fator = _arg(argv, "--fator", 2.0, float)
    minimo = _arg(argv, "--minimo", 15.0, float)
    somente = _arg(argv, "--somente-acima", 1.3, float)
    janela = _arg(argv, "--janela-km", 5.0, float)
    poupar = _arg(argv, "--poupar", "Itajai_Mirim:10")

    poupados = {}
    for p in poupar.split(","):
        if ":" in p:
            r, km = p.split(":")
            poupados[r.strip()] = float(km)

    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    rios = sorted({d["rio"] for d in S})
    medidas = {r: serie_medida(pasta, r, janela) for r in rios}
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"alvo   : lamina SIG-SC (mediana movel {janela:.0f} km) "
          f"x {fator}, piso {minimo:.0f} m\n")

    novos = {}
    por_rio = {}
    for i, d in enumerate(S):
        rio = d["rio"]
        m = medidas.get(rio)
        if m is None:
            continue
        km = d["rs"] / 1000.0
        if rio in poupados and km <= poupados[rio]:
            continue
        sta = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        lb, rb = float(d["lb"]), float(d["rb"])
        atual = rb - lb
        alvo = max(minimo, float(np.interp(km, m[0], m[1])) * fator)
        alvo = min(alvo, atual)
        if atual < somente * alvo:
            continue
        st2, z2, nb_l, nb_r = encolher(sta, z, lb, rb, alvo)
        novos[i] = {"sta": st2, "z": z2, "lb": nb_l, "rb": nb_r}
        e = por_rio.setdefault(rio, [0, [], []])
        e[0] += 1
        e[1].append(atual)
        e[2].append(nb_r - nb_l)

    print(f"{'rio':16s} {'secoes':>6s} {'antes med':>9s} {'depois med':>10s}")
    for rio in sorted(por_rio):
        n, a, b = por_rio[rio]
        print(f"{rio:16s} {n:6d} {np.median(a):8.0f}m {np.median(b):9.0f}m")
    if not novos:
        print("nada a encolher")
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
            if l.startswith("#Mann="):
                partes = l.split("=")[1].split(",")
                cnt = int(partes[0])
                j += 1
                vals = []
                while j < len(linhas) and len(vals) < 3 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    vals += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                             if x[c:c + 8].strip()]
                    j += 1
                if cnt == 3 and len(vals) == 9:
                    vals[3] = nv["lb"]
                    vals[6] = nv["rb"]
                saida.append("#Mann= %d , %s , 0 " %
                             (cnt, partes[1].strip() if len(partes) > 1
                              else "0"))
                saida += _col(vals, 8, 3)
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%.2f,%.2f" % (nv["lb"], nv["rb"]))
                j += 1
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    ok_tal, ok_sobe, ok_larg, ok_anc = 0, 0, 0, 0
    for i, (a, b) in enumerate(zip(S, B)):
        if i not in novos:
            continue
        za = np.asarray(a["z"], float)
        zb = np.asarray(b["z"], float)
        if abs(za.min() - zb.min()) < 1e-3:
            ok_tal += 1
        sa = np.asarray(a["sta"], float)
        sb = np.asarray(b["sta"], float)
        # 5 cm: a gravacao tem 2 casas e em barranco ingreme o
        # arredondamento da ESTACAO desloca a cota interpolada uns mm
        if np.all(zb - np.interp(sb, sa, za) > -0.05):
            ok_sobe += 1
        if b["rb"] - b["lb"] <= a["rb"] - a["lb"] + 1e-3:
            ok_larg += 1
        if np.min(np.abs(sb - b["lb"])) < 1e-6 \
                and np.min(np.abs(sb - b["rb"])) < 1e-6:
            ok_anc += 1
    n = len(novos)
    print(f"   talvegue intacto : {ok_tal}/{n}   (tem de ser {n})")
    print(f"   nada rebaixado   : {ok_sobe}/{n}   (tem de ser {n})")
    print(f"   canal nao cresceu: {ok_larg}/{n}   (tem de ser {n})")
    print(f"   banco em estacao : {ok_anc}/{n}   (tem de ser {n} -- "
          f"exigencia do RAS)")


if __name__ == "__main__":
    main(sys.argv[1:])
