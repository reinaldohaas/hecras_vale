# -*- coding: utf-8 -*-
"""Poe rugosidade de MONTANHA onde o leito e ingreme, secao a secao.

    python scripts/rugosidade_da_serra.py taha_ai.g95 --saida h01

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  O gerador escreveu n=0,035 de canal em toda parte -- ate nas rampas de
  3% da serra, que na realidade sao cascatas de matacao com n=0,06-0,10
  (Jarrett 1984; Chow). Com 0,035 a agua dispara, Froude passa de 1 e o
  solver 1D oscila exatamente nesses trechos (Acu R1 km 139-166, Norte
  km 10-30...). Corrigir o n e fisica, nao estabilizacao artificial.

O QUE SE FAZ

  Declividade local do talvegue (mediana entre a secao e as vizinhas):

     > 0,5%  ->  n_canal >= 0,045
     > 1%    ->  n_canal >= 0,060
     > 2%    ->  n_canal >= 0,080
     > 4%    ->  n_canal >= 0,100

  So se AUMENTA o n (planicie intacta); margens ganham +0,01 sobre o
  canal onde o canal subiu. So mexe em #Mann de 3 quebras (o padrao do
  gerador). CONFERENCIA relendo o gravado.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402

DEGRAUS = [(0.04, 0.100), (0.02, 0.080), (0.01, 0.060), (0.005, 0.045)]


def n_da_declividade(decl):
    for lim, n in DEGRAUS:
        if decl > lim:
            return n
    return None


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "h01")
    # --forcar rio,reach,rs0,rs1,n : n minimo explicito numa janela
    # (ex.: Salto dos Piloes, corredeira de matacoes, n 0,12)
    forcados = []
    for k, a in enumerate(argv):
        if a == "--forcar":
            rio_f, reach_f, r0, r1, nf = argv[k + 1].split(",")
            forcados.append((rio_f.strip(), reach_f.strip(),
                             float(r0), float(r1), float(nf)))
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    # declividade local por secao (mediana das diferencas com vizinhas)
    por_reach = {}
    for i, d in enumerate(S):
        por_reach.setdefault((d["rio"], d["reach"]), []).append(i)
    decl = {}
    for chave, idx in por_reach.items():
        idx.sort(key=lambda i: -S[i]["rs"])
        tal = [float(np.asarray(S[i]["z"], float).min()) for i in idx]
        rs = [S[i]["rs"] for i in idx]
        for j, i in enumerate(idx):
            ds = []
            if j > 0:
                ds.append((tal[j - 1] - tal[j]) / max(rs[j - 1] - rs[j], 1))
            if j + 1 < len(idx):
                ds.append((tal[j] - tal[j + 1]) / max(rs[j] - rs[j + 1], 1))
            decl[i] = float(np.median(ds)) if ds else 0.0

    novos = {}
    por_rio = {}
    for i, d in enumerate(S):
        alvo = n_da_declividade(max(decl.get(i, 0.0), 0.0))
        for rio_f, reach_f, r0, r1, nf in forcados:
            if d["rio"] == rio_f and d["reach"] == reach_f \
                    and r0 <= d["rs"] <= r1:
                alvo = max(alvo or 0.0, nf)
        if alvo is not None:
            # chave por (rio, reach, RS): indice de ler_secoes NAO segue
            # o arquivo quando ha estrutura (barragem) no meio
            novos[(d["rio"], d["reach"], round(d["rs"], 2))] = alvo
            e = por_rio.setdefault(d["rio"], [0, 0.0])
            e[0] += 1
            e[1] = max(e[1], alvo)

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")
    print(f"{'rio':16s} {'secoes':>6s} {'n max':>6s}")
    for rio in sorted(por_rio):
        n, nmax = por_rio[rio]
        print(f"{rio:16s} {n:6d} {nmax:6.3f}")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    saida, j, mudadas = [], 0, 0
    rio_c = reach_c = None
    chave = None
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio_c = p[0].strip()
            reach_c = p[1].strip() if len(p) > 1 else ""
        if l.startswith("Type RM Length L Ch R"):
            p = l.split("=", 1)[1].split(",")
            try:
                chave = (rio_c, reach_c, round(float(p[1]), 2))
            except ValueError:
                chave = None
        alvo = novos.get(chave)
        if alvo is not None and l.startswith("#Mann="):
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
                if vals[4] < alvo:
                    vals[4] = alvo
                    vals[1] = max(vals[1], alvo + 0.01)
                    vals[7] = max(vals[7], alvo + 0.01)
                    mudadas += 1
            saida.append("#Mann= %d , %s , 0 "
                         % (cnt, partes[1].strip() if len(partes) > 1
                            else "0"))
            saida += _col(vals, 8, 3)
            continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    print(f"\nsecoes com n aumentado: {mudadas}")
    print("CONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    t = open(novo, encoding="latin-1", errors="replace").read()
    for n in ("0.100", "0.080", "0.060"):
        print(f"   ocorrencias de n={n}: {t.count(n)}")


if __name__ == "__main__":
    main(sys.argv[1:])
