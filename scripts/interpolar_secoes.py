# -*- coding: utf-8 -*-
"""Insere secoes INTERPOLADAS onde o desnivel entre vizinhas e grande.

    python scripts/interpolar_secoes.py taha_ai.g95 --saida h02 \
        --trecho Itajai_Acu,R1,139000,166000 [--trecho ...] \
        --queda-max 2.0

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Nas rampas da serra as secoes estao a ~1 km com 20-30 m de desnivel
  entre vizinhas. Diferencas finitas 1D nao resolvem um degrau desses:
  o solver oscila ali NAO importa o n nem o passo de tempo (testado:
  n=0,10 ja gravado, 30 s piora, LPI mais duro piora). A regra classica
  e secao a cada 100-200 m em trecho ingreme.

O QUE SE FAZ, so DENTRO dos trechos pedidos

  Entre cada par de secoes vizinhas com queda de talvegue maior que
  `--queda-max` m, insere as secoes intermediarias necessarias para que
  nenhum salto passe do limite. Cada interpolada e a MISTURA linear das
  duas vizinhas (perfil reamostrado em 80 estacoes normalizadas, bancos
  e n interpolados), com RS proporcional -- o mesmo que o RAS faz no "XS
  Interpolation", gravado em texto.

  CONFERENCIA relendo o gravado: nenhum salto acima do limite nos
  trechos, contagem de secoes bate com o previsto.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402

N_PONTOS = 160


def reamostrar(sta, z):
    s0, s1 = float(sta[0]), float(sta[-1])
    u = np.linspace(0.0, 1.0, N_PONTOS)
    su = s0 + u * (s1 - s0)
    return u, np.interp(su, sta, z), s0, s1


def misturar(A, B, w):
    """Secao interpolada: w*A + (1-w)*B (w=1 na montante A)."""
    uA, zA, a0, a1 = reamostrar(np.asarray(A["sta"], float),
                                np.asarray(A["z"], float))
    uB, zB, b0, b1 = reamostrar(np.asarray(B["sta"], float),
                                np.asarray(B["z"], float))
    s0 = w * a0 + (1 - w) * b0
    s1 = w * a1 + (1 - w) * b1
    u = np.linspace(0.0, 1.0, N_PONTOS)
    # 2 casas ja aqui: Sta/Elev, #Mann e Bank Sta gravam esta mesma grade
    # (com casas diferentes o RAS acusa "Manning's n value not set")
    sta = np.round(s0 + u * (s1 - s0), 2)
    z = w * zA + (1 - w) * zB
    # a reamostragem perde o fundo do canal estreito: crava o talvegue
    # misturado no ponto mais baixo do perfil misturado
    alvo_tal = w * float(np.asarray(A["z"], float).min()) \
        + (1 - w) * float(np.asarray(B["z"], float).min())
    z[int(np.argmin(z))] = min(float(z.min()), alvo_tal)

    def frac(x, lo, hi):
        return (x - lo) / max(hi - lo, 1e-9)

    lb = s0 + (s1 - s0) * (w * frac(A["lb"], a0, a1)
                           + (1 - w) * frac(B["lb"], b0, b1))
    rb = s0 + (s1 - s0) * (w * frac(A["rb"], a0, a1)
                           + (1 - w) * frac(B["rb"], b0, b1))
    # bancos ancorados em estacoes existentes (exigencia do RAS)
    lb = float(sta[int(np.argmin(np.abs(sta - lb)))])
    rb = float(sta[int(np.argmin(np.abs(sta - rb)))])
    if rb <= lb:
        rb = float(sta[min(int(np.argmin(np.abs(sta - lb))) + 1,
                           N_PONTOS - 1)])
    return sta, z, lb, rb


def bloco_secao(rs, comp, sta, z, lb, rb, mann):
    b = ["Type RM Length L Ch R = 1 ,%-8s,%g,%g,%g"
         % (f"{rs:.1f}", comp, comp, comp),
         "BEGIN DESCRIPTION:", "interpolada (interpolar_secoes.py)",
         "END DESCRIPTION:",
         "#Sta/Elev= %d " % len(sta)]
    v = []
    for a, c in zip(sta, z):
        v += [a, c]
    b += _col(v, 8, 2)
    b.append("#Mann= 3 , 0 , 0 ")
    b += _col([sta[0], mann[0], 0, lb, mann[1], 0, rb, mann[2], 0], 8, 3)
    b.append("Bank Sta=%.2f,%.2f" % (lb, rb))
    b.append("XS HTab Starting El and Incr=%.2f,0.3, 60 "
             % (float(np.min(z)) + 0.15))
    b.append("Exp/Cntr=0.3,0.1")
    return b


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "h02")
    queda = _arg(argv, "--queda-max", 2.0, float)
    trechos = []
    for k, a in enumerate(argv):
        if a == "--trecho":
            rio, reach, r0, r1 = argv[k + 1].split(",")
            trechos.append((rio.strip(), reach.strip(),
                            float(r0), float(r1)))
    todos = "--todos" in argv
    if not trechos and not todos:
        raise SystemExit("faltou --trecho rio,reach,rs_min,rs_max "
                         "(ou --todos)")

    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    if todos:
        vistos = sorted({(d["rio"], d["reach"]) for d in S})
        trechos = [(rio, reach, 0.0, 1e12) for rio, reach in vistos]
    # mann por secao (le do arquivo na passada de escrita; aqui so padrao)
    inserir = {}          # indice da secao de MONTANTE -> [blocos]
    previsto = 0
    for rio, reach, r0, r1 in trechos:
        # so secoes de verdade (tipo 1); estruturas (barragem = tipo 5)
        # nem entram no par nem podem ficar ENTRE um par interpolado
        idx = [i for i, d in enumerate(S)
               if d["rio"] == rio and d["reach"] == reach
               and r0 <= d["rs"] <= r1 and d["tipo"] == "1"]
        estruturas = [d["rs"] for d in S
                      if d["rio"] == rio and d["reach"] == reach
                      and d["tipo"] != "1"]
        idx.sort(key=lambda i: -S[i]["rs"])
        n_ins = 0
        for a, b in zip(idx, idx[1:]):
            A, B = S[a], S[b]
            if any(B["rs"] < rs_e < A["rs"] for rs_e in estruturas):
                continue
            zA = float(np.asarray(A["z"], float).min())
            zB = float(np.asarray(B["z"], float).min())
            salto = zA - zB
            if salto <= queda:
                continue
            n = int(np.ceil(salto / queda)) - 1
            n = min(n, 20)
            comp = (A["rs"] - B["rs"]) / (n + 1)
            blocos = []
            for k in range(1, n + 1):
                w = 1 - k / (n + 1)
                sta, z, lb, rb = misturar(A, B, w)
                rs_k = A["rs"] - k * comp
                blocos.append((rs_k, comp, sta, z, lb, rb))
            # chave por (rio, reach, RS): o indice de ler_secoes NAO
            # acompanha o arquivo quando ha estruturas (sem Sta/Elev
            # elas nem entram na lista) -- foi a dessincronia que
            # reescreveu o cabecalho da barragem
            inserir[(rio, reach, round(A["rs"], 2))] = blocos
            n_ins += n
        print(f"   {rio} {reach} RS {r0:.0f}-{r1:.0f}: {len(idx)} secoes, "
              f"+{n_ins} interpoladas")
        previsto += n_ins

    if not inserir:
        print("nada a interpolar")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    # mann da secao corrente para herdar
    saida, j = [], 0
    rio_c = reach_c = None
    mann_atual = (0.06, 0.045, 0.06)
    pendente = None
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio_c = p[0].strip()
            reach_c = p[1].strip() if len(p) > 1 else ""
        if l.startswith("Type RM Length L Ch R"):
            # antes de abrir a proxima secao, descarrega interpoladas
            if pendente is not None:
                for rs_k, comp, sta, z, lb, rb in pendente:
                    saida += bloco_secao(rs_k, comp, sta, z, lb, rb,
                                         mann_atual) + [""]
                pendente = None
            p = l.split("=", 1)[1].split(",")
            try:
                chave = (rio_c, reach_c, round(float(p[1]), 2))
            except ValueError:
                chave = None
            if chave in inserir and p[0].strip() == "1":
                pendente = inserir.pop(chave)
                # encurta os comprimentos da secao de montante
                comp = pendente[0][1]
                saida.append("Type RM Length L Ch R = 1 ,%s,%g,%g,%g"
                             % (p[1].strip(), comp, comp, comp))
                j += 1
                continue
        if l.startswith("#Mann=") and not pendente:
            try:
                cnt = int(l.split("=")[1].split(",")[0])
                vals = []
                jj = j + 1
                while jj < len(linhas) and len(vals) < 3 * cnt:
                    x = linhas[jj]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    vals += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                             if x[c:c + 8].strip()]
                    jj += 1
                if cnt == 3 and len(vals) == 9:
                    mann_atual = (vals[1], vals[4], vals[7])
            except (ValueError, IndexError):
                pass
        saida.append(l)
        j += 1
    if pendente is not None:
        for rs_k, comp, sta, z, lb, rb in pendente:
            saida += bloco_secao(rs_k, comp, sta, z, lb, rb,
                                 mann_atual) + [""]
    escrever(novo, "\n".join(saida))

    print(f"\nCONFERENCIA (relendo o arquivo gravado)")
    if inserir:
        print(f"   AVISO: {len(inserir)} chaves nao encontradas no "
              f"arquivo: {list(inserir)[:3]}")
    t_novo = open(novo, encoding="latin-1", errors="replace").read()
    t_vel = open(entrada, encoding="latin-1", errors="replace").read()
    n5v = t_vel.count("Type RM Length L Ch R = 5")
    n5n = t_novo.count("Type RM Length L Ch R = 5")
    print(f"   estruturas (Type 5): {n5v} -> {n5n}   (nao pode mudar)")
    B2 = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B2)}   "
          f"(previsto {len(S) + previsto})")
    for rio, reach, r0, r1 in trechos:
        idx = [i for i, d in enumerate(B2)
               if d["rio"] == rio and d["reach"] == reach
               and r0 <= d["rs"] <= r1]
        idx.sort(key=lambda i: -B2[i]["rs"])
        tal = [float(np.asarray(B2[i]["z"], float).min()) for i in idx]
        pior = max((a - b for a, b in zip(tal, tal[1:])), default=0.0)
        print(f"   {rio} {reach}: pior salto {pior:.2f} m "
              f"(limite {queda:.1f})")


if __name__ == "__main__":
    main(sys.argv[1:])
