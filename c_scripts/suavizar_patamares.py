# -*- coding: utf-8 -*-
"""Troca os patamares planos do leito por rampa entre as cotas que existem.

    python scripts/suavizar_patamares.py modelo/mirim_t30/mirim_t30.g19 \
        --saida g24

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O QUE E UM PATAMAR, E POR QUE ELE NAO E RIO

  Medido no `so_mirim.g01` original: 509 dos 1.417 pares de secoes vizinhas
  (36%) tem cota de leito IDENTICA, e os patamares sao longos -- 38 secoes em
  exatamente 126,00 m ao longo de 875 m, 38 em 100,50 m ao longo de 1.426 m,
  32 em 162,99 m, 31 em 150,50 m. As cotas sao redondas.

  Isso e assinatura de perfil extraido de MDE grosseiro, com o fundo do rio
  achatado no degrau do dado. Rio natural nao tem 1,4 km de leito na mesma
  cota ao centimetro. E os espelhos entre patamares viram declividade local de
  ate 5,32% em 25 m, que e o que trava o solver: os quatro estouros da rodada
  (RS 118252, 118399, 119284, 121628) estao todos na borda de um patamar.

NAO SE INVENTA COTA

  A rampa vai da cota da secao IMEDIATAMENTE ANTES do patamar ate a da
  IMEDIATAMENTE DEPOIS -- dois valores que ja estavam no arquivo -- e e
  distribuida pela distancia real entre elas (`Length Ch`). Nenhuma cota nova
  aparece fora desse intervalo, e patamar sem vizinha dos dois lados (no comeco
  ou no fim do reach) NAO e tocado.

COMO A COTA E APLICADA AO PERFIL

  Mover so o ponto mais fundo deixaria um bico. O deslocamento e aplicado a
  TODOS os pontos entre as margens, com peso proporcional a profundidade:
  o ponto mais fundo recebe o deslocamento inteiro e os pontos na cota das
  margens recebem zero. Assim a secao desce ou sobe o fundo sem descolar da
  margem, e a forma da calha e preservada.

  O `XS HTab Starting El and Incr` acompanha, 2 cm acima do novo talvegue --
  a mesma convencao do resto do modelo.

FORA DO CANAL NADA MUDA. A planicie, as margens e a largura ficam iguais.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from ras_io import escrever            # noqa: E402

TOL_IGUAL = 0.005    # m; abaixo disto duas cotas sao "a mesma"
MIN_PATAMAR = 3      # secoes seguidas para valer suavizar
FOLGA_HTAB = 0.02


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def _col(v, larg, dec):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%*.*f" % (larg, dec, x)
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g24")
    minp = int(_arg(argv, "--minimo", MIN_PATAMAR, float))
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S_arq = ler_secoes(entrada)
    ordem = sorted(range(len(S_arq)), key=lambda i: -S_arq[i]["rs"])
    S = [S_arq[i] for i in ordem]
    z = np.array([float(d["z"].min()) for d in S])
    rs = np.array([d["rs"] for d in S])
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(S)}")

    # ---- acha os patamares
    igual = np.abs(np.diff(z)) < TOL_IGUAL
    patamares, i = [], 0
    while i < len(igual):
        if igual[i]:
            j = i
            while j < len(igual) and igual[j]:
                j += 1
            if (j - i + 1) >= minp:
                patamares.append((i, j))       # secoes i..j na mesma cota
            i = j
        else:
            i += 1
    print(f"\npatamares com {minp}+ secoes: {len(patamares)}")
    n_sec = sum(b - a + 1 for a, b in patamares)
    print(f"   secoes neles: {n_sec}  ({100*n_sec/len(S):.0f}% do modelo)")
    if patamares:
        comp = [rs[a] - rs[b] for a, b in patamares]
        print(f"   extensao: mediana {np.median(comp):.0f} m   "
              f"max {max(comp):.0f} m")

    # ---- a rampa, ancorada SO em quem nao e patamar
    # Rampar cada patamar contra as suas vizinhas ORIGINAIS nao serve: a
    # vizinha pode ser o extremo de OUTRO patamar, que esta sendo rampado ao
    # mesmo tempo, e entre dois patamares rampados de forma independente a
    # monotonicidade se perde. Medido assim: 3 contradeclives criados, um
    # deles de 0,85 m. A interpolacao aqui e UMA SO, sobre o rio inteiro,
    # usando como ancora apenas as secoes que nao pertencem a patamar algum --
    # e essas nao mudam de cota.
    e_patamar = np.zeros(len(z), bool)
    for a, b in patamares:
        e_patamar[a:b + 1] = True
    tocadas, sem_borda = {}, 0
    if e_patamar[0] or e_patamar[-1]:
        # os extremos nao tem ancora de um dos lados: ficam como estao
        i = 0
        while i < len(z) and e_patamar[i]:
            e_patamar[i] = False
            i += 1
            sem_borda += 1
        i = len(z) - 1
        while i >= 0 and e_patamar[i]:
            e_patamar[i] = False
            i -= 1
            sem_borda += 1
    anc = ~e_patamar
    novo_z = z.copy()
    x = rs.max() - rs                      # cresce para jusante
    novo_z[e_patamar] = np.interp(x[e_patamar], x[anc], z[anc])
    # o arquivo grava cota em %8.2f; depois do arredondamento um par pode
    # subir 1 cm. So os pontos rampados sao rebaixados para manter o perfil
    # nao-crescente -- ancora nenhuma e tocada.
    novo_z = np.round(novo_z, 2)
    for k in range(1, len(novo_z)):
        if e_patamar[k] and novo_z[k] > novo_z[k - 1]:
            novo_z[k] = novo_z[k - 1]
    mud = np.abs(novo_z - z)
    for k in np.flatnonzero(mud > 1e-9):
        tocadas[ordem[int(k)]] = float(novo_z[k] - z[k])
    print(f"   sem vizinha dos dois lados (nao tocados): {sem_borda}")
    print(f"\nsecoes com leito ajustado: {len(tocadas)}")
    if tocadas:
        v = np.array(list(tocadas.values()))
        print(f"   ajuste: mediana {np.median(np.abs(v)):.3f} m   "
              f"p90 {np.percentile(np.abs(v),90):.3f}   "
              f"max {np.abs(v).max():.3f} m")
        print(f"   desceram {int((v<0).sum())}   subiram {int((v>0).sum())}")

    # ---- reescreve
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j, n_htab = -1, [], 0, 0
    zmin_novo = {}
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        dz = tocadas.get(i_sec)
        if dz is not None:
            d = S_arq[i_sec]
            if l.startswith("#Sta/Elev"):
                st = np.asarray(d["sta"], float)
                zz = np.asarray(d["z"], float).copy()
                lb, rb = float(d["lb"]), float(d["rb"])
                dentro = (st >= lb - 1e-9) & (st <= rb + 1e-9)
                if dentro.any():
                    zc = zz[dentro]
                    z_marg = zc.max()
                    prof = z_marg - zc
                    p = prof.max()
                    peso = prof / p if p > 1e-9 else np.zeros_like(prof)
                    zz[dentro] = zc + dz * peso
                zmin_novo[i_sec] = float(zz.min())
                v = []
                for aa, bb in zip(st, zz):
                    v += [aa, bb]
                saida.append("#Sta/Elev= %d " % len(st))
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
            if l.startswith("XS HTab Starting El and Incr") \
                    and i_sec in zmin_novo:
                q = [x.strip() for x in l.split("=", 1)[1].split(",")]
                saida.append("XS HTab Starting El and Incr="
                             f"{zmin_novo[i_sec]+FOLGA_HTAB:.2f},"
                             f"{float(q[1]):.3f}, {int(q[2])} ")
                n_htab += 1
                j += 1
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))
    print(f"HTab reancorado em {n_htab} secoes")

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    A2.sort(key=lambda d: -d["rs"])
    B2 = ler_secoes(novo)
    B2.sort(key=lambda d: -d["rs"])
    za = np.array([float(x["z"].min()) for x in A2])
    zb = np.array([float(x["z"].min()) for x in B2])
    ca = np.array([float(x["rb"] - x["lb"]) for x in A2])
    cb = np.array([float(x["rb"] - x["lb"]) for x in B2])
    la = np.array([float(x["sta"][-1]) for x in A2])
    lb2 = np.array([float(x["sta"][-1]) for x in B2])
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    print(f"   largura do canal mudou : max {np.abs(cb-ca).max():.6f} m "
          "(tem de ser zero)")
    print(f"   largura da secao mudou : max {np.abs(lb2-la).max():.6f} m "
          "(tem de ser zero)")
    print(f"   leito mudou            : max {np.abs(zb-za).max():.3f} m")
    print(f"   leito fora do intervalo original [{za.min():.2f}, "
          f"{za.max():.2f}]: "
          f"{int(((zb < za.min()-1e-6) | (zb > za.max()+1e-6)).sum())}"
          "  (tem de ser zero -- nao se inventa cota)")
    for rot, v in (("antes", za), ("depois", zb)):
        d = np.diff(v)
        ig = int((np.abs(d) < TOL_IGUAL).sum())
        adv = int((d > 1e-9).sum())
        ch = np.array([float(x["len_ch"]) for x in
                       (A2 if rot == "antes" else B2)])[:-1]
        s = np.abs(d) / np.maximum(ch, 1e-9)
        print(f"   {rot:<6}: pares iguais {ig:4d}   sobem para jusante {adv:3d}"
              f"   declividade >2% {int((s>0.02).sum()):3d}   "
              f"max {100*s.max():.2f}%")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
