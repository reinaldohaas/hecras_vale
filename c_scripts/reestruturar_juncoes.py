# -*- coding: utf-8 -*-
"""Parte o reach unico em tres, com duas juncoes onde o canal vai entrar.

    python scripts/reestruturar_juncoes.py modelo/mirim_t30/mirim_t30.g02 \
        --rs-montante 20359.1 --rs-jusante 1221.0 --saida g04

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

POR QUE PARTIR ANTES DE TER O CANAL

  O Canal Retificado e o meandro antigo sao dois caminhos entre os MESMOS dois
  pontos. No HEC-RAS isso e um reach por caminho, ligados por duas juncoes --
  e o modelo de hoje tem UM reach so, sem juncao nenhuma. Enquanto ele for um
  reach unico nao ha onde encaixar o canal.

  Partir primeiro tem uma vantagem que vale o passo: com um caminho so, as
  duas juncoes ficam com uma entrada e uma saida, e a hidraulica TEM de sair
  identica a de antes. Isso torna o passo VERIFICAVEL -- se o resultado mudar,
  o defeito esta na partição, e nao no canal que ainda nem existe.

O QUE MUDA E O QUE NAO

  Nenhuma secao e criada, apagada ou alterada: as 1.418 continuam iguais, com
  as mesmas cotas, estacas, margens e Manning. O que muda e a MOLDURA:

    - `River Reach=` passa de um para tres (R1, R2, R3)
    - `Reach XY` (a polilinha do eixo) e cortada nos dois pontos de juncao
    - a ULTIMA secao de cada reach recebe comprimento zero, porque o vao ate
      o reach seguinte passa a ser declarado na juncao (`Junc L&A`) -- e assim
      que o HEC-RAS conta essa distancia, e deixar o comprimento antigo a
      contaria duas vezes
    - dois blocos `Junct Name=`, que no arquivo vem ANTES dos reaches

ONDE AS JUNCOES FICAM

  Nos dois pontos em que o canal encosta no rio modelado: RS 20359,1 a
  montante e RS 1221,0 a jusante. Medido antes, as pontas do canal caem a
  2,4 m e 11,7 m do eixo do modelo -- ou seja, ele comeca e termina em cima
  do rio, e nao ha ambiguidade sobre onde parte.

QUANDO O CANAL CHEGAR

  Ele entra como um quarto reach, ligado as duas juncoes ja existentes: vira
  `Up River,Reach` na de jusante e `Dn River,Reach` na de montante, com o seu
  proprio `Junc L&A`. Nada do que esta aqui precisa ser refeito.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever    # noqa: E402

L16 = 16


def _pad(s):
    return f"{s:<{L16}}"[:L16]


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def ler(entrada):
    t = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_rr = next(i for i, l in enumerate(t) if l.startswith("River Reach="))
    cab = t[:i_rr]
    p = t[i_rr].split("=", 1)[1].split(",")
    rio = p[0].strip()
    n = int(t[i_rr + 1].split("=")[1])
    xy, j = [], i_rr + 2
    while len(xy) < 2 * n:
        xy += [float(t[j][c:c + L16]) for c in range(0, len(t[j]), L16)
               if t[j][c:c + L16].strip()]
        j += 1
    xy = np.array(xy).reshape(-1, 2)
    txt = next((l for l in t[j:j + 4] if l.startswith("Rch Text X Y")),
               "Rch Text X Y=0,0,0,0")
    ini = [i for i, l in enumerate(t) if l.startswith("Type RM Length L Ch R")]
    blocos, rss = [], []
    for a, b in zip(ini, ini[1:] + [len(t)]):
        blocos.append(list(t[a:b]))
        rss.append(float(t[a].split(",")[1]))
    return cab, rio, xy, txt, blocos, np.array(rss)


def _xy(coords):
    saida, linha = [], ""
    for k, (x, y) in enumerate(coords):
        linha += "%16.4f%16.4f" % (x, y)
        if (k + 1) % 2 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString, Point
    from shapely.ops import substring
    entrada = argv[0]
    rsa = _arg(argv, "--rs-montante", 20359.1, float)
    rsb = _arg(argv, "--rs-jusante", 1221.0, float)
    ext = _arg(argv, "--saida", "g04")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    cab, rio, xy, txt, blocos, rss = ler(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"rio '{rio}'   {len(blocos)} secoes   "
          f"RS {rss.max():.1f} a {rss.min():.1f}")

    ka = int(np.argmin(np.abs(rss - rsa)))
    kb = int(np.argmin(np.abs(rss - rsb)))
    if not (0 < ka < kb < len(blocos) - 1):
        raise SystemExit("os pontos de juncao nao partem o rio em tres")
    grupos = [("R1", range(0, ka + 1)),
              ("R2", range(ka + 1, kb + 1)),
              ("R3", range(kb + 1, len(blocos)))]
    print(f"\ncorte a montante: RS {rss[ka]:.1f}   "
          f"a jusante: RS {rss[kb]:.1f}")
    for nome, r in grupos:
        r = list(r)
        print(f"   {nome}: {len(r):5d} secoes   "
              f"RS {rss[r[0]]:9.1f} a {rss[r[-1]]:8.1f}")

    # ---- onde cortar a polilinha do eixo
    eixo = LineString(xy)
    cen = []
    for b in blocos:
        i = next(i for i, l in enumerate(b) if l.startswith("XS GIS Cut Line"))
        v = [float(b[i + 1][c:c + L16]) for c in range(0, len(b[i + 1]), L16)
             if b[i + 1][c:c + L16].strip()]
        cen.append(0.5 * (np.array(v[:2]) + np.array(v[-2:])))
    cen = np.array(cen)
    sj = []
    for k in (ka, kb):
        s0 = eixo.project(Point(*cen[k]))
        s1 = eixo.project(Point(*cen[k + 1]))
        sj.append(0.5 * (s0 + s1))
    if not (0 < sj[0] < sj[1] < eixo.length):
        raise SystemExit(f"cortes do eixo fora de ordem: {sj}")
    trechos = [substring(eixo, 0, sj[0]),
               substring(eixo, sj[0], sj[1]),
               substring(eixo, sj[1], eixo.length)]
    print(f"\neixo partido: {eixo.length/1000:.2f} km -> "
          + " + ".join(f"{t.length/1000:.2f}" for t in trechos) + " km")

    # ---- juncoes
    juncs = []
    for n, (k, nome, up, dn) in enumerate(
            ((ka, "Bifurcacao", "R1", "R2"),
             (kb, "Reencontro", "R2", "R3"))):
        P = np.array(eixo.interpolate(sj[n]).coords[0])
        comp = float(rss[k] - rss[k + 1])
        juncs.append([
            f"Junct Name={_pad(nome)}",
            f"Junct Desc=Canal Retificado, 0 , 0 , 0 ,0",
            f"Junct X Y & Text X Y={P[0]:.2f},{P[1]:.2f},"
            f"{P[0]+800:.2f},{P[1]+800:.2f}",
            f"Up River,Reach={_pad(rio)},{_pad(up)}",
            f"Dn River,Reach={_pad(rio)},{_pad(dn)}",
            f"Junc L&A={comp:.2f},0",
            "",
        ])
        print(f"   juncao '{nome}' em ({P[0]:.0f}, {P[1]:.0f})   "
              f"vao de {comp:.2f} m  ({up} -> {dn})")

    # ---- monta o arquivo
    saida = list(cab)
    for j in juncs:
        saida += j
    for n, (nome, r) in enumerate(grupos):
        r = list(r)
        saida.append(f"River Reach={_pad(rio)},{_pad(nome)}")
        c = np.asarray(trechos[n].coords)
        saida.append(f"Reach XY= {len(c)} ")
        saida += _xy(c)
        saida.append(txt)
        saida.append("")
        for i, k in enumerate(r):
            b = list(blocos[k])
            if i == len(r) - 1:      # ultimo do reach: comprimento zero
                b[0] = re.sub(
                    r"^(Type RM Length L Ch R\s*=\s*[^,]+,[^,]+),.*$",
                    r"\1,%8.2f,%8.2f,%8.2f" % (0.0, 0.0, 0.0), b[0])
            saida += b
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    t2 = open(novo, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    nrr = len(re.findall(r"(?m)^River Reach=", t2))
    njn = len(re.findall(r"(?m)^Junct Name=", t2))
    nxs = len(re.findall(r"(?m)^Type RM Length L Ch R", t2))
    print(f"   reaches {nrr} (era 1)   juncoes {njn} (era 0)   "
          f"secoes {nxs} (era {len(blocos)})")
    # as secoes nao podem ter mudado
    a = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    for chave in ("#Sta/Elev=", "Bank Sta=", "#Mann=",
                  "XS HTab Starting El and Incr=", "XS GIS Cut Line="):
        na = a.count(chave)
        nb = t2.count(chave)
        print(f"   {chave:<32} {na} -> {nb}   "
              f"{'ok' if na == nb else 'DIVERGIU'}")
    # comprimento total do rio: soma dos Ch + os vaos das juncoes
    def soma(txt_):
        v = [float(m.group(1)) for m in re.finditer(
            r"(?m)^Type RM Length L Ch R\s*=\s*[^,]+,[^,]+,[^,]+,\s*([-\d.]+)",
            txt_)]
        return sum(v)
    sa, sb = soma(a), soma(t2)
    vaos = sum(float(m.group(1)) for m in
               re.finditer(r"(?m)^Junc L&A=([-\d.]+)", t2))
    print(f"   soma dos comprimentos de trecho: {sa:.2f} -> "
          f"{sb:.2f} + {vaos:.2f} de juncao = {sb+vaos:.2f} m")
    print(f"   diferenca: {sa-(sb+vaos):+.2f} m  (tem de ser zero)")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
