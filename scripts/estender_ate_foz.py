# -*- coding: utf-8 -*-
"""Estende o reach ate a foz, trazendo as secoes que ja existem no original.

    python scripts/estender_ate_foz.py modelo/mirim_novo/mirim_novo.g01 \
        --original modelo/so_mirim.g01 --saida g02

A ENTRADA NAO E TOCADA. Sai um .gXX novo, ao lado dela.

O DEFEITO

  O `mirim_novo` termina em (729585, 7023251), que e o ULTIMO VERTICE do Canal
  Retificado no OpenStreetMap -- distancia zero. Isso fica 1.058 m antes da
  confluencia com o Itajai-Acu, e o vetor do OSM provavelmente acaba ali porque
  e onde a BR-101 cruza. Medido no MDT, o corredor que falta tem um aterro a
  805 m: o fundo sobe de 0,05 m para 2,93 m e nao volta.

  A consequencia nao e de desenho. O contorno de jusante e um HIDROGRAMA DE
  MARE do Itajai-Acu, e hoje ele e imposto 1 km rio acima de onde a mare entra,
  do lado de dentro do aterro.

O QUE SE FAZ

  As duas secoes que faltam JA EXISTEM no `so_mirim.g01` -- RS 746,07 e 75,00,
  esta ultima na foz. Elas sao copiadas como estao, sem uma cota alterada.

  O estaqueamento e refeito somando a distancia que passou a existir, para
  manter a identidade deste modelo: RS igual a distancia ate a foz, com a foz
  em zero. Conferido no fim que `Length Ch` continua batendo com a diferenca
  de RS.

  O `Reach XY` e prolongado pelo trecho correspondente do eixo ORIGINAL, e nao
  por uma reta: entre o fim do canal e a foz o rio curva, e ligar em linha reta
  poria o eixo em cima do aterro.

O QUE ISTO TRAZ JUNTO, e precisa ser dito

  A secao da foz declara 640 m de CANAL onde a vizinha tem 76 m. Esse salto de
  8,4 vezes e a origem do degrau hidraulico que ja se via no modelo original --
  a lamina cai 5 m nos ultimos 671 m. Trazer a secao traz o degrau. Ele e do
  dado original, nao desta operacao, e some quando a secao da foz for revista.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from ras_io import escrever            # noqa: E402

L16 = 16


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def blocos(g):
    """(cabecalho, River Reach, xy do eixo, Rch Text, lista de blocos)."""
    t = open(g, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i = next(k for k, l in enumerate(t) if l.startswith("River Reach="))
    cab = t[:i]
    rr = t[i]
    n = int(t[i + 1].split("=")[1])
    v, j = [], i + 2
    while len(v) < 2 * n:
        v += [float(t[j][c:c + L16]) for c in range(0, len(t[j]), L16)
              if t[j][c:c + L16].strip()]
        j += 1
    xy = np.array(v).reshape(-1, 2)
    txt = next((l for l in t[j:j + 4] if l.startswith("Rch Text X Y")),
               "Rch Text X Y=0,0,0,0")
    ini = [k for k, l in enumerate(t) if l.startswith("Type RM Length L Ch R")]
    bl = [list(t[a:b]) for a, b in zip(ini, ini[1:] + [len(t)])]
    rs = [float(b[0].split(",")[1]) for b in bl]
    return cab, rr, xy, txt, bl, rs


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
    from shapely.geometry import Point, LineString
    from shapely.ops import substring
    entrada = argv[0]
    orig = _arg(argv, "--original", "modelo/so_mirim.g01")
    ext = _arg(argv, "--saida", "g02")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    cab, rr, xy, txt, bl, rs = blocos(entrada)
    ordem = np.argsort(-np.array(rs))
    bl = [bl[i] for i in ordem]
    rs = [rs[i] for i in ordem]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(bl)}   RS {rs[0]:.2f} a {rs[-1]:.2f}")

    # ---- as duas que faltam, do original
    So = ler_secoes(orig)
    So.sort(key=lambda d: -d["rs"])
    _, _, xyo, _, blo, rso = blocos(orig)
    oo = np.argsort(-np.array(rso))
    blo = [blo[i] for i in oo]
    rso = [rso[i] for i in oo]
    faltam = [(blo[i], rso[i]) for i in range(len(rso)) if rso[i] < 1000.0]
    print(f"secoes a trazer do original: {len(faltam)}   "
          f"RS {[round(r,2) for _, r in faltam]}")
    if not faltam:
        raise SystemExit("nao achei as secoes de jusante no original")

    eixo_o = LineString(xyo)
    fim = np.array(xy[-1])
    cen = []
    for b, r in faltam:
        i = next(k for k, l in enumerate(b) if l.startswith("XS GIS Cut Line"))
        v = [float(b[i + 1][c:c + L16]) for c in range(0, len(b[i + 1]), L16)
             if b[i + 1][c:c + L16].strip()]
        cen.append(0.5 * (np.array(v[:2]) + np.array(v[-2:])))
    # NAO PROJETAR O FIM DO CANAL NO EIXO ORIGINAL. Aquele eixo desce pelos
    # MEANDROS, e o fim do canal cai perto de um trecho que, ao longo da
    # polilinha, vem DEPOIS da secao de jusante -- a projecao devolveu vao
    # negativo de -786 m. O vao aqui e a distancia direta ate a secao
    # seguinte, que sao 329 m de rio quase reto.
    vao = float(np.hypot(*(cen[0] - fim)))
    # e o comprimento seguinte e o DECLARADO no original (`Length Ch`), que e
    # o que o solver usa -- nao uma medida minha sobre a polilinha
    m_ch = re.match(r"^Type RM Length L Ch R\s*=\s*[^,]+,[^,]+,[^,]+,\s*"
                    r"([-\d.]+)", faltam[0][0][0])
    comp2 = float(m_ch.group(1))
    if vao <= 0 or comp2 <= 0:
        raise SystemExit(f"vaos invalidos: vao={vao}, comp={comp2}")
    print(f"\nvaos:")
    print(f"   fim do canal -> RS {faltam[0][1]:.2f} : {vao:.2f} m  (direta)")
    print(f"   RS {faltam[0][1]:.2f} -> RS {faltam[1][1]:.2f} (foz) : "
          f"{comp2:.2f} m  (Length Ch do original)")
    s_sec = [eixo_o.project(Point(*p)) for p in cen]
    desl = vao + comp2
    print(f"   deslocamento do estaqueamento: +{desl:.2f} m  "
          f"(a foz passa a ser RS 0)")

    # ---- reescreve os RS e os comprimentos
    saida = list(cab)
    saida.append(rr)
    # o eixo novo termina com a CAUDA do eixo original, de RS 746 ate a foz;
    # o pedaco entre o fim do canal e ali fica como um segmento reto, que sao
    # 329 m onde o rio e quase reto
    cauda = np.asarray(substring(eixo_o, min(s_sec), max(s_sec)).coords)
    if np.hypot(*(cauda[0] - cen[1])) < np.hypot(*(cauda[0] - cen[0])):
        cauda = cauda[::-1]
    # O EIXO ORIGINAL NAO CHEGA NA SECAO DA FOZ. Ele termina 230 m antes do
    # centro dela, e uma secao que nao cruza o proprio reach e recusada:
    # "XS doesn't intersect the associated Reach". A cauda passa pelo centro
    # de cada secao trazida e segue um pouco alem da ultima, para que ela
    # CRUZE em vez de tocar na ponta -- o mesmo cuidado de `prolongar_eixo.py`.
    if np.hypot(*(cauda[-1] - cen[1])) > 1.0:
        cauda = np.vstack([cauda, cen[1]])
    d = cauda[-1] - cauda[-2]
    d = d / max(float(np.hypot(*d)), 1e-9)
    cauda = np.vstack([cauda, cauda[-1] + 60.0 * d])
    eixo_novo = np.vstack([xy, cauda])
    saida.append(f"Reach XY= {len(eixo_novo)} ")
    saida += _xy(eixo_novo)
    saida.append(txt)
    saida.append("")
    novos_rs = [r + desl for r in rs] + [comp2, 0.0]
    todos = bl + [faltam[0][0], faltam[1][0]]
    for k, b in enumerate(todos):
        c = (novos_rs[k] - novos_rs[k + 1]) if k < len(todos) - 1 else 0.0
        b = list(b)
        b[0] = re.sub(
            r"^(Type RM Length L Ch R\s*=\s*[^,]+),[^,]+,.*$",
            r"\1,%.2f,%8.2f,%8.2f,%8.2f" % (novos_rs[k], c, c, c), b[0])
        saida += b
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    B.sort(key=lambda d: -d["rs"])
    rb = np.array([d["rs"] for d in B])
    ch = np.array([float(d["len_ch"]) for d in B])
    A = ler_secoes(entrada)
    print(f"   secoes                 : {len(A)} -> {len(B)}   "
          f"(esperado {len(A)}+2)")
    print(f"   RS                     : {rb.max():.2f} a {rb.min():.2f}")
    print(f"   RS estritamente decrescente: "
          f"{bool((np.diff(rb) < 0).all())}")
    err = np.abs(ch[:-1] - (-np.diff(rb))).max()
    print(f"   Ch == diferenca de RS  : erro maximo {err:.6f} m")
    print(f"   ultima secao com Ch=0  : {ch[-1] == 0.0}")
    print(f"   extensao               : {rb.max()-rb.min():.2f} m   "
          f"(era {max(rs)-min(rs):.2f} m)")
    ea = LineString(xy)
    eb = list(ler_eixos(novo).values())[0]
    print(f"   eixo                   : {ea.length/1000:.3f} -> "
          f"{eb.length/1000:.3f} km")
    C = np.asarray(B[-1]["cut"], float)
    p = 0.5 * (C[0] + C[-1])
    print(f"   ultima secao agora em  : ({p[0]:.0f}, {p[1]:.0f})")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
