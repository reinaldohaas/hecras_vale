# -*- coding: utf-8 -*-
"""Faz o rio descer pelo canal e abandona o meandro: um reach so, sem juncao.

    python scripts/unificar_canal.py modelo/mirim_t30/mirim_t30.g11 --saida g12

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O QUE MUDA

  A bifurcacao deixa de existir. O rio passa a ser um unico reach:

      R1 (1375 secoes)  +  canal (51)  +  R3 (2)  =  1428 secoes

  As 41 secoes do curso meandrico saem do roteamento. Elas nao sao apagadas do
  mundo -- o leito antigo continua la, no terreno --, mas deixam de conduzir
  vazao, que e o que "vira planicie de inundacao" significa em modelo 1D.

  E COM ELAS SAEM AS JUNCOES. Com um caminho so, cada juncao ficaria com uma
  entrada e uma saida, e o HEC-RAS recusa isso ("Junctions are for flow
  confluences and splits"). Um reach unico e a forma correta -- e leva junto
  todos os erros que a bifurcacao trouxe.

O ESTACIONAMENTO CONTINUA SENDO DISTANCIA

  Neste modelo o RS e exatamente a distancia acumulada: conferido, `Length Ch`
  bate com a diferenca de RS com erro 0,000000 m. Quebrar isso seria legal para
  o HEC-RAS e pessimo para quem le. Entao o canal e o trecho de jusante sao
  RENUMERADOS para manter a identidade, consumindo os vaos das juncoes:

      R1     ate RS 20359,10          (intocado)
      vao da bifurcacao               301,67 m
      canal  de 20057,43 a 12557,43   (7.500 m)
      vao do reencontro               474,94 m
      R3     de 12082,49 a 11411,42   (671,07 m)

  O rio encurta de 141,35 km para 130,01 km -- os 11,34 km que o canal corta.
  Como o RS de jusante muda, o `.uNN` tem de acompanhar: ver
  `reestruturar_fluxo.py --unico`.

O QUE ISTO NAO FAZ, E PRECISA SER DITO

  Nao poe o meandro abandonado como armazenamento. Medido: o curso antigo esta
  a 1.498 m do eixo do canal na mediana e a 2.535 m no extremo, e NENHUMA das
  51 secoes do canal o alcanca. Para engoli-lo como planicie as secoes
  precisariam de cerca de 5 km de largura, atravessando varzea plana -- onde a
  agua nao corre na direcao da secao e a secao 1D nao descreve nada.

  A forma certa de guardar essa varzea e area de armazenamento ligada
  lateralmente ao canal, e isso e outra etapa. Enquanto nao existir, o modelo
  perde o volume que o meandro segura na cheia.
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
    i0 = next(i for i, l in enumerate(t)
              if l.startswith(("Junct Name=", "River Reach=")))
    cab = t[:i0]
    juncs, reaches = {}, []
    i = i0
    atual = None
    while i < len(t):
        l = t[i]
        if l.startswith("Junct Name="):
            nome = l.split("=", 1)[1].strip()
            j = i + 1
            la = []
            while j < len(t) and not t[j].startswith(("Junct Name=",
                                                      "River Reach=")):
                if t[j].startswith("Junc L&A="):
                    la.append(float(t[j].split("=")[1].split(",")[0]))
                j += 1
            juncs[nome] = la
            i = j
            continue
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            n = int(t[i + 1].split("=")[1])
            v, j = [], i + 2
            while len(v) < 2 * n:
                v += [float(t[j][c:c + L16]) for c in range(0, len(t[j]), L16)
                      if t[j][c:c + L16].strip()]
                j += 1
            txt = next((x for x in t[j:j + 3]
                        if x.startswith("Rch Text X Y")),
                       "Rch Text X Y=0,0,0,0")
            atual = {"rio": p[0].strip(), "reach": p[1].strip(),
                     "xy": np.array(v).reshape(-1, 2), "txt": txt,
                     "blocos": [], "rs": []}
            reaches.append(atual)
            i = j
            continue
        if l.startswith("Type RM Length L Ch R"):
            k = i + 1
            while k < len(t) and not t[k].startswith(
                    ("Type RM Length L Ch R", "River Reach=", "Junct Name=")):
                k += 1
            atual["blocos"].append(list(t[i:k]))
            atual["rs"].append(float(t[i].split(",")[1]))
            i = k
            continue
        i += 1
    return cab, juncs, reaches


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g12")
    fora = _arg(argv, "--descartar", "R2")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    cab, juncs, reaches = ler(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")
    for r in reaches:
        print(f"   {r['rio']:<14} {r['reach']:<4} {len(r['blocos']):5d} secoes"
              f"   RS {max(r['rs']):9.1f} a {min(r['rs']):8.1f}")
    print(f"   juncoes: {', '.join(f'{k} {v}' for k, v in juncs.items())}")

    def pega(rio, rch):
        for r in reaches:
            if (r["rio"], r["reach"]) == (rio, rch):
                return r
        raise SystemExit(f"nao achei {rio},{rch}")

    rio = reaches[0]["rio"]
    R1 = pega(rio, "R1")
    CAN = next(r for r in reaches if r["rio"] != rio)
    R3 = pega(rio, "R3")
    descartado = pega(rio, fora)
    # OS VAOS SAO MEDIDOS, e nao lidos de `Junc L&A`: aquele campo foi escrito
    # por `inserir_canal.py` com 100 m de reserva, e nao com a distancia real
    # -- usa-lo encurtaria o rio em 580 m sem que nada avisasse. A distancia
    # que importa e a que separa as secoes extremas de cada trecho.
    def centro(bloco):
        i = next(k for k, l in enumerate(bloco)
                 if l.startswith("XS GIS Cut Line"))
        v = [float(bloco[i + 1][c:c + L16])
             for c in range(0, len(bloco[i + 1]), L16)
             if bloco[i + 1][c:c + L16].strip()]
        return 0.5 * (np.array(v[:2]) + np.array(v[-2:]))

    def extremo(r, maior):
        k = int(np.argmax(r["rs"]) if maior else np.argmin(r["rs"]))
        return centro(r["blocos"][k])

    g1 = float(np.hypot(*(extremo(CAN, True) - extremo(R1, False))))
    g2 = float(np.hypot(*(extremo(R3, True) - extremo(CAN, False))))
    dec = [j[0] for j in juncs.values() if j]
    print(f"\ndescartado: {rio},{fora}  ({len(descartado['blocos'])} secoes)")
    print(f"vaos MEDIDOS entre as secoes extremas: {g1:.2f} m e {g2:.2f} m")
    print(f"   (o que estava declarado em Junc L&A: "
          f"{', '.join('%.2f' % x for x in dec)} m)")

    # ---- novo estacionamento, mantendo RS == distancia
    ordem = []
    rs0 = min(R1["rs"])
    for b, rs in sorted(zip(R1["blocos"], R1["rs"]), key=lambda x: -x[1]):
        ordem.append((b, rs))
    base_can = rs0 - g1
    can = sorted(zip(CAN["blocos"], CAN["rs"]), key=lambda x: -x[1])
    rs_can_max = max(CAN["rs"])
    for b, rs in can:
        ordem.append((b, base_can - (rs_can_max - rs)))
    base_r3 = ordem[-1][1] - g2
    r3 = sorted(zip(R3["blocos"], R3["rs"]), key=lambda x: -x[1])
    rs_r3_max = max(R3["rs"])
    for b, rs in r3:
        ordem.append((b, base_r3 - (rs_r3_max - rs)))

    rs_novo = np.array([x[1] for x in ordem])
    if (np.diff(rs_novo) >= 0).any():
        raise SystemExit("o RS novo nao ficou estritamente decrescente")
    print(f"\nreach unico: {len(ordem)} secoes   "
          f"RS {rs_novo[0]:.2f} a {rs_novo[-1]:.2f}")
    print(f"   extensao: {(rs_novo[0]-rs_novo[-1])/1000:.2f} km   "
          f"(era {(max(R1['rs'])-min(R3['rs']))/1000:.2f} km)")
    print(f"   canal   : RS {ordem[len(R1['blocos'])][1]:.2f} a "
          f"{ordem[len(R1['blocos'])+len(can)-1][1]:.2f}")
    print(f"   jusante : RS {base_r3:.2f} a {rs_novo[-1]:.2f}")

    # ---- reescreve
    saida = list(cab)
    saida.append(f"River Reach={_pad(rio)},{_pad('R1')}")
    xy = np.vstack([R1["xy"], CAN["xy"], R3["xy"]])
    saida.append(f"Reach XY= {len(xy)} ")
    lin = ""
    for k, (x, y) in enumerate(xy):
        lin += "%16.4f%16.4f" % (x, y)
        if (k + 1) % 2 == 0:
            saida.append(lin)
            lin = ""
    if lin:
        saida.append(lin)
    saida.append(R1["txt"])
    saida.append("")
    for k, (b, rs) in enumerate(ordem):
        comp = (rs - ordem[k + 1][1]) if k < len(ordem) - 1 else 0.0
        b = list(b)
        b[0] = re.sub(r"^(Type RM Length L Ch R\s*=\s*[^,]+),.*$",
                      r"\1,%.2f,%8.2f,%8.2f,%8.2f" % (rs, comp, comp, comp),
                      b[0])
        saida += b
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    from qc_secoes import ler_secoes
    B = ler_secoes(novo)
    B.sort(key=lambda d: -d["rs"])
    rs = np.array([d["rs"] for d in B])
    ch = np.array([float(d["len_ch"]) for d in B])
    t2 = open(novo, encoding="latin-1", errors="replace").read()
    print(f"   secoes           : {len(B)}   "
          f"(esperado {len(R1['blocos'])}+{len(can)}+{len(r3)}="
          f"{len(R1['blocos'])+len(can)+len(r3)})")
    print(f"   River Reach=     : {t2.count('River Reach=')}   (esperado 1)")
    print(f"   Junct Name=      : {t2.count('Junct Name=')}   (esperado 0)")
    print(f"   RS decrescente   : {bool((np.diff(rs) < 0).all())}")
    err = np.abs(ch[:-1] - (-np.diff(rs))).max() if len(B) > 1 else 0.0
    print(f"   Ch == diferenca de RS: erro maximo {err:.6f} m")
    print(f"   ultima secao com Ch=0: {ch[-1] == 0.0}")
    a = open(entrada, encoding="latin-1", errors="replace").read()
    for chave in ("#Sta/Elev=", "Bank Sta=", "XS GIS Cut Line="):
        print(f"   {chave:<20} {a.count(chave)} -> {t2.count(chave)}  "
              f"(saem {len(descartado['blocos'])})")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
