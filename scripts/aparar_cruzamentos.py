# -*- coding: utf-8 -*-
"""Apara pares de cutlines que se CRUZAM -- inclusive as nao vizinhas.

    python scripts/aparar_cruzamentos.py taha_ai_novo/taha_ai.g01 --saida g15

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Duas cutlines que se tocam representam a mesma agua duas vezes, e sao a
  materia-prima das dobras: a edge line liga as PONTAS das secoes vizinhas e
  a bank line liga as margens -- onde cutlines se cruzam, essas linhas
  derivadas dobram sobre si mesmas ("edge line self intersection",
  "Polyline has self intersections", "XS intersects > 2 banklines").

  O passe de vizinhas do corrigir_cutlines so compara pares CONSECUTIVOS no
  mesmo reach. Medido no g01 reparado: 67 pares cruzados, e a maioria NAO e
  vizinha -- e a secao comprida de um lado do meandro alcancando a volta
  seguinte do rio, ate de outro reach.

O QUE SE FAZ

  Para cada par cruzado, apara-se a secao em que o cruzamento cai mais LONGE
  do canal -- na ponta, alem de lb-FOLGA..rb+FOLGA -- recuando RECUO m antes
  do ponto de cruzamento. O CANAL NUNCA E INVADIDO: se nas duas secoes o
  cruzamento cair dentro da faixa protegida, o par fica como esta e sai no
  relatorio. Cotas sobreviventes, talvegue e largura de canal nao mudam.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import (mapa_reaches, _col, _fmt,
                               _arg)                   # noqa: E402
from ras_io import escrever                            # noqa: E402

FOLGA = 8.0     # m preservados alem da margem (o piso do corrigir_cutlines)
RECUO = 3.0     # m antes do ponto de cruzamento


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    from shapely.strtree import STRtree
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g15")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    mapa = mapa_reaches(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    A = [np.asarray(d["cut"][0], float) for d in S]
    B = [np.asarray(d["cut"][-1], float) for d in S]
    L0 = [float(np.hypot(*(b - a))) for a, b in zip(A, B)]
    U = [(b - a) / max(l, 1e-9) for a, b, l in zip(A, B, L0)]
    lim = [[0.0, l] for l in L0]

    def linha(i):
        return LineString([A[i] + lim[i][0] * U[i], A[i] + lim[i][1] * U[i]])

    def cruzam(i, j):
        li, lj = linha(i), linha(j)
        if not li.intersects(lj) or li.touches(lj):
            return None
        x = li.intersection(lj)
        p = x if x.geom_type == "Point" else x.centroid
        return np.array([p.x, p.y])

    def apara(i, p):
        """Tenta aparar a secao i no ponto p. True se conseguiu."""
        s = float(np.dot(p - A[i], U[i]))
        lb, rb = float(S[i]["lb"]), float(S[i]["rb"])
        if s > rb + FOLGA:
            lim[i][1] = min(lim[i][1], s - RECUO)
            return True
        if s < lb - FOLGA:
            lim[i][0] = max(lim[i][0], s + RECUO)
            return True
        return False

    tree = STRtree([LineString([a, b]) for a, b in zip(A, B)])
    pares = set()
    for i in range(len(S)):
        for j in tree.query(linha(i)):
            j = int(j)
            if j > i:
                pares.add((i, j))

    presos, aparados = [], set()
    for _ in range(4):          # aparar um par pode desfazer outro; itera
        mexeu = False
        for i, j in sorted(pares):
            p = cruzam(i, j)
            if p is None:
                continue
            si = float(np.dot(p - A[i], U[i]))
            sj = float(np.dot(p - A[j], U[j]))
            di = min(abs(si - S[i]["lb"]), abs(si - S[i]["rb"]))
            dj = min(abs(sj - S[j]["lb"]), abs(sj - S[j]["rb"]))
            # a vitima preferida: onde o cruzamento cai mais longe do canal
            ordem = (i, j) if di >= dj else (j, i)
            ok = apara(ordem[0], p) or apara(ordem[1], p)
            if ok:
                aparados |= {k for k in ordem if lim[k] != [0.0, L0[k]]}
                mexeu = True
        if not mexeu:
            break
    restantes = [(i, j) for i, j in pares if cruzam(i, j) is not None]
    for i, j in restantes:
        presos.append((mapa[i], S[i]["rs"], mapa[j], S[j]["rs"]))

    print(f"pares cruzados tratados: {len(pares and aparados)} secoes "
          f"aparadas")
    print(f"pares ainda cruzando (canal protegido dos dois lados): "
          f"{len(restantes)}")
    for a, ra, b, rb_ in presos[:12]:
        print(f"   fica: {a[0]} RS {ra:.1f}  x  {b[0]} RS {rb_:.1f}")

    # ------------------------------------------------ reescreve (padrao
    # identico ao corrigir_cutlines: so saem pontos das pontas)
    TOL = 0.005
    novos = {}
    for i in sorted(aparados):
        s0, s1 = lim[i]
        if s0 <= TOL and s1 >= L0[i] - TOL:
            continue
        st = np.asarray(S[i]["sta"], float)
        z = np.asarray(S[i]["z"], float)
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
        keep = [0]
        for k in range(1, len(ns)):
            if round(ns[k], 2) > round(ns[keep[-1]], 2):
                keep.append(k)
        ns, nz = ns[keep], nz[keep]
        novos[i] = {"sta": ns, "z": nz, "lb": float(S[i]["lb"]) - s0,
                    "rb": float(S[i]["rb"]) - s0, "desl": s0,
                    "cut": (A[i] + s0 * U[i], A[i] + s1 * U[i])}

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
                val[0] = 0.0
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
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print(f"   secoes: {len(A2)} -> {len(B2)}   (nao pode mudar)")
    tal = max(abs(float(a['z'].min()) - float(b['z'].min()))
              for a, b in zip(A2, B2))
    print(f"   talvegue mudou no maximo {tal:.6f}  (tem de ser zero)")
    Lb = [LineString(np.asarray(d['cut'], float)) for d in B2]
    tb = STRtree(Lb)
    n = 0
    for i, ln in enumerate(Lb):
        for jx in tb.query(ln):
            jx = int(jx)
            if jx > i and Lb[i].intersects(Lb[jx]) \
                    and not Lb[i].touches(Lb[jx]):
                n += 1
    print(f"   pares de cutlines cruzando: {n}")


if __name__ == "__main__":
    main(sys.argv[1:])
