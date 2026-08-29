# -*- coding: utf-8 -*-
"""Aplica o teto de declividade do gerador (decl_teto) por PREENCHIMENTO.

    python scripts/limitar_declividade.py taha_ai_novo/taha_ai.g01 \
        --saida g25 --teto 0.05

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  O contrato do gerador (opcoes.json) diz `decl_teto: 0.05` -- e o leito
  gravado tem trechos de 7 a 9% (Benedito 38832->38682: 8,73% em 150 m).
  Numa cascata dessas o unsteady 1D nao para em pe na vazao baixa: a lamina
  afina a centimetros, o passo seguinte poca, e o solver oscila ate abortar
  (medido em quatro rodadas seguidas, sempre com epicentro nas cascatas).

O QUE SE FAZ

  Andando rio abaixo, secao cuja queda ate a vizinha exceda `teto * dx`
  tem a VIZINHA DE JUSANTE preenchida ate o limite: o fundo do poco de
  queda sobe, nunca mais que o necessario, e o preenchimento se propaga
  ate a declividade natural voltar a ficar abaixo do teto. E o aterro de
  poco de queda -- o 1D nao resolve cachoeira de qualquer forma; o que se
  perde e profundidade morta de poco, o que se ganha e um perfil que o
  solver aguenta. So o canal (lb..rb) sobe; planicie intacta; HTab
  acompanha o talvegue novo.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g25")
    teto = _arg(argv, "--teto", 0.05, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    por = {}
    for i, d in enumerate(S):
        por.setdefault((d["rio"], d["reach"]), []).append((i, d))
    for k in por:
        por[k].sort(key=lambda t: -t[1]["rs"])

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"teto   : {teto:.1%}\n")

    novos = {}
    for k, secs in por.items():
        z_ant = None
        rs_ant = None
        for i, d in secs:
            zt = float(np.asarray(d["z"], float).min())
            zt = novos[i]["z"].min() if i in novos else zt
            if z_ant is not None:
                dx = rs_ant - d["rs"]
                piso = z_ant - teto * dx
                if zt < piso - 0.01:
                    delta = float(piso - zt)
                    st = np.asarray(d["sta"], float)
                    z = np.asarray(d["z"], float).copy()
                    # o piso vale para a SECAO INTEIRA: na garganta o fundo
                    # do vale aparece fora das margens tambem, e aterrar so
                    # o canal deixaria o talvegue global (que e o que o
                    # solver ve) abaixo do teto -- foi medido: 8,13% de
                    # declividade sobraram com o aterro so no canal
                    z = np.maximum(z, piso)
                    novos[i] = {"sta": st, "z": z,
                                "htab": float(z.min()) + 0.15}
                    print(f"   {k[0]:13s} {k[1]:3s} RS {d['rs']:9.1f}  "
                          f"poco preenchido +{delta:5.2f} m "
                          f"(decl era {(z_ant-zt)/max(dx,1):.1%})")
                    zt = float(z.min())
            z_ant, rs_ant = zt, d["rs"]

    if not novos:
        print("nenhuma queda acima do teto")
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
            if l.startswith("XS HTab Starting El and Incr="):
                resto = l.split("=", 1)[1].split(",")
                saida.append("XS HTab Starting El and Incr=%.2f,%s,%s"
                             % (nv["htab"], resto[1], resto[2]))
                j += 1
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    por2 = {}
    for d in B:
        por2.setdefault((d["rio"], d["reach"]), []).append(d)
    pior = 0.0
    for k, secs in por2.items():
        secs.sort(key=lambda d: -d["rs"])
        for a, b in zip(secs, secs[1:]):
            za = float(np.asarray(a["z"], float).min())
            zb = float(np.asarray(b["z"], float).min())
            dx = a["rs"] - b["rs"]
            if dx > 0:
                pior = max(pior, (za - zb) / dx)
    print(f"   declividade maxima restante: {pior:.2%}   "
          f"(tem de ser <= {teto:.0%} + folga)")


if __name__ == "__main__":
    main(sys.argv[1:])
