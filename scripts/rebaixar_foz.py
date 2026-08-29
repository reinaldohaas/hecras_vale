# -*- coding: utf-8 -*-
"""Rampa o leito do reach de montante ate ENCONTRAR o leito da juncao.

    python scripts/rebaixar_foz.py taha_ai_novo/taha_ai.g01 --saida g21 \
        --dist 1500 --minimo 0.3

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO (medido no taha_ai_novo, .bco de 26/08/2026)

  A escavacao sintetica de cada rio foi calculada por rio, sem coordenar as
  juncoes: o Iraputa chega na foz com leito 3,33 m ACIMA da cabeca do Norte
  R2; o Sul, 4,29 m acima do Acu. Na inicializacao do unsteady o espelho do
  receptor fica ABAIXO do leito do tributario -- foz seca sobre degrau -- e
  o solver de juncao (igual espelho) explode ali no primeiro passo:
  Q passa de -617 para 499.152 e depois 23.442.830 m3/s em tres iteracoes,
  cota vai a 2.605 m. E dai a instabilidade contamina a rede inteira.

O QUE SE FAZ

  Para cada juncao, e cada reach de montante cuja ULTIMA secao tenha
  talvegue mais de `--minimo` m acima do leito da PRIMEIRA secao de jusante:
  as secoes dos ultimos `--dist` m recebem um rebaixamento LINEAR do canal
  (so entre lb e rb), de zero no inicio da rampa ate o degrau inteiro na
  foz. So se REBAIXA -- nunca se sobe leito. O `XS HTab Starting El` de
  cada secao mexida acompanha o talvegue novo (leito + 0,15 m, a convencao
  do arquivo). A planicie (fora de lb..rb) nao muda.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def juncoes_do_g01(linhas):
    """[(nome, [ups], dn)] lidos do texto."""
    out, nome, ups, dn = [], None, [], None
    for l in linhas:
        if l.startswith("Junct Name="):
            if nome:
                out.append((nome, ups, dn))
            nome, ups, dn = l.split("=", 1)[1].strip(), [], None
        elif l.startswith("Up River,Reach="):
            p = l.split("=", 1)[1].split(",")
            ups.append((p[0].strip(), p[1].strip()))
        elif l.startswith("Dn River,Reach="):
            p = l.split("=", 1)[1].split(",")
            dn = (p[0].strip(), p[1].strip())
    if nome:
        out.append((nome, ups, dn))
    return out


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g21")
    dist = _arg(argv, "--dist", 1500.0, float)
    minimo = _arg(argv, "--minimo", 0.3, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    S = ler_secoes(entrada)
    por = {}
    for i, d in enumerate(S):
        por.setdefault((d["rio"], d["reach"]), []).append((i, d))
    for k in por:
        por[k].sort(key=lambda t: -t[1]["rs"])

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"rampa  : {dist:.0f} m   degrau minimo: {minimo} m\n")

    rebaixo = {}          # indice global -> quanto abaixar o canal
    for nome, ups, dn in juncoes_do_g01(linhas):
        if dn not in por:
            continue
        zdn = float(np.asarray(por[dn][0][1]["z"], float).min())
        for u in ups:
            if u not in por:
                continue
            i_last, d_last = por[u][-1]
            zu = float(np.asarray(d_last["z"], float).min())
            degrau = zu - zdn
            if degrau <= minimo:
                continue
            rs_last = d_last["rs"]
            n_mex = 0
            for i, d in por[u]:
                s = d["rs"] - rs_last          # distancia ate a foz
                if s > dist:
                    continue
                delta = degrau * (1.0 - s / dist)
                atual = rebaixo.get(i, 0.0)
                rebaixo[i] = max(atual, delta)
                n_mex += 1
            print(f"   {nome:18s} {u[0]:13s} {u[1]:3s}: degrau "
                  f"{degrau:+5.2f} m, rampa em {n_mex} secoes")

    if not rebaixo:
        print("nenhum degrau de juncao acima do minimo")
        return

    # ------------------------------------------------------- reescreve
    # chave por (rio, reach, RS): o indice de ler_secoes NAO acompanha o
    # arquivo quando ha ESTRUTURA (barragem Type 5) no meio -- o gravador
    # por indice escrevia o perfil rebaixado na secao ERRADA dali em
    # diante (pego na foz do Luis Alves, 28/08)
    novos = {}
    for i, delta in rebaixo.items():
        d = S[i]
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float).copy()
        m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
        z[m] = z[m] - delta
        novos[(d["rio"], d["reach"], round(d["rs"], 2))] = \
            {"sta": st, "z": z, "htab": float(z.min()) + 0.15}

    saida, j = [], 0
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
            except (ValueError, IndexError):
                chave = None
        nv = novos.get(chave)
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

    # ------------------------------------------------------ conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    por2 = {}
    for d in B:
        por2.setdefault((d["rio"], d["reach"]), []).append(d)
    for k in por2:
        por2[k].sort(key=lambda d: -d["rs"])
    pior = 0.0
    for nome, ups, dn in juncoes_do_g01(linhas):
        if dn not in por2:
            continue
        zdn = float(np.asarray(por2[dn][0]["z"], float).min())
        for u in ups:
            if u not in por2:
                continue
            zu = float(np.asarray(por2[u][-1]["z"], float).min())
            pior = max(pior, zu - zdn)
    print(f"   maior degrau de juncao restante: {pior:+.2f} m "
          f"(era ate +4.29)")
    sobe = 0
    for k, secs in por2.items():
        zz = [float(np.asarray(d['z'], float).min()) for d in secs]
        sobe += sum(1 for a, b in zip(zz, zz[1:]) if b > a + 0.01)
    print(f"   pares vizinhos com leito SUBINDO rio abaixo: {sobe}")


if __name__ == "__main__":
    main(sys.argv[1:])
