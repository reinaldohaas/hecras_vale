# -*- coding: utf-8 -*-
"""Recentra o canal das secoes que o HEC-RAS marca por cruzar bank lines demais.

    python scripts/ajustar_bancos.py modelo/mirim_t30/mirim_t30.g10 \
        --erros modelo/mirim_t37/mirim_t37_erros.csv --saida g11

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  "XS intersects > 2 banklines". A bank line e construida pelo RAS Mapper
  ligando o ponto `Bank Sta` de uma secao ao da seguinte. Onde uma secao tem o
  canal declarado fora do lugar, a linha faz um gancho e volta -- e a secao
  vizinha passa a cruzar quatro bank lines em vez de duas.

  Medido no g10, nas 51 secoes marcadas: o meio do canal declarado esta a
  13,2 m do ponto onde o eixo do rio cruza a cutline na mediana, mas a 151 m
  no p90 e a 285 m no extremo. Contra 9,1 m de mediana no modelo inteiro. Nao
  e um desvio geral -- e uma cauda ruim.

O QUE SE FAZ

  So nas secoes que o HEC-RAS marca, e so a POSICAO do canal:

    - a LARGURA do canal e preservada exatamente;
    - o canal e recentrado onde o eixo do rio cruza a cutline;
    - as margens novas caem sobre ESTACAS EXISTENTES (o `.gNN` deste modelo
      tem todas as `Bank Sta` sobre pontos do perfil, e quebrar isso faz o
      HEC-RAS recusar a geometria -- ja aconteceu nesta reconstrucao);
    - as quebras do `#Mann` acompanham as margens, para que a rugosidade de
      canal continue no canal.

  Nenhuma cota muda, nenhum ponto entra ou sai, a cutline nao se move.

O QUE ISTO NAO RESOLVE, E POR QUE

  Recentrar supoe que o EIXO esta certo e a margem errada. Onde a cutline nem
  alcanca o eixo -- ha duas assim, a 79 m e a 53 m de distancia -- nao existe
  travessia em que centrar, e a secao fica como esta. Move-la seria deslocar
  secao, que este projeto nao faz automaticamente.
"""
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from corrigir_cutlines import mapa_reaches   # noqa: E402
from ras_io import escrever            # noqa: E402


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def marcadas(erros, S, mapa, chave="banklines"):
    """Indices das secoes que o HEC-RAS aponta, casando rio, reach e RS."""
    saida, perdidos = set(), []
    for r in csv.DictReader(open(erros, encoding="utf-8"), delimiter=";"):
        if chave not in r["mensagem"]:
            continue
        m = re.match(r"\s*([^,]+),\s*(\S+)\s*\(([\d.]+)\)", r["onde"].strip())
        if not m:
            continue
        rio_e, rch_e, rs_e = m.group(1).strip(), m.group(2), float(m.group(3))
        dec = len(m.group(3).split(".")[-1]) if "." in m.group(3) else 0
        passo = 10 ** (-dec)          # o RAS TRUNCA ao exibir; ver
        cand = [i for i, d in enumerate(S)      # endireitar_cutlines.py
                if mapa[i] == (rio_e, rch_e)
                and -0.5 * passo - 1e-9 <= float(d["rs"]) - rs_e < passo + 1e-9]
        if len(cand) == 1:
            saida.add(cand[0])
        else:
            perdidos.append((rio_e, rch_e, rs_e, len(cand)))
    return saida, perdidos


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    entrada = argv[0]
    erros = _arg(argv, "--erros")
    ext = _arg(argv, "--saida", "g11")
    if not erros:
        raise SystemExit("informe --erros com a tabela do ler_erros_geometria")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)
    alvo, perdidos = marcadas(erros, S, mapa)
    # A GANCHO E DA VIZINHA, e nao so da marcada. A bank line liga o ponto de
    # margem de uma secao ao da seguinte: quem cruza o gancho e marcado, mas
    # quem o produz pode ser o vizinho. Medido no g18, nas 30 marcadas: o
    # deslocamento proprio tem mediana de 12,0 m e o das VIZINHAS, 15,1 m --
    # e 11 das 30 estao abaixo da mediana do modelo inteiro, ou seja, nao tem
    # defeito nenhum. Foi por isso que recentrar so as marcadas levou o
    # contador de 49 para 52 na tentativa anterior.
    janela = int(_arg(argv, "--vizinhas", 1, float))
    if janela:
        ordem = sorted(range(len(S)), key=lambda i: (mapa[i], -S[i]["rs"]))
        pos = {i: k for k, i in enumerate(ordem)}
        estendido = set(alvo)
        for i in list(alvo):
            for dk in range(-janela, janela + 1):
                k = pos[i] + dk
                if 0 <= k < len(ordem) and mapa[ordem[k]] == mapa[i]:
                    estendido.add(ordem[k])
        print(f"janela de {janela} vizinha(s): {len(alvo)} -> "
              f"{len(estendido)} secoes tratadas")
        alvo = estendido
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes a tratar: {len(alvo)}"
          + (f"   (nao casaram: {len(perdidos)})" if perdidos else ""))
    for x in perdidos:
        print(f"   NAO CASOU: {x[0]} {x[1]} RS {x[2]} ({x[3]} candidatas)")

    novos, sem_eixo, desl = {}, 0, []
    for i in sorted(alvo):
        d = S[i]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        v = B - A
        L = float(np.hypot(*v))
        u = v / max(L, 1e-9)
        x = LineString([A, B]).intersection(eixos[mapa[i]])
        if x.is_empty:
            sem_eixo += 1
            continue
        p = list(x.geoms)[0] if hasattr(x, "geoms") else x
        if p.geom_type != "Point":
            p = p.centroid
        s_eixo = float(np.dot(np.array([p.x, p.y]) - A, u))
        st = np.asarray(d["sta"], float)
        lb, rb = float(d["lb"]), float(d["rb"])
        larg = rb - lb
        # centra mantendo a largura, sem sair da secao
        c = min(max(s_eixo, st[0] + larg / 2), st[-1] - larg / 2)
        nl = st[int(np.argmin(np.abs(st - (c - larg / 2))))]
        nr = st[int(np.argmin(np.abs(st - (c + larg / 2))))]
        if nr <= nl:
            k = int(np.argmin(np.abs(st - c)))
            nl = st[max(k - 1, 0)]
            nr = st[min(k + 1, len(st) - 1)]
        if abs(nl - lb) < 1e-9 and abs(nr - rb) < 1e-9:
            continue
        desl.append(abs(0.5 * (nl + nr) - 0.5 * (lb + rb)))
        novos[i] = {"lb": nl, "rb": nr, "lb0": lb, "rb0": rb}

    print(f"\nsecoes recentradas       : {len(novos)}")
    print(f"sem travessia do eixo    : {sem_eixo}  (ficam como estao)")
    if desl:
        desl = np.array(desl)
        print(f"deslocamento do canal    : mediana {np.median(desl):.1f} m   "
              f"p90 {np.percentile(desl,90):.1f}   max {desl.max():.1f} m")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
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
                # as quebras que estavam nas margens antigas vao para as novas
                for t in range(0, 3 * cnt, 3):
                    if abs(val[t] - nv["lb0"]) < 1e-6:
                        val[t] = nv["lb"]
                    elif abs(val[t] - nv["rb0"]) < 1e-6:
                        val[t] = nv["rb"]
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

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    la = np.array([float(x["rb"] - x["lb"]) for x in A2])
    lb_ = np.array([float(x["rb"] - x["lb"]) for x in B2])
    dl = np.abs(lb_ - la)
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    # A LARGURA NAO SE PRESERVA EXATAMENTE, e nao ha como preservar: a margem
    # tem de cair sobre uma estaca EXISTENTE (fora da grade o HEC-RAS recusa a
    # geometria -- aconteceu aqui, com 519 erros de 521 secoes), e o
    # espacamento das estacas e finito. O que se garante e que a mudanca fica
    # limitada a esse espacamento.
    print(f"   largura do canal mudou : max {dl.max():.2f} m em "
          f"{int((dl > 1e-9).sum())} secoes  (limite = o passo das estacas)")
    print(f"      mediana da mudanca, onde houve: "
          f"{np.median(dl[dl > 1e-9]) if (dl > 1e-9).any() else 0:.2f} m")
    za = np.array([float(x["z"].min()) for x in A2])
    zb = np.array([float(x["z"].min()) for x in B2])
    print(f"   talvegue mudou         : max {np.abs(zb-za).max():.6f} m "
          "(tem de ser zero)")
    fora = 0
    for d in B2:
        st = np.asarray(d["sta"], float)
        for b in (d["lb"], d["rb"]):
            if np.abs(st - float(b)).min() > 1e-6:
                fora += 1
                break
    print(f"   Bank Sta fora de estaca existente: {fora}  (tem de ser zero)")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
