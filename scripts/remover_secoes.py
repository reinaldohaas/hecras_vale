# -*- coding: utf-8 -*-
"""Remove secoes NOMEADAS, somando o comprimento delas a vizinha de montante.

    python scripts/remover_secoes.py taha_ai_novo/taha_ai.g06 --saida g08 \
        --alvo "Itajai_Mirim,R1,142991.95" --alvo "Taio,R1,44088.99"

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

QUANDO REMOVER E O CONSERTO CERTO

  Num gancho de meandro ha cutlines INTERCALADAS: a de RS maior e cruzada
  pelo eixo DEPOIS da de RS menor, e o Validate Geometry acusa "River
  Station out of order". Nao existe eixo que cruze as duas na ordem -- foi
  tentado, por emenda pelo centro do canal e por aglomerado, e cada emenda
  empurrava o defeito para a vizinha. A secao intercalada nao tem como
  ficar; o conjunto minimo a sair e a MAIOR SUBSEQUENCIA CRESCENTE das
  estacoes de cruzamento (o que nao esta nela, sai).

O QUE MUDA, E O QUE NAO

  O bloco inteiro da secao sai do arquivo. O `Type RM Length L Ch R` da
  secao DE MONTANTE (a anterior no mesmo reach) recebe a soma dos proprios
  comprimentos com os da removida -- o rio nao encurta. Nada mais muda.

  RECUSA-SE a remover: a primeira ou a ultima secao de um reach (sao a
  conexao com contorno ou juncao) e secao citada no `.u01` (enderecos de
  condicao de contorno quebrariam).
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever            # noqa: E402


def _alvos(argv):
    out = []
    for i, a in enumerate(argv):
        if a == "--alvo":
            rio, rch, rs = argv[i + 1].split(",")
            out.append((rio.strip(), rch.strip(), float(rs)))
    return out


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g08"
    alvos = _alvos(argv)
    if not alvos:
        raise SystemExit("preciso de ao menos um --alvo Rio,Reach,RS")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    # protegidos: RS citados no .u01
    u01 = os.path.join(raiz, base + ".u01")
    protegidos = set()
    if os.path.exists(u01):
        for l in open(u01, encoding="latin-1", errors="replace"):
            if re.match(r"^(Boundary Location|Lateral Inflow Location)", l):
                protegidos |= {x for x in re.findall(r"[-\d.]+",
                                                     l.split("=", 1)[1])}

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")

    # mapa: bloco de cada secao, com reach e RS
    secoes = []          # (i_ini, reach, rs)
    ch = None
    for i, l in enumerate(linhas):
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            ch = (p[0].strip(), p[1].strip())
        elif l.startswith("Type RM Length L Ch R"):
            corpo = l.split("=", 1)[1].split(",")
            secoes.append({"i": i, "ch": ch, "rs": float(corpo[1]),
                           "len": [float(x) for x in corpo[2:5]]})
    fim = [s["i"] for s in secoes[1:]] + [len(linhas)]
    # o bloco da ultima secao de um reach termina antes do proximo
    # "River Reach="; o corte em `fim` pelo inicio da secao seguinte
    # incluiria o cabecalho do reach seguinte -- encolhe ate ele
    for k, s in enumerate(secoes):
        j = s["i"]
        f = fim[k]
        for t in range(j + 1, f):
            if linhas[t].startswith("River Reach="):
                f = t
                break
        s["fim"] = f

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    remover = {}
    for rio, rch, rs in alvos:
        cand = [k for k, s in enumerate(secoes)
                if s["ch"] == (rio, rch) and abs(s["rs"] - rs) < 0.05]
        if len(cand) != 1:
            raise SystemExit(f"alvo {rio},{rch},{rs}: "
                             f"{len(cand)} candidatas -- recusado")
        k = cand[0]
        s = secoes[k]
        no_reach = [t for t, x in enumerate(secoes) if x["ch"] == s["ch"]]
        if k == no_reach[0] or k == no_reach[-1]:
            raise SystemExit(f"{rio} RS {rs} e extremo do reach -- recusado")
        if any(abs(float(p) - s["rs"]) < 0.05 for p in protegidos):
            raise SystemExit(f"{rio} RS {rs} esta no .u01 -- recusado")
        remover[k] = s
        print(f"   sai {rio:14s} {rch:3s} RS {s['rs']:10.2f}   "
              f"(comprimentos {s['len'][1]:.2f} m somados a montante)")

    # soma comprimento na secao de montante (a anterior no mesmo reach)
    for k, s in sorted(remover.items()):
        ant = max(t for t in range(k) if secoes[t]["ch"] == s["ch"]
                  and t not in remover)
        a = secoes[ant]
        soma = [x + y for x, y in zip(a["len"], s["len"])]
        a["len"] = soma
        l = linhas[a["i"]]
        cab = l.split("=", 1)[1].split(",")
        linhas[a["i"]] = ("Type RM Length L Ch R = %s,%s,%.2f,%.2f,%.2f"
                          % (cab[0].strip(), cab[1].strip(), *soma))

    saida = []
    i = 0
    blocos = {s["i"]: s["fim"] for s in remover.values()}
    while i < len(linhas):
        if i in blocos:
            i = blocos[i]
            continue
        saida.append(linhas[i])
        i += 1
    escrever(novo, "\n".join(saida))

    # -------------------------------------------------------- conferencia
    from qc_secoes import ler_secoes
    import numpy as np
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    print(f"   secoes: {len(A2)} -> {len(B2)}   "
          f"(esperado {len(A2) - len(remover)})")
    print("   comprimento total (soma dos Ch): "
          f"{sum(d['len_ch'] for d in A2):.2f} -> "
          f"{sum(d['len_ch'] for d in B2):.2f}   (nao pode mudar)")


if __name__ == "__main__":
    main(sys.argv[1:])
