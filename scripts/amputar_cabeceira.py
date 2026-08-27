# -*- coding: utf-8 -*-
"""Amputa a cabeceira torrencial de um reach; o contorno desce junto.

    python scripts/amputar_cabeceira.py taha_ai_novo/taha_ai.g01 \
        --reach "Benedito,R1" --rs-corte 28266.24 --saida g31

A GEOMETRIA DE ENTRADA NAO E TOCADA (sai .gXX novo). O .uNN do projeto E
EDITADO (com backup .antes_da_amputacao): o Boundary Location e o Initial
Flow Loc do reach passam para a nova secao de topo.

POR QUE AMPUTAR

  O alto Benedito e uma garganta com 3 a 9% de declividade por 17 km. Em
  19 rodadas medidas, toda configuracao numerica morreu ali: o 1D
  unsteady nao sustenta vazao baixa em corredeira sustentada, e o poco
  que se forma contamina a rede. A garganta e desabitada e so acrescenta
  minutos de defasagem ao hidrograma; a mancha de cheia que importa e a
  do vale. O contorno de vazao passa a entrar na frente da montanha --
  decisao de projeto padrao em modelo de varzea.

O QUE SE FAZ

  Todas as secoes do reach ACIMA de --rs-corte saem; a polilinha do eixo
  fica (so cartografia). No .uNN, o RS do Boundary Location e do Initial
  Flow Loc do reach viram o RS da nova secao de topo. Lateral que comece
  acima do corte tem o inicio movido para a nova secao de topo.
"""
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from ras_io import escrever            # noqa: E402


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    rio, rch = [x.strip() for x in _arg(argv, "--reach").split(",")]
    corte = _arg(argv, "--rs-corte", None, float)
    ext = _arg(argv, "--saida", "g31")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    # blocos de secao com reach
    secoes, ch = [], None
    for i, l in enumerate(linhas):
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            ch = (p[0].strip(), p[1].strip())
        elif l.startswith("Type RM Length L Ch R"):
            rs = float(l.split("=", 1)[1].split(",")[1])
            secoes.append({"i": i, "ch": ch, "rs": rs})
    fim = [s["i"] for s in secoes[1:]] + [len(linhas)]
    for k, s in enumerate(secoes):
        f = fim[k]
        for t in range(s["i"] + 1, f):
            if linhas[t].startswith("River Reach="):
                f = t
                break
        s["fim"] = f

    alvo = [s for s in secoes if s["ch"] == (rio, rch) and s["rs"] > corte]
    resto = [s for s in secoes if s["ch"] == (rio, rch) and s["rs"] <= corte]
    if not alvo:
        raise SystemExit("nada acima do corte")
    if len(resto) < 2:
        raise SystemExit("sobrariam menos de 2 secoes -- recusado")
    novo_topo = max(s["rs"] for s in resto)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"   {rio} {rch}: saem {len(alvo)} secoes acima de RS {corte}")
    print(f"   novo topo: RS {novo_topo}")

    blocos = {s["i"]: s["fim"] for s in alvo}
    saida, i = [], 0
    while i < len(linhas):
        if i in blocos:
            i = blocos[i]
            continue
        saida.append(linhas[i])
        i += 1
    escrever(novo, "\n".join(saida))

    # ------------------------------------------------------------- u01
    u01 = os.path.join(raiz, base + ".u01")
    shutil.copy2(u01, u01 + ".antes_da_amputacao")
    t = open(u01, encoding="latin-1", errors="replace").read()
    rs_txt = ("%.2f" % novo_topo).rstrip("0").rstrip(".")

    def bl(m):
        return (m.group(0) if float(m.group(3)) <= corte else
                "Boundary Location=%s,%s,%-8s,%s" %
                (m.group(1), m.group(2), rs_txt, m.group(4)))
    t = re.sub(r"Boundary Location=(%s\s*),(%s\s*),([\d.]+)\s*,(.*)"
               % (re.escape(rio), re.escape(rch)), bl, t)

    def ifl(m):
        return ("Initial Flow Loc=%s,%s,%-8s,%s" %
                (m.group(1), m.group(2), rs_txt, m.group(4))
                if float(m.group(3)) > corte else m.group(0))
    t = re.sub(r"Initial Flow Loc=(%s\s*),(%s\s*),([\d.]+)\s*,(.*)"
               % (re.escape(rio), re.escape(rch)), ifl, t)

    # o RAS PROIBE lateral comecando na secao de topo do reach: lateral
    # movida para o topo desce para a SEGUNDA secao
    segundo = sorted((s["rs"] for s in resto), reverse=True)[1]
    seg_txt = ("%.2f" % segundo).rstrip("0").rstrip(".")
    partes = re.split(r"(?=^Boundary Location=)", t, flags=re.M)
    for k, b in enumerate(partes):
        if ("Uniform Lateral Inflow" in b
                and re.match(r"Boundary Location=%s\s*,%s\s*,%s\s*,"
                             % (re.escape(rio), re.escape(rch),
                                re.escape(rs_txt)), b)):
            partes[k] = b.replace(
                b.split("\n")[0],
                re.sub(r"^(Boundary Location=[^,]+,[^,]+,)[^,]+"
                       % (), r"\g<1>%-8s" % seg_txt,
                       b.split("\n")[0]), 1)
            print(f"   lateral no topo movida para RS {seg_txt}")
    t = "".join(partes)
    escrever(u01, t)
    print(f"   u01: contorno e inicial de {rio} movidos para RS {rs_txt}"
          f"   (backup .antes_da_amputacao)")

    B = ler_secoes(novo)
    n_r = sum(1 for d in B if d["rio"] == rio and d["reach"] == rch)
    print("\nCONFERENCIA (relendo o gravado)")
    print(f"   secoes de {rio} {rch}: {len(alvo)+len(resto)} -> {n_r}")
    print(f"   secoes totais: {len(secoes)} -> {len(B)}")


if __name__ == "__main__":
    main(sys.argv[1:])
