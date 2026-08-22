# -*- coding: utf-8 -*-
"""Ajusta o .uNN depois que o reach unico virou tres.

    python scripts/reestruturar_fluxo.py modelo/mirim_t30/mirim_t30.u01 \
        --geom modelo/mirim_t30/mirim_t30.g04 --saida u02

A ENTRADA NAO E TOCADA. Sai um .uNN novo.

POR QUE O FLUXO PRECISA MUDAR JUNTO

  As condicoes de contorno sao endereçadas por `rio, reach, RS`. Partir o
  reach nao move nenhuma secao, mas muda o ENDERECO de duas delas:

  1. O contorno de jusante estava em `R1, RS 75.00`. A RS 75 agora pertence a
     R3. Sem trocar o nome do reach o HEC-RAS nao acha o contorno -- e um
     modelo sem contorno de jusante nao roda.

  2. A ENTRADA LATERAL UNIFORME ia de RS 141397,3 a 746,07, ou seja, cobria o
     rio INTEIRO. Um intervalo unico nao pode mais atravessar tres reaches:
     ele e partido em um por reach, com a vazao rateada PELO COMPRIMENTO de
     cada pedaco.

O RATEIO PRESERVA O TOTAL, E ISSO E UMA ESCOLHA

  Os vaos das juncoes (776,61 m somados) nao pertencem a reach nenhum, entao
  nao podem receber entrada lateral. Ha duas maneiras de ratear:

    - preservar a VAZAO TOTAL: cada reach leva L_reach / (soma dos L_reach).
      O total entra inteiro; a taxa por metro sobe 0,55%.
    - preservar a TAXA POR METRO: cada reach leva L_reach / L_original.
      A taxa fica igual; 0,55% da vazao se perde.

  Aqui se preserva o TOTAL, porque o proximo passo e comparar o modelo
  partido com o inteiro, e balanco de volume que nao fecha por 0,55% embaralha
  justamente a medida que se quer. A escolha esta impressa no relatorio.

CONDICAO INICIAL POR REACH

  `Initial RS` existia so para R1. Cada reach precisa da sua, e a vazao
  inicial cresce para jusante junto com a entrada lateral: R1 comeca com a
  vazao de montante, e os de baixo somam a parcela lateral que ja passou.
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


def _col8(v):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%8.2f" % x
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def reaches_da_geom(g):
    """(rio, reach, rs_montante, rs_jusante, [todos os RS]) por reach."""
    t = open(g, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    saida, atual, rss = [], None, []
    for l in t:
        if l.startswith("River Reach="):
            if atual:
                saida.append((atual[0], atual[1], max(rss), min(rss),
                              sorted(rss, reverse=True)))
            p = l.split("=", 1)[1].split(",")
            atual = (p[0].strip(), p[1].strip())
            rss = []
        elif l.startswith("Type RM Length L Ch R"):
            rss.append(float(l.split(",")[1]))
    if atual:
        saida.append((atual[0], atual[1], max(rss), min(rss),
                      sorted(rss, reverse=True)))
    return saida


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    geom = _arg(argv, "--geom")
    ext = _arg(argv, "--saida", "u02")
    if not geom:
        raise SystemExit("informe --geom com a geometria ja partida")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    R = reaches_da_geom(geom)
    rio = R[0][0]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"geometria: {len(R)} reaches")
    for r in R:
        print(f"   {r[1]:<4} RS {r[2]:9.1f} a {r[3]:8.1f}   "
              f"({r[2]-r[3]:9.1f} m)")

    t = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")

    # ---- localiza os blocos de contorno
    ini = [i for i, l in enumerate(t) if l.startswith("Boundary Location=")]
    fim = ini[1:] + [len(t)]
    blocos = [(a, b) for a, b in zip(ini, fim)]

    def rs_do(l):
        p = l.split("=", 1)[1].split(",")
        a = p[2].strip()
        b = p[3].strip()
        return (float(a) if a else None), (float(b) if b else None)

    def reach_de(rs):
        for r in R:
            if r[3] - 1e-6 <= rs <= r[2] + 1e-6:
                return r[1], rs
        # RS FORA DE TODOS OS REACHES. Acontece quando a geometria foi
        # RENUMERADA -- ao trocar 19 km de meandro por 7,5 km de canal, o
        # contorno de jusante deixa de estar em RS 75,00. Sem tratar isso o
        # contorno some e o modelo nao roda. Vai para a secao mais proxima,
        # e o ajuste e impresso, porque mover contorno nao pode ser silencioso.
        melhor = min(((abs(x - rs), r[1], x) for r in R for x in r[4]),
                     key=lambda t: t[0])
        print(f"      RS {rs:.2f} nao existe mais -- levado para "
              f"{melhor[2]:.2f} em {melhor[1]} (a {melhor[0]:.0f} m)")
        return melhor[1], melhor[2]

    saida = list(t[:ini[0]])
    # ---- condicao inicial por reach
    saida = [l for l in saida if not l.startswith("Initial RS=")]

    lat_tot = None
    novos_bc = []
    for a, b in blocos:
        l = t[a]
        rs0, rs1 = rs_do(l)
        corpo = t[a + 1:b]
        tipo = next((x.split("=")[0] for x in corpo
                     if "Hydrograph=" in x), "?")
        if rs1 is None:
            # contorno pontual: corrige o reach e, se houve renumeracao, o RS
            print(f"\n   {tipo:<34} RS {rs0:9.2f}")
            rc, rs_novo = reach_de(rs0)
            p = l.split("=", 1)[1].split(",")
            p[1] = _pad(rc)
            p[2] = f"{rs_novo:<8.2f}"
            novo_l = "Boundary Location=" + ",".join(p)
            novos_bc.append([novo_l] + corpo)
            print(f"      -> reach {rc}, RS {rs_novo:.2f}")
            continue

        # ---- intervalo: rateia por reach
        i_h = next(i for i, x in enumerate(corpo) if "Hydrograph=" in x)
        n = int(corpo[i_h].split("=")[1])
        vals, j = [], i_h + 1
        while len(vals) < n and j < len(corpo):
            x = corpo[j]
            if not x.strip() or re.match(r"^[A-Za-z]", x):
                break
            vals += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                     if x[c:c + 8].strip()]
            j += 1
        vals = np.array(vals[:n])
        cauda = corpo[j:]
        lat_tot = vals.sum()
        print(f"\n   {tipo:<34} RS {rs0:9.2f} a {rs1:8.2f}  "
              f"-- cobre o rio inteiro")
        # A ENTRADA LATERAL NAO PODE TERMINAR NA ULTIMA SECAO DO REACH.
        # O HEC-RAS recusa com "Downstream RS is at the bottom of a reach.
        # Uniform lateral inflows cannot end on the downstream cross section
        # of a reach" -- e o modelo nem chega a computar. O limite de jusante
        # de cada pedaco recua uma secao quando cai na ultima.
        pedacos = []
        for r in R:
            lo = max(r[3], min(rs0, rs1))
            hi = min(r[2], max(rs0, rs1))
            if abs(lo - r[3]) < 1e-6 and len(r[4]) > 1:
                lo = r[4][-2]
            if hi - lo > 1.0:
                pedacos.append((r[1], lo, hi, hi - lo))
        soma = sum(p[3] for p in pedacos)
        print(f"      comprimento util somado: {soma:.1f} m "
              f"(o vao das juncoes nao recebe entrada)")
        for rc, lo, hi, L in pedacos:
            f = L / soma
            v = vals * f
            novos_bc.append(
                [f"Boundary Location={_pad(rio)},{_pad(rc)},"
                 f"{hi:<8.1f},{lo:<8.2f},{_pad('')},{_pad('')}",
                 corpo[i_h - 1] if i_h else "Interval=1HOUR",
                 f"{tipo}= {n} "] + _col8(v) + list(cauda))
            print(f"      {rc}: RS {hi:9.1f} a {lo:8.1f}  "
                  f"{L:9.1f} m  ->  {100*f:5.1f}% da vazao "
                  f"(pico {v.max():.2f} m3/s)")

    # ---- condicoes iniciais
    q0_mont = 2.19
    ic = []
    acum = 0.0
    for r in R:
        rsi = r[2]
        ic.append(f"Initial RS={_pad(rio)},{_pad(r[1])},"
                  f"{rsi:.0f},{q0_mont + acum:.2f}")
        for bc in novos_bc:
            if f",{_pad(r[1])}," in bc[0] and "Uniform Lateral" in " ".join(bc[:3]):
                i_h = next(i for i, x in enumerate(bc) if "Hydrograph=" in x)
                acum += float(bc[i_h + 1][:8])
    print("\n   condicoes iniciais:")
    for x in ic:
        print("      " + x)

    saida += ic + [""]
    for bc in novos_bc:
        saida += bc
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    t2 = open(novo, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    print(f"   contornos: {len(blocos)} -> "
          f"{len(re.findall(r'(?m)^Boundary Location=', t2))}")
    print(f"   Initial RS: {len(re.findall(r'(?m)^Initial RS=', t2))} "
          f"(um por reach)")
    tot = 0.0
    for m in re.finditer(r"(?m)^Uniform Lateral Inflow Hydrograph=\s*(\d+)",
                         t2):
        n = int(m.group(1))
        resto = t2[m.end():].split("\n")[1:]
        v = []
        for x in resto:
            if len(v) >= n or not x.strip() or re.match(r"^[A-Za-z]", x):
                break
            v += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                  if x[c:c + 8].strip()]
        tot += sum(v[:n])
    print(f"   entrada lateral somada: {lat_tot:.2f} -> {tot:.2f} m3/s"
          f"   diferenca {tot-lat_tot:+.4f} (tem de ser zero)")
    # todo contorno aponta para um reach que existe?
    nomes = {r[1] for r in R}
    ruim = [m.group(1).strip() for m in
            re.finditer(r"(?m)^Boundary Location=[^,]+,([^,]+),", t2)
            if m.group(1).strip() not in nomes]
    print(f"   contornos apontando para reach inexistente: "
          f"{len(ruim)}  {ruim if ruim else ''}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
