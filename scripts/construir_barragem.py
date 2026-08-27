# -*- coding: utf-8 -*-
"""Constroi uma BARRAGEM (estrutura inline) no g01, com vertedouro e fenda.

    python scripts/construir_barragem.py taha_ai.g01 --saida g99 \
        --rio Itajai_Sul --reach R1 --rs 33500 \
        --crista 399.0 --topo 402.0 --larg-vertedouro 100 \
        --fenda 1.3 --nome "Barragem Sul (Ituporanga)"

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

COMO A BARRAGEM E ESCRITA

  Bloco `Type RM Length L Ch R = 3` entre as duas secoes vizinhas de
  `--rs`, com perfil de crista (#Inline Weir SE):

    ombreiras na cota `--topo` ate as bordas da secao de montante,
    vertedouro na cota `--crista` com `--larg-vertedouro` m centrado no
    talvegue, e uma FENDA de `--fenda` m ate o leito -- aproximacao fixa
    dos condutos de fundo sempre abertos das barragens de contencao do
    Alto Vale (Oeste: 7 condutos, 163 m3/s; Sul: 194 m3/s -- JICA 2011,
    Anexo A, figs. 7.5.5-7.5.9). Largura da fenda calibra a capacidade.

  Sem comportas: barragem seca de soleira livre, como operavam em 1983.

  CONFERENCIA relendo o gravado: bloco presente, perfil fecha com as
  cotas pedidas, secoes inalteradas em numero.
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
    ext = _arg(argv, "--saida", "g99")
    rio = _arg(argv, "--rio")
    reach = _arg(argv, "--reach")
    rs = _arg(argv, "--rs", None, float)
    crista = _arg(argv, "--crista", None, float)
    topo = _arg(argv, "--topo", None, float)
    larg_v = _arg(argv, "--larg-vertedouro", 100.0, float)
    fenda = _arg(argv, "--fenda", 1.3, float)
    rampa = _arg(argv, "--rampa", 0.0, float)
    coef = _arg(argv, "--coef", 1.7, float)
    nome = _arg(argv, "--nome", "Barragem")
    if None in (rio, reach, rs, crista, topo):
        raise SystemExit("faltou --rio/--reach/--rs/--crista/--topo")

    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    monte = jusante = None
    for d in S:
        if d["rio"] != rio or d["reach"] != reach:
            continue
        if d["rs"] > rs and (monte is None or d["rs"] < monte["rs"]):
            monte = d
        if d["rs"] < rs and (jusante is None or d["rs"] > jusante["rs"]):
            jusante = d
    if monte is None or jusante is None:
        raise SystemExit(f"nao ha secoes vizinhas de RS {rs} em {rio} {reach}")

    st = np.asarray(monte["sta"], float)
    z = np.asarray(monte["z"], float)
    m = (st >= monte["lb"] - 1e-6) & (st <= monte["rb"] + 1e-6)
    ct = float(st[m][int(np.argmin(z[m]))])      # estacao do talvegue
    leito = float(z[m].min())
    s0, s1 = float(st[0]), float(st[-1])
    v0 = max(s0 + 1.0, ct - larg_v / 2)
    v1 = min(s1 - 1.0, ct + larg_v / 2)
    f0, f1 = ct - fenda / 2, ct + fenda / 2
    dist = (monte["rs"] - rs) * 0.9              # m ate a secao de montante

    # degrau vertical = estacao repetida, como no BaldEagle.g01 oficial.
    # --rampa > 0: vertedouro em V raso (centro `rampa` m abaixo das
    # pontas) para o vertimento ENGAJAR gradualmente -- crista plana de
    # 100 m engajando de uma vez e um choque numerico no pe
    if rampa > 0:
        perfil = [(s0, topo), (v0, topo), (v0, crista),
                  (ct - 10, crista - rampa), (f0, crista - rampa),
                  (f0, leito), (f1, leito), (f1, crista - rampa),
                  (ct + 10, crista - rampa), (v1, crista), (v1, topo),
                  (s1, topo)]
    else:
        perfil = [(s0, topo), (v0, topo), (v0, crista),
                  (f0, crista), (f0, leito), (f1, leito),
                  (f1, crista), (v1, crista), (v1, topo),
                  (s1, topo)]
    print(f"{nome}: {rio} {reach} RS {rs:.0f}")
    print(f"   entre RS {monte['rs']:.1f} e {jusante['rs']:.1f}")
    print(f"   leito {leito:.2f}  fenda {fenda:.2f} m  vertedouro "
          f"{crista:.2f} m x {v1-v0:.0f} m  topo {topo:.2f} m")

    # ordem e grafia copiadas do BaldEagle.g01 (exemplo oficial HEC):
    # Type 5, comprimentos em branco, SE antes da linha de parametros,
    # cabecalho "IW Dist,..." SEM '=' e valores na linha seguinte
    bloco = ["Type RM Length L Ch R = 5 ,%-8s,,," % (f"{rs:.0f}"),
             "BEGIN DESCRIPTION:",
             nome + " -- JICA 2011 Anexo A",
             "END DESCRIPTION:",
             "IW Pilot Flow=0",
             "#Inline Weir SE= %d " % len(perfil)]
    v = []
    for a, b in perfil:
        v += [a, b]
    bloco += _col(v, 8, 2)
    bloco += ["IW Dist,WD,Coef,Skew,MaxSub,Min_El,Is_Ogee,SpillHt,DesHd",
              "%g,%g,%g,0,0.95,, 0 ,,,," % (max(dist, 10.0), 10.0, coef)]

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    saida, i = [], 0
    rio_c = reach_c = None
    inserido = False
    while i < len(linhas):
        l = linhas[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio_c, reach_c = p[0].strip(), p[1].strip()
        if (not inserido and rio_c == rio and reach_c == reach
                and l.startswith("Type RM Length L Ch R =")):
            try:
                rs_l = float(l.split("=", 1)[1].split(",")[1])
            except (ValueError, IndexError):
                rs_l = None
            if rs_l is not None and rs_l < rs:
                saida += bloco + [""]
                inserido = True
        saida.append(l)
        i += 1
    if not inserido:
        raise SystemExit("nao achei onde inserir -- nada gravado")
    escrever(novo, "\n".join(saida))

    print(f"\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    t = open(novo, encoding="latin-1", errors="replace").read()
    tem = "#Inline Weir SE=" in t and nome in t
    print(f"   bloco da barragem presente: {tem}   (tem de ser True)")
    print(f"   saida: {novo}")


if __name__ == "__main__":
    main(sys.argv[1:])
