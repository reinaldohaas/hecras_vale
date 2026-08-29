# -*- coding: utf-8 -*-
"""Afina a geometria: espacamento minimo em LARGURAS DE CANAL, nao em metros.

    python scripts/afinar_secoes.py modelo/so_mirim.g08 --saida g09 --minimo 1.5

A ENTRADA NAO E TOCADA. Sai um .gXX novo, com menos secoes.

POR QUE EM LARGURAS, E NAO EM METROS

  O criterio util nao e uma distancia fixa: e a distancia comparada a escala do
  rio. A guia usual do HEC-RAS e no minimo UMA largura de canal entre secoes, e
  cerca de duas em curva. Medido no g08:

      dx entre secoes : p10 25,0   mediana 47,3   p90 251,3 m
      largura do canal: mediana 52,6 m
      dx / largura    : mediana 0,89  --  54% dos vaos abaixo de UMA largura

  E essa densidade nao vem de exigencia numerica. Com o dt de 15 s do plano e
  as velocidades da propria rodada (0,31 a 0,92 m/s), o Courant fica entre 0,09
  e 0,39: a onda anda 4,8 m por passo e as secoes estao a 47 m. Ela veio do
  criterio de Samuels, que ja foi desligado no padrao -- mas o g08 herdou a
  densidade de quando ele estava ligado.

  O custo esta medido: no g08 o solver bate no teto de 40 iteracoes em 90% dos
  passos, gastando 36,3 iteracoes por passo contra as 2 de um passo que
  converge. Cada iteracao varre as 1418 secoes.

O QUE SE PRESERVA

  AS SECOES CORTADAS DO TERRENO NUNCA SAO APAGADAS. Elas sao amostra real; as
  interpoladas sao preenchimento. So interpolada perto demais sai. A primeira e
  a ultima tambem ficam sempre -- sao a conexao com o contorno.

  E OS COMPRIMENTOS DE TRECHO SAO REFEITOS. `Type RM Length L Ch R` guarda a
  distancia desta secao ate a PROXIMA. Apagar uma sem somar o vao deixaria o
  rio mais curto do que e. Conferido neste arquivo: LOB == Ch == ROB e
  Ch == diferenca de RS com erro maximo de 0,000000 m, entao os novos
  comprimentos saem do RS das que ficam, e a ultima recebe zero.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes   # noqa: E402
from ras_io import escrever        # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g09"
    K = float(argv[argv.index("--minimo") + 1]) if "--minimo" in argv else 1.5
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    ordem = sorted(range(len(S)), key=lambda i: -S[i]["rs"])
    import pickle
    est = pickle.load(open(os.path.join(raiz, f"estado_{base}.pkl"), "rb"))
    chave = next(iter(est["xs_pronto"]))
    real = {round(float(x["rs"]), 2) for x in est["xs_pronto"][chave]
            if not x.get("interpolada")}

    rs = np.array([S[i]["rs"] for i in ordem])
    lc = np.array([float(S[i]["rb"] - S[i]["lb"]) for i in ordem])
    ereal = np.array([round(r, 2) in real for r in rs])

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"minimo : {K:.2f} x a largura do canal")
    print(f"secoes : {len(S)}   reais {int(ereal.sum())}   "
          f"interpoladas {int((~ereal).sum())}")

    # ---- quem fica
    fica = np.zeros(len(ordem), bool)
    fica[0] = fica[-1] = True
    ult = 0
    for k in range(1, len(ordem) - 1):
        if ereal[k] or (rs[ult] - rs[k]) >= K * 0.5 * (lc[k] + lc[ult]):
            fica[k] = True
            ult = k
    mantidos = [ordem[k] for k in np.flatnonzero(fica)]
    rs_f = rs[fica]
    dx = -np.diff(rs_f)
    lcf = lc[fica]
    r = dx / np.maximum(0.5 * (lcf[:-1] + lcf[1:]), 1e-9)
    print(f"\nficam {len(mantidos)} secoes  (saem {len(S)-len(mantidos)})   "
          f"reais preservadas: {int(ereal[fica].sum())} de {int(ereal.sum())}")
    print(f"   dx      : p10 {np.percentile(dx,10):6.1f}  mediana "
          f"{np.median(dx):6.1f}  p90 {np.percentile(dx,90):6.1f} m")
    print(f"   dx/larg : p10 {np.percentile(r,10):6.2f}  mediana "
          f"{np.median(r):6.2f}   abaixo de 1: {int((r<1).sum())}")

    # ---- comprimentos de trecho novos, do RS de quem ficou
    comp = {}
    for a, b in zip(mantidos, mantidos[1:]):
        comp[a] = S[a]["rs"] - S[b]["rs"]
    comp[mantidos[-1]] = 0.0

    # ---- reescreve mantendo so os blocos das secoes que ficam
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    ini = [i for i, l in enumerate(linhas) if l.startswith("Type RM Length L Ch R")]
    fim = ini[1:] + [len(linhas)]
    manter = set(mantidos)
    saida = linhas[:ini[0]]
    for k, (a, b) in enumerate(zip(ini, fim)):
        if k not in manter:
            continue
        bloco = list(linhas[a:b])
        c = comp[k]
        bloco[0] = re.sub(
            r"^(Type RM Length L Ch R\s*=\s*[^,]+,[^,]+),.*$",
            r"\1,%8.2f,%8.2f,%8.2f" % (c, c, c), bloco[0])
        saida += bloco
    txt = "\n".join(saida)
    t0 = linhas[0].split("=", 1)[1] if "=" in linhas[0] else ""
    if t0:
        txt = txt.replace("Geom Title=" + t0,
                          "Geom Title=" + t0 + " + afinado", 1)
    escrever(novo, txt)

    # ---- conferencia
    B = ler_secoes(novo)
    B.sort(key=lambda d: -d["rs"])
    rsb = np.array([d["rs"] for d in B])
    chb = np.array([d["len_ch"] for d in B])
    err = np.abs(chb[:-1] - (-np.diff(rsb))).max() if len(B) > 1 else 0.0
    print("\nCONFERENCIA")
    print(f"   secoes no arquivo novo      : {len(B)}")
    print(f"   Ch == diferenca de RS       : erro maximo {err:.6f} m")
    print(f"   ultima secao com Ch=0       : {chb[-1] == 0.0}")
    print(f"   extensao do rio             : {rs[0]-rs[-1]:.1f} m -> "
          f"{rsb[0]-rsb[-1]:.1f} m  (tem de ser igual)")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
