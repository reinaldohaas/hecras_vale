# -*- coding: utf-8 -*-
"""Cria a juncao que liga reaches cujos eixos JA se encontram num ponto.

    python scripts/criar_juncao.py taha_ai_novo/taha_ai.g01 --saida g19 \
        --nome Foz_Rio_do_Sul \
        --up "Itajai_Oeste,R4" --up "Itajai_Sul,R1" --dn "Itajai_Acu,R1"

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO QUE ISTO CORRIGE

  O construtor da rede fez o snap dos eixos -- Oeste R4, Sul R1 e Acu R1
  terminam/nascem no MESMO ponto, a 0,0 m um do outro -- e as ultimas
  secoes ja tem comprimento zero, mas o bloco `Junct Name=` nunca foi
  escrito. Sem ele o HEC-RAS ve dois reaches pendurados e recusa o unsteady:
  "River Itajai_Oeste Reach R4 needs an downstream boundary condition."
  Sem a juncao, alem do erro, TODA A AGUA do Sul e do Oeste (~1.900 m3/s de
  picos somados) simplesmente nao entrava no Acu.

O QUE SE FAZ

  So se INSERE o bloco da juncao, no padrao das nove irmas do arquivo
  (Junc L&A=1.00 para reach que termina em RS 75, como Foz_Testo e
  Foz_Itajai_Norte fazem). RECUSA-SE a criar se os extremos dos eixos nao
  coincidirem a menos de `TOL_M` metros -- juncao por nome com eixos
  separados e mapa mentiroso.

  O `.u01` NAO e tocado aqui: se o reach de jusante tinha hidrograma de
  vazao na cabeceira, ele passa a ser invalido (abaixo de juncao) e tem de
  ser tratado a parte -- o script avisa.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_geometria import ler_eixos     # noqa: E402
from ras_io import escrever            # noqa: E402

TOL_M = 1.0


def _lista(argv, chave):
    out = []
    for i, a in enumerate(argv):
        if a == chave:
            rio, rch = argv[i + 1].split(",")
            out.append((rio.strip(), rch.strip()))
    return out


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g19"
    nome = argv[argv.index("--nome") + 1] if "--nome" in argv else None
    ups = _lista(argv, "--up")
    dns = _lista(argv, "--dn")
    if not (nome and ups and len(dns) == 1):
        raise SystemExit("preciso de --nome, ao menos um --up e um --dn")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    eixos = ler_eixos(entrada)
    fim = {}
    for ch in ups:
        if ch not in eixos:
            raise SystemExit(f"reach {ch} nao existe")
        fim[ch] = np.asarray(eixos[ch].coords, float)[-1]
    dn = dns[0]
    if dn not in eixos:
        raise SystemExit(f"reach {dn} nao existe")
    p0 = np.asarray(eixos[dn].coords, float)[0]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")
    for ch, p in fim.items():
        d = float(np.hypot(*(p - p0)))
        print(f"   {ch[0]:14s} {ch[1]:3s} termina a {d:6.2f} m do inicio "
              f"de {dn[0]} {dn[1]}")
        if d > TOL_M:
            raise SystemExit(f"   extremo de {ch} a {d:.2f} m > {TOL_M} m "
                             "-- snap primeiro, juncao depois")

    t = open(entrada, encoding="latin-1", errors="replace").read()
    fim_de_linha = "\r\n" if "\r\n" in t[:2000] else "\n"
    linhas = t.replace("\r\n", "\n").replace("\r", "\n").split("\n")

    if any(l.startswith("Junct Name=") and nome in l for l in linhas):
        raise SystemExit(f"juncao {nome} ja existe")

    bloco = [f"Junct Name={nome:<16s}",
             "Junct Desc=Confluencia, 0 , 0 , 0 ,0",
             "Junct X Y & Text X Y=%.2f,%.2f,%.2f,%.2f"
             % (p0[0], p0[1], p0[0] + 800.0, p0[1] + 800.0)]
    for rio, rch in ups:
        bloco.append(f"Up River,Reach={rio:<16s},{rch:<16s}")
    bloco.append(f"Dn River,Reach={dn[0]:<16s},{dn[1]:<16s}")
    for _ in ups:
        bloco.append("Junc L&A=1.00,0")
    bloco.append("")

    # insere antes do primeiro "River Reach=" (depois das juncoes irmas)
    idx = next(i for i, l in enumerate(linhas)
               if l.startswith("River Reach="))
    saida = linhas[:idx] + bloco + linhas[idx:]
    escrever(novo, "\n".join(saida))

    # -------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    t2 = open(novo, encoding="latin-1", errors="replace").read()
    n_j = t2.count("Junct Name=")
    print(f"   juncoes: {t.count('Junct Name=')} -> {n_j}")
    from qc_secoes import ler_secoes
    print(f"   secoes : {len(ler_secoes(entrada))} -> "
          f"{len(ler_secoes(novo))}   (nao pode mudar)")
    print(f"\nATENCAO: se {dn[0]} {dn[1]} tem 'Flow Hydrograph' de cabeceira "
          "no .u01, ele ficou invalido (reach abaixo de juncao) -- trate no "
          "fluxo antes de rodar.")


if __name__ == "__main__":
    main(sys.argv[1:])
