# -*- coding: utf-8 -*-
"""Abre canal piloto em secao de canal MORTO-PLANO.

    python scripts/abrir_piloto.py taha_ai_novo/taha_ai.g01 --saida g27

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

POR QUE

  Secao cujo canal (lb..rb) e uma chapa plana nao conduz vazao baixa: a
  lamina de centimetros espalhada por centenas de metros tem condutancia
  quase nula, a secao vira barragem e o solver poca dezenas de metros a
  montante (medido: 45 m de poco atras do Benedito RS 37482, aplainado
  pelo aterro do limitar_declividade). O gerador do modelo especifica
  canal piloto de 25 m x 1,5 m (opcoes: pilot_largura, pilot_prof) --
  aqui ele e devolvido a quem o perdeu.

O QUE SE FAZ

  Secao com desnivel de canal < `--relevo` m recebe um entalhe retangular
  de `--largura` x `--prof` no CENTRO do canal, com paredes nas estacas
  existentes mais proximas. HTab acompanha o fundo novo. Nada mais muda.
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
    ext = _arg(argv, "--saida", "g27")
    so_rios = _arg(argv, "--rios", None)
    so_rios = set(so_rios.split(",")) if so_rios else None
    relevo = _arg(argv, "--relevo", 0.5, float)
    largura = _arg(argv, "--largura", 25.0, float)
    prof = _arg(argv, "--prof", 1.5, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    novos = {}
    for i, d in enumerate(S):
        if so_rios and d["rio"] not in so_rios:
            continue
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
        if m.sum() < 2 or (z[m].max() - z[m].min()) >= relevo:
            continue
        centro = 0.5 * (d["lb"] + d["rb"])
        s0, s1 = centro - largura / 2, centro + largura / 2
        fundo = float(z[m].min()) - prof
        z2 = z.copy()
        dentro = (st >= s0) & (st <= s1) & m
        if dentro.sum() < 2:
            # garante ao menos os pontos das bordas do entalhe
            dentro = m & (np.abs(st - centro) <= largura)
        z2[dentro] = fundo
        novos[(d["rio"], d["reach"], round(d["rs"], 2))] = {"sta": st, "z": z2, "htab": fundo + 0.15}
        print(f"   {d['rio']:13s} {d['reach']:3s} RS {d['rs']:9.1f}  "
              f"piloto {largura:.0f}x{prof:.1f} m aberto "
              f"(fundo {fundo:.2f})")

    if not novos:
        print("nenhum canal morto-plano")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
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

    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    resta = 0
    for d in B:
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
        if m.sum() >= 2 and (z[m].max() - z[m].min()) < relevo:
            resta += 1
    print(f"   canais mortos-planos restantes: {resta}   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
