# -*- coding: utf-8 -*-
"""Encaixa a estaca de margem na grade de estacas do proprio perfil.

    python scripts/encaixar_margens.py modelo/so_mirim.g07 --saida g08

A ENTRADA NAO E TOCADA. Muda EXCLUSIVAMENTE `Bank Sta=` e as quebras do
`#Mann=`, e apenas onde a margem nao caia sobre uma estaca existente.

POR QUE

  O HEC-RAS exige que o valor de `Bank Sta=` COINCIDA com um valor presente no
  bloco `#Sta/Elev`. Nao basta estar dentro da faixa. Quando nao coincide, ele
  recusa a rodar -- sem computar nada -- e escreve em <plano>.data_errors.txt:

      "Left bank station not in station elevation data."
      "Right bank station not in station elevation data."

  Medido no g07: 519 margens esquerdas e 521 direitas fora da grade, com erro
  de ate 2,31 m; no g01 e no g06, zero. O defeito nasceu ao recortar as 541
  secoes: as estacas novas sairam de `linspace(0, L, n)` e as margens foram
  postas em `L/2 +- largura_do_canal/2`, valores que quase nunca caem na grade.

  Esta e a TERCEIRA regra do formato .g01 da mesma familia -- as outras duas
  sao a quebra do Manning ser a estaca de margem, e o HTab ter de acompanhar
  o leito.

O QUE MUDA, E O QUE NAO

  A margem anda para a estaca EXISTENTE mais proxima -- no maximo meio passo
  da grade. Nao se cria estaca nova, nao se mexe em cota, largura, cutline,
  RS, Manning (os valores de n) nem HTab. As quebras do Manning acompanham a
  margem, com a precisao por coluna: estaca %8.2f, n %8.3f, terceiro inteiro.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes  # noqa: E402


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g08"
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    ib = [i for i, l in enumerate(linhas) if l.startswith("Bank Sta=")]
    if len(ib) != len(S):
        raise SystemExit(f"Bank Sta ({len(ib)}) e secoes ({len(S)}) nao batem")

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")

    mud, desloc = {}, []
    for k, d in enumerate(S):
        st = np.asarray(d["sta"], float)
        lb, rb = float(d["lb"]), float(d["rb"])
        il = int(np.argmin(np.abs(st - lb)))
        ir = int(np.argmin(np.abs(st - rb)))
        if ir <= il:                      # nao deixar o canal colapsar
            ir = min(len(st) - 1, il + 1)
        nlb, nrb = float(st[il]), float(st[ir])
        if abs(nlb - lb) < 1e-6 and abs(nrb - rb) < 1e-6:
            continue
        mud[ib[k]] = (nlb, nrb, lb, rb)
        desloc += [abs(nlb - lb), abs(nrb - rb)]

    print(f"secoes com margem encaixada: {len(mud)} de {len(S)}")
    if desloc:
        v = np.array(desloc)
        print(f"   deslocamento da margem: mediana {np.median(v):.3f}  "
              f"p90 {np.percentile(v, 90):.3f}  max {v.max():.3f} m")

    # ---- reescreve Bank Sta e, na mesma secao, as quebras do Manning
    saida, i = [], 0
    while i < len(linhas):
        l = linhas[i]
        if i in mud:
            nlb, nrb, lb, rb = mud[i]
            saida.append(f"Bank Sta={_fmt(nlb)},{_fmt(nrb)}")
            i += 1
            # procura o #Mann= desta secao
            j = i
            while j < len(linhas) and not linhas[j].startswith("#Mann="):
                if linhas[j].startswith("Type RM Length"):
                    j = -1
                    break
                j += 1
            if j > 0 and j < len(linhas):
                while i < j:
                    saida.append(linhas[i]); i += 1
                n = int(linhas[j].split("=")[1].split(",")[0])
                bruto, k2 = [], j + 1
                while k2 < len(linhas) and len(bruto) < 3 * n:
                    x = linhas[k2]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                              if x[c:c + 8].strip()]
                    k2 += 1
                vv = [float(x) for x in bruto[:3 * n]]
                for t in range(0, 3 * n, 3):
                    if abs(vv[t] - lb) < 1e-6:
                        vv[t] = nlb
                    elif abs(vv[t] - rb) < 1e-6:
                        vv[t] = nrb
                saida.append(linhas[j])
                lin, corpo = "", []
                for t, x in enumerate(vv):
                    lin += ("%8.2f" % x if t % 3 == 0 else
                            "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                    if (t + 1) % 9 == 0:
                        corpo.append(lin); lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                i = k2
            continue
        saida.append(l)
        i += 1

    txt = "\n".join(saida)
    t0 = linhas[0].split("=", 1)[1]
    txt = txt.replace("Geom Title=" + t0, "Geom Title=" + t0 + " + margens na grade", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)
    print("cotas, largura, cutline, RS, valores de n e HTab: intocados")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
