# -*- coding: utf-8 -*-
"""Reancora a tabela hidraulica (HTab) no talvegue de cada secao.

    python scripts/ajustar_htab.py modelo/so_mirim.g05 --saida g06

A ENTRADA NAO E TOCADA. Muda EXCLUSIVAMENTE a linha
`XS HTab Starting El and Incr=`, e dela apenas o PRIMEIRO campo. O incremento
e a contagem de pontos ficam como estao.

POR QUE

  A tabela hidraulica de uma secao e uma lista de N cotas, do fundo para cima,
  onde o HEC-RAS pre-calcula area, perimetro e conducao. Se ela comeca ACIMA
  do leito, o solver nao tem valores para as laminas baixas: extrapola, e
  extrapolar area molhada perto do fundo diverge no primeiro passo.

  O g01 respeitava isso -- a tabela comecava 2 cm acima do talvegue em TODAS
  as 1418 secoes. Ao recortar o perfil no MDT eu baixei o leito em 4,30 m na
  mediana (24,55 m no pior caso) e nao baixei a tabela junto. Resultado
  medido no g05: a tabela passou a comecar 2,225 m acima do leito na mediana,
  17,02 m no pior caso, e ficou acima do fundo em 1.179 das 1.418 secoes.

  O HEC-RAS caiu no primeiro passo, e disse onde:

      *****  Warning!  Extrapolated above Cross Section Table at:  *****
          Itajai_Mirim R1    R.S.   125829.6
      Minimum error exceeds allowable tolerance at 01AUG2026 00:00:30
      Simulation went unstable at: 01AUG2026 00:00:45

  Naquela secao o talvegue caiu de 139,18 para 130,01 m e a tabela continuou
  comecando em 139,20 -- 9,19 m acima do fundo novo.

O INCREMENTO NAO MUDA, e isso foi conferido antes: com a contagem atual a
faixa coberta e de 54,9 m na mediana contra 30,2 m de altura de secao, e
NENHUMA das 1418 exigiria incremento maior. Baixar o inicio basta.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes  # noqa: E402

FOLGA = 0.02      # m acima do talvegue -- a mesma do g01


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g06"
    raiz = os.path.dirname(entrada) or "."
    nome = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{nome}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    secoes = ler_secoes(entrada)
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    idx = [i for i, l in enumerate(linhas)
           if l.startswith("XS HTab Starting El and Incr")]
    if len(idx) != len(secoes):
        raise SystemExit(f"blocos HTab ({len(idx)}) e secoes ({len(secoes)}) "
                         f"nao batem")

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")

    antes, depois, curto = [], [], []
    for k, d in enumerate(secoes):
        p = [x.strip() for x in linhas[idx[k]].split("=", 1)[1].split(",")]
        el, inc, n = float(p[0]), float(p[1]), int(p[2])
        zmin = float(np.min(d["z"])); zmax = float(np.max(d["z"]))
        antes.append(el - zmin)
        novo_el = zmin + FOLGA
        teto = novo_el + inc * (n - 1)
        if teto < zmax:
            curto.append((d["rs"], zmax - teto))
        linhas[idx[k]] = (f"XS HTab Starting El and Incr={novo_el:.2f},"
                          f"{inc:.3f}, {n} ")
        depois.append(novo_el - zmin)

    txt = "\n".join(linhas)
    t0 = linhas[0].split("=", 1)[1]
    txt = txt.replace("Geom Title=" + t0, "Geom Title=" + t0 + " + htab", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)

    a = np.array(antes); b = np.array(depois)
    print(f"\ninicio do HTab menos o talvegue:")
    print(f"   antes : mediana {np.median(a):+7.3f}  p90 {np.percentile(a,90):+7.2f}"
          f"  max {a.max():+7.2f} m   acima do leito em {(a>0.05).sum()} secoes")
    print(f"   depois: mediana {np.median(b):+7.3f}  p90 {np.percentile(b,90):+7.2f}"
          f"  max {b.max():+7.2f} m   acima do leito em {(b>0.05).sum()} secoes")
    if curto:
        print(f"\nATENCAO: em {len(curto)} secoes o teto da tabela ficaria ABAIXO "
              f"do topo da secao (faltaria ate {max(c[1] for c in curto):.2f} m)")
        for rs, f in sorted(curto, key=lambda c: -c[1])[:8]:
            print(f"   RS {rs:>10.2f}  faltam {f:.2f} m")
    else:
        print("\nteto da tabela cobre o topo em todas as secoes")
    print("incremento e contagem de pontos: inalterados")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
