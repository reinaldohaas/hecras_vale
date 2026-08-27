# -*- coding: utf-8 -*-
"""Recorta o perfil das cross sections no MDT do SIG-SC a 1 m.

    python scripts/recortar_do_mdt.py modelo/so_mirim.prj --saida g03

O ORIGINAL NAO E TOCADO, nem o .prj. Sai um .gXX novo.

A REGRA, e por que ela e esta

    z_novo = min(z_HEC-RAS, z_MDT)      onde o MDT tem dado
    z_novo = z_HEC-RAS                  onde o MDT nao tem

  O perfil atual foi cortado do Copernicus GLO-30, que e modelo de SUPERFICIE:
  traz dossel e lamina d'agua. Medido contra o SIG-SC a 1 m, ele esta 7,15 m
  ACIMA do chao no overbank (mediana; p90 +17,7 m; maximo +41,6 m), em 92% das
  secoes. Nao e batimetria -- e um vale inteiro erguido.

  O `min` resolve as duas exigencias ao mesmo tempo:

    NAO INVENTA FUNDO. O MDT so pode BAIXAR o perfil, e so ate o chao que ele
    mediu. Nunca levanta nada.

    PRESERVA A BATIMETRIA. Onde o perfil ja esta abaixo do MDT -- 468 das
    1418 secoes -- ele fica como esta. O MDT nao enxerga sob a agua, entao
    nao tem autoridade para levantar leito.

  Medido antes de gravar: 72,0% dos pontos vem do MDT, 27,6% ficam com o
  HEC-RAS por ja estarem mais baixos, 0,4% ficam por falta de dado.

O QUE NAO MUDA: estacas, contagem de pontos, cutline, estacas de margem,
Manning, htab, comprimentos de trecho. Nenhuma secao e deslocada, estendida
ou reorientada -- as categorias 3, 4 e 5 da auditoria ficam intactas, de
proposito, porque sao outro problema.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import geometria_do_projeto, ler_secoes  # noqa: E402
from mdt_sigsc import MosaicoSigsc  # noqa: E402


def _formatar(sta, z, por_linha=10):
    v = []
    for a, b in zip(sta, z):
        v.append(a); v.append(b)
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%8.2f" % x
        if (i + 1) % por_linha == 0:
            saida.append(linha); linha = ""
    if linha:
        saida.append(linha)
    return saida


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    prj = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g03"
    g01 = geometria_do_projeto(prj)
    raiz = os.path.dirname(prj) or "."
    nome = os.path.splitext(os.path.basename(prj))[0]
    novo = os.path.join(raiz, f"{nome}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(g01):
        raise SystemExit("a saida coincide com a original -- recusado")
    lista = os.path.join(raiz, f"sigsc_tiles_{nome}.txt")
    if not os.path.exists(lista):
        raise SystemExit(f"falta {lista} -- rode antes o scripts/qc_geometria.py")

    print(f"origem : {g01}   (intocado)")
    print(f"MDT    : SIG-SC 1 m, {len(open(lista).read().split(chr(10)))} folhas")
    print(f"saida  : {novo}")
    secoes = ler_secoes(g01)
    mdt = MosaicoSigsc(tiles=open(lista).read().split("\n"))

    linhas = open(g01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    idx = [i for i, l in enumerate(linhas) if l.startswith("#Sta/Elev")]
    assert len(idx) == len(secoes)

    subst = {}
    n_mdt = n_hec = n_vazio = 0
    baixou = []
    thal_saiu = 0
    for k, d in enumerate(secoes):
        st, z = d["sta"], d["z"]
        cut = LineString(d["cut"])
        f = np.clip((st - st[0]) / max(st[-1] - st[0], 1e-9), 0.0, 1.0)
        Q = [cut.interpolate(float(x), normalized=True) for x in f]
        zm = mdt.cota([p.x for p in Q], [p.y for p in Q])
        ok = np.isfinite(zm)
        usa = ok & (zm < z)
        z2 = z.copy()
        z2[usa] = zm[usa]
        n_mdt += int(usa.sum())
        n_hec += int((ok & ~usa).sum())
        n_vazio += int((~ok).sum())
        baixou.append(float(np.median(z - z2)))
        lb, rb = d["lb"], d["rb"]
        if lb is not None:
            antes = lb <= st[int(np.argmin(z))] <= rb
            depois = lb <= st[int(np.argmin(z2))] <= rb
            if antes and not depois:
                thal_saiu += 1
        subst[idx[k]] = _formatar(st, z2)

    fora, i = [], 0
    while i < len(linhas):
        fora.append(linhas[i])
        if i in subst:
            n = int(linhas[i].split("=")[1])
            j, lidos = i + 1, 0
            while j < len(linhas) and lidos < n * 2:
                l = linhas[j]
                if not l.strip() or l[:1].isalpha() or l[:1] == "#":
                    break
                lidos += len([1 for c in range(0, len(l), 8) if l[c:c + 8].strip()])
                j += 1
            fora += subst[i]
            i = j
            continue
        i += 1
    txt = "\n".join(fora).replace(f"Geom Title={nome}",
                                  f"Geom Title={nome} (perfil do MDT SIG-SC 1 m)", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)

    b = np.array(baixou)
    tot = n_mdt + n_hec + n_vazio
    print()
    print(f"pontos do MDT              : {n_mdt:7d}  {100*n_mdt/tot:5.1f}%")
    print(f"pontos mantidos (ja abaixo): {n_hec:7d}  {100*n_hec/tot:5.1f}%")
    print(f"pontos mantidos (sem MDT)  : {n_vazio:7d}  {100*n_vazio/tot:5.1f}%")
    print(f"\nquanto cada secao baixou (mediana da secao):")
    print(f"   p10 {np.percentile(b,10):+.2f}   mediana {np.median(b):+.2f}   "
          f"p90 {np.percentile(b,90):+.2f}   maximo {b.max():+.2f} m")
    print(f"\ntalvegue que saiu do canal : {thal_saiu} secoes")
    print("estacas / cutline / margens / Manning / htab : intocados")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
