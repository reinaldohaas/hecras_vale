# -*- coding: utf-8 -*-
"""Corrige o overbank das cross sections contra o DEM, gravando uma geometria NOVA.

    python scripts/corrigir_secoes.py modelo/so_mirim.prj --saida g02 --so-cortadas
    python scripts/corrigir_secoes.py modelo/so_mirim.prj --saida g02 --rs-min 60000 --rs-max 70000

O ARQUIVO ORIGINAL NAO E TOCADO, nem o .prj. Sai um .gXX novo, que so passa a
valer quando alguem apontar um plano para ele.

O que muda, e o que nao muda:

  MUDAM APENAS AS COTAS, E APENAS FORA DAS MARGENS. Para cada ponto com
  estaca fora de [lb, rb]:
      - onde o DEM existe, a cota passa a ser a do DEM   (item 2)
      - onde o DEM falta e o ponto e spike, interpola-se
        linearmente entre os vizinhos validos             (item 1)

  NAO MUDAM: as estacas, a contagem de pontos, a cutline, as estacas de
  margem, o Manning, o htab, os comprimentos de trecho, nem uma linha entre
  `lb` e `rb`. A BATIMETRIA FICA BYTE A BYTE COMO ESTAVA -- e ela e dado, nao
  erro: o DEM de superficie ve a lamina d'agua, nao o fundo.

  NAO SE DESLOCA NEM SE ESTENDE SECAO NENHUMA. O talvegue do DEM estar perto
  da extremidade e diagnostico, nunca gatilho.

Ao final relata, entre outras coisas, em quantas secoes o ponto mais baixo
passou a cair FORA do canal -- que e o efeito colateral a vigiar quando o
overbank do DEM desce abaixo do leito.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import (Dem, amostrar_dem, geometria_do_projeto,  # noqa: E402
                       ler_secoes, _mediana_movel, SPIKE)


def _formatar(sta, z, por_linha=10):
    """Coluna fixa de 8 caracteres, 10 valores por linha -- como o RAS grava."""
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


def corrigir_perfil(d, zd):
    """Novo vetor de cotas. Devolve (z_novo, n_dem, n_spike)."""
    st, z = d["sta"], d["z"].copy()
    lb, rb = d.get("lb"), d.get("rb")
    if lb is None or rb is None or rb <= lb:
        return z, 0, 0
    fora = (st < lb) | (st > rb)
    tem = np.isfinite(zd)

    # item 2 -- overbank do DEM
    m2 = fora & tem
    n_dem = int(m2.sum())
    z[m2] = zd[m2]

    # item 1 -- spike onde o DEM nao alcanca
    m1 = fora & ~tem
    n_spike = 0
    if m1.any():
        mm = _mediana_movel(z)
        pico = m1 & (np.abs(z - mm) > SPIKE)
        for i in np.flatnonzero(pico):
            a = i - 1
            while a >= 0 and pico[a]:
                a -= 1
            b = i + 1
            while b < len(z) and pico[b]:
                b += 1
            if a < 0 or b >= len(z):
                continue
            z[i] = float(np.interp(st[i], [st[a], st[b]], [z[a], z[b]]))
            n_spike += 1
    return z, n_dem, n_spike


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    prj = argv[0]
    saida_ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g02"
    g01 = geometria_do_projeto(prj)
    raiz = os.path.dirname(prj) or "."
    nome = os.path.splitext(os.path.basename(prj))[0]
    import glob
    tif = glob.glob(os.path.join(raiz, "Terrain", f"{nome}_Terreno*.tif"))[0]
    novo = os.path.join(raiz, f"{nome}.{saida_ext}")
    if os.path.abspath(novo) == os.path.abspath(g01):
        raise SystemExit("a saida coincide com a geometria original -- recusado")

    print(f"origem : {g01}   (intocado)")
    print(f"DEM    : {tif}")
    print(f"saida  : {novo}")

    secoes = ler_secoes(g01)
    dem = Dem(tif)

    # ---- quais secoes entram
    alvo = np.ones(len(secoes), bool)
    if "--so-cortadas" in argv:
        import pickle
        est = pickle.load(open(os.path.join(raiz, f"estado_{nome}.pkl"), "rb"))
        chave = next(iter(est["xs_pronto"]))
        flag = {round(float(x["rs"]), 2): bool(x.get("interpolada"))
                for x in est["xs_pronto"][chave]}
        alvo = np.array([flag.get(round(d["rs"], 2)) is False for d in secoes])
        print(f"escopo : SO as cortadas do terreno -- {alvo.sum()} de {len(secoes)}")
    for k, cmp_ in (("--rs-min", np.greater_equal), ("--rs-max", np.less_equal)):
        if k in argv:
            v = float(argv[argv.index(k) + 1])
            alvo &= cmp_(np.array([d["rs"] for d in secoes]), v)
    if not alvo.any():
        raise SystemExit("nenhuma secao no escopo")

    # ---- reescreve os blocos
    linhas = open(g01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    idx_bloco = [i for i, l in enumerate(linhas) if l.startswith("#Sta/Elev")]
    assert len(idx_bloco) == len(secoes), "blocos e secoes nao batem"

    tot_dem = tot_spike = n_sec = 0
    fugiu = []
    subst = {}
    for k, d in enumerate(secoes):
        if not alvo[k]:
            continue
        zd = amostrar_dem(d, dem)
        z2, nd, ns = corrigir_perfil(d, zd)
        if nd == 0 and ns == 0:
            continue
        lb, rb = d["lb"], d["rb"]
        antes_fora = not (lb <= d["sta"][int(np.argmin(d["z"]))] <= rb)
        depois_fora = not (lb <= d["sta"][int(np.argmin(z2))] <= rb)
        if depois_fora and not antes_fora:
            fugiu.append(d["rs"])
        n_sec += 1; tot_dem += nd; tot_spike += ns
        subst[idx_bloco[k]] = _formatar(d["sta"], z2)

    fora = []
    i = 0
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
    txt = "\n".join(fora)
    txt = txt.replace(f"Geom Title={nome}",
                      f"Geom Title={nome} (overbank do DEM)", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)

    print()
    print(f"secoes alteradas          : {n_sec}")
    print(f"pontos vindos do DEM      : {tot_dem}")
    print(f"spikes interpolados       : {tot_spike}")
    print(f"canal (entre lb e rb)     : intocado")
    print(f"estacas / cutline / margens/ Manning / htab : intocados")
    if fugiu:
        print(f"\nATENCAO: em {len(fugiu)} secoes o ponto mais baixo passou a "
              f"cair FORA do canal (o overbank do DEM desceu abaixo do leito).")
        print("   RS:", ", ".join(f"{r:.0f}" for r in fugiu[:12]),
              "..." if len(fugiu) > 12 else "")
    else:
        print("\nem nenhuma secao o ponto mais baixo saiu do canal")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
