# -*- coding: utf-8 -*-
"""Alarga as secoes de um trecho, estendendo a planicie com o MDT.

    python scripts/alargar_trecho.py modelo/so_mirim.g08 --saida g10 \
        --rs-min 124000 --rs-max 128000 --largura 800

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

POR QUE, E DE ONDE VEM O NUMERO

  Comparado ao modelo de referencia do mesmo rio (itajaim_hecras), que roda:

                       referencia g04    este modelo
      largura da secao      800 m           132 m
      planicie por lado     360 m            41 m
      canal ocupa            10%             37%
      altura da secao      82,9 m          34,4 m

  Com 41 m de planicie de cada lado, a cheia que extravasa sai pela borda da
  secao, e o HEC-RAS avisa "Extrapolated above Cross Section Table" -- o que
  acontece em todas as rodadas deste modelo e em nenhuma da referencia.

O QUE SE FAZ, E O QUE NAO SE TOCA

  A secao e estendida PARA FORA, nas duas pontas, ao longo da propria cutline,
  ate a largura alvo. Os pontos novos recebem a cota do MDT SIG-SC 1 m.

  NADA do que ja existe muda: nem uma cota, nem o canal, nem a posicao relativa
  de nada. As estacas sao rebaseadas para comecar em zero -- e como `lb`, `rb`
  e as quebras do Manning andam todas pela MESMA constante (a extensao da
  esquerda), elas continuam caindo exatamente sobre estacas existentes.

  O HTab e reancorado SE a extensao encontrar chao mais baixo que o leito
  atual; do contrario fica como esta.

  Isto ESTENDE secoes, o que e proibido como regra automatica neste projeto.
  Aqui e um experimento delimitado por intervalo de RS, pedido explicitamente.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                      # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402
from ras_io import escrever                            # noqa: E402

PASSO = 10.0        # m entre pontos novos -- a referencia usa 8 m
FOLGA_HTAB = 0.02


def _col8(v):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%8.2f" % x
        if (i + 1) % 10 == 0:
            saida.append(linha); linha = ""
    if linha:
        saida.append(linha)
    return saida


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g10"
    rmin = float(argv[argv.index("--rs-min") + 1])
    rmax = float(argv[argv.index("--rs-max") + 1])
    alvo = float(argv[argv.index("--largura") + 1]) if "--largura" in argv else 800.0
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    alvos = [i for i, d in enumerate(S) if rmin <= d["rs"] <= rmax
             and float(d["sta"][-1] - d["sta"][0]) < alvo]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"trecho : RS {rmin:.0f} a {rmax:.0f}   largura alvo {alvo:.0f} m")
    print(f"secoes no trecho a alargar: {len(alvos)} de {len(S)}")
    if not alvos:
        raise SystemExit("nada a alargar")

    # tiles: o dominio cresce com a extensao
    P = np.vstack([np.asarray(S[i]["cut"], float) for i in alvos])
    folga = alvo
    bbox = (P[:, 0].min() - folga, P[:, 1].min() - folga,
            P[:, 0].max() + folga, P[:, 1].max() + folga)
    lista = os.path.join(raiz, f"sigsc_tiles_{base}_alargado.txt")
    if os.path.exists(lista):
        tiles = open(lista).read().split("\n")
    else:
        tiles = tiles_do_dominio(bbox)
        open(lista, "w").write("\n".join(tiles))
    print(f"MDT: {len(tiles)} folhas sobre o trecho alargado")
    mdt = MosaicoSigsc(tiles=tiles)

    novos, sem_dado, htab_mudou = {}, 0, 0
    larg_antes, larg_depois = [], []
    for i in alvos:
        d = S[i]
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        L = float(st[-1] - st[0])
        A = np.array(d["cut"][0], float)
        B = np.array(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        add = 0.5 * (alvo - L)
        n = max(1, int(round(add / PASSO)))
        pas = add / n

        esq_s = np.array([-(k * pas) for k in range(n, 0, -1)])
        dir_s = np.array([L + k * pas for k in range(1, n + 1)])
        Pe = [A + s_ * u for s_ in esq_s]
        Pd = [A + s_ * u for s_ in dir_s]
        ze = mdt.cota([p[0] for p in Pe], [p[1] for p in Pe])
        zd = mdt.cota([p[0] for p in Pd], [p[1] for p in Pd])
        # onde falta MDT, repete a cota da ponta -- nao inventa relevo
        f = ~np.isfinite(ze); sem_dado += int(f.sum()); ze[f] = z[0]
        f = ~np.isfinite(zd); sem_dado += int(f.sum()); zd[f] = z[-1]

        ns = np.concatenate([esq_s, st, dir_s]) + add      # rebaseia em zero
        nz = np.concatenate([ze, z, zd])
        nlb = float(d["lb"]) + add
        nrb = float(d["rb"]) + add
        novos[i] = {"sta": ns, "z": nz, "lb": nlb, "rb": nrb, "add": add,
                    "cut": np.array([A - add * u, B + add * u]),
                    "baixou": float(nz.min()) < float(z.min()) - 1e-9}
        larg_antes.append(L); larg_depois.append(float(ns[-1] - ns[0]))

    print(f"pontos sem MDT (repetiram a cota da ponta): {sem_dado}")
    print(f"largura: mediana {np.median(larg_antes):.0f} -> "
          f"{np.median(larg_depois):.0f} m")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec = -1
    saida, j = [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        if i_sec in novos:
            nv = novos[i_sec]
            if l.startswith("XS GIS Cut Line"):
                saida.append("XS GIS Cut Line= 2 ")
                saida.append("".join("%16.4f" % x for x in
                                     (nv["cut"][0][0], nv["cut"][0][1],
                                      nv["cut"][1][0], nv["cut"][1][1])))
                j += 1
                while j < len(linhas) and linhas[j].strip() and \
                        linhas[j][:1] in " -0123456789":
                    j += 1
                continue
            if l.startswith("#Sta/Elev"):
                v = []
                for a, b in zip(nv["sta"], nv["z"]):
                    v.append(a); v.append(b)
                saida.append("#Sta/Elev= %d " % len(nv["sta"]))
                saida += _col8(v)
                cnt = int(l.split("=")[1]); j += 1; lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(x), 8)
                                  if x[c:c + 8].strip()])
                    j += 1
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%s,%s" % (_fmt(nv["lb"]), _fmt(nv["rb"])))
                j += 1
                continue
            if l.startswith("#Mann="):
                cnt = int(l.split("=")[1].split(",")[0])
                bruto, k2 = [], j + 1
                while k2 < len(linhas) and len(bruto) < 3 * cnt:
                    x = linhas[k2]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                              if x[c:c + 8].strip()]
                    k2 += 1
                vv = [float(x) for x in bruto[:3 * cnt]]
                for t in range(0, 3 * cnt, 3):
                    vv[t] = 0.0 if t == 0 else vv[t] + nv["add"]
                saida.append(l)
                lin, corpo = "", []
                for t, x in enumerate(vv):
                    lin += ("%8.2f" % x if t % 3 == 0 else
                            "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                    if (t + 1) % 9 == 0:
                        corpo.append(lin); lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                j = k2
                continue
            if l.startswith("XS HTab Starting El and Incr") and nv["baixou"]:
                p_ = [x.strip() for x in l.split("=", 1)[1].split(",")]
                el = float(np.min(nv["z"])) + FOLGA_HTAB
                saida.append(f"XS HTab Starting El and Incr={el:.2f},"
                             f"{float(p_[1]):.3f}, {int(p_[2])} ")
                htab_mudou += 1
                j += 1
                continue
        saida.append(l)
        j += 1

    txt = "\n".join(saida)
    t0 = linhas[0].split("=", 1)[1] if "=" in linhas[0] else ""
    if t0:
        txt = txt.replace("Geom Title=" + t0, "Geom Title=" + t0 + " + alargado", 1)
    escrever(novo, txt)
    print(f"HTab reancorado em {htab_mudou} secoes (a extensao achou chao mais baixo)")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
