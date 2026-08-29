# -*- coding: utf-8 -*-
"""Confere a geometria de um projeto do jeito que o RAS Mapper confere.

Roda sobre o `<projeto>.g01.hdf` e responde, com numero, as perguntas que a
Errors Layer do RAS Mapper responde com uma caixa de dialogo que so diz
quantos erros existem.

    python scripts/auditar_geometria.py                    todos os projetos
    python scripts/auditar_geometria.py modelo/so_mirim    um so
    python scripts/auditar_geometria.py --detalhe          lista as secoes

As checagens, e de onde vem cada uma. As frases entre aspas foram extraidas
do proprio RasMapperLib.dll, e sao as mensagens que o RAS Mapper usa:

  1. cutline x estacas   "The polyline length must match the last station
                          minus the first station."  O comprimento da linha de
                          corte tem de ser igual a ultima estaca menos a
                          primeira. Sem isso o RAS nao consegue mapear estaca
                          -> posicao no mapa, e as bank lines e edge lines que
                          saem dai nao coincidem com a secao.
  2. cutline x eixo      "Some cross-sections do not cross a river line -
                          invalid geometry."
  3. estacas             fora de ordem, repetidas, ou secao com menos de dois
                          pontos.
  4. margens             margem esquerda >= direita, fora da faixa de estacas,
                          ou canal ocupando menos de 2% da secao.
  5. RS                  repetido dentro do trecho, ou fora de ordem.
  6. espacamento         vaos abaixo do piso de 25 m, onde o 1D deixa de valer,
                          e o salto entre vaos vizinhos.
  7. preprocessador      se ha bank lines, edge lines e tabelas -- isto e, se o
                          projeto chegou a computar alguma vez.

Nao conserta nada. So mede, e diz onde.
"""
import os
import sys
import glob

import numpy as np

TOL_CUTLINE = 0.01      # m; abaixo disto e arredondamento do float, nao erro
PISO_DX = 25.0          # m; o mesmo espacamento_piso do config
CANAL_MIN = 0.02        # fracao da largura da secao


def _ler(hdf):
    import h5py
    from shapely.geometry import LineString
    f = h5py.File(hdf, "r")
    G = f["Geometry"]
    grupos = set(G) | set(f)
    if "Cross Sections" not in G:
        return None
    xa = G["Cross Sections/Attributes"][:]
    xi = G["Cross Sections/Polyline Info"][:]
    xp = G["Cross Sections/Polyline Points"][:].astype(float)
    si = G["Cross Sections/Station Elevation Info"][:]
    sv = G["Cross Sections/Station Elevation Values"][:].astype(float)
    eixos = {}
    if "River Centerlines" in G:
        ra = G["River Centerlines/Attributes"][:]
        ri = G["River Centerlines/Polyline Info"][:]
        rp = G["River Centerlines/Polyline Points"][:].astype(float)
        for k in range(len(ra)):
            o, n = int(ri[k, 0]), int(ri[k, 1])
            eixos[(ra[k]["River Name"].decode().strip(),
                   ra[k]["Reach Name"].decode().strip())] = LineString(rp[o:o + n])
    sec = []
    for k in range(len(xa)):
        o, n = int(xi[k, 0]), int(xi[k, 1])
        so, sn = int(si[k, 0]), int(si[k, 1])
        sec.append(dict(
            rio=xa[k]["River"].decode().strip(),
            trecho=xa[k]["Reach"].decode().strip(),
            rs=float(xa[k]["RS"]),
            cut=LineString(xp[o:o + n]),
            sta=sv[so:so + sn, 0],
            z=sv[so:so + sn, 1],
            lb=float(xa[k]["Left Bank"]),
            rb=float(xa[k]["Right Bank"]),
        ))
    return sec, eixos, grupos


def auditar(hdf, detalhe=False):
    lido = _ler(hdf)
    if lido is None:
        return None
    sec, eixos, grupos = lido
    n = len(sec)
    p = {}

    # 1. o comprimento da cutline TEM de ser a faixa de estacas
    d = np.array([s["cut"].length - (s["sta"][-1] - s["sta"][0]) for s in sec])
    p["cutline != estacas"] = np.flatnonzero(np.abs(d) > TOL_CUTLINE)
    p["_desvio"] = d

    # 2. a cutline tem de cruzar o eixo do proprio trecho, uma vez
    fora, multi = [], []
    for i, s in enumerate(sec):
        e = eixos.get((s["rio"], s["trecho"]))
        if e is None or not s["cut"].intersects(e):
            fora.append(i); continue
        x = s["cut"].intersection(e)
        if x.geom_type.startswith("Multi") and len(x.geoms) > 1:
            multi.append(i)
    p["cutline nao cruza o eixo"] = np.array(fora, int)
    p["cutline cruza o eixo 2+ vezes"] = np.array(multi, int)

    # 3. estacas
    p["estacas fora de ordem"] = np.array(
        [i for i, s in enumerate(sec) if np.any(np.diff(s["sta"]) < 0)], int)
    p["estacas repetidas"] = np.array(
        [i for i, s in enumerate(sec) if np.any(np.diff(s["sta"]) == 0)], int)
    p["menos de 2 pontos"] = np.array(
        [i for i, s in enumerate(sec) if len(s["sta"]) < 2], int)

    # 4. margens
    larg = np.array([s["sta"][-1] - s["sta"][0] for s in sec])
    lb = np.array([s["lb"] for s in sec]); rb = np.array([s["rb"] for s in sec])
    s0 = np.array([s["sta"][0] for s in sec]); s1 = np.array([s["sta"][-1] for s in sec])
    p["margem direita <= esquerda"] = np.flatnonzero(rb <= lb)
    p["margem fora das estacas"] = np.flatnonzero((lb < s0 - 1e-6) | (rb > s1 + 1e-6))
    p["canal < 2% da secao"] = np.flatnonzero((rb - lb) / np.maximum(larg, 1e-9) < CANAL_MIN)
    centro = ((lb + rb) / 2 - s0) / np.maximum(larg, 1e-9)
    p["canal fora do terco central"] = np.flatnonzero((centro < 1 / 3) | (centro > 2 / 3))

    # 5. RS e espacamento, por trecho
    rsdup, rsord, curto, salto = [], [], [], []
    for chave in {(s["rio"], s["trecho"]) for s in sec}:
        idx = [i for i, s in enumerate(sec) if (s["rio"], s["trecho"]) == chave]
        idx.sort(key=lambda i: -sec[i]["rs"])
        rs = np.array([sec[i]["rs"] for i in idx])
        vis = {}
        for j, i in enumerate(idx):
            r = round(rs[j], 2)
            if r in vis: rsdup.append(i)
            vis[r] = i
        dx = -np.diff(rs)
        curto += [idx[j] for j in np.flatnonzero(dx < PISO_DX)]
        rsord += [idx[j + 1] for j in np.flatnonzero(dx <= 0)]
        if len(dx) > 2:
            razao = np.maximum(dx[1:] / np.maximum(dx[:-1], 1e-9),
                               dx[:-1] / np.maximum(dx[1:], 1e-9))
            salto += [idx[j + 1] for j in np.flatnonzero(razao > 4)]
    p["RS repetido no trecho"] = np.array(sorted(set(rsdup)), int)
    p["RS fora de ordem"] = np.array(sorted(set(rsord)), int)
    p["vao abaixo de %.0f m" % PISO_DX] = np.array(sorted(set(curto)), int)
    p["salto de espacamento > 4x"] = np.array(sorted(set(salto)), int)

    p["_n"] = n
    p["_grupos"] = grupos
    p["_sec"] = sec
    return p


def _linha(nome, idx, n):
    q = len(idx)
    marca = "   " if q == 0 else ("<<<" if q == n else " ! ")
    return "   %s %-32s %5d  %5.1f%%" % (marca, nome, q, 100.0 * q / max(n, 1))


def relatorio(hdf, detalhe=False):
    p = auditar(hdf, detalhe)
    if p is None:
        print("%s -- sem secoes" % hdf); return
    n = p["_n"]
    print("=" * 62)
    print("%s   %d secoes" % (hdf, n))
    print("=" * 62)
    for k in [k for k in p if not k.startswith("_")]:
        print(_linha(k, p[k], n))
    d = np.abs(p["_desvio"])
    print("   desvio cutline-estacas: mediana %.4f m, p90 %.3f m, maximo %.2f m"
          % (np.median(d), np.percentile(d, 90), d.max()))
    pre = "SIM" if "GeomPreprocess" in p["_grupos"] else "NAO"
    bl = "SIM" if "River Bank Lines" in p["_grupos"] else "NAO"
    print("   ja computou (GeomPreprocess=%s, River Bank Lines=%s)" % (pre, bl))
    total = len({i for k, v in p.items() if not k.startswith("_") for i in v})
    print("   SECOES COM AO MENOS UM PROBLEMA: %d de %d (%.0f%%)"
          % (total, n, 100.0 * total / max(n, 1)))
    if detalhe:
        sec = p["_sec"]
        for k in [k for k in p if not k.startswith("_")]:
            if len(p[k]) == 0 or len(p[k]) == n:
                continue
            amostra = ", ".join("%s RS %.0f" % (sec[i]["rio"], sec[i]["rs"])
                                for i in p[k][:8])
            print("      %s: %s%s" % (k, amostra, " ..." if len(p[k]) > 8 else ""))
    print()


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    det = "--detalhe" in sys.argv
    if args:
        alvos = []
        for a in args:
            alvos += [a] if a.endswith(".hdf") else glob.glob(a + ".g01.hdf")
    else:
        alvos = sorted(glob.glob("modelo/*.g01.hdf") + glob.glob("modelo/*/*.g01.hdf"))
    if not alvos:
        print("nada para auditar"); sys.exit(1)
    for h in alvos:
        if os.path.getsize(h) == 0:
            continue
        try:
            relatorio(h, det)
        except Exception as e:                                   # noqa: BLE001
            print("%s -- falhou: %s %s\n" % (h, type(e).__name__, e))
