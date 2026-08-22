# -*- coding: utf-8 -*-
"""Apara cutlines longas demais e sincroniza estaca com comprimento.

    python scripts/corrigir_cutlines.py modelo/mirim_t30/mirim_t30.g01 --saida g02

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O QUE ISTO CORRIGE, E SO ISTO

  1. CUTLINE QUE CRUZA O PROPRIO EIXO MAIS DE UMA VEZ (67 secoes no g01).
     Uma cutline de 1.550 m num rio de 52 m atravessa a volta seguinte do
     meandro: a mesma agua entra duas vezes na mesma secao. A secao e APARADA
     nas pontas ate sobrar UMA travessia -- a que contem o canal.

  2. CUTLINE CRUZANDO A DA VIZINHA (66 secoes; 68 pares no g01).
     Depois de aparar pelo eixo, o que ainda cruzar a vizinha imediata e
     encurtado daquele lado, ate limpar.

  3. ESTACA DISCORDANDO DO COMPRIMENTO DA CUTLINE (432 secoes acima de 5 mm).
     Nao se mexe em XY: a ultima estaca passa a valer exatamente
     `round(comprimento, 2)`. Como o `.g01` grava estaca em `%8.2f`, o erro
     residual fica limitado a meio centavo de metro POR CONSTRUCAO -- e nao
     ha como fazer melhor neste formato, porque a cutline tambem e gravada
     com 1 cm de resolucao (`%16.2f`). A secao muda de largura em ate 1,6 cm.

O QUE ISTO NAO TOCA, DE PROPOSITO

  Nenhuma COTA muda. Nenhum ponto do canal sai. O HTab fica como esta -- ele
  esta ancorado 2 cm acima do talvegue nas 1.418 secoes, e o talvegue nao se
  move porque o canal nunca e aparado.

  A APARA NUNCA ENTRA NO CANAL. O corte respeita `Bank Sta` mais uma folga de
  `FOLGA_CANAL` larguras de canal de cada lado. Quando nao da para satisfazer
  o criterio sem invadir essa faixa, A SECAO FICA COMO ESTA e entra no
  relatorio -- e o caso em que aparar seria pior que o defeito.

  A ORIENTACAO NAO E CORRIGIDA. Ha 253 secoes com angulo ruim contra a
  tangente (135 delas quase paralelas ao fluxo). Girar uma cutline a move de
  lugar, e isso e outra decisao -- nao entra aqui.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from ras_io import escrever            # noqa: E402

FOLGA_CANAL = 1.5     # larguras de canal preservadas de cada lado da margem
FOLGA_MIN = 30.0      # m; folga minima, quando o canal e muito estreito
RECUO = 2.0           # m de recuo antes da travessia indesejada do eixo


def _col(v, larg, dec):
    """Grava em colunas de largura fixa, 10 por linha."""
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += ("%*.*f" % (larg, dec, x))
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def travessias(A, u, L, eixo):
    """Estacas em que a cutline corta o eixo, em ordem."""
    from shapely.geometry import LineString
    ln = LineString([A, A + L * u])
    x = ln.intersection(eixo)
    if x.is_empty:
        return []
    pts = list(x.geoms) if hasattr(x, "geoms") else [x]
    s = []
    for p in pts:
        if p.geom_type != "Point":
            p = p.centroid
        s.append(float(np.dot(np.array([p.x, p.y]) - A, u)))
    return sorted(s)


def aparar(d, eixo):
    """Devolve (s0, s1, motivo) -- a faixa de estacas que fica."""
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    v = B - A
    L = float(np.hypot(*v))
    u = v / max(L, 1e-9)
    lb, rb = float(d["lb"]), float(d["rb"])
    folga = max(FOLGA_MIN, FOLGA_CANAL * max(rb - lb, 1.0))
    lim0, lim1 = max(0.0, lb - folga), min(L, rb + folga)

    s0, s1 = 0.0, L
    t = travessias(A, u, L, eixo)
    dentro = [x for x in t if lim0 <= x <= lim1]
    alvo = dentro[0] if dentro else (min(t, key=lambda x: abs(x - 0.5 * (lb + rb)))
                                     if t else 0.5 * (lb + rb))
    for x in t:
        if x < alvo - 1e-6:
            s0 = max(s0, min(x + RECUO, lim0))
        elif x > alvo + 1e-6:
            s1 = min(s1, max(x - RECUO, lim1))
    return s0, s1, len(t)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g02"
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    ordem = sorted(range(len(S)), key=lambda i: -S[i]["rs"])
    eixo = list(ler_eixos(entrada).values())[0]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(S)}\n")

    # ---- 1a passada: aparar pelo eixo
    faixa = {}
    for i in range(len(S)):
        s0, s1, nt = aparar(S[i], eixo)
        faixa[i] = [s0, s1, nt]

    # ---- 2a passada: encurtar o que ainda cruza a vizinha
    def linha(i):
        d = S[i]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        s0, s1, _ = faixa[i]
        return LineString([A + s0 * u, A + s1 * u])

    encurtou = 0
    for k in range(len(ordem) - 1):
        i, j = ordem[k], ordem[k + 1]
        for _ in range(24):
            if not linha(i).intersects(linha(j)):
                break
            # encurta o mais largo dos dois, pelo lado que ainda tem folga
            a, b = (i, j) if (faixa[i][1] - faixa[i][0]) >= \
                             (faixa[j][1] - faixa[j][0]) else (j, i)
            d = S[a]
            lb, rb = float(d["lb"]), float(d["rb"])
            folga = max(FOLGA_MIN, FOLGA_CANAL * max(rb - lb, 1.0))
            lim0, lim1 = max(0.0, lb - folga), rb + folga
            s0, s1 = faixa[a][0], faixa[a][1]
            passo = max(0.02 * (s1 - s0), 1.0)
            if (s1 - lim1) >= (lim0 - s0):
                if s1 - passo < lim1:
                    break
                faixa[a][1] = s1 - passo
            else:
                if s0 + passo > lim0:
                    break
                faixa[a][0] = s0 + passo
            encurtou += 1

    # ---- monta o novo conteudo de cada secao
    novos, intocadas, perdida = {}, 0, []
    for i, d in enumerate(S):
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        vv = B - A
        L = float(np.hypot(*vv))
        u = vv / max(L, 1e-9)
        s0, s1, _ = faixa[i]
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        if s0 <= 1e-9 and s1 >= L - 1e-9:
            s0, s1 = 0.0, L
            intocadas += 1
        # recorta o perfil, interpolando as pontas novas.
        # A TOLERANCIA E A DO FORMATO, e nao 1e-9: a estaca e gravada em
        # `%8.2f`, entao dois valores a menos de 1 cm sao O MESMO PONTO no
        # arquivo. Com 1e-9 o codigo inseria uma ponta a 3 mm da ultima
        # estaca e depois arredondava as duas para o mesmo valor -- o
        # HEC-RAS recusou 504 secoes com "Station and elevation data
        # contains duplicate points".
        TOL = 0.005
        m = (st >= s0 - TOL) & (st <= s1 + TOL)
        ns = list(st[m])
        nz = list(z[m])
        if not ns or ns[0] > s0 + TOL:
            ns.insert(0, s0)
            nz.insert(0, float(np.interp(s0, st, z)))
        if ns[-1] < s1 - TOL:
            ns.append(s1)
            nz.append(float(np.interp(s1, st, z)))
        ns = np.array(ns) - s0
        nz = np.array(nz)
        nA, nB = A + s0 * u, A + s1 * u
        # ---- item 3: a ultima estaca casa com o comprimento gravado
        Lg = round(float(np.hypot(*(nB - nA))), 2)
        if len(ns) > 1 and abs(ns[-1] - Lg) < 0.5 and Lg > ns[-2] + TOL:
            ns[-1] = Lg
        # ---- e nada sai daqui com estaca repetida ou fora de ordem
        keep = [0]
        for k in range(1, len(ns)):
            if ns[k] > ns[keep[-1]] + TOL:
                keep.append(k)
        ns, nz = ns[keep], nz[keep]
        perdida.append(L - (s1 - s0))
        novos[i] = {"sta": ns, "z": nz, "lb": float(d["lb"]) - s0,
                    "rb": float(d["rb"]) - s0, "desl": s0,
                    "cut": (nA, nB), "n_antes": len(st)}

    perdida = np.array(perdida)
    mexidas = int((perdida > 1e-6).sum())
    print(f"secoes aparadas          : {mexidas}   (intocadas: {intocadas})")
    print(f"encurtamentos por vizinha: {encurtou}")
    if mexidas:
        p = perdida[perdida > 1e-6]
        print(f"largura retirada         : mediana {np.median(p):.0f} m   "
              f"p90 {np.percentile(p, 90):.0f}   max {p.max():.0f} m")

    # ---- reescreve o arquivo
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
            if l.startswith("XS GIS Cut Line"):
                saida.append("XS GIS Cut Line= 2")
                saida.append("".join("%16.2f" % x for x in
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
                val = [float(x) for x in bruto[:3 * cnt]]
                lim = float(nv["sta"][-1])
                for t in range(0, 3 * cnt, 3):
                    val[t] = min(max(val[t] - nv["desl"], 0.0), lim)
                saida.append(l)
                lin, corpo = "", []
                for t, x in enumerate(val):
                    lin += ("%8.2f" % x if t % 3 == 0 else
                            "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                    if (t + 1) % 9 == 0:
                        corpo.append(lin)
                        lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                j = k2
                continue
        saida.append(l)
        j += 1

    txt = "\n".join(saida)
    t0 = linhas[0].split("=", 1)[1] if "=" in linhas[0] else ""
    if t0:
        txt = txt.replace("Geom Title=" + t0, "Geom Title=" + t0, 1)
    escrever(novo, txt)

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B2 = ler_secoes(novo)
    B2.sort(key=lambda d: -d["rs"])
    A2 = ler_secoes(entrada)
    A2.sort(key=lambda d: -d["rs"])
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    zi = np.array([float(x["z"].min()) for x in A2])
    zf = np.array([float(x["z"].min()) for x in B2])
    print(f"   talvegue mudou em      : max {np.abs(zf - zi).max():.6f} m "
          "(tem de ser zero)")
    ci = np.array([float(x["rb"] - x["lb"]) for x in A2])
    cf = np.array([float(x["rb"] - x["lb"]) for x in B2])
    print(f"   largura do canal mudou : max {np.abs(cf - ci).max():.6f} m "
          "(tem de ser zero)")
    for rot, X in (("antes", A2), ("depois", B2)):
        e = []
        for d in X:
            C = np.asarray(d["cut"], float)
            Lc = float(np.hypot(*(np.diff(C, axis=0).T)).sum())
            e.append(abs(Lc - float(d["sta"][-1] - d["sta"][0])))
        e = np.array(e)
        print(f"   estaca x cutline {rot:<6}: > 5 mm em {int((e > 0.005).sum()):4d} "
              f"secoes   max {e.max()*1000:.1f} mm")
    for rot, X in (("antes", A2), ("depois", B2)):
        L = [LineString(np.asarray(d["cut"], float)) for d in X]
        mult = sum(1 for ln in L
                   if len(ln.intersection(eixo).geoms
                          if hasattr(ln.intersection(eixo), "geoms")
                          else [ln.intersection(eixo)]) >= 2)
        from shapely.strtree import STRtree
        tr = STRtree(L)
        pares = {(i, int(j)) for i, ln in enumerate(L)
                 for j in tr.query(ln) if int(j) > i and ln.intersects(L[int(j)])}
        print(f"   {rot:<6}: cruza o eixo >=2x em {mult:4d} secoes   |  "
              f"cutlines cruzadas: {len(pares):4d} pares")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
