# -*- coding: utf-8 -*-
"""Apara cutline larga demais para a curva em que ela esta.

    python scripts/corrigir_cutlines.py modelo/mirim_t30/mirim_t30.g07 --saida g08

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O QUE ISTO ATACA

  Duas familias do Validate Geometry, que sao o mesmo defeito visto de dois
  angulos -- medido no g07, com quatro reaches:

    446  "This is a edge line self intersection point"
     64  "XS must intersect exactly one Reach"
      5  "XS naming error - River Station out of order"  (sintoma: com a
         cutline cortando o eixo duas vezes, a projecao fica ambigua)

  A causa e uma so: secao larga demais para o raio da curva. Medido no g07,
  57 secoes cortam o proprio eixo mais de uma vez (41 cortam 2x, 12 cortam 3x,
  3 cortam 4x e uma corta 5x), e a largura mediana no reach do meandro e de
  623 m -- com maximo de 1.550 m -- num rio cujo canal tem 52 a 100 m.

O CRITERIO E O RAIO DA CURVA, E NAO UM NUMERO FIXO

  Numa curva de raio R, a secao que se estende alem de R pelo lado de DENTRO
  passa do centro de curvatura: a partir dali ela anda para tras, cruza as
  vizinhas e dobra sobre o proprio rio. E o mecanismo que produz tanto a
  travessia repetida do eixo quanto o emaranhado das edge lines.

  Entao a meia-largura do lado interno fica limitada a `K_CURV * R`, com R
  medido no eixo do PROPRIO reach numa janela proporcional a largura do canal.
  Onde a curva e aberta, R e grande e nada e aparado.

A FOLGA CEDE EM VEZ DE TRAVAR

  A versao anterior protegia 1,5 largura de canal de cada lado como constante,
  e por isso conseguiu aparar so 31 secoes: a travessia indesejada fica, na
  mediana, a 18 m da margem -- as vezes a 0 m --, dentro da faixa protegida.
  Proteger 78 m de planicie onde a volta seguinte do meandro passa a 18 m nao
  descreve terreno nenhum: ali nao HA planicie, ha um pescoco.

  Aqui a folga comeca em 1,5 largura e vai cedendo por tentativas ate
  `FOLGA_PISO` metros alem da margem. O CANAL NUNCA E INVADIDO: entre `lb` e
  `rb` nada e cortado, em nenhuma hipotese. Quando nem com a folga minima o
  criterio se satisfaz, A SECAO FICA COMO ESTA e entra no relatorio.

CIENTE DE REACHES

  Cada secao e conferida contra o eixo do SEU reach -- com quatro reaches, usar
  um eixo so mediria a coisa errada -- e tambem contra os eixos dos OUTROS,
  porque cutline que atravessa o ramo vizinho e ambigua para o HEC-RAS.

NADA DE COTA MUDA

  So saem pontos das PONTAS. Talvegue, largura de canal e as cotas que
  sobrevivem sao identicos; a conferencia no fim mede isso.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from ras_io import escrever            # noqa: E402

K_CURV = 0.80        # fracao do raio de curvatura admitida do lado de dentro
FOLGA_MULT = 1.5     # folga desejada, em larguras de canal
FOLGA_PISO = 8.0     # m; o minimo que se preserva alem da margem
RECUO = 3.0          # m de recuo antes da travessia indesejada
TOL = 0.005          # m; resolucao da estaca no formato do .gNN


def _col(v, larg, dec):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%*.*f" % (larg, dec, x)
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def mapa_reaches(g):
    t = open(g, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    ch, saida = None, []
    for l in t:
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            ch = (p[0].strip(), p[1].strip())
        elif l.startswith("Type RM Length L Ch R"):
            saida.append(ch)
    return saida


def curvatura(eixo, s, janela):
    """Raio de curvatura e sinal do giro, no ponto `s` do eixo."""
    a = max(s - janela, 0.0)
    b = min(s + janela, eixo.length)
    P0 = np.array(eixo.interpolate(a).coords[0])
    P1 = np.array(eixo.interpolate(0.5 * (a + b)).coords[0])
    P2 = np.array(eixo.interpolate(b).coords[0])
    v1, v2 = P1 - P0, P2 - P1
    cr = v1[0] * v2[1] - v1[1] * v2[0]
    d01 = np.hypot(*(P1 - P0))
    d12 = np.hypot(*(P2 - P1))
    d02 = np.hypot(*(P2 - P0))
    area2 = abs(cr)
    if area2 < 1e-9 or d01 * d12 * d02 == 0:
        return np.inf, 0.0
    R = d01 * d12 * d02 / (2.0 * area2)      # raio do circulo pelos 3 pontos
    return R, np.sign(cr)


def travessias(A, u, s0, s1, eixo):
    from shapely.geometry import LineString
    ln = LineString([A + s0 * u, A + s1 * u])
    x = ln.intersection(eixo)
    if x.is_empty:
        return []
    pts = list(x.geoms) if hasattr(x, "geoms") else [x]
    out = []
    for p in pts:
        if p.geom_type != "Point":
            p = p.centroid
        out.append(float(np.dot(np.array([p.x, p.y]) - A, u)))
    return sorted(out)


def faixa(d, eixo, outros):
    """(s0, s1, motivo) -- o quanto da secao fica, e por que foi aparado."""
    from shapely.geometry import Point
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    v = B - A
    L = float(np.hypot(*v))
    u = v / max(L, 1e-9)
    lb, rb = float(d["lb"]), float(d["rb"])
    larg = max(rb - lb, 1.0)

    # ---- limite de curvatura, do lado de dentro da curva
    M = 0.5 * (A + B)
    s_eixo = eixo.project(Point(*M))
    R, giro = curvatura(eixo, s_eixo, max(2.0 * larg, 40.0))
    lim_curv = K_CURV * R if np.isfinite(R) else np.inf
    # `giro` > 0: o eixo vira para a esquerda, e o lado INTERNO e o esquerdo.
    # A cutline vai da esquerda para a direita, entao o lado interno e o de
    # estaca BAIXA quando giro > 0.
    meio = 0.5 * (lb + rb)
    c_esq = lim_curv if giro > 0 else np.inf
    c_dir = lim_curv if giro < 0 else np.inf

    for k in range(9):
        folga = max(FOLGA_PISO, FOLGA_MULT * larg * (1.0 - k / 8.0))
        lim0 = max(0.0, lb - folga)
        lim1 = min(L, rb + folga)
        s0 = max(0.0, min(lim0, meio - c_esq))
        s1 = min(L, max(lim1, meio + c_dir))
        t = travessias(A, u, s0, s1, eixo)
        dentro = [x for x in t if lim0 <= x <= lim1]
        alvo = (dentro[0] if dentro else
                (min(t, key=lambda x: abs(x - meio)) if t else meio))
        for x in t:
            if x < alvo - 1e-6:
                s0 = max(s0, min(x + RECUO, lim0))
            elif x > alvo + 1e-6:
                s1 = min(s1, max(x - RECUO, lim1))
        n_prop = len(travessias(A, u, s0, s1, eixo))
        n_out = sum(len(travessias(A, u, s0, s1, e)) for e in outros)
        if n_prop == 1 and n_out == 0:
            return s0, s1, "ok"
        if n_out:
            # tambem atravessa outro ramo: encolhe pelas duas pontas
            for e in outros:
                for x in travessias(A, u, s0, s1, e):
                    if x < meio:
                        s0 = max(s0, min(x + RECUO, lim0))
                    else:
                        s1 = min(s1, max(x - RECUO, lim1))
            if (len(travessias(A, u, s0, s1, eixo)) == 1
                    and not any(travessias(A, u, s0, s1, e) for e in outros)):
                return s0, s1, "ok (outro ramo)"
    return s0, s1, "insuficiente"


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g08")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(S)}   reaches: {len(eixos)}\n")

    lim = {}
    motivos = {}
    for i, d in enumerate(S):
        ch = mapa[i]
        outros = [e for k, e in eixos.items() if k != ch]
        s0, s1, mot = faixa(d, eixos[ch], outros)
        lim[i] = [s0, s1]
        motivos[i] = mot

    # ---- vizinhas que ainda se cruzam, dentro do mesmo reach
    def linha(i):
        d = S[i]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        return LineString([A + lim[i][0] * u, A + lim[i][1] * u])

    ordem = sorted(range(len(S)), key=lambda i: (mapa[i], -S[i]["rs"]))
    encurtou = 0
    for k in range(len(ordem) - 1):
        i, j = ordem[k], ordem[k + 1]
        if mapa[i] != mapa[j]:
            continue
        for _ in range(20):
            if not linha(i).intersects(linha(j)):
                break
            a = i if (lim[i][1] - lim[i][0]) >= (lim[j][1] - lim[j][0]) else j
            d = S[a]
            lb, rb = float(d["lb"]), float(d["rb"])
            p0 = max(0.0, lb - FOLGA_PISO)
            p1 = rb + FOLGA_PISO
            s0, s1 = lim[a]
            passo = max(0.03 * (s1 - s0), 1.0)
            if (s1 - p1) >= (p0 - s0):
                if s1 - passo < p1:
                    break
                lim[a][1] = s1 - passo
            else:
                if s0 + passo > p0:
                    break
                lim[a][0] = s0 + passo
            encurtou += 1

    # ---- reescreve
    novos, tirado = {}, []
    for i, d in enumerate(S):
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        vv = B - A
        L = float(np.hypot(*vv))
        u = vv / max(L, 1e-9)
        s0, s1 = lim[i]
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= s0 - TOL) & (st <= s1 + TOL)
        ns, nz = list(st[m]), list(z[m])
        if not ns or ns[0] > s0 + TOL:
            ns.insert(0, s0)
            nz.insert(0, float(np.interp(s0, st, z)))
        if ns[-1] < s1 - TOL:
            ns.append(s1)
            nz.append(float(np.interp(s1, st, z)))
        ns = np.array(ns) - s0
        nz = np.array(nz)
        nA, nB = A + s0 * u, A + s1 * u
        Lg = round(float(np.hypot(*(nB - nA))), 2)
        if len(ns) > 1 and abs(ns[-1] - Lg) < 0.5 and Lg > ns[-2] + TOL:
            ns[-1] = Lg
        keep = [0]
        for k in range(1, len(ns)):
            if ns[k] > ns[keep[-1]] + TOL:
                keep.append(k)
        ns, nz = ns[keep], nz[keep]
        tirado.append(L - (s1 - s0))
        novos[i] = {"sta": ns, "z": nz, "lb": float(d["lb"]) - s0,
                    "rb": float(d["rb"]) - s0, "desl": s0, "cut": (nA, nB)}

    tirado = np.array(tirado)
    mex = int((tirado > 1e-6).sum())
    print(f"secoes aparadas          : {mex}   "
          f"(intocadas: {len(S)-mex})")
    print(f"encurtamentos por vizinha: {encurtou}")
    if mex:
        p = tirado[tirado > 1e-6]
        print(f"largura retirada         : mediana {np.median(p):.0f} m   "
              f"p90 {np.percentile(p,90):.0f}   max {p.max():.0f} m")
    ruins = [i for i, m in motivos.items() if m == "insuficiente"]
    print(f"secoes em que o criterio NAO foi satisfeito: {len(ruins)}")

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
                saida.append("Bank Sta=%s,%s"
                             % (_fmt(nv["lb"]), _fmt(nv["rb"])))
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
                topo = float(nv["sta"][-1])
                for t in range(0, 3 * cnt, 3):
                    val[t] = min(max(val[t] - nv["desl"], 0.0), topo)
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
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    for rot, f in (("talvegue", lambda x: float(x["z"].min())),
                   ("largura do canal", lambda x: float(x["rb"] - x["lb"]))):
        a = np.array([f(x) for x in A2])
        b = np.array([f(x) for x in B2])
        print(f"   {rot:<22} mudou no maximo "
              f"{np.abs(b-a).max():.6f}  (tem de ser zero)")
    rep = sum(1 for d in B2
              if (np.diff(np.round(np.asarray(d['sta'], float), 2)) <= 0).any())
    print(f"   secoes com estaca repetida ou fora de ordem: {rep}")
    eixos2 = ler_eixos(novo)
    mapa2 = mapa_reaches(novo)
    mult = 0
    for i, d in enumerate(B2):
        ln = LineString(np.asarray(d["cut"], float))
        x = ln.intersection(eixos2[mapa2[i]])
        n = 0 if x.is_empty else (len(x.geoms) if hasattr(x, "geoms") else 1)
        if n != 1:
            mult += 1
    print(f"   secoes que nao cruzam o proprio eixo exatamente 1x: {mult}"
          f"   (era 58)")
    la = np.array([float(x["sta"][-1]) for x in A2])
    lbb = np.array([float(x["sta"][-1]) for x in B2])
    print(f"   largura da secao: mediana {np.median(la):.0f} -> "
          f"{np.median(lbb):.0f} m   max {la.max():.0f} -> {lbb.max():.0f} m")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
