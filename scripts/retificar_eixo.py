# -*- coding: utf-8 -*-
"""Retifica o eixo onde ele serpenteia DENTRO do proprio canal.

    python scripts/retificar_eixo.py taha_ai_novo/taha_ai.g03 --saida g04

A ENTRADA NAO E TOCADA. Sai um .gXX novo. So o `Reach XY` muda: nenhuma
secao, cota, margem, Manning ou comprimento (`Type RM Length`) e tocado.

O DEFEITO

  O eixo veio do traçado fino do relevo e guarda meandros de amplitude menor
  que a LARGURA DO CANAL ESCAVADO (medido: sinuosidade ate 2,8 numa janela de
  400 m, com canal de 66 a 486 m). Uma cutline perpendicular atravessa esse
  zigue-zague 2 ou 3 vezes -- e o Validate Geometry acusa "XS must intersect
  exactly one Reach" e, com a projecao ambigua, "River Station out of order".

POR QUE PODE

  A hidraulica 1D nao le a polilinha do eixo: as distancias entre secoes sao
  as do `Type RM Length L Ch R`, que nao mudam aqui. O eixo e a moldura
  cartografica -- projecao, ordenacao e interpolacao do RAS Mapper -- e por
  isso tem de cruzar cada cutline UMA vez. Um canal de 100+ m de largura nao
  faz meandro de 50 m de raio; o eixo que o descreve tambem nao deve fazer.

O QUE SE FAZ

  Nas cutlines cruzadas != 1 vez, o trecho de eixo entre os cruzamentos
  LIMPOS vizinhos e substituido pela poligonal que liga: cruzamento limpo
  anterior -> centro do canal de cada cutline problematica -> cruzamento
  limpo seguinte. Se a emenda quebrar uma vizinha antes limpa, a janela
  cresce e a vizinha entra na emenda -- ate 10 rodadas. Os EXTREMOS do eixo
  (juncoes) nunca mudam.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                        # noqa: E402
from qc_geometria import ler_eixos                      # noqa: E402
from corrigir_cutlines import mapa_reaches, _arg        # noqa: E402
from ras_io import escrever                             # noqa: E402


def n_cruz(cut, eixo):
    from shapely.geometry import LineString
    x = LineString(cut).intersection(eixo)
    if x.is_empty:
        return 0
    return len(x.geoms) if hasattr(x, "geoms") else 1


def ponto_meio(d):
    """Centro do canal, em coordenadas do mapa."""
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
    return A + 0.5 * (float(d["lb"]) + float(d["rb"])) * u


FURO = 12.0     # m de cada lado ao furar a cutline
NUDGE = 15.0    # m de recuo da emenda em relacao ao cruzamento do vizinho


def furo(d, rumo, dist=FURO):
    """Dois pontos que FURAM a cutline pelo centro do canal.

    Passar por um ponto exatamente SOBRE a cutline e tangencia: o numero de
    cruzamentos vira loteria de arredondamento (foi medido -- a conferencia
    interna dizia 1 e o arquivo relido dizia 2 ou 0). Aqui o eixo atravessa
    na NORMAL da cutline, com `dist` m de folga de cada lado, no sentido de
    `rumo` (montante -> jusante). Quando os dois vizinhos estao do MESMO
    lado da cutline (gancho de meandro), o segmento de volta re-corta a
    cutline perto da ponta; `dist` maior empurra esse retorno para alem da
    ponta, onde a cutline ja acabou.
    """
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
    n = np.array([-u[1], u[0]])
    if np.dot(n, rumo) < 0:
        n = -n
    m = ponto_meio(d)
    return [m - dist * n, m + dist * n]


def cruzamento_unico(d, eixo):
    """O ponto do cruzamento, quando ele e um so."""
    from shapely.geometry import LineString
    x = LineString(np.asarray(d["cut"], float)).intersection(eixo)
    if x.is_empty:
        return None
    pts = list(x.geoms) if hasattr(x, "geoms") else [x]
    if len(pts) != 1 or pts[0].geom_type != "Point":
        return None
    return np.array([pts[0].x, pts[0].y])


def inversoes(secs, eixo):
    """Pares de secoes cruzadas FORA DA ORDEM do RS ao longo do eixo.

    O RAS ordena as secoes projetando o cruzamento no eixo; quando duas
    cutlines se intercalam num gancho de meandro, a de RS maior e cruzada
    DEPOIS da de RS menor e o validador acusa "River Station out of order".
    Marcar o par como ruim poe os dois na mesma emenda, que os fura na
    sequencia certa.
    """
    from shapely.geometry import LineString, Point
    est = []
    for d in secs:
        x = LineString(np.asarray(d["cut"], float)).intersection(eixo)
        if x.is_empty:
            est.append(np.nan)
            continue
        pts = list(x.geoms) if hasattr(x, "geoms") else [x]
        p = pts[0] if pts[0].geom_type == "Point" else pts[0].centroid
        est.append(eixo.project(p))
    fora = set()
    for k in range(1, len(est)):
        if not (np.isnan(est[k]) or np.isnan(est[k - 1])) \
                and est[k] <= est[k - 1]:
            fora |= {k - 1, k}
    return fora


def retificar_reach(V, secs):
    """V: vertices (n,2) do eixo.  secs: secoes do reach, na ordem do arquivo
    (RS decrescente = montante -> jusante, o sentido do eixo).
    Devolve (V_novo, n_emendas) ou (V, 0)."""
    from shapely.geometry import LineString
    eixo = LineString(V)
    ruins = {i for i, d in enumerate(secs)
             if n_cruz(np.asarray(d["cut"], float), eixo) != 1}
    ruins |= inversoes(secs, eixo)
    if not ruins:
        return V, 0

    for _ in range(10):
        eixo = LineString(V)
        # ancoras: para cada secao LIMPA, a estaca (ao longo do eixo) do seu
        # cruzamento; para as ruins, a do centro do canal projetado
        anc = np.empty(len(secs))
        pts = [None] * len(secs)
        for i, d in enumerate(secs):
            if i in ruins:
                pts[i] = ponto_meio(d)
                from shapely.geometry import Point
                anc[i] = eixo.project(Point(*pts[i]))
            else:
                p = cruzamento_unico(d, eixo)
                if p is None:          # degenerou: trata como ruim
                    ruins.add(i)
                    pts[i] = ponto_meio(d)
                    from shapely.geometry import Point
                    anc[i] = eixo.project(Point(*pts[i]))
                else:
                    pts[i] = p
                    from shapely.geometry import Point
                    anc[i] = eixo.project(Point(*p))

        # janelas: runs consecutivos de ruins
        ordem = sorted(ruins)
        runs, atual = [], [ordem[0]]
        for i in ordem[1:]:
            if i == atual[-1] + 1:
                atual.append(i)
            else:
                runs.append(atual)
                atual = [i]
        runs.append(atual)

        dist = np.hypot(*np.diff(V, axis=0).T)
        s_acum = np.concatenate([[0.0], np.cumsum(dist)])
        cortes = []      # (s_ini, s_fim, cadeia de pontos)
        for run in runs:
            i0, i1 = run[0], run[-1]
            # NUDGE: a emenda entra ANTES do cruzamento do vizinho limpo e
            # sai DEPOIS dele. Ancorar exatamente no cruzamento poe um
            # vertice SOBRE a cutline do vizinho, e o motor do RAS conta
            # vertice-na-linha como dois cruzamentos (foi medido: 9 secoes
            # vizinhas de emenda reprovaram so por isso).
            s_ini = max(anc[i0 - 1] - NUDGE, 0.0) if i0 > 0 else 0.0
            s_fim = (min(anc[i1 + 1] + NUDGE, s_acum[-1])
                     if i1 + 1 < len(secs) else s_acum[-1])
            if s_fim <= s_ini:      # ordem invertida: alarga na marra
                s_ini, s_fim = min(s_ini, s_fim), max(s_ini, s_fim)
            p_a = np.array(eixo.interpolate(s_ini).coords[0])
            p_b = np.array(eixo.interpolate(s_fim).coords[0])
            rumo = p_b - p_a
            cadeia = [q for i in run for q in furo(secs[i], rumo)]
            cortes.append((s_ini, s_fim, cadeia))
        cortes.sort(key=lambda c: c[0])

        novo = []
        s_j = 0.0
        for (s_ini, s_fim, cadeia) in cortes:
            m = (s_acum >= s_j) & (s_acum < s_ini)
            novo.extend(V[m])
            eixo_l = LineString(V)
            p_ini = np.array(eixo_l.interpolate(s_ini).coords[0])
            p_fim = np.array(eixo_l.interpolate(s_fim).coords[0])
            novo.append(p_ini)
            novo.extend(cadeia)
            novo.append(p_fim)
            s_j = s_fim
        m = s_acum >= s_j
        novo.extend(V[m])
        # extremos originais preservados
        novo[0] = V[0]
        novo[-1] = V[-1]
        NV = np.array(novo)
        # tira vertice repetido
        keep = [0]
        for t in range(1, len(NV)):
            if np.hypot(*(NV[t] - NV[keep[-1]])) > 0.01:
                keep.append(t)
        NV = NV[keep]

        eixo2 = LineString(NV)
        ruins2 = {i for i, d in enumerate(secs)
                  if n_cruz(np.asarray(d["cut"], float), eixo2) != 1}
        ruins2 |= inversoes(secs, eixo2)
        if not ruins2:
            return NV, len(runs)
        if ruins2 <= ruins:
            V = NV          # melhorou ou igual: consolida e tenta de novo
            ruins = ruins2
        else:
            ruins |= ruins2  # vizinha quebrou: entra na emenda

    # ---- fallback: janela pelos VIZINHOS, chao a chao pelo centro do canal.
    # Quando o eixo antigo tem lacos maiores que o espacamento das secoes, a
    # projecao que ancora a emenda cai no laco errado e a rodada nao fecha.
    # Aqui a ancora deixa de ser projecao: e o proprio centro do canal das
    # secoes vizinhas, e o trecho entre elas vira a poligonal
    # mid(i-1) -> mid(i) -> mid(i+1).
    from shapely.geometry import LineString, Point
    tentativas = {}
    for _ in range(12):
        eixo = LineString(V)
        ruins = sorted({i for i, d in enumerate(secs)
                        if n_cruz(np.asarray(d["cut"], float), eixo) != 1}
                       | inversoes(secs, eixo))
        if not ruins:
            return V, 1

        # o AGLOMERADO: do primeiro ruim, engole vizinhos ruins com ate 2
        # secoes limpas de intervalo -- num gancho de meandro as cutlines se
        # intercalam e consertar uma por vez so muda o lugar do defeito
        # (foi medido: 6 fora de ordem viraram 8 cruzamentos errados).
        # A emenda unica fura TODAS as secoes do aglomerado na ordem do RS.
        i0 = i1 = ruins[0]
        for r in ruins[1:]:
            if r - i1 <= 3:
                i1 = r
        grupo = list(range(i0, i1 + 1))
        tentativas[i0] = tentativas.get(i0, -1) + 1
        dist = FURO * (3 ** min(tentativas[i0], 3))

        def ancora(k):
            """Ponto de apoio na secao k: o cruzamento limpo, se existir."""
            p = cruzamento_unico(secs[k], eixo)
            return p if p is not None else ponto_meio(secs[k])

        mids = [ponto_meio(secs[k]) for k in grupo]
        cadeia = []
        for t, k in enumerate(grupo):
            ref_a = mids[t - 1] if t > 0 else None
            ref_b = mids[t + 1] if t + 1 < len(mids) else None
            if ref_a is None and ref_b is None:
                rumo = np.array([1.0, 0.0])
            elif ref_a is None:
                rumo = ref_b - mids[t]
            elif ref_b is None:
                rumo = mids[t] - ref_a
            else:
                rumo = ref_b - ref_a
            cadeia += furo(secs[k], rumo, dist)

        if i0 == 0:
            # comeca na cabeceira: o extremo de montante pode andar, e o
            # eixo passa a NASCER 40 m acima do canal da primeira secao
            u = cadeia[0] - cadeia[1]
            u = u / max(np.hypot(*u), 1e-9)
            inicio = [cadeia[0] + 40.0 * u]
            antes = []
            s_a = None
        else:
            m0 = ancora(i0 - 1)
            s_a = max(eixo.project(Point(*m0)) - NUDGE, 0.0)
            antes = [np.asarray(p, float) for p in eixo.coords
                     if eixo.project(Point(*p)) < s_a - 0.01]
            inicio = [np.array(eixo.interpolate(s_a).coords[0])]
        if i1 + 1 < len(secs):
            m1 = ancora(i1 + 1)
            s_b = min(eixo.project(Point(*m1)) + NUDGE, eixo.length)
            if s_a is not None and s_b <= s_a:
                s_b = min(s_a + 2 * NUDGE, eixo.length)
            depois = [np.asarray(p, float) for p in eixo.coords
                      if eixo.project(Point(*p)) > s_b + 0.01]
            fim = [np.array(eixo.interpolate(s_b).coords[0])]
        else:
            depois = []
            fim = [np.asarray(eixo.coords[-1], float)]
        V = np.array(antes + inicio + cadeia + fim + depois)
        V[-1] = np.asarray(eixo.coords[-1], float)
        keep = [0]
        for t in range(1, len(V)):
            if np.hypot(*(V[t] - V[keep[-1]])) > 0.01:
                keep.append(t)
        V = V[keep]

    eixo = LineString(V)
    resta = ({i for i, d in enumerate(secs)
              if n_cruz(np.asarray(d["cut"], float), eixo) != 1}
             | inversoes(secs, eixo))
    return V, (-1 if resta else 1)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g04")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    por_reach = {}
    for d, ch in zip(S, mapa):
        por_reach.setdefault(ch, []).append(d)

    novos_eixos = {}
    for ch, eixo in eixos.items():
        V = np.asarray(eixo.coords, float)
        secs = por_reach.get(ch, [])
        if not secs:
            continue
        NV, n = retificar_reach(V.copy(), secs)
        if n == 0:
            continue
        rot = "NAO convergiu" if n < 0 else f"{n} emenda(s)"
        print(f"   {ch[0]:14s} {ch[1]:3s}  vertices {len(V)} -> {len(NV)}   "
              f"{rot}")
        novos_eixos[ch] = NV

    if not novos_eixos:
        print("nenhum eixo precisou de emenda")
        return

    # ------------------------------------------------ reescreve os Reach XY
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    saida, j, ch = [], 0, None
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            ch = (p[0].strip(), p[1].strip())
        if l.startswith("Reach XY=") and ch in novos_eixos:
            NV = novos_eixos[ch]
            saida.append(f"Reach XY= {len(NV)} ")
            plano = NV.reshape(-1)
            lin = ""
            for t, x in enumerate(plano):
                lin += "%16.4f" % x
                if (t + 1) % 4 == 0:
                    saida.append(lin)
                    lin = ""
            if lin:
                saida.append(lin)
            cnt = int(l.split("=")[1])
            j += 1
            lidos = 0
            while j < len(linhas) and lidos < 2 * cnt:
                x = linhas[j]
                if not x.strip() or x[:1].isalpha():
                    break
                lidos += len(x) // 16
                j += 1
            continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    # -------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    S2 = ler_secoes(novo)
    E2 = ler_eixos(novo)
    M2 = mapa_reaches(novo)
    mult = 0
    for d, ch2 in zip(S2, M2):
        n = n_cruz(np.asarray(d["cut"], float), E2[ch2])
        if n != 1:
            mult += 1
            print(f"   ainda != 1 cruzamento: {ch2[0]} RS {d['rs']:.2f} ({n}x)")
    print(f"   secoes cruzando o proprio eixo != 1 vez: {mult}   "
          "(tem de ser zero)")
    por2 = {}
    for d, ch2 in zip(S2, M2):
        por2.setdefault(ch2, []).append(d)
    fora = sum(len(inversoes(v, E2[k])) for k, v in por2.items())
    print(f"   secoes cruzadas fora da ordem do RS: {fora}   "
          "(tem de ser zero)")
    for chq, eixo in E2.items():
        Va = np.asarray(eixos[chq].coords, float)
        Vb = np.asarray(eixo.coords, float)
        # o extremo de MONTANTE pode andar (cabeceira); o de JUSANTE e
        # juncao ou foz, e nao pode
        if not np.allclose(Va[-1], Vb[-1], atol=0.01):
            print(f"   EXTREMO DE JUSANTE MUDOU em {chq} -- nao podia")


if __name__ == "__main__":
    main(sys.argv[1:])
