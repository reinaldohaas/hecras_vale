# -*- coding: utf-8 -*-
"""Liga os rios avulsos numa REDE, com a topologia do legado.

    python scripts/construir_rede.py
    python scripts/construir_rede.py --saida modelo/itajai_rede/itajai_rede.g01

Cada rio foi construido sozinho (`construir_rio.py`), com contorno provisorio.
Este script os costura numa unica geometria de rede, do jeito que o modelo
legado que roda -- `legado/Itajai_Rede_1983.g01` -- os liga:

    Rio_do_Sul   Itajai_Sul  + Itajai_Oeste  -> Itajai_Acu (cabeceira)
    Ibirama      Itajai_Norte                -> Itajai_Acu
    Indaial      Rio_Benedito                -> Itajai_Acu
    Itajai       Itajai_Mirim                -> Itajai_Acu -> mar

O TRONCO E O ACU, PARTIDO NAS JUNCOES

  No legado o Acu e um so rio em 4 reaches (R1..R4), cortado nos 3 pontos onde
  um afluente entra. Aqui o corte e o MESMO -- lido das faixas de RS de cada
  reach do Acu no legado, sem inventar estaca --, mas so nas juncoes ATIVAS.

  Uma juncao so entra se o afluente dela tem geometria pronta. O Benedito
  ainda nao tem g02 (o eixo de montante sobe a encosta -- ver
  `diagnostico_benedito.py`), entao a juncao Indaial ficaria com UMA entrada e
  UMA saida, o que o HEC-RAS recusa. Nesse caso o tronco NAO se parte ali: o
  Acu passa reto por Indaial, e a rede sai com os rios que fecham. Quando o
  Benedito ganhar eixo e g02, ele entra sozinho -- o script se ajusta.

O QUE NAO E INVENTADO

  As coordenadas das juncoes, os comprimentos `Junc L&A` e os pontos de corte
  do Acu vem todos do legado. As secoes e os eixos sao os que ja estavam nos
  g02/g01 de cada rio. A ultima secao de cada reach que morre numa juncao tem
  L Ch R = 0, como no legado -- a juncao e que conduz dali para a frente.
"""
import argparse
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from ras_io import escrever    # noqa: E402

LEGADO = "legado/Itajai_Rede_1983.g01"
TRONCO = "Itajai_Acu"

# rio (nome no legado) -> geometria pronta. .g02 (com batimetria) de preferencia.
FONTES = {
    "Itajai_Acu":   ["modelo/itajai_acu/itajai_acu.g02",
                     "modelo/itajai_acu/itajai_acu.g01"],
    "Itajai_Norte": ["modelo/itajai_norte/itajai_norte.g02",
                     "modelo/itajai_norte/itajai_norte.g01"],
    "Itajai_Oeste": ["modelo/itajai_oeste/itajai_oeste.g02",
                     "modelo/itajai_oeste/itajai_oeste.g01"],
    "Itajai_Sul":   ["modelo/itajai_sul/itajai_sul.g02",
                     "modelo/itajai_sul/itajai_sul.g01"],
    "Itajai_Mirim": ["modelo/itajai_mirim/itajai_mirim.g02",
                     "modelo/itajai_mirim/itajai_mirim.g01"],
    # Benedito SO com g02: o g01 cru sobe a encosta a montante (ver
    # diagnostico_benedito.py) e derrubaria a rede inteira. Enquanto nao houver
    # g02, ele nao entra, e a juncao Indaial fica inativa.
    "Rio_Benedito": ["modelo/rio_benedito/rio_benedito.g02"],
}


def campo16(nome):
    return f"{nome:<16.16}"


# ------------------------------------------------------------------ leitura
def ler_excluir(path):
    """Le um CSV river;reach;rs de secoes a NAO incluir na rede.

    Sao as secoes que o proprio validador do HEC-RAS rejeita na costura -- as
    que, na confluencia, cruzam o reach vizinho. Removê-las nao inventa nada:
    so tira a sobreposicao que a juncao ja conduz. Casa por (river, rs), que e
    unico dentro de um rio mesmo depois de partir o tronco em R1..Rk.
    """
    # river -> lista de RS a descartar. O casamento e por TOLERANCIA: o
    # validador do RAS mostra a RS com 1 decimal (140150.5) e a geometria tem 2
    # (140150.51); um casamento exato erraria por 0,01 e a secao ficaria. Como
    # as secoes distam ~150 m, 0,6 m de folga identifica sem ambiguidade.
    excl = {}
    if path and os.path.exists(path):
        import csv
        for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";"):
            excl.setdefault(r["river"].strip(), []).append(float(r["rs"]))
    return excl


def _excluida(river, rs, excl):
    return any(abs(rs - e) < 0.6 for e in excl.get(river, ()))


def ler_geom(path, excl=None):
    """Devolve (preambulo, [reach]). Cada reach: dict com river, reach, xy_head
    (linhas do bloco 'Reach XY='), xy (Nx2), secoes [ (rs, [linhas]) ].

    `excl` = dict river -> [rs] a descartar (casamento por tolerancia)."""
    excl = excl or {}
    t = open(path, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i0 = next(i for i, l in enumerate(t) if l.startswith("River Reach="))
    preamb = t[:i0]
    reaches = []
    i = i0
    while i < len(t):
        if not t[i].startswith("River Reach="):
            i += 1
            continue
        river, reach = [x.strip() for x in t[i].split("=", 1)[1].split(",")]
        j = i + 1
        xy_head = []
        while j < len(t) and not t[j].startswith("Type RM"):
            xy_head.append(t[j])
            j += 1
        xy = _parse_xy(xy_head)
        secoes = []
        while j < len(t) and not t[j].startswith("River Reach="):
            if t[j].startswith("Type RM"):
                rs = float(t[j].split(",")[1])
                bloco = [t[j]]
                j += 1
                while (j < len(t) and not t[j].startswith("Type RM")
                       and not t[j].startswith("River Reach=")):
                    bloco.append(t[j])
                    j += 1
                # tira brancos do fim do bloco (ficam entre secoes)
                while bloco and not bloco[-1].strip():
                    bloco.pop()
                if not _excluida(river, rs, excl):
                    secoes.append((rs, bloco))
            else:
                j += 1
        secoes.sort(key=lambda s: -s[0])
        reaches.append({"river": river, "reach": reach, "xy_head": xy_head,
                        "xy": xy, "secoes": secoes})
        i = j
    return preamb, reaches


def _parse_xy(xy_head):
    if not xy_head or not xy_head[0].startswith("Reach XY"):
        return np.zeros((0, 2))
    v = []
    for L in xy_head[1:]:
        if "=" in L or not L.strip():      # 'Rch Text X Y=0,0' e brancos
            continue
        v += [float(L[c:c + 16]) for c in range(0, len(L), 16)
              if L[c:c + 16].strip()]
    a = np.array(v)
    return a.reshape(-1, 2) if len(a) >= 2 else np.zeros((0, 2))


def emit_xy(pts):
    """Bloco 'Reach XY= N' + coordenadas em colunas de 16, 4 por linha,
    fechado com a linha 'Rch Text X Y=0,0' que o HEC-RAS espera."""
    out = [f"Reach XY= {len(pts)} "]
    flat = pts.reshape(-1)
    linha = ""
    for i, x in enumerate(flat):
        linha += "%16.4f" % x
        if (i + 1) % 4 == 0:
            out.append(linha)
            linha = ""
    if linha:
        out.append(linha)
    out.append("Rch Text X Y=0,0")
    return out


# ------------------------------------------------------------- topologia
def ler_juncoes(path):
    """[ {name, xy, ups:[(river,reach)], dn:(river,reach), la:[float]} ]."""
    t = open(path, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    juncts, i = [], 0
    while i < len(t):
        if t[i].startswith("Junct Name="):
            j = {"name": t[i].split("=", 1)[1].strip(), "xy": None,
                 "ups": [], "dn": None, "la": []}
            k = i + 1
            while k < len(t) and not t[k].startswith("Junct Name=") \
                    and not t[k].startswith("River Reach="):
                l = t[k]
                if l.startswith("Junct X Y"):
                    v = l.split("=", 1)[1].split(",")
                    j["xy"] = (float(v[0]), float(v[1]))
                elif l.startswith("Up River,Reach="):
                    j["ups"].append(tuple(x.strip() for x in
                                          l.split("=", 1)[1].split(",")))
                elif l.startswith("Dn River,Reach="):
                    j["dn"] = tuple(x.strip() for x in
                                    l.split("=", 1)[1].split(","))
                elif l.startswith("Junc L&A="):
                    j["la"].append(float(l.split("=", 1)[1].split(",")[0]))
                k += 1
            juncts.append(j)
            i = k
        else:
            i += 1
    return juncts


def faixas_tronco(path, tronco):
    """RS (max,min) de cada reach do tronco no legado, de montante p/ jusante."""
    t = open(path, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    cur, rs = None, {}
    for l in t:
        if l.startswith("River Reach="):
            cur = tuple(x.strip() for x in l.split("=", 1)[1].split(","))
        elif l.startswith("Type RM") and cur and cur[0] == tronco:
            rs.setdefault(cur[1], []).append(float(l.split(",")[1]))
    faixas = [(r, max(v), min(v)) for r, v in rs.items()]
    faixas.sort(key=lambda f: -f[1])
    return faixas    # [(reach, rsmax, rsmin), ...]


# ------------------------------------------------------------- montagem
def zera_ultima(secoes):
    """L Ch R = 0 na ultima secao (ela morre numa juncao ou na foz)."""
    if not secoes:
        return
    rs, bloco = secoes[-1]
    h = bloco[0]
    pre = h.split("=", 1)[0]
    campos = [x.strip() for x in h.split("=", 1)[1].split(",")]
    campos[2] = campos[3] = campos[4] = "0.00"
    bloco[0] = f"{pre}= {campos[0]} ,{campos[1]},0.00,0.00,0.00"


def reach_block(river, reach, xy_pts, secoes):
    out = [f"River Reach={campo16(river)},{campo16(reach)}"]
    out += emit_xy(xy_pts) if len(xy_pts) else ["Reach XY= 0 "]
    for _, bloco in secoes:
        out += bloco
    out.append("")
    return out


def parte_tronco(acu, splits):
    """acu (reach unico) -> lista de reaches R1..Rk cortados nas RS de `splits`.

    `splits` = RS crescentes de corte (fronteiras internas). Cada sub-reach
    leva as secoes na sua faixa e um Reach XY recortado no ponto de corte.
    """
    sec = acu["secoes"]                       # ja ordenadas RS desc
    xy = acu["xy"]
    bordas = [np.inf] + sorted(splits, reverse=True) + [-np.inf]
    subs = []
    for a, b in zip(bordas[:-1], bordas[1:]):
        s = [(rs, bl) for rs, bl in sec if b < rs <= a] if a != np.inf \
            else [(rs, bl) for rs, bl in sec if rs > b]
        if s:
            subs.append(s)
    # Reach XY de cada sub: recorta a polilinha nos pontos mais proximos do
    # centro da primeira e da ultima secao do sub-reach.
    reaches = []
    for k, s in enumerate(subs, 1):
        c0 = _centro(s[0][1])
        c1 = _centro(s[-1][1])
        i0 = int(np.argmin(np.hypot(xy[:, 0] - c0[0], xy[:, 1] - c0[1]))) \
            if len(xy) else 0
        i1 = int(np.argmin(np.hypot(xy[:, 0] - c1[0], xy[:, 1] - c1[1]))) \
            if len(xy) else 0
        lo, hi = min(i0, i1), max(i0, i1)
        sub_xy = xy[lo:hi + 1] if len(xy) and hi > lo else \
            np.array([c0, c1])
        reaches.append({"reach": f"R{k}", "secoes": s, "xy": sub_xy})
    return reaches


def _centro(bloco):
    for i, l in enumerate(bloco):
        if l.startswith("XS GIS Cut Line"):
            pts = []
            j = i + 1
            while j < len(bloco) and not bloco[j].lstrip()[:1].isalpha() \
                    and not bloco[j].startswith("#"):
                L = bloco[j]
                pts += [float(L[c:c + 16]) for c in range(0, len(L), 16)
                        if L[c:c + 16].strip()]
                j += 1
            a = np.array(pts).reshape(-1, 2)
            return a.mean(0) if len(a) else np.array([0.0, 0.0])
    return np.array([0.0, 0.0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--saida", default="modelo/itajai_rede/itajai_rede.g01")
    ap.add_argument("--legado", default=LEGADO)
    ap.add_argument("--excluir", default=None,
                    help="CSV river;reach;rs de secoes a descartar na costura")
    a = ap.parse_args()
    excl = ler_excluir(a.excluir)
    if excl:
        n = sum(len(v) for v in excl.values())
        print(f"excluindo {n} secao(oes) da costura (validador do RAS)")

    disp = {}
    for rio, cs in FONTES.items():
        for c in cs:
            if os.path.exists(c):
                disp[rio] = c
                break
    print("geometrias disponiveis:")
    for rio, c in disp.items():
        print(f"   {rio:14} {c}  ({'batimetria' if c.endswith('g02') else 'MDT cru'})")
    if TRONCO not in disp:
        raise SystemExit(f"sem geometria do tronco {TRONCO}")

    juncts = ler_juncoes(a.legado)
    faixas = faixas_tronco(a.legado, TRONCO)
    print(f"\ntronco {TRONCO} no legado: {len(faixas)} reaches")
    for r, mx, mn in faixas:
        print(f"   {r}: RS {mx:.0f} -> {mn:.0f}")

    # ---- quais juncoes estao ativas (afluente nao-tronco disponivel)
    # e em que RS do tronco elas cortam.
    #
    # fronteira interna = entre reach do tronco e o seguinte, RS media.
    fronteira = {}                      # dn-reach do tronco -> RS de corte
    for (r, mx, mn), (r2, mx2, mn2) in zip(faixas, faixas[1:]):
        fronteira[r2] = (mn + mx2) / 2.0
    ativas, splits = [], []
    for j in juncts:
        trib = [u for u in j["ups"] if u[0] != TRONCO]
        up_tronco = [u for u in j["ups"] if u[0] == TRONCO]
        falta = [u for u in trib if u[0] not in disp]
        if falta:
            print(f"\njuncao {j['name']}: afluente sem geometria "
                  f"({', '.join(u[0] for u in falta)}) -- o tronco passa reto")
            continue
        ativas.append(j)
        if up_tronco and j["dn"] and j["dn"][0] == TRONCO:
            rs = fronteira.get(j["dn"][1])
            if rs is not None:
                splits.append(rs)
    splits = sorted(set(splits), reverse=True)
    print(f"\njuncoes ativas: {', '.join(j['name'] for j in ativas)}")
    print(f"cortes do tronco em RS: "
          f"{', '.join('%.0f' % s for s in splits) or '(nenhum)'}")

    # ---- monta reaches
    preamb0, rA = ler_geom(disp[TRONCO], excl)
    acu = rA[0]
    sub_acu = parte_tronco(acu, splits)
    # renumera juncoes -> nomes de reach do tronco (R1..Rk de montante p/ jus.)
    # mapa: dn-reach original do legado -> novo Rk. Vamos so usar ordem.
    print(f"\n{TRONCO} partido em {len(sub_acu)} reach(es): "
          + ", ".join(f"{s['reach']}({len(s['secoes'])} sec)" for s in sub_acu))

    corpo = []
    # tronco
    for s in sub_acu[:-1]:
        zera_ultima(s["secoes"])
    for s in sub_acu:
        corpo += reach_block(TRONCO, s["reach"], s["xy"], s["secoes"])
    # afluentes
    trib_reaches = {}      # (river,reach) -> presente
    for rio, c in disp.items():
        if rio == TRONCO:
            continue
        _, rr = ler_geom(c, excl)
        for r in rr:
            zera_ultima(r["secoes"])
            corpo += reach_block(rio, "R1", r["xy"], r["secoes"])
            trib_reaches[(rio, "R1")] = True

    # ---- juncoes: reescreve up/dn do tronco para os R1..Rk novos
    # constroi um mapa RS-> reach do tronco novo
    def tronco_reach_por_rs(rs):
        # devolve o Rk cujo intervalo contem rs
        bordas = [np.inf] + splits + [-np.inf]
        for k, (a2, b2) in enumerate(zip(bordas[:-1], bordas[1:]), 1):
            if b2 < rs <= a2 or (a2 == np.inf and rs > b2):
                return f"R{k}"
        return f"R{len(sub_acu)}"

    jblocks = []
    for j in ativas:
        trib = [u for u in j["ups"] if u[0] != TRONCO]
        up_tronco = [u for u in j["ups"] if u[0] == TRONCO]
        # RS de corte desta juncao (se ha up do tronco)
        rs_corte = fronteira.get(j["dn"][1]) if j["dn"][0] == TRONCO else None
        ups_novo = []
        for u in trib:
            ups_novo.append((u[0], "R1"))
        for u in up_tronco:
            # o reach do tronco a MONTANTE do corte
            ups_novo.append((TRONCO, tronco_reach_por_rs(rs_corte + 1)
                             if rs_corte else "R1"))
        if j["dn"][0] == TRONCO:
            dn = (TRONCO, tronco_reach_por_rs((rs_corte - 1) if rs_corte
                                              else faixas[0][2]))
        else:
            dn = j["dn"]
        jblocks.append((j, ups_novo, dn))

    # ---- preambulo: usa o do tronco, recomputa Viewing Rectangle sobre tudo
    todos_xy = [s["xy"] for s in sub_acu if len(s["xy"])]
    for rio, c in disp.items():
        if rio == TRONCO:
            continue
        _, rr = ler_geom(c, excl)
        for r in rr:
            if len(r["xy"]):
                todos_xy.append(r["xy"])
    P = np.vstack(todos_xy)
    preamb = []
    for l in preamb0:
        if l.startswith("Viewing Rectangle"):
            preamb.append(f"Viewing Rectangle= {P[:,0].min():.6f} , "
                          f"{P[:,0].max():.6f} , {P[:,1].max():.6f} , "
                          f"{P[:,1].min():.6f} ")
        elif l.startswith("Geom Title"):
            preamb.append("Geom Title=Itajai_Rede (MDT + batimetria 1983)")
        else:
            preamb.append(l)

    out = list(preamb)
    for j, ups, dn in jblocks:
        out.append(f"Junct Name={campo16(j['name'])}")
        out.append("Junct Desc=Confluencia, 0 , 0 , 0 ,0")
        if j["xy"]:
            x, y = j["xy"]
            out.append(f"Junct X Y & Text X Y={x:.2f},{y:.2f},"
                       f"{x+800:.2f},{y+800:.2f}")
        for u in ups:
            out.append(f"Up River,Reach={campo16(u[0])},{campo16(u[1])}")
        out.append(f"Dn River,Reach={campo16(dn[0])},{campo16(dn[1])}")
        for _ in ups:
            out.append("Junc L&A=150.00,0")
        out.append("")
    out += corpo

    os.makedirs(os.path.dirname(a.saida) or ".", exist_ok=True)
    escrever(a.saida, "\n".join(out))
    n_sec = sum(len(s["secoes"]) for s in sub_acu) + \
        sum(1 for _ in ())  # tronco
    print(f"\nrede -> {a.saida}")
    print("juncoes escritas:")
    for j, ups, dn in jblocks:
        print(f"   {j['name']:12} {' + '.join(f'{u[0]},{u[1]}' for u in ups)}"
              f"  ->  {dn[0]},{dn[1]}")
    return a.saida


if __name__ == "__main__":
    main()
