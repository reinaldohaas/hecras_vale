# -*- coding: utf-8 -*-
"""Encaixa o reach do Canal Retificado na geometria ja partida em juncoes.

    python scripts/inserir_canal.py modelo/mirim_t30/mirim_t30.g04 --saida g05

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

=====================================================================
  O LEITO DESTE ARQUIVO E PROVISORIO. NAO USE OS RESULTADOS.
=====================================================================

  A batimetria do canal nao existe. Aqui o leito e uma INTERPOLACAO LINEAR
  entre as duas unicas cotas medidas -- -0,76 m na ponta de montante e
  -2,68 m na de jusante, ambas vindas das secoes do proprio modelo -- posta
  plana ao longo dos 45 m de canal.

  Isso NAO e uma estimativa da batimetria. E o minimo necessario para que o
  HEC-RAS aceite o reach e para que se possa VERIFICAR O ENCANAMENTO: se as
  duas juncoes estao bem formadas, se a vazao se divide e se reencontra, se o
  balanco de volume fecha. Nivel, velocidade e divisao de vazao produzidos por
  este arquivo nao valem como resultado.

  Existe porque o HEC-RAS recusa juncao com um reach entrando e um saindo
  ("Junctions are for flow confluences and splits"). Sem o canal, a partição
  em tres reaches nao computa, e nao havia como testa-la antes do
  levantamento chegar. Quando a batimetria vier, so as cotas dentro do canal
  mudam -- o resto deste arquivo continua valendo.

O QUE VEM DE ONDE

  tracado e cutlines   `doc/canal/canal_secoes.csv`, gerado de
                       `canal_itajai_mirim.geojson` (OpenStreetMap)
  planicie e margens   MDT SIG-SC 1 m -- dado real
  leito                INTERPOLADO entre duas cotas medidas -- provisorio
  Manning              copiado da secao do modelo na bifurcacao
                       (0,063 nas planicies, 0,035 no canal)
  HTab                 2 cm acima do leito, incremento 0,100 e 500 pontos,
                       a mesma convencao do resto do modelo
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever    # noqa: E402

L16 = 16
CSV = "doc/canal/canal_secoes.csv"
EIXO = "doc/canal/canal_eixo.geojson"
Z_MONT, Z_JUS = -0.76, -2.68
N_PLANICIE = 34        # pontos por planicie depois do afinamento
N_CANAL = 8            # pontos dentro do canal
N_MANN, N_CANAL_MANN = 0.063, 0.035
HTAB_INCR, HTAB_N = 0.100, 500


def _pad(s):
    return f"{s:<{L16}}"[:L16]


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


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


def secoes_do_canal():
    """Le o CSV e devolve as secoes prontas, com o leito interpolado."""
    por_rs = defaultdict(list)
    for r in csv.DictReader(open(CSV, encoding="utf-8"), delimiter=";"):
        por_rs[float(r["rs"])].append(r)
    rss = sorted(por_rs, reverse=True)
    rmax, rmin = max(rss), min(rss)
    saida = []
    for rs in rss:
        L = por_rs[rs]
        st = np.array([float(x["estaca"]) for x in L])
        z = np.array([float(x["z"]) if x["z"] else np.nan for x in L])
        canal = np.array([x["origem"] == "A LEVANTAR" for x in L])
        xy = np.array([[float(x["x"]), float(x["y"])] for x in L])
        # leito provisorio: reta entre as duas cotas medidas
        f = (rmax - rs) / max(rmax - rmin, 1e-9)
        zb = Z_MONT + (Z_JUS - Z_MONT) * f
        z[canal] = zb
        # afinamento: a planicie nao precisa de um ponto a cada 2 m
        esq = np.flatnonzero(~canal & (st < st[canal].min()))
        dir_ = np.flatnonzero(~canal & (st > st[canal].max()))
        ic = np.flatnonzero(canal)
        sel = np.unique(np.r_[
            esq[np.linspace(0, len(esq) - 1, N_PLANICIE).astype(int)],
            ic[np.linspace(0, len(ic) - 1, N_CANAL).astype(int)],
            dir_[np.linspace(0, len(dir_) - 1, N_PLANICIE).astype(int)]])
        st2, z2 = st[sel], z[sel]
        ok = np.isfinite(z2)
        st2, z2 = st2[ok], z2[ok]
        base = st2[0]
        saida.append({
            "rs": rs, "sta": st2 - base, "z": z2,
            "lb": st[canal].min() - base, "rb": st[canal].max() - base,
            "cut": (xy[0], xy[-1]), "zb": zb})
    return saida


def bloco_secao(d, comp):
    b = [f"Type RM Length L Ch R = 1 ,{d['rs']:.2f},{comp:8.2f},"
         f"{comp:8.2f},{comp:8.2f}",
         f"Bank Sta={_fmt(d['lb'])},{_fmt(d['rb'])}",
         "XS GIS Cut Line= 2",
         "".join("%16.2f" % x for x in (d["cut"][0][0], d["cut"][0][1],
                                        d["cut"][1][0], d["cut"][1][1]))]
    v = []
    for a, c in zip(d["sta"], d["z"]):
        v += [a, c]
    b.append("#Sta/Elev= %d " % len(d["sta"]))
    b += _col(v, 8, 2)
    b.append("#Mann= 3 , 0 , 0 ")
    mv = [0.0, N_MANN, 0, float(d["lb"]), N_CANAL_MANN, 0,
          float(d["rb"]), N_MANN, 0]
    lin = ""
    for t, x in enumerate(mv):
        lin += ("%8.2f" % x if t % 3 == 0 else
                "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
    b.append(lin)
    b.append(f"XS HTab Starting El and Incr={d['zb']+0.02:.2f},"
             f"{HTAB_INCR:.3f}, {HTAB_N} ")
    b.append("XS HTab Horizontal Distribution=-1,-1,-1")
    b.append("Exp/Cntr=0.3,0.1")
    b.append("")
    return b


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g05")
    nome = _arg(argv, "--reach", "R1")
    # O CANAL PRECISA DE NOME DE RIO PROPRIO, e nao so de reach proprio.
    # Com os dois ramos sob 'Itajai_Mirim' o RAS Mapper costura as bank lines
    # POR NOME DE RIO e gera DOIS jogos completos, cada um varrendo os 141 km
    # -- medido: 4 polilinhas de 146, 153, 133 e 135 km onde deviam existir
    # duas. Toda secao passa a cruzar 4 bank lines em vez de 2, e o Validate
    # Geometry acusa "XS intersects > 2 banklines" em 1.401 das 1.469 secoes,
    # mais "Multiple upstream/downstream connections to Rivers with the same
    # name". Com nome proprio a bifurcacao fica sem ambiguidade.
    rio_canal = _arg(argv, "--rio", "Canal_Retif")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = secoes_do_canal()
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"canal: {len(S)} secoes   RS {S[0]['rs']:.0f} a {S[-1]['rs']:.0f}")
    print(f"   pontos por secao: {np.median([len(d['sta']) for d in S]):.0f}"
          f"   (limite do HEC-RAS: 500)")
    print(f"   leito PROVISORIO de {S[0]['zb']:+.2f} a {S[-1]['zb']:+.2f} m "
          "-- interpolado entre as duas cotas medidas")
    larg = np.array([float(d["rb"] - d["lb"]) for d in S])
    print(f"   canal: {np.median(larg):.1f} m   "
          f"secao: {np.median([d['sta'][-1] for d in S]):.0f} m")

    t = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    rio = next(l.split("=", 1)[1].split(",")[0].strip()
               for l in t if l.startswith("River Reach="))

    # ---- as juncoes ganham o canal
    saida, i = [], 0
    while i < len(t):
        l = t[i]
        if l.startswith("Junct Name="):
            j = i
            blk = []
            while j < len(t) and (t[j].strip() or not blk):
                blk.append(t[j])
                j += 1
                if j < len(t) and t[j].startswith(("Junct Name=",
                                                   "River Reach=")):
                    break
            jn = blk[0].split("=", 1)[1].strip()
            ups = [x for x in blk if x.startswith("Up River,Reach=")]
            dns = [x for x in blk if x.startswith("Dn River,Reach=")]
            cab = [x for x in blk if not x.startswith(("Up River,Reach=",
                                                       "Dn River,Reach=",
                                                       "Junc L&A="))]
            if jn.startswith("Bifurcacao"):
                dns.append(f"Dn River,Reach={_pad(rio_canal)},{_pad(nome)}")
            else:
                ups.append(f"Up River,Reach={_pad(rio_canal)},{_pad(nome)}")
            # uma linha por par (montante, jusante)
            la = [f"Junc L&A={comp},0" for comp in
                  ["100.00"] * (len(ups) * len(dns))]
            cab = [x for x in cab if x.strip()]
            saida += cab + ups + dns + la + [""]
            print(f"   juncao '{jn.strip()}': {len(ups)} entra(m), "
                  f"{len(dns)} sai(em), {len(la)} par(es) de comprimento")
            i = j
            continue
        saida.append(l)
        i += 1

    # ---- o reach do canal, no fim
    eixo = np.asarray(json.load(open(EIXO))["features"][0]
                      ["geometry"]["coordinates"], float)
    saida.append(f"River Reach={_pad(rio_canal)},{_pad(nome)}")
    saida.append(f"Reach XY= {len(eixo)} ")
    lin = ""
    for k, (x, y) in enumerate(eixo):
        lin += "%16.4f%16.4f" % (x, y)
        if (k + 1) % 2 == 0:
            saida.append(lin)
            lin = ""
    if lin:
        saida.append(lin)
    saida.append("Rch Text X Y=0,0,0,0")
    saida.append("")
    for k, d in enumerate(S):
        comp = (S[k]["rs"] - S[k + 1]["rs"]) if k < len(S) - 1 else 0.0
        saida += bloco_secao(d, comp)

    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    t2 = open(novo, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    a = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    for chave, esp in (("River Reach=", 4), ("Junct Name=", 2)):
        print(f"   {chave:<16} {a.count(chave)} -> {t2.count(chave)} "
              f"(esperado {esp})")
    nxs_a = len(re.findall(r"(?m)^Type RM Length L Ch R", a))
    nxs_b = len(re.findall(r"(?m)^Type RM Length L Ch R", t2))
    print(f"   secoes           {nxs_a} -> {nxs_b} "
          f"(esperado {nxs_a}+{len(S)}={nxs_a+len(S)})")
    npt = [len(d["sta"]) for d in S]
    print(f"   pontos por secao do canal: min {min(npt)} max {max(npt)}")
    mono = all((np.diff(d["sta"]) > 0).all() for d in S)
    print(f"   estacas estritamente crescentes em todas: {mono}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
