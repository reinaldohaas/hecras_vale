# -*- coding: utf-8 -*-
"""Reamostra o perfil do canal e poe as margens no topo do encaixe.

    python scripts/refazer_perfil_canal.py modelo/mirim_t30/mirim_t30.g19 \
        --rs-min 12773 --rs-max 20274 --saida g23

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

POR QUE O PERFIL PRECISA SER REFEITO ANTES DAS MARGENS

  As secoes do canal foram afinadas para 76 pontos, com ~14 m entre pontos na
  planicie. A margem do canal tem uns 10 m de projecao horizontal, entao ela
  CAIU FORA do perfil: numa secao tipica ele vai de 2,41 m na estaca 460 para
  0,06 m na 476 e cai direto para o leito. Nao ha margem no arquivo -- ha um
  degrau. Nao adianta mover `Bank Sta` sobre um perfil que nao tem o que
  marcar.

  Aqui a faixa central e reamostrada a 2 m no MDT SIG-SC 1 m, e a margem passa
  a existir como geometria. Fora da faixa o perfil antigo e mantido.

AS MARGENS SAO MEDIDAS, E NAO ARBITRADAS

  Medido a 2 m em seis secoes do canal, o encaixe NAO e uniforme: a lamina tem
  46 a 64 m (os 45 m informados descrevem bem a agua), mas o encaixe ate o topo
  da margem vai de 60 a 186 m. Fixar 144 m para todas seria trocar um numero
  errado por outro.

  Cada margem e posta no TOPO DO ENCAIXE daquela secao: andando para fora a
  partir da agua, o primeiro ponto que alcanca a cota da planicie menos
  `FOLGA_TOPO`. A planicie e a mediana do proprio perfil alem de +-100 m.

O LEITO CONTINUA PROVISORIO, E SO ELE

  Onde o MDT ve AGUA -- cota ate `LAMINA_TOL` acima do minimo central -- ele
  esta medindo a superficie livre, e nao o fundo. Ali o perfil mantem a cota
  provisoria que ja existia. Todo o resto passa a ser terreno medido.

  Ou seja: a margem e o encaixe entram como dado real; o fundo continua sendo
  a interpolacao entre as duas cotas medidas, esperando levantamento.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qc_secoes import ler_secoes                        # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio     # noqa: E402
from ras_io import escrever                             # noqa: E402

FAIXA = 120.0       # m para cada lado do centro, reamostrados
PASSO = 2.0         # m entre pontos novos
LAMINA_TOL = 0.25   # m acima do minimo central ainda e agua
FOLGA_TOPO = 0.50   # m abaixo da planicie ja conta como topo de margem
TOL = 0.005


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


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    rmin = _arg(argv, "--rs-min", 12773.0, float)
    rmax = _arg(argv, "--rs-max", 20274.0, float)
    ext = _arg(argv, "--saida", "g23")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S_arq = ler_secoes(entrada)
    alvo = [k for k, d in enumerate(S_arq) if rmin <= float(d["rs"]) <= rmax]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes do canal (RS {rmin:.0f}..{rmax:.0f}): {len(alvo)}")
    if not alvo:
        raise SystemExit("nenhuma secao no intervalo")

    off = np.arange(-FAIXA, FAIXA + PASSO / 2, PASSO)
    pts = []
    for k in alvo:
        d = S_arq[k]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        c = A + 0.5 * (float(d["lb"]) + float(d["rb"])) * u
        for o in off:
            pts.append(c + o * u)
    pts = np.array(pts)
    bb = (pts[:, 0].min() - 40, pts[:, 1].min() - 40,
          pts[:, 0].max() + 40, pts[:, 1].max() + 40)
    print(f"MDT: {len(pts)} pontos a {PASSO:g} m sobre "
          f"{len(tiles_do_dominio(bb))} folhas")
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(bb))
    Z = mdt.cota(pts[:, 0], pts[:, 1]).reshape(len(alvo), len(off))

    novos, larg_ant, larg_novo, npt = {}, [], [], []
    sem = 0
    for q, k in enumerate(alvo):
        d = S_arq[k]
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        lb, rb = float(d["lb"]), float(d["rb"])
        sc = 0.5 * (lb + rb)
        zn = Z[q]
        if not np.isfinite(zn).any():
            sem += 1
            continue
        ref = np.nanmedian(np.r_[zn[off < -100], zn[off > 100]])
        agua = np.nanmin(zn[np.abs(off) <= 40])
        e_agua = np.isfinite(zn) & (zn <= agua + LAMINA_TOL)
        # o leito provisorio que ja estava la, sob a agua
        z_leito = float(z[(st >= lb - TOL) & (st <= rb + TOL)].min())
        # topo do encaixe: primeiro ponto, para fora, que alcanca a planicie
        i0 = int(np.flatnonzero(e_agua).min())
        i1 = int(np.flatnonzero(e_agua).max())
        alvo_topo = ref - FOLGA_TOPO
        le = i0
        while le > 0 and not (np.isfinite(zn[le]) and zn[le] >= alvo_topo):
            le -= 1
        ri = i1
        while ri < len(off) - 1 and not (np.isfinite(zn[ri])
                                         and zn[ri] >= alvo_topo):
            ri += 1
        # perfil novo: fora da faixa fica o antigo; dentro, o MDT, com o
        # leito provisorio preservado sob a agua
        s_novo = sc + off
        z_novo = zn.copy()
        z_novo[e_agua] = z_leito
        bom = np.isfinite(z_novo)
        s_novo, z_novo = s_novo[bom], z_novo[bom]
        fora = (st < sc - FAIXA - TOL) | (st > sc + FAIXA + TOL)
        ns = np.concatenate([st[fora & (st < sc)], s_novo,
                             st[fora & (st > sc)]])
        nz = np.concatenate([z[fora & (st < sc)], z_novo,
                             z[fora & (st > sc)]])
        o = np.argsort(ns)
        ns, nz = ns[o], nz[o]
        keep = [0]
        for j in range(1, len(ns)):
            if ns[j] > ns[keep[-1]] + TOL:
                keep.append(j)
        ns, nz = ns[keep], nz[keep]
        # margens sobre estacas EXISTENTES do perfil novo
        nl = ns[int(np.argmin(np.abs(ns - (sc + off[le]))))]
        nr = ns[int(np.argmin(np.abs(ns - (sc + off[ri]))))]
        if nr <= nl:
            continue
        larg_ant.append(rb - lb)
        larg_novo.append(nr - nl)
        npt.append(len(ns))
        novos[k] = {"sta": ns, "z": nz, "lb": nl, "rb": nr,
                    "lb0": lb, "rb0": rb}

    larg_ant = np.array(larg_ant)
    larg_novo = np.array(larg_novo)
    print(f"\nsecoes refeitas          : {len(novos)}"
          + (f"   (sem MDT: {sem})" if sem else ""))
    print(f"pontos por secao         : mediana {np.median(npt):.0f}   "
          f"max {max(npt)}   (limite do HEC-RAS: 500)")
    print(f"largura do CANAL         : {np.median(larg_ant):.0f} -> "
          f"{np.median(larg_novo):.0f} m   "
          f"(p10 {np.percentile(larg_novo,10):.0f}, "
          f"p90 {np.percentile(larg_novo,90):.0f}, max {larg_novo.max():.0f})")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
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
                for t in range(0, 3 * cnt, 3):
                    if abs(val[t] - nv["lb0"]) < 1e-6:
                        val[t] = nv["lb"]
                    elif abs(val[t] - nv["rb0"]) < 1e-6:
                        val[t] = nv["rb"]
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
    za = np.array([float(x["z"].min()) for x in A2])
    zb = np.array([float(x["z"].min()) for x in B2])
    print(f"   talvegue mudou         : max {np.abs(zb-za).max():.6f} m "
          "(tem de ser zero -- o leito nao foi tocado)")
    la = np.array([float(x["sta"][-1]) for x in A2])
    lb2 = np.array([float(x["sta"][-1]) for x in B2])
    print(f"   largura da secao mudou : max {np.abs(lb2-la).max():.6f} m "
          "(tem de ser zero)")
    fora = 0
    for d in B2:
        st = np.asarray(d["sta"], float)
        for b in (d["lb"], d["rb"]):
            if np.abs(st - float(b)).min() > 1e-6:
                fora += 1
                break
    print(f"   Bank Sta fora de estaca existente: {fora}  (tem de ser zero)")
    rep = sum(1 for d in B2
              if (np.diff(np.round(np.asarray(d['sta'], float), 2)) <= 0).any())
    print(f"   estaca repetida ou fora de ordem : {rep}")
    n = max(len(x["sta"]) for x in B2)
    print(f"   maior contagem de pontos         : {n}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
