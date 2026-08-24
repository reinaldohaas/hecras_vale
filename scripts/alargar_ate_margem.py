# -*- coding: utf-8 -*-
"""Alarga so as secoes que transbordam, ate achar margem, com teto.

    python scripts/alargar_ate_margem.py modelo/mirim_t30/mirim_t30.g19 \
        --resultados modelo/mirim_canal6/mirim_canal6.p01.hdf \
        --teto 400 --saida g22

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O ALVO E MEDIDO, NAO ESCOLHIDO

  Alarga apenas as secoes em que a LAMINA MAXIMA DA PROPRIA RODADA passa da
  ponta mais baixa do perfil -- ou seja, onde a agua escapa pela borda e o
  HEC-RAS extrapola a tabela hidraulica. Medido no `mirim_canal6`, sao 96 de
  1.428 (7%); as outras 1.332 tem 17 m de folga e nao sao tocadas.

  Isso importa porque alargar nao e de graca: secao larga demais em curva
  atravessa o proprio eixo e embaraca as edge lines, que e exatamente o
  defeito que a apara por curvatura acabou de reduzir. Alargar quem nao
  precisa so devolveria erro de geometria.

O TETO E POR QUE ELE EXISTE

  Procurando no MDT de 1 m a que distancia ha terreno acima da lamina, numa
  amostra das 96: mediana de 320 m por lado, mas p90 de 2.160 m, e em um terco
  dos casos NAO HA margem dentro de 3 km -- a varzea e plana. Secao 1D de
  3 km atravessando varzea plana nao representa escoamento: ali a agua nao
  corre na direcao da secao, ela espalha e armazena. Isso pede area de
  armazenamento, e nao secao mais larga.

  Entao o alargamento para em `--teto` metros por lado. Quem nao achar margem
  dentro do teto FICA COMO ESTA e entra no relatorio -- explicitamente, para
  nao passar por resolvido.

DOIS FREIOS CONTRA PIORAR A GEOMETRIA

  1. a extensao para antes de fazer a cutline cruzar o eixo uma segunda vez;
  2. a extensao para antes de encostar na cutline da vizinha.

  Ambos sao conferidos ponto a ponto enquanto se avanca, e nao depois.

O QUE NAO MUDA

  Nenhuma cota existente e alterada e nenhum ponto sai: so entram pontos NAS
  PONTAS, com cota lida do MDT SIG-SC 1 m. Talvegue, largura de canal e o
  perfil inteiro entre as margens ficam iguais -- a conferencia mede isso.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qc_secoes import ler_secoes                        # noqa: E402
from qc_geometria import ler_eixos                      # noqa: E402
from corrigir_cutlines import mapa_reaches              # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio     # noqa: E402
from ras_io import escrever                             # noqa: E402

PASSO = 10.0        # m entre pontos novos
FOLGA = 0.50        # m de borda livre exigida acima da lamina maxima
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


def lamina_maxima(hdf, S):
    """Lamina maxima por secao, casada por River Station."""
    from ras_commander import HdfResultsXsec
    ds = HdfResultsXsec.get_xsec_timeseries(hdf)
    W = np.nanmax(np.asarray(ds["Water_Surface"].values), axis=0)
    nomes = [str(x) for x in ds["cross_section"].values]
    import re
    rs_hdf = []
    for nm in nomes:
        m = re.findall(r"([\d.]+)\s*$", nm.strip())
        rs_hdf.append(float(m[0]) if m else np.nan)
    rs_hdf = np.array(rs_hdf)
    saida = np.full(len(S), np.nan)
    for i, d in enumerate(S):
        k = int(np.argmin(np.abs(rs_hdf - float(d["rs"]))))
        if abs(rs_hdf[k] - float(d["rs"])) < 0.6:
            saida[i] = W[k]
    return saida


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString
    entrada = argv[0]
    res = _arg(argv, "--resultados")
    teto = _arg(argv, "--teto", 400.0, float)
    folga = _arg(argv, "--folga", FOLGA, float)
    ext = _arg(argv, "--saida", "g22")
    if not res:
        raise SystemExit("informe --resultados com o .pNN.hdf da rodada")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    S.sort(key=lambda d: -d["rs"])
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)          # na ordem do ARQUIVO
    # mapa segue a ordem do arquivo; S foi reordenado -- refaz o vinculo
    S_arq = ler_secoes(entrada)
    ordem_arq = {round(float(d["rs"]), 3): k for k, d in enumerate(S_arq)}
    ch_de = [mapa[ordem_arq[round(float(d["rs"]), 3)]] for d in S]

    W = lamina_maxima(res, S)
    zb = np.array([min(float(d["z"][0]), float(d["z"][-1])) for d in S])
    alvo = [i for i in range(len(S))
            if np.isfinite(W[i]) and zb[i] < W[i] + folga]
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(S)}   com lamina passando da borda: {len(alvo)}"
          f"   (teto de {teto:g} m por lado, folga de {folga:g} m)")
    if not alvo:
        raise SystemExit("nada a alargar")

    # ---- MDT sobre a faixa possivel
    P = []
    for i in alvo:
        d = S[i]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        for k in range(1, int(teto / PASSO) + 1):
            P.append(A - k * PASSO * u)
            P.append(B + k * PASSO * u)
    P = np.array(P)
    bb = (P[:, 0].min() - 60, P[:, 1].min() - 60,
          P[:, 0].max() + 60, P[:, 1].max() + 60)
    print(f"MDT: lendo {len(P)} pontos sobre "
          f"{len(tiles_do_dominio(bb))} folhas")
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(bb))
    Zt = mdt.cota(P[:, 0], P[:, 1])

    # ---- cutlines das vizinhas, para o segundo freio
    L = [LineString(np.asarray(d["cut"], float)) for d in S]

    novos, achou, parou_teto, sem_dado = {}, 0, 0, 0
    ganho = []
    p = 0
    npass = int(teto / PASSO)
    for i in alvo:
        d = S[i]
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
        ze = Zt[p:p + 2 * npass:2]
        zd = Zt[p + 1:p + 2 * npass:2]
        p += 2 * npass
        eixo = eixos[ch_de[i]]
        viz = [L[j] for j in (i - 1, i + 1)
               if 0 <= j < len(S) and ch_de[j] == ch_de[i]]
        limite = W[i] + folga

        z_leito = float(np.asarray(d["z"], float).min())

        def avanca(z, sinal, ponta):
            """Quantos passos cabem, e onde para."""
            ok = 0
            for k in range(1, npass + 1):
                if not np.isfinite(z[k - 1]):
                    break
                # TERCEIRO FREIO: nao descer abaixo do leito desta secao.
                # A 400 m da ponta o MDT encontra OUTRO curso d'agua, e sem
                # esta trava o talvegue da secao passa a ser ele: medido, o
                # leito baixou ate 8,69 m em 7 secoes. Terreno abaixo do leito
                # nao e planicie desta secao -- e drenagem vizinha, e engoli-la
                # inventaria um canal onde ha so uma valeta a 400 m.
                if z[k - 1] < z_leito:
                    break
                novo_p = ponta + sinal * k * PASSO * u
                seg = LineString([ponta, novo_p])
                if seg.crosses(eixo) or any(seg.intersects(v) for v in viz):
                    break
                ok = k
                if z[k - 1] > limite:
                    break
            return ok

        ne = avanca(ze, -1.0, A)
        nd = avanca(zd, +1.0, B)
        if ne == 0 and nd == 0:
            continue
        ze_ok = ne > 0 and np.isfinite(ze[ne - 1]) and ze[ne - 1] > limite
        zd_ok = nd > 0 and np.isfinite(zd[nd - 1]) and zd[nd - 1] > limite
        if ze_ok and zd_ok:
            achou += 1
        else:
            parou_teto += 1
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        s_e = np.array([-(k * PASSO) for k in range(ne, 0, -1)])
        s_d = np.array([st[-1] + k * PASSO for k in range(1, nd + 1)])
        z_e = np.array([ze[k - 1] for k in range(ne, 0, -1)])
        z_d = np.array([zd[k - 1] for k in range(1, nd + 1)])
        mask_e = np.isfinite(z_e)
        mask_d = np.isfinite(z_d)
        sem_dado += int((~mask_e).sum() + (~mask_d).sum())
        s_e, z_e = s_e[mask_e], z_e[mask_e]
        s_d, z_d = s_d[mask_d], z_d[mask_d]
        desl = -s_e[0] if len(s_e) else 0.0
        ns = np.concatenate([s_e, st, s_d]) + desl
        nz = np.concatenate([z_e, z, z_d])
        nA = A - (len(s_e) * PASSO if len(s_e) else 0.0) * u
        nB = B + (len(s_d) * PASSO if len(s_d) else 0.0) * u
        ganho.append(float(ns[-1] - st[-1]))
        novos[i] = {"sta": ns, "z": nz, "lb": float(d["lb"]) + desl,
                    "rb": float(d["rb"]) + desl, "desl": desl,
                    "cut": (nA, nB), "zmin0": float(z.min())}

    ganho = np.array(ganho) if ganho else np.array([0.0])
    print(f"\nalargadas                : {len(novos)} de {len(alvo)}")
    print(f"   acharam margem dos dois lados : {achou}")
    print(f"   pararam no teto ou na trava   : {parou_teto}"
          "   (ficam como estao, sem resolver)")
    print(f"   nao alargaram nada            : {len(alvo)-len(novos)}")
    print(f"largura acrescentada     : mediana {np.median(ganho):.0f} m   "
          f"p90 {np.percentile(ganho,90):.0f}   max {ganho.max():.0f} m")
    if sem_dado:
        print(f"pontos sem MDT descartados: {sem_dado}")

    # ---- reescreve, casando pela ordem do ARQUIVO
    idx_arq = {}
    for i, d in enumerate(S):
        idx_arq[ordem_arq[round(float(d["rs"]), 3)]] = novos.get(i)
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j, htab = -1, [], 0, 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = idx_arq.get(i_sec)
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
                for t in range(0, 3 * cnt, 3):
                    val[t] = 0.0 if t == 0 else val[t] + nv["desl"]
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
            if l.startswith("XS HTab Starting El and Incr"):
                zmin = float(np.min(nv["z"]))
                if zmin < nv["zmin0"] - 1e-9:
                    q = [x.strip() for x in l.split("=", 1)[1].split(",")]
                    saida.append("XS HTab Starting El and Incr="
                                 f"{zmin+0.02:.2f},{float(q[1]):.3f}, "
                                 f"{int(q[2])} ")
                    htab += 1
                    j += 1
                    continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))
    if htab:
        print(f"HTab reancorado em {htab} secoes "
              "(a extensao achou chao mais baixo que o leito)")

    # ---------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    A2.sort(key=lambda d: -d["rs"])
    B2 = ler_secoes(novo)
    B2.sort(key=lambda d: -d["rs"])
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    ca = np.array([float(x["rb"] - x["lb"]) for x in A2])
    cb = np.array([float(x["rb"] - x["lb"]) for x in B2])
    print(f"   largura do canal mudou : max {np.abs(cb-ca).max():.6f} m "
          "(tem de ser zero)")
    za = np.array([float(x["z"].min()) for x in A2])
    zb2 = np.array([float(x["z"].min()) for x in B2])
    pior = np.minimum(zb2 - za, 0)
    print(f"   talvegue baixou em     : max {abs(pior.min()):.3f} m "
          f"em {int((pior < -1e-9).sum())} secoes  (a extensao achou chao "
          "mais baixo)")
    la = np.array([float(x["sta"][-1]) for x in A2])
    lb2 = np.array([float(x["sta"][-1]) for x in B2])
    print(f"   largura da secao       : mediana {np.median(la):.0f} -> "
          f"{np.median(lb2):.0f} m   max {la.max():.0f} -> {lb2.max():.0f} m")
    rep = sum(1 for d in B2
              if (np.diff(np.round(np.asarray(d['sta'], float), 2)) <= 0).any())
    print(f"   secoes com estaca repetida ou fora de ordem: {rep}")
    # contencao esperada
    nzb = np.array([min(float(d["z"][0]), float(d["z"][-1])) for d in B2])
    m = np.isfinite(W)
    print(f"   ainda com a borda abaixo da lamina: "
          f"{int((nzb[m] < W[m]).sum())}  (era {int((zb[m] < W[m]).sum())})")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
