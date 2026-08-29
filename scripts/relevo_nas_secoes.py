# -*- coding: utf-8 -*-
"""Troca o RELEVO das secoes: perfil reamostrado do MDT 1 m do SIG-SC.

    python scripts/relevo_nas_secoes.py taha_ai.g01.antes_do_reparo_1983 \
        --saida r00 [--figura doc/figuras/relevo_trocado.png]

A ENTRADA NAO E TOCADA. Sai um .gXX novo -- a ETAPA 0 do construir_1983:
o esqueleto (eixos, cutlines, RS, juncoes, u01) fica; so as COTAS mudam,
de Copernicus 30 m para SIG-SC 1 m.

COMO

  Cada secao e reamostrada ao longo da PROPRIA `XS GIS Cut Line`
  (passo ~2 m, reduzido a <=450 pontos por mediana em blocos). O zero do
  SIG-SC e vazio ([[curar-zero]]): ponto sem dado herda a cota do perfil
  ANTIGO interpolada na mesma estacao -- nunca zero.

  O CANAL E REESCAVADO: o laser ve a lamina d'agua, nao o fundo. Dentro
  do vao `Bank Sta` original a cota nova nao pode ficar ACIMA da lamina
  do MDT nem o talvegue subir: o fundo do canal desce ao talvegue
  ORIGINAL da secao (a batimetria sintetica calibrada), rampando das
  margens. Fora dos bancos vale o MDT puro.

  CONFERENCIA: contagem de secoes intacta, talvegue de cada secao <=
  talvegue antigo + 0,01, nenhum ponto em zero absoluto, e uma figura
  com 6 transectos sorteados (antigo x novo) para o olho conferir.
"""
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402

MAX_PONTOS = 450
PASSO = 2.0
LARGURAS = "doc/larguras_g95"   # lamina congelada (a mesma da campanha)
LARG_MIN = 12.0
BORDA = 20.0        # m de esponja em cada emenda SIG-SC <-> Copernicus:
                    # sem ela cada vazio do SIG-SC (colarinho de folha)
                    # vira DEGRAU vertical entre os dois datums
COLAR = 0.20        # m: quase-zero onde o antigo diz terra alta = colarinho


_larg_cache = {}


def larg_lamina(rio, km):
    """lamina mediana-movel (5 km) do CSV congelado; LARG_MIN sem dado."""
    if rio not in _larg_cache:
        import csv
        arq = os.path.join(LARGURAS, f"{rio}.csv")
        pares = []
        if os.path.exists(arq):
            for r in csv.reader(open(arq, encoding="utf-8"),
                                delimiter=";"):
                if r[0] == "dist_foz_km" or not r[1]:
                    continue
                pares.append((float(r[0]), float(r[1])))
        pares.sort()
        if len(pares) >= 3:
            dd = np.array([p[0] for p in pares])
            ww = np.array([p[1] for p in pares])
            suave = np.array([np.median(ww[max(0, i - 5):i + 6])
                              for i in range(len(ww))])
            _larg_cache[rio] = (dd, suave)
        else:
            _larg_cache[rio] = None
    m = _larg_cache[rio]
    if m is None:
        return LARG_MIN
    return float(max(LARG_MIN, np.interp(km, m[0], m[1])))


def ler_cutlines(g01):
    """indice da secao (ordem de Type RM) -> [(x,y), ...]"""
    linhas = open(g01, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out, i_sec, j = {}, -1, 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        if l.startswith("XS GIS Cut Line="):
            n = int(l.split("=")[1])
            vals, j2 = [], j + 1
            while j2 < len(linhas) and len(vals) < 2 * n:
                x = linhas[j2]
                if not x.strip() or x[:1].isalpha():
                    break
                vals += [float(x[c:c + 16]) for c in range(0, len(x), 16)
                         if x[c:c + 16].strip()]
                j2 += 1
            out[i_sec] = [(vals[k], vals[k + 1])
                          for k in range(0, len(vals) - 1, 2)]
            j = j2
            continue
        j += 1
    return out


def amostrar_cutline(pontos, comprimento_alvo):
    """estacoes (0..L) e coordenadas a cada PASSO, reescalonadas para o
    comprimento da secao original (cutline e secao podem divergir ~1%)."""
    P = np.asarray(pontos, float)
    seg = np.hypot(*np.diff(P, axis=0).T)
    L = float(seg.sum())
    s_acum = np.concatenate([[0.0], np.cumsum(seg)])
    s = np.arange(0.0, L + PASSO / 2, PASSO)
    xs = np.interp(s, s_acum, P[:, 0])
    ys = np.interp(s, s_acum, P[:, 1])
    est = s * (comprimento_alvo / max(L, 1e-9))
    return est, xs, ys


def reduzir(est, z, max_p=MAX_PONTOS):
    if len(est) <= max_p:
        return est, z
    blocos = np.array_split(np.arange(len(est)), max_p)
    e2 = np.array([est[b].mean() for b in blocos])
    z2 = np.array([np.median(z[b]) for b in blocos])
    # o ponto mais baixo do bloco nao pode sumir na mediana
    k = int(np.argmin([z[b].min() for b in blocos]))
    z2[k] = z[blocos[k]].min()
    return e2, z2


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "r00")
    fig_arq = _arg(argv, "--figura", "doc/figuras/relevo_trocado.png")
    # --fundo talvegue: cava ao talvegue ANTIGO (sintetico calibrado).
    # --fundo lamina: cava a `--profundidade` m abaixo do ESPELHO do
    #   MDT -- o espelho e agua real e desce monotonico; e o unico fundo
    #   sao sem batimetria nos trechos onde o sintetico e escadaria
    #   (Cedros restaurado explodia com 20 m de erro pelo leito antigo)
    fundo = _arg(argv, "--fundo", "talvegue")
    prof = _arg(argv, "--profundidade", 2.0, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    cut = ler_cutlines(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes: {len(S)}   cutlines: {len(cut)}")

    # dominio -> tiles
    xs_all = [p[0] for c in cut.values() for p in c]
    ys_all = [p[1] for c in cut.values() for p in c]
    bbox = (min(xs_all) - 100, min(ys_all) - 100,
            max(xs_all) + 100, max(ys_all) + 100)
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(bbox))
    print(f"tiles do MDT no dominio: {len(mdt.caminhos)}")

    novos, sem_dado, exemplos = {}, 0, []
    for i, d in enumerate(S):
        if d["tipo"] != "1" or i not in cut:
            continue
        sta0 = np.asarray(d["sta"], float)
        z0 = np.asarray(d["z"], float)
        est, cx, cy = amostrar_cutline(cut[i], float(sta0[-1] - sta0[0]))
        est = est + float(sta0[0])
        znovo = mdt.cota(cx, cy)
        z_ant = np.interp(est, sta0, z0)
        vazio = ~np.isfinite(znovo)
        # colarinho de folha: quase-zero onde o antigo diz terra alta
        vazio |= (np.abs(np.nan_to_num(znovo)) < COLAR) & (z_ant > 2.0)
        if vazio.all():
            sem_dado += 1
            continue
        # vazio herda o perfil antigo (Copernicus) na mesma estacao...
        znovo[vazio] = z_ant[vazio]
        # ...e a EMENDA leva esponja de BORDA m: peso do SIG-SC cai
        # linearmente a zero ao se aproximar de um trecho vazio
        if vazio.any() and not vazio.all():
            dist = np.full(len(est), np.inf)
            iv = np.flatnonzero(vazio)
            for k in np.flatnonzero(~vazio):
                dist[k] = np.min(np.abs(est[k] - est[iv]))
            w = np.clip(dist / BORDA, 0.0, 1.0)
            znovo = w * znovo + (1.0 - w) * z_ant
        est, znovo = reduzir(est, znovo)
        # reescavar o canal: dentro dos bancos originais, fundo nao sobe.
        # NAO e um espigao de 1 ponto (r00/r01: canal de ~10 m carregando
        # o Acu inteiro, explosao no canion) -- cava-se um fundo PLANO ao
        # talvegue antigo na LARGURA DA LAMINA medida do rio, centrado no
        # ponto mais fundo do MDT, com paredes em rampa 1:2
        m = (est >= d["lb"] - 1e-6) & (est <= d["rb"] + 1e-6)
        tal0 = float(z0[(sta0 >= d["lb"]) & (sta0 <= d["rb"])].min()) \
            if ((sta0 >= d["lb"]) & (sta0 <= d["rb"])).any() \
            else float(z0.min())
        if m.any():
            i_min = np.flatnonzero(m)[int(np.argmin(znovo[m]))]
            c = est[i_min]
            W = larg_lamina(d["rio"], d["rs"] / 1000.0)
            if fundo == "lamina":
                alvo_f = float(znovo[i_min]) - prof
            else:
                alvo_f = tal0
            dcen = np.abs(est - c)
            piso = alvo_f + np.maximum(0.0, (dcen - W / 2.0)) * 0.5
            cava = m & (dcen <= W / 2.0 + 8.0)
            znovo[cava] = np.minimum(znovo[cava], piso[cava])
        est = np.round(est, 2)
        novos[i] = {"sta": est, "z": znovo,
                    "htab": float(znovo.min()) + 0.15,
                    "lb": float(est[int(np.argmin(np.abs(est - d['lb'])))]),
                    "rb": float(est[int(np.argmin(np.abs(est - d['rb'])))])}
        if len(exemplos) < 6 and i % max(1, len(S) // 7) == 0:
            exemplos.append((d, sta0, z0, est, znovo))
        if i % 200 == 0:
            print(f"   ... {i}/{len(S)}")
    mdt.fechar()
    print(f"reamostradas: {len(novos)}   sem dado no MDT: {sem_dado}")

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
            if l.startswith("#Mann="):
                partes = l.split("=")[1].split(",")
                cnt = int(partes[0])
                j += 1
                vals = []
                while j < len(linhas) and len(vals) < 3 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    vals += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                             if x[c:c + 8].strip()]
                    j += 1
                if cnt == 3 and len(vals) == 9:
                    vals[0] = float(nv["sta"][0])
                    vals[3] = nv["lb"]
                    vals[6] = nv["rb"]
                saida.append("#Mann= %d , %s , 0 "
                             % (cnt, partes[1].strip()
                                if len(partes) > 1 else "0"))
                saida += _col(vals, 8, 3)
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%.2f,%.2f" % (nv["lb"], nv["rb"]))
                j += 1
                continue
            if l.startswith("XS HTab Starting El and Incr="):
                resto = l.split("=", 1)[1].split(",")
                saida.append("XS HTab Starting El and Incr=%.2f,%s,%s"
                             % (nv["htab"], resto[1], resto[2]))
                j += 1
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    sobe, zerado = 0, 0
    for i, (a, b) in enumerate(zip(S, B)):
        if i not in novos:
            continue
        za = np.asarray(a["z"], float)
        zb = np.asarray(b["z"], float)
        if zb.min() > za.min() + 0.01:
            sobe += 1
        if (np.abs(zb) < 1e-6).any():
            zerado += 1
    if fundo == "lamina":
        print(f"   talvegue subiu vs sintetico: {sobe}   (esperado em "
              f"--fundo lamina: o leito novo e espelho-{prof:.1f} m)")
    else:
        print(f"   talvegue subiu: {sobe}   (tem de ser 0)")
    print(f"   pontos em zero absoluto: {zerado}   (tem de ser 0)")

    if exemplos:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, eixos = plt.subplots(3, 2, figsize=(14, 10))
        for ax, (d, s0, z0, s1, z1) in zip(eixos.flat, exemplos):
            ax.plot(s0, z0, "-", color="0.6", lw=1,
                    label="Copernicus (antigo)")
            ax.plot(s1, z1, "-", color="crimson", lw=1.2,
                    label="SIG-SC 1 m (novo)")
            ax.set_title(f"{d['rio']} {d['reach']} RS {d['rs']:.0f}",
                         fontsize=9)
            ax.grid(alpha=0.3)
            ax.tick_params(labelsize=8)
        eixos.flat[0].legend(fontsize=8)
        fig.suptitle("Relevo trocado: 6 secoes sorteadas, antigo x novo")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        os.makedirs(os.path.dirname(fig_arq) or ".", exist_ok=True)
        fig.savefig(fig_arq, dpi=130)
        print(f"   figura de inspecao: {fig_arq}")


if __name__ == "__main__":
    main(sys.argv[1:])
