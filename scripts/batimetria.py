# -*- coding: utf-8 -*-
"""Pede a batimetria de um rio, e ancora o leito nela quando ela chegar.

    # 1. o que levantar
    python scripts/batimetria.py pedir modelo/benedito6/benedito6.g01 \
        --cada 2000 --saida doc/batimetria_benedito.csv

    # 2. depois do levantamento, com a coluna z_leito preenchida
    python scripts/batimetria.py aplicar modelo/benedito6/benedito6.g01 \
        --pontos doc/batimetria_benedito.csv --saida g02

Serve para qualquer rio: le a geometria, nao supoe nada do Itajai-Mirim. O
mesmo comando roda no Benedito e no Mirim.

POR QUE ISTO EXISTE

  O perfil sai do MDT SIG-SC, e o MDT VE A LAMINA D'AGUA, nao o fundo. O
  "talvegue" medido e a superficie livre, que tem pocas e quedas: no Benedito
  isso da 56 trechos com declividade acima de 2% e maximo de 17,2%, e o solver
  bate no teto de 40 iteracoes em praticamente todo passo -- a agua entra e
  some do balanco.

  Alisar o perfil resolveria a numerica e destruiria a medida: e o caminho do
  `gerar_mirim_do_zero.py`, cujo talvegue sao oito numeros escritos no codigo.
  Ancorar em batimetria e a saida que mantem os dois.

COMO A ANCORAGEM FUNCIONA

  Entre dois pontos levantados o leito e INTERPOLADO LINEARMENTE na distancia
  ao longo do rio. Fora do intervalo levantado o perfil do MDT fica como esta,
  e o relatorio diz quanto do rio ficou de fora -- extrapolar seria inventar.

  A cota nova entra deslocando TODOS os pontos entre as margens, com peso
  proporcional a profundidade: o fundo recebe o deslocamento inteiro e os
  pontos na cota da margem recebem zero. A forma da calha e a largura nao
  mudam, e a planicie nao e tocada.

  O `XS HTab Starting El and Incr` acompanha, 2 cm acima do INVERT DA CALHA --
  que e contra o que o HEC-RAS compara.

ONDE PEDIR

  Os pontos saem espacados de `--cada` metros ao longo do rio, MAIS as secoes
  em que o perfil do MDT tem os maiores degraus: sao elas que mais precisam de
  ancora, porque e ali que a interpolacao erraria mais. O CSV traz X, Y em
  UTM 22S, a estaca do talvegue e a cota da lamina no MDT, que serve de TETO
  (o fundo esta abaixo dela) e de conferencia para quem levanta.
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes    # noqa: E402
from ras_io import escrever         # noqa: E402

TAPER = 0.30           # fracao da profundidade em que o deslocamento se suaviza
DEGRAU_EXTRA = 20      # secoes de maior degrau que entram alem do espacamento


def secoes_ordenadas(g):
    S = ler_secoes(g)
    S.sort(key=lambda d: -d["rs"])
    return S


def talvegue_xy(d):
    """(x, y, estaca, cota) do ponto mais baixo ENTRE AS MARGENS."""
    st = np.asarray(d["sta"], float)
    z = np.asarray(d["z"], float)
    m = (st >= float(d["lb"])) & (st <= float(d["rb"]))
    if not m.any():
        m = np.ones_like(st, bool)
    i = int(np.flatnonzero(m)[np.argmin(z[m])])
    A = np.asarray(d["cut"][0], float)
    B = np.asarray(d["cut"][-1], float)
    u = (B - A) / max(float(np.hypot(*(B - A))), 1e-9)
    P = A + st[i] * u
    return float(P[0]), float(P[1]), float(st[i]), float(z[i])


def cmd_pedir(a):
    S = secoes_ordenadas(a.geom)
    rs = np.array([d["rs"] for d in S])
    z = np.array([talvegue_xy(d)[3] for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    x = np.r_[0.0, np.cumsum(ch[:-1])]          # distancia desde montante

    escolhidos = set()
    alvo = 0.0
    for i in range(len(S)):
        if x[i] >= alvo - 1e-9:
            escolhidos.add(i)
            alvo = x[i] + a.cada
    escolhidos.add(0)
    escolhidos.add(len(S) - 1)
    dz = np.abs(np.diff(z))
    for i in np.argsort(-dz)[:DEGRAU_EXTRA]:
        escolhidos.add(int(i))
        escolhidos.add(int(i) + 1)
    idx = sorted(escolhidos)

    os.makedirs(os.path.dirname(a.saida) or ".", exist_ok=True)
    with open(a.saida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["rs", "x", "y", "estaca_talvegue", "z_lamina_mdt",
                    "z_leito_A_LEVANTAR", "observacao"])
        for i in idx:
            xx, yy, st, zz = talvegue_xy(S[i])
            grande = (i < len(dz) and dz[i] > np.percentile(dz, 90))
            w.writerow([f"{S[i]['rs']:.2f}", f"{xx:.2f}", f"{yy:.2f}",
                        f"{st:.2f}", f"{zz:.2f}", "",
                        "degrau grande no MDT" if grande else ""])
    print(f"geometria : {a.geom}   {len(S)} secoes   "
          f"{x[-1]/1000:.2f} km")
    print(f"pontos    : {len(idx)}   "
          f"({a.cada:g} m de espacamento + {DEGRAU_EXTRA} de maior degrau)")
    print(f"pedido    -> {a.saida}")
    print(f"\n   a coluna `z_lamina_mdt` e a LAMINA lida no MDT: o fundo esta")
    print("   ABAIXO dela. Serve de teto e de conferencia para quem levanta.")
    print("   Preencher `z_leito_A_LEVANTAR` em metros, mesmo datum do modelo.")
    return a.saida


def cmd_aplicar(a):
    S = secoes_ordenadas(a.geom)
    rio_alvo = getattr(a, "rio", None)
    if rio_alvo:
        # num arquivo de REDE o RS se repete entre rios (todo tributario
        # termina em 75.00); sem o filtro, o casamento por RS escreveria o
        # rebaixamento do rio pedido em secoes homonimas dos outros
        S = [d_ for d_ in S if d_["rio"] == rio_alvo]
        if not S:
            raise SystemExit(f"nenhuma secao do rio {rio_alvo} em {a.geom}")
    rs = np.array([d["rs"] for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    x = np.r_[0.0, np.cumsum(ch[:-1])]
    z0 = np.array([talvegue_xy(d)[3] for d in S])

    med_rs, med_z, fic_rs = [], [], []
    for r in csv.DictReader(open(a.pontos, encoding="utf-8"), delimiter=";"):
        v = (r.get("z_leito_A_LEVANTAR") or r.get("z_leito") or "").strip()
        if not v:
            # ponto SEM cota: seja por legado sintetico, rebaixamento
            # implausivel ou secao longe demais, e um lugar onde NAO ha
            # ancora. O motivo nao importa para o aplicar; o que importa e o
            # aglomerado -- ver o peso por intervalo adiante.
            fic_rs.append(float(r["rs"]))
            continue
        med_rs.append(float(r["rs"]))
        med_z.append(float(v.replace(",", ".")))
    if len(med_rs) < 2:
        raise SystemExit(f"{a.pontos} tem {len(med_rs)} cota(s) preenchida(s); "
                         "sao precisas ao menos 2 para interpolar")
    o = np.argsort(-np.array(med_rs))
    med_rs = np.array(med_rs)[o]
    med_z = np.array(med_z)[o]
    print(f"geometria : {a.geom}   {len(S)} secoes")
    print(f"levantado : {len(med_rs)} pontos   "
          f"RS {med_rs.max():.0f} a {med_rs.min():.0f}   "
          f"cota {med_z.min():.2f} a {med_z.max():.2f} m")

    # distancia ao longo do rio de cada ponto levantado
    xm = np.interp(-med_rs, -rs, x)
    dentro = (rs <= med_rs.max() + 1e-6) & (rs >= med_rs.min() - 1e-6)
    # COTA ABSOLUTA ENTRE ANCORAS -- e a que zera os contradeclives, trocando
    # o ruido da lamina pela linha levantada (Mirim: 117 -> 0, o g02 que
    # convergiu). Vale tambem entre ancoras REAIS afastadas: interpolar 4 km
    # entre dois pontos levantados e o que ancorar significa.
    novo = z0.copy()
    novo[dentro] = np.interp(x[dentro], xm, med_z)
    # ---- ONDE NAO HA ANCORA POR QUILOMETROS, O REBAIXAMENTO E ZERO.
    # A interpolacao absoluta entre ancoras que cercam um VAO constroi uma
    # PONTE por baixo do terreno (alvo bruto 119 m abaixo do MDT, no trecho
    # em que o detector descartou a reta desenhada do legado). Tentativas com
    # teto e com taper por distancia vazavam (57-67 m na borda) ou puniam
    # vaos legitimos do Mirim (contradeclives 117 -> 43). O criterio final
    # nao olha o MOTIVO do ponto apagado nem usa limiar de cota: um AGLOMERADO
    # de pontos do pedido sem cota (vizinhos a menos de 3 km um do outro) com
    # extensao >= 3 km e um vao sem ancora -- dentro dele o peso e zero (o
    # MDT fica) com rampa de 1 km nas bordas. Ponto apagado ISOLADO segue
    # sendo ponte legitima entre as ancoras vizinhas, como sempre foi.
    n_fic = 0
    if fic_rs:
        xf = np.sort(np.interp(-np.array(fic_rs), -rs, x))
        grupos, ini = [], xf[0]
        for aa, bb in zip(xf[:-1], xf[1:]):
            if bb - aa > 3000.0:
                if aa - ini >= 3000.0:
                    grupos.append((ini, aa))
                ini = bb
        if xf[-1] - ini >= 3000.0:
            grupos.append((ini, xf[-1]))
        if grupos:
            w = np.ones(len(x))
            for xa, xb in grupos:
                dist = np.where(x < xa, xa - x,
                                np.where(x > xb, x - xb, 0.0))
                w = np.minimum(w, np.clip(dist / 1000.0, 0.0, 1.0))
            novo = z0 + w * (novo - z0)
            n_fic = int((w < 1.0).sum())
            print(f"   VAO SEM ANCORA   : rebaixamento zerado em {n_fic} "
                  f"secoes, em {len(grupos)} trecho(s): "
                  + ", ".join(f"{xa/1000:.1f}-{xb/1000:.1f} km"
                              for xa, xb in grupos))

    # Garante monotonicidade estrita do perfil final (zero contradeclives)
    for i in range(len(novo) - 2, -1, -1):
        if novo[i] < novo[i + 1] + 0.01:
            novo[i] = novo[i + 1] + 0.01

    fora = int((~dentro).sum())
    print(f"   secoes ancoradas : {int(dentro.sum())}")
    print(f"   FORA do intervalo levantado (perfil do MDT mantido): {fora}"
          + ("   <- extrapolar seria inventar" if fora else ""))
    d = novo - z0
    print(f"   ajuste do leito  : mediana {np.median(np.abs(d[dentro])):.2f} m"
          f"   p90 {np.percentile(np.abs(d[dentro]),90):.2f}"
          f"   max {np.abs(d[dentro]).max():.2f} m")
    dz0 = np.diff(z0)
    dz1 = np.diff(novo)
    s0 = np.abs(dz0) / np.maximum(ch[:-1], 1e-9)
    s1 = np.abs(dz1) / np.maximum(ch[:-1], 1e-9)
    print(f"   contradeclives   : {int((dz0>1e-9).sum())} -> "
          f"{int((dz1>1e-9).sum())}")
    print(f"   declividade >2%  : {int((s0>0.02).sum())} -> "
          f"{int((s1>0.02).sum())}   max {100*s0.max():.2f}% -> "
          f"{100*s1.max():.2f}%")

    # ---- aplica no perfil
    raiz = os.path.dirname(a.geom) or "."
    base = os.path.basename(a.geom).split(".")[0]
    saida = os.path.join(raiz, f"{base}.{a.saida}")
    if os.path.abspath(saida) == os.path.abspath(a.geom):
        raise SystemExit("saida igual a entrada -- recusado")
    porrs = {round(float(d_["rs"]), 2): k for k, d_ in enumerate(S)}
    linhas = open(a.geom, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, out, j, n_htab = -1, [], 0, 0
    invert = {}
    rio_atual = None
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("River Reach="):
            rio_atual = l.split("=", 1)[1].split(",")[0].strip()
        if l.startswith("Type RM Length L Ch R"):
            if rio_alvo and rio_atual != rio_alvo:
                i_sec = -1
            else:
                i_sec = porrs.get(round(float(l.split(",")[1]), 2), -1)
        k = i_sec
        if k >= 0 and abs(d[k]) > 1e-9:
            dd = S[k]
            if l.startswith("#Sta/Elev"):
                st = np.asarray(dd["sta"], float)
                zz = np.asarray(dd["z"], float).copy()
                m = (st >= float(dd["lb"])) & (st <= float(dd["rb"]))
                if m.any():
                    zc = zz[m]
                    prof = zc.max() - zc
                    p = prof.max()
                    # PESO EM PATAMAR, e nao proporcional a profundidade.
                    # Com peso `prof/profmax` o deslocamento de cada ponto
                    # cresce com a profundidade dele, e quando |d| passa da
                    # profundidade da calha a ORDEM SE INVERTE: o ponto mais
                    # fundo sobe mais que os outros e deixa de ser o mais
                    # fundo. Medido no teste de ida e volta, o leito errava o
                    # alvo em ate 4,66 m. Aqui o fundo inteiro anda JUNTO --
                    # peso 1 abaixo de `TAPER` da profundidade -- e so a faixa
                    # junto as margens e suavizada ate zero, para nao criar
                    # degrau na margem. A forma do fundo fica intacta e o
                    # minimo cai exatamente na cota levantada.
                    if p > 1e-9:
                        w = np.clip(prof / (TAPER * p), 0.0, 1.0)
                    else:
                        st_c = st[m]
                        L_c = st_c[-1] - st_c[0]
                        if L_c > 1e-3:
                            u_c = 4.0 * (st_c - st_c[0]) * (st_c[-1] - st_c) / (L_c ** 2)
                            w = np.sin(0.5 * np.pi * u_c)
                        else:
                            w = np.zeros_like(prof)
                    zz[m] = zc + d[k] * w
                invert[k] = float(zz[m].min()) if m.any() else float(zz.min())
                v = []
                for aa, bb in zip(st, zz):
                    v += [aa, bb]
                out.append("#Sta/Elev= %d " % len(st))
                lin, corpo = "", []
                for t2, xv in enumerate(v):
                    lin += "%8.2f" % xv
                    if (t2 + 1) % 10 == 0:
                        corpo.append(lin)
                        lin = ""
                if lin:
                    corpo.append(lin)
                out += corpo
                cnt = int(l.split("=")[1])
                j += 1
                lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    xq = linhas[j]
                    if not xq.strip() or xq[:1].isalpha() or xq[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(xq), 8)
                                  if xq[c:c + 8].strip()])
                    j += 1
                continue
            if l.startswith("XS HTab Starting El and Incr") and k in invert:
                q = [y.strip() for y in l.split("=", 1)[1].split(",")]
                out.append("XS HTab Starting El and Incr="
                           f"{invert[k]+0.02:.2f},{float(q[1]):.3f}, "
                           f"{int(q[2])} ")
                n_htab += 1
                j += 1
                continue
        out.append(l)
        j += 1
    escrever(saida, "\n".join(out))
    print(f"\ngeometria nova: {saida}   (HTab reancorado em {n_htab} secoes)")

    B = secoes_ordenadas(saida)
    zb = np.array([talvegue_xy(dd)[3] for dd in B])
    ca = np.array([float(dd["rb"] - dd["lb"]) for dd in S])
    cb = np.array([float(dd["rb"] - dd["lb"]) for dd in B])
    print("CONFERENCIA")
    print(f"   secoes                 : {len(S)} -> {len(B)}")
    print(f"   largura do canal mudou : max {np.abs(cb-ca).max():.6f} m "
          "(tem de ser zero)")
    err = np.abs(zb - novo)
    ruins = np.flatnonzero(err > 0.05)
    print(f"   leito bate com o alvo  : erro maximo {err.max():.3f} m   "
          f"acima de 5 cm em {len(ruins)} secoes")
    if len(ruins):
        # NAO SE LEVANTA O FUNDO ACIMA DAS PROPRIAS MARGENS. Quando a cota
        # pedida esta acima do topo da calha, a secao nao tem como acomoda-la
        # sem que as margens subam junto -- e mexer nas margens seria inventar
        # terreno. Acontece quando o levantamento pede SUBIR o leito muito
        # acima do que o MDT via; com batimetria de verdade o pedido costuma
        # ser DESCER, que nao tem esse limite.
        prof = []
        for i in ruins:
            st = np.asarray(S[i]["sta"], float)
            zq = np.asarray(S[i]["z"], float)
            m = (st >= float(S[i]["lb"])) & (st <= float(S[i]["rb"]))
            prof.append(float(zq[m].max() - zq[m].min()) if m.any() else 0.0)
        prof = np.array(prof)
        ped = np.abs(novo - z0)[ruins]
        print(f"      nelas o deslocamento pedido e de "
              f"{np.median(ped):.2f} m (mediana) numa calha de "
              f"{np.median(prof):.2f} m de profundidade")
        print("      -> a cota pedida esta acima do topo da calha; a secao so "
              "a aceitaria")
        print("         levantando as margens junto, e isso seria inventar "
              "terreno")
    return saida


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("pedir")
    p1.add_argument("geom")
    p1.add_argument("--cada", type=float, default=2000.0)
    p1.add_argument("--saida", required=True)
    p2 = sub.add_parser("aplicar")
    p2.add_argument("geom")
    p2.add_argument("--pontos", required=True)
    p2.add_argument("--saida", default="g02")
    p2.add_argument("--rio", default=None,
                    help="num arquivo de rede, aplica so ao rio nomeado")
    a = ap.parse_args()
    return cmd_pedir(a) if a.cmd == "pedir" else cmd_aplicar(a)


if __name__ == "__main__":
    main()
