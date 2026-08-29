# -*- coding: utf-8 -*-
"""Preenche o pedido de batimetria com o LEVANTAMENTO que ja existe no repo.

    python scripts/batimetria_do_legado.py doc/batimetria_itajai_mirim.csv \
        --rio Itajai_Mirim

A batimetria estava no repositorio o tempo todo. `legado/Itajai_Rede_1983.g01`
traz 1.240 secoes levantadas dos seis rios, com calha de verdade:

    Itajai_Mirim   258 secoes   calha mediana  8,73 m
    Itajai_Acu     285          10,91
    Itajai_Norte   366          10,43
    Rio_Benedito   144           8,76
    Itajai_Sul     110           7,51
    Itajai_Oeste    77           8,28

Contra 0,02 m de mediana na geometria tirada do MDT -- porque o MDT ve a
LAMINA e nao o fundo. Uma calha de 44 m por 2 cm nao conduz cheia nenhuma: o
Mirim batia nas 40 iteracoes em todo instante, terminava com "Solution Solver
Failed" e 92,38% de erro de volume.

ISTO NAO E INVENTAR COTA

  E o unico levantamento de fundo que existe para esta bacia, e ele ja estava
  aqui. O eixo tambem e o mesmo -- conferido: a distancia entre o eixo de
  `eixos_do_relevo.geojson` e o `Reach XY` do legado e ZERO em cinco dos seis
  rios, e os comprimentos batem ate a segunda casa. O `eixos_do_relevo` veio
  do legado, e nao do terreno.

  O que este script faz e so casar, por posicao, cada ponto do pedido com a
  secao levantada mais proxima, e escrever a cota de fundo DELA na coluna
  `z_leito_A_LEVANTAR`. Quem aplica e o `batimetria.py aplicar`, que interpola
  ao longo do rio e rebaixa a calha DENTRO DAS MARGENS, sem tocar na planicie
  -- que continua vindo do MDT medido.

  Fica gravado na coluna `observacao` de onde veio cada cota e a que distancia
  estava a secao levantada. Ponto que so ache secao a mais de `--limite`
  metros fica em branco: melhor faltar do que casar com o rio errado.

O QUE ELE NAO RESOLVE

  O levantamento e de 1983 e o leito move. Onde a distancia sai grande, ou
  onde o pedido cai fora do trecho levantado, a coluna fica vazia e o
  `aplicar` mantem o perfil do MDT ali -- extrapolar seria inventar.
"""
import argparse
import csv
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))

LEGADO = "legado/Itajai_Rede_1983.g01"


def secoes_levantadas(g, rio, completas=False):
    """(rs, x, y, invert, topo) de cada secao daquele rio no `.gNN` legado.

    Le o texto direto: `qc_secoes` nao devolve a que rio cada secao pertence,
    e num arquivo de rede isso e o que importa.

    Com `completas=True` devolve, em vez do array-resumo, uma lista de dicts
    com o PERFIL inteiro (rs, x, y, sta, z, lb, rb) -- e o que o
    `rio_do_relevo.py` usa para adotar a secao levantada onde o MDT e vazio
    (o estuario e lamina d'agua, e no SIG-SC lamina = 0.0 = nodata).
    """
    t = open(g, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    cur, out, i = None, [], 0
    while i < len(t):
        l = t[i]
        if l.startswith("River Reach="):
            cur = l.split("=", 1)[1].split(",")[0].strip()
        elif l.startswith("Type RM Length") and cur == rio:
            rs = float(l.split(",")[1])
            lb = rb = sta = z = cut = None
            j = i + 1
            while (j < len(t) and not t[j].startswith("Type RM Length")
                   and not t[j].startswith("River Reach=")):
                if t[j].startswith("Bank Sta="):
                    lb, rb = [float(x) for x in t[j].split("=", 1)[1].split(",")]
                elif t[j].startswith("XS GIS Cut Line="):
                    cut, j = _bloco(t, j, 16)
                    cut = cut.reshape(-1, 2)
                    continue
                elif t[j].startswith("#Sta/Elev="):
                    v, j = _bloco(t, j, 8)
                    sta, z = v[0::2], v[1::2]
                    continue
                j += 1
            if sta is not None and cut is not None and lb is not None:
                m = (sta >= lb) & (sta <= rb)
                inv = float(z[m].min()) if m.sum() >= 2 else float(z.min())
                topo = float(min(np.interp(lb, sta, z), np.interp(rb, sta, z)))
                c = cut.mean(0)
                if completas:
                    out.append({"rs": rs, "x": float(c[0]), "y": float(c[1]),
                                "sta": sta, "z": z, "lb": lb, "rb": rb,
                                "cut": (tuple(cut[0]), tuple(cut[-1]))})
                else:
                    out.append((rs, c[0], c[1], inv, topo))
            i = j - 1
        i += 1
    if completas:
        return out or None
    return np.array(out) if out else None


def _bloco(t, j, larg):
    """Le `N` valores em colunas de largura fixa a partir da linha `j`."""
    n = int(t[j].split("=")[1])
    v, k = [], j + 1
    while k < len(t) and len(v) < n * (2 if larg == 16 else 2):
        L = t[k]
        v += [float(L[c:c + larg]) for c in range(0, len(L), larg)
              if L[c:c + larg].strip()]
        k += 1
        if not L.strip():
            break
    return np.array(v[:2 * n]), k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pedido")
    ap.add_argument("--rio", required=True)
    ap.add_argument("--legado", default=LEGADO)
    ap.add_argument("--limite", type=float, default=1500.0,
                    help="m; alem disso nao casa e deixa em branco")
    ap.add_argument("--saida", default=None)
    a = ap.parse_args()

    L = secoes_levantadas(a.legado, a.rio)
    if L is None:
        raise SystemExit(f"'{a.rio}' nao esta em {a.legado}")
    prof = L[:, 4] - L[:, 3]
    print(f"levantamento: {a.legado}")
    print(f"   {a.rio}: {len(L)} secoes   RS {L[:,0].max():.0f} a "
          f"{L[:,0].min():.0f}")
    print(f"   fundo {L[:,3].min():.2f} a {L[:,3].max():.2f} m   "
          f"calha mediana {np.median(prof):.2f} m")

    # ---- DETECTOR DE FICCAO. Nas cabeceiras do Acu e do Benedito o "fundo
    # levantado" do legado e uma RETA DESENHADA: declive exatamente 8,00 m/km
    # com residuo rms de 1-2 MILIMETROS ao longo de 73-119 secoes -- nos dois
    # rios, a mesma constante. Rio de verdade nao faz isso (controles: Mirim
    # rms 15,4 m; Benedito de jusante 1,15 m). O titulo do proprio legado
    # avisa: "rede real ANA + relevo DEM". Ancorar nisso pediria cavar ate
    # 284 m ("eixo alto" foi o diagnostico errado; o eixo esta no vale que o
    # lidar mede -- procurado um caminho mais baixo num corredor de +-1500 m,
    # nao ha). Secao cujo entorno ajusta uma reta com rms < 5 cm e sintetica
    # e NAO vira ancora: o MDT fica valendo ali.
    o = np.argsort(-L[:, 0])
    rs_s, inv_s = L[o, 0], L[o, 3]
    sint_s = np.zeros(len(L), bool)
    W = 7
    for i in range(len(rs_s)):
        j0, j1 = max(0, i - W), min(len(rs_s), i + W + 1)
        if j1 - j0 < 5:
            continue
        c = np.polyfit(rs_s[j0:j1], inv_s[j0:j1], 1)
        # A reta PLANA e real: o canal dragado do estuario e reto POR PROJETO
        # (0,10 m/km, rms 0, fundo -10,85 -- profundidade de porto mantida).
        # A ficcao tem DECLIVE: 8,00 m/km desenhados morro acima. Sem a
        # condicao de declive o detector descartava o estuario junto, e o
        # aplicar do Acu piorava (ajuste max 67 m, declividades 45 -> 59).
        if abs(c[0]) >= 0.001 and \
                (inv_s[j0:j1] - np.polyval(c, rs_s[j0:j1])).std() < 0.05:
            sint_s[i] = True
    sintetica = np.zeros(len(L), bool)
    sintetica[o] = sint_s
    if sintetica.any():
        rr = rs_s[sint_s]
        print(f"   LEGADO SINTETICO: {int(sintetica.sum())} de {len(L)} "
              f"secoes sao reta desenhada (residuo local < 5 cm), "
              f"RS {rr.max():.0f} a {rr.min():.0f} -- nao viram ancora")

    linhas = list(csv.DictReader(open(a.pedido, encoding="utf-8"),
                                 delimiter=";"))
    if not linhas:
        raise SystemExit(f"{a.pedido} esta vazio")
    campos = list(linhas[0].keys())
    n_ok, n_fic, dists, quedas = 0, 0, [], []
    for r in linhas:
        x, y = float(r["x"]), float(r["y"])
        rs_p = float(r["rs"])
        # O CASAMENTO E POR RS, nao por distancia 2D. O centro da secao do
        # legado e esquematico; no cinturao de meandros ele fica longe do
        # ponto do pedido e o limite de 300 m em XY deixava trechos REAIS
        # orfaos de ancora (Oeste: 27 em branco, e o aplicar mantinha o ruido
        k = int(np.argmin(np.abs(L[:, 0] - rs_p)))
        d_rs = abs(float(L[k, 0]) - rs_p)
        d = float(np.hypot(L[k, 1] - x, L[k, 2] - y))
        if d_rs > 800.0 or d > a.limite:
            r["observacao"] = (f"sem secao levantada (RS mais proximo a "
                               f"{d_rs:.0f} m, centro a {d:.0f} m)")
            continue
        if sintetica[k]:
            r["z_leito_A_LEVANTAR"] = ""
            r["observacao"] = (f"legado sintetico em RS {L[k,0]:.0f} (reta "
                               "desenhada, residuo < 5 cm) -- MDT vale aqui")
            n_fic += 1
            continue
        # A ZONA DE TRANSICAO da reta desenhada foi "suavizada" (rms 8 m --
        # escapa do teste de linearidade) mas ainda pede cortar 55-71 m numa
        # calha levantada de ~11 m: implausivel do mesmo jeito. Este filtro
        # por ponto ja existiu e FALHOU quando o aplicar interpolava por cima
        # do buraco; agora o aplicar zera POR INTERVALO os pontos marcados, e
        # o filtro volta a ser seguro.
        zm_ = float(r.get("z_lamina_mdt") or "nan")
        if np.isfinite(zm_) and zm_ - L[k, 3] > 25.0:
            r["z_leito_A_LEVANTAR"] = ""
            r["observacao"] = (f"rebaixamento implausivel ({zm_-L[k,3]:.0f} m "
                               "> 25 numa calha de ~11) -- MDT vale aqui")
            n_fic += 1
            continue
        r["z_leito_A_LEVANTAR"] = f"{L[k,3]:.2f}"
        zm = float(r.get("z_lamina_mdt") or "nan")
        quedas.append(zm - L[k, 3])
        dists.append(d)
        n_ok += 1
        r["observacao"] = (f"1983 RS {L[k,0]:.2f} a {d:.0f} m; "
                           f"calha levantada {L[k,4]-L[k,3]:.2f} m")

    saida = a.saida or a.pedido
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, campos, delimiter=";")
        w.writeheader()
        w.writerows(linhas)
    print(f"\npedido    : {a.pedido}   {len(linhas)} pontos")
    print(f"   casados : {n_ok}   em branco: {len(linhas)-n_ok}"
          + (f"   (destes, {n_fic} por legado sintetico)" if n_fic else ""))
    if dists:
        dists = np.array(dists)
        quedas = np.array(quedas)
        print(f"   distancia ate a secao levantada: mediana "
              f"{np.median(dists):.0f} m   max {dists.max():.0f} m")
        print(f"   REBAIXAMENTO do leito (lamina do MDT - fundo levantado): "
              f"mediana {np.median(quedas):.2f} m   "
              f"p10 {np.percentile(quedas,10):.2f}   "
              f"p90 {np.percentile(quedas,90):.2f}")
        if (quedas < 0).any():
            print(f"   ATENCAO: {int((quedas<0).sum())} ponto(s) com fundo "
                  "levantado ACIMA da lamina do MDT -- conferir antes de "
                  "aplicar")
    print(f"\npreenchido -> {saida}")
    print("   aplicar com:")
    print(f"   python scripts/batimetria.py aplicar <geom.g01> "
          f"--pontos {saida} --saida g02")
    return saida


if __name__ == "__main__":
    main()
