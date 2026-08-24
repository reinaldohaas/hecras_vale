# -*- coding: utf-8 -*-
"""Geometria de um rio a partir do RELEVO, sem nada esculpido.

    python scripts/rio_do_relevo.py --rio Rio_Benedito --saida modelo/benedito

NAO IMPORTA NADA DA CADEIA DE CORRECOES. Usa so `mdt_sigsc` (acesso ao MDT) e
`ras_io` (gravar em CRLF, que o HEC-RAS exige). Tudo o mais e calculado aqui, a
partir do terreno.

O QUE ESTE GERADOR NAO FAZ, e por que ele existe

  `gerar_mirim_do_zero.py` produz um modelo que converge em 2 iteracoes e tem 2
  erros de geometria -- mas o talvegue dele sao OITO NUMEROS escritos no
  codigo, interpolados por PCHIP, e a calha e uma parabola:

      z[calha] = z_alvo + (cota_margem - z_alvo) * (dist_norm ** 2)
      z_lob    = max(z_lob, z_alvo + 2.5)
      z[esq]   = np.maximum(z[esq], z_lob)
      z[0]     = max(z[0], 4.50)

  Medido contra o MDT em 64 secoes: dentro da calha a mediana e -0,50 m mas o
  p90 chega a +39 m; na planicie a mediana e zero e o p90 e +25 m, que e o
  `np.maximum` levantando o terreno onde ele desce. Converge porque e liso por
  construcao, e nao porque descreve o rio.

  Aqui NADA e esculpido: o perfil e o que o MDT da, ponto a ponto.

O QUE E MEDIDO, E COMO

  talvegue    o ponto mais baixo do MDT perto do eixo. RESSALVA QUE NAO SE
              resolve: o MDT ve a LAMINA D'AGUA, nao o fundo. Onde ha agua, o
              "talvegue" e a superficie livre. Sem batimetria nao ha como
              saber o fundo, e inventa-lo e o que este gerador recusa fazer.

  margens     andando para fora do talvegue, o primeiro ponto de cada lado que
              sobe `FOLGA_CALHA` acima dele. E o topo do encaixe, medido.

  meia-largura  continua para fora ate subir `ALVO_SECAO` acima do talvegue,
              com teto em `MEIA_MAX`. Onde a varzea e plana o teto manda, e o
              relatorio diz em quantas secoes isso aconteceu -- ali a secao
              1D nao contem a cheia, e isso pede armazenamento ou 2D.

  Manning     `N_CALHA` e `N_PLANICIE`, constantes e declaradas. Nao ha dado de
              rugosidade nesta bacia; fingir que ha seria o mesmo erro.

O PERFIL LONGITUDINAL SAI CRU

  Sem suavizacao, sem declividade minima forcada, sem monotonicidade imposta.
  O relatorio mede quantos contradeclives e quantos degraus o terreno traz, e
  a decisao de tratar isso e de quem le -- com `--monotono` disponivel, que
  aplica regressao isotonica: ela ajusta uma curva nao-crescente aos valores
  MEDIDOS, sem sair da faixa deles.
"""
import argparse
import json
import os
import sys

import numpy as np
from shapely.geometry import LineString

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402
from ras_io import escrever                            # noqa: E402

EIXOS = "eixos_do_relevo.geojson"
DX = 150.0            # m entre secoes
PASSO = 4.0           # m entre pontos amostrados na cutline
JANELA = 60.0         # m para cada lado, ao medir a tangente do eixo
MEIA_MAX = 400.0      # m; teto da meia-largura
MEIA_MIN = 60.0       # m; piso, para a secao nunca degenerar
FOLGA_CALHA = 1.5     # m acima do talvegue = topo da margem
ALVO_SECAO = 8.0      # m acima do talvegue = onde a secao pode parar
BUSCA = 500.0         # m; ate onde se procura terreno alto
N_CALHA, N_PLANICIE = 0.032, 0.055
WKT = ('PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",'
       'DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],'
       'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
       'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
       'PARAMETER["False_Northing",10000000.0],'
       'PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],'
       'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')


def eixo_do_rio(nome, caminho=EIXOS):
    d = json.load(open(caminho, encoding="utf-8"))
    for f in d["features"]:
        if f["properties"].get("nome") == nome:
            return LineString(np.asarray(f["geometry"]["coordinates"], float))
    raise SystemExit(f"'{nome}' nao esta em {caminho}. Ha: "
                     f"{[f['properties'].get('nome') for f in d['features']]}")


def isotonica(z):
    """Maior curva NAO-CRESCENTE que melhor ajusta z (pool adjacent violators).

    Nao inventa valor fora da faixa medida: cada patamar da saida e a MEDIA de
    um bloco de valores de entrada.
    """
    v = [float(z[0])]
    w = [1.0]
    for x in z[1:]:
        v.append(float(x))
        w.append(1.0)
        while len(v) > 1 and v[-2] < v[-1]:      # violou o nao-crescente
            x2 = v.pop()
            w2 = w.pop()
            x1 = v.pop()
            w1 = w.pop()
            v.append((x1 * w1 + x2 * w2) / (w1 + w2))
            w.append(w1 + w2)
    saida = []
    for x, k in zip(v, w):
        saida += [x] * int(k)
    return np.array(saida)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rio", required=True)
    ap.add_argument("--saida", required=True)
    ap.add_argument("--reach", default="R1")
    ap.add_argument("--dx", type=float, default=DX)
    ap.add_argument("--monotono", action="store_true",
                    help="ajusta o talvegue por regressao isotonica")
    a = ap.parse_args()

    eixo = eixo_do_rio(a.rio)
    L = eixo.length
    est = np.arange(0.0, L, a.dx)
    if L - est[-1] > 20.0:
        est = np.append(est, L)
    print(f"rio    : {a.rio}   eixo {L/1000:.2f} km   {len(est)} secoes "
          f"a cada {a.dx:g} m")

    # ---- geometria das cutlines e amostragem do MDT
    off = np.arange(-BUSCA, BUSCA + PASSO / 2, PASSO)
    base, normais, pts = [], [], []
    for s in est:
        p = np.array(eixo.interpolate(s).coords[0])
        q0 = np.array(eixo.interpolate(max(s - JANELA, 0.0)).coords[0])
        q1 = np.array(eixo.interpolate(min(s + JANELA, L)).coords[0])
        t = q1 - q0
        nt = float(np.hypot(*t))
        if nt < 1e-6:
            continue
        t /= nt
        n = np.array([-t[1], t[0]])
        base.append((float(s), p))
        normais.append(n)
        for o in off:
            pts.append(p + o * n)
    pts = np.array(pts)
    bb = (pts[:, 0].min() - 60, pts[:, 1].min() - 60,
          pts[:, 0].max() + 60, pts[:, 1].max() + 60)
    tiles = tiles_do_dominio(bb)
    print(f"MDT    : {len(pts)} pontos a {PASSO:g} m sobre {len(tiles)} folhas")
    Z = MosaicoSigsc(tiles=tiles).cota(pts[:, 0], pts[:, 1]) \
        .reshape(len(base), len(off))

    # ---- cada secao, medida
    secoes, no_teto, sem_dado = [], 0, 0
    for k, ((s, p), n) in enumerate(zip(base, normais)):
        z = Z[k]
        if not np.isfinite(z).any():
            sem_dado += 1
            continue
        centro = np.abs(off) <= 30.0
        if not np.isfinite(z[centro]).any():
            sem_dado += 1
            continue
        i0 = int(np.nanargmin(np.where(centro, z, np.nan)))
        zt = float(z[i0])

        def anda(sinal, alvo):
            i = i0
            while 0 < i < len(off) - 1:
                i += sinal
                if not np.isfinite(z[i]):
                    continue
                if z[i] >= zt + alvo:
                    return i
            return None

        ie = anda(-1, FOLGA_CALHA)
        idr = anda(+1, FOLGA_CALHA)
        se = anda(-1, ALVO_SECAO)
        sd = anda(+1, ALVO_SECAO)
        lim_e = abs(off[se]) if se is not None else MEIA_MAX
        lim_d = abs(off[sd]) if sd is not None else MEIA_MAX
        if se is None or sd is None:
            no_teto += 1
        me = float(np.clip(lim_e, MEIA_MIN, MEIA_MAX))
        md = float(np.clip(lim_d, MEIA_MIN, MEIA_MAX))
        dentro = (off >= -me) & (off <= md) & np.isfinite(z)
        if dentro.sum() < 8:
            sem_dado += 1
            continue
        sta = off[dentro] + me
        zz = z[dentro]
        lb = (off[ie] + me) if ie is not None else float(sta[0])
        rb = (off[idr] + me) if idr is not None else float(sta[-1])
        lb = float(sta[np.argmin(np.abs(sta - lb))])
        rb = float(sta[np.argmin(np.abs(sta - rb))])
        if rb <= lb:
            j = int(np.argmin(np.abs(sta - (off[i0] + me))))
            lb = float(sta[max(j - 1, 0)])
            rb = float(sta[min(j + 1, len(sta) - 1)])
        A = p - me * n
        B = p + md * n
        secoes.append({"s": s, "rs": round(float(L - s), 2),
                       "cut": (A, B), "sta": np.round(sta, 2),
                       "z": np.round(zz, 2), "lb": lb, "rb": rb,
                       "zt": float(zz.min())})

    print(f"secoes : {len(secoes)}   sem MDT utilizavel: {sem_dado}   "
          f"pararam no teto de {MEIA_MAX:g} m: {no_teto}")

    # ---- talvegue: cru, ou isotonico
    zt = np.array([s["zt"] for s in secoes])
    cru = zt.copy()
    if a.monotono:
        novo = isotonica(zt)
        for s, z0, z1 in zip(secoes, cru, novo):
            d = z1 - z0
            st = s["sta"]
            m = (st >= s["lb"]) & (st <= s["rb"])
            if m.any():
                prof = s["z"][m].max() - s["z"][m]
                pmax = prof.max()
                peso = prof / pmax if pmax > 1e-9 else np.zeros_like(prof)
                s["z"][m] = np.round(s["z"][m] + d * peso, 2)
            s["zt"] = float(s["z"].min())
        zt = np.array([s["zt"] for s in secoes])
        print(f"talvegue: regressao isotonica aplicada   "
              f"ajuste mediano {np.median(np.abs(novo-cru)):.2f} m   "
              f"max {np.abs(novo-cru).max():.2f} m")

    # ---- escreve
    os.makedirs(a.saida, exist_ok=True)
    nome = os.path.basename(a.saida.rstrip("/\\"))
    g = os.path.join(a.saida, f"{nome}.g01")
    l = [f"Geom Title={nome}", "Program Version=7.01"]
    P = np.vstack([np.vstack(s["cut"]) for s in secoes])
    l.append("Viewing Rectangle= %.2f , %.2f , %.2f , %.2f "
             % (P[:, 0].min(), P[:, 0].max(), P[:, 1].max(), P[:, 1].min()))
    l.append("Spatial Reference System=" + WKT)
    l.append("")
    l.append(f"River Reach={a.rio:<16.16},{a.reach:<16.16}")
    c = list(eixo.coords)
    l.append(f"Reach XY= {len(c)} ")
    ss = [f"{x:16.4f}{y:16.4f}" for x, y in c]
    for k in range(0, len(ss), 2):
        l.append("".join(ss[k:k + 2]))
    l.append("Rch Text X Y=0,0,0,0")
    l.append("")
    for i, s in enumerate(secoes):
        d = (round(float(secoes[i + 1]["s"] - s["s"]), 2)
             if i + 1 < len(secoes) else 0.0)
        l.append(f"Type RM Length L Ch R = 1 ,{s['rs']:.2f},"
                 f"{d:8.2f},{d:8.2f},{d:8.2f}")
        l.append(f"Bank Sta={s['lb']:.2f},{s['rb']:.2f}")
        l.append("XS GIS Cut Line= 2")
        l.append("".join(f"{q[0]:16.2f}{q[1]:16.2f}" for q in s["cut"]))
        l.append(f"#Sta/Elev= {len(s['sta'])} ")
        pf = [f"{x:8.2f}{y:8.2f}" for x, y in zip(s["sta"], s["z"])]
        for k in range(0, len(pf), 5):
            l.append("".join(pf[k:k + 5]))
        l.append("#Mann= 3 , 0 , 0 ")
        l.append(f"{s['sta'][0]:8.2f}{N_PLANICIE:8.3f}{0:8d}"
                 f"{s['lb']:8.2f}{N_CALHA:8.3f}{0:8d}"
                 f"{s['rb']:8.2f}{N_PLANICIE:8.3f}{0:8d}")
        l.append(f"XS HTab Starting El and Incr={s['zt']+0.02:.2f},0.100, 500 ")
        l.append("XS HTab Horizontal Distribution=-1,-1,-1")
        l.append("XS Rating Curve= 0 ,0")
        l.append("Exp/Cntr=0.3,0.1")
        l.append("")
    escrever(g, "\n".join(l))

    # ---- o que o terreno entregou, sem maquiagem
    lc = np.array([s["rb"] - s["lb"] for s in secoes])
    ls = np.array([s["sta"][-1] for s in secoes])
    npt = np.array([len(s["sta"]) for s in secoes])
    dz = np.diff(zt)
    ch = np.array([secoes[i + 1]["s"] - secoes[i]["s"]
                   for i in range(len(secoes) - 1)])
    decl = np.abs(dz) / np.maximum(ch, 1e-9)
    print(f"\ngeometria: {g}")
    print(f"   talvegue     : {zt.min():.2f} a {zt.max():.2f} m")
    print(f"   sobem p/ jusante: {int((dz > 1e-9).sum())} de {len(dz)}"
          f"   (o terreno traz isto; --monotono trata)")
    print(f"   pares de cota igual: {int((np.abs(dz) < 0.005).sum())}")
    print(f"   declividade  : mediana {np.median(decl):.5f}   "
          f">2% em {int((decl > 0.02).sum())} trechos   max {decl.max():.3f}")
    print(f"   calha        : mediana {np.median(lc):.0f} m   "
          f"p10 {np.percentile(lc,10):.0f}   p90 {np.percentile(lc,90):.0f}")
    print(f"   secao        : mediana {np.median(ls):.0f} m   "
          f"max {ls.max():.0f} m")
    print(f"   pontos/secao : mediana {np.median(npt):.0f}   max {npt.max()}"
          "   (limite do HEC-RAS: 500)")
    return g


if __name__ == "__main__":
    main()
