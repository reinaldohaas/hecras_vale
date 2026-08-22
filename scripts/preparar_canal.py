# -*- coding: utf-8 -*-
"""Monta as secoes do Canal Retificado -- tudo menos o fundo.

    python scripts/preparar_canal.py --dx 150 --largura 45 --meia 500

O QUE SAI, E O QUE FICA EM ABERTO

  Sai a geometria completa do canal EXCETO as cotas dentro do canal:

      eixo, estacas, linhas de corte, orientacao        prontos
      planicie e margens, cortadas do MDT SIG-SC 1 m    prontas
      posicao das margens (Bank Sta)                    pronta
      COTA DE LEITO entre as margens                    EM ABERTO

  O leito fica vazio de proposito. O MDT enxerga a LAMINA D'AGUA, nao o
  fundo: ele da 0,06 m no canal enquanto a batimetria do modelo nas duas
  pontas da -0,76 m e -2,68 m -- erro de pelo menos 2,74 m. Preencher com o
  MDT produziria um canal raso e falso, com cara de dado medido.

  A cota do MDT dentro do canal NAO e jogada fora: vai na coluna
  `z_mdt_lamina`, identificada como lamina. Serve de teto para o leito
  (o fundo esta abaixo dela), nunca de leito.

DE ONDE VEM CADA COISA

  eixo      `canal_itajai_mirim.geojson` (OpenStreetMap, `waterway=canal`,
            vias 290409480 e 290766868), WGS84 -> EPSG:31982
  planicie  MDT SIG-SC 1 m, folhas lidas direto
  largura   45 m, informada pelo usuario. Confere com o MDT: a largura da
            lamina sai em 52 m de mediana e 44 m no p10 -- o MDT inclui um
            pedaco da margem acima d'agua, entao 45 m e coerente.

DUAS DECISOES QUE NAO SAO MINHAS, E POR QUE

  EXTENSAO LATERAL. A planicie e plana em torno de 2,7 a 3,4 m ate 500 m dos
  dois lados, e a lamina maxima do modelo nesse trecho e 6,05 a 6,13 m. Para
  CONTER essa cheia a secao precisaria de cerca de 3,3 km de largura, e em 21
  de 31 perfis o terreno nem alcanca 6,13 m dentro de 2,5 km. Secao 1D de
  3 km atravessando varzea plana nao representa escoamento -- ali a agua nao
  corre na direcao da secao. Isto pede area de armazenamento ou 2D, e nao
  secao mais larga. O `--meia` fica em 500 m por padrao, explicitamente como
  ESCOLHA, e nao como resultado de medida.

  ESPACAMENTO. `--dx` 150 m, cerca de 3,3 larguras de canal. O canal e reto e
  prismatico, e com dt de 15 s e velocidade de ordem 1 m/s o Courant fica em
  0,1. Nao ha exigencia numerica de adensar.

O ESTACIONAMENTO E PROVISORIO

  RS vai de `dx*(n-1)` ate 0, medido do extremo de JUSANTE. Ao juntar isto ao
  modelo, o canal vira um reach entre duas juncoes e o RS tera de ser
  renumerado junto com a divisao do reach existente -- o que ainda nao foi
  feito e nao esta neste arquivo.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402
from qc_secoes import ler_secoes                        # noqa: E402
from qc_geometria import ler_eixos                      # noqa: E402

GEOJSON = r"C:\Users\haas\Downloads\canal_itajai_mirim.geojson"
MODELO = "modelo/mirim_t30/mirim_t30.g01"
SAIDA = "doc/canal"
PASSO_MDT = 2.0        # m entre pontos amostrados na cutline


def _arg(argv, chave, padrao, tipo=float):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def eixo_do_canal(p=GEOJSON):
    from pyproj import Transformer
    from shapely.geometry import LineString
    from shapely.ops import linemerge
    tr = Transformer.from_crs("EPSG:4326", "EPSG:31982", always_xy=True)
    d = json.load(open(p, encoding="utf-8"))
    segs = []
    for f in d["features"]:
        a = np.asarray(f["geometry"]["coordinates"], float)
        x, y = tr.transform(a[:, 0], a[:, 1])
        segs.append(LineString(np.c_[x, y]))
    c = linemerge(segs)
    return max(c.geoms, key=lambda g: g.length) if hasattr(c, "geoms") else c


def main(argv):
    from shapely.geometry import LineString, Point, mapping
    dx = _arg(argv, "--dx", 150.0)
    larg = _arg(argv, "--largura", 45.0)
    meia = _arg(argv, "--meia", 500.0)
    os.makedirs(SAIDA, exist_ok=True)

    canal = eixo_do_canal()
    # ---- sentido: de montante para jusante, conferido contra o modelo
    eixo_m = list(ler_eixos(MODELO).values())[0]
    S = ler_secoes(MODELO)
    S.sort(key=lambda d: -d["rs"])
    cen = np.array([0.5 * (np.asarray(d["cut"][0], float)
                           + np.asarray(d["cut"][-1], float)) for d in S])
    rs = np.array([d["rs"] for d in S])

    def rs_de(P):
        k = int(np.argmin(np.hypot(cen[:, 0] - P[0], cen[:, 1] - P[1])))
        return rs[k]
    A, B = np.array(canal.coords[0]), np.array(canal.coords[-1])
    if rs_de(A) < rs_de(B):                 # RS cresce para MONTANTE
        canal = LineString(list(canal.coords)[::-1])
        A, B = B, A
    print(f"eixo do canal: {canal.length/1000:.3f} km   "
          f"de RS ~{rs_de(A):.0f} (montante) a RS ~{rs_de(B):.0f} (jusante)")
    print(f"canal de {larg:g} m   secoes a cada {dx:g} m   "
          f"meia-largura {meia:g} m")

    s = np.arange(0.0, canal.length, dx)
    n = len(s)
    print(f"secoes: {n}")

    # ---- amostragem do MDT
    off = np.arange(-meia, meia + PASSO_MDT / 2, PASSO_MDT)
    pts, eixos_n, centros = [], [], []
    for si in s:
        p0 = np.array(canal.interpolate(max(si - 25.0, 0.0)).coords[0])
        p1 = np.array(canal.interpolate(
            min(si + 25.0, canal.length)).coords[0])
        t = p1 - p0
        t /= max(float(np.hypot(*t)), 1e-9)
        nv = np.array([-t[1], t[0]])       # esquerda olhando para jusante
        c = np.array(canal.interpolate(si).coords[0])
        centros.append(c)
        eixos_n.append(nv)
        for o in off:
            pts.append(c - o * nv)         # -meia = margem esquerda
    pts = np.array(pts)
    bb = (pts[:, 0].min() - 100, pts[:, 1].min() - 100,
          pts[:, 0].max() + 100, pts[:, 1].max() + 100)
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(bb))
    Z = mdt.cota(pts[:, 0], pts[:, 1]).reshape(n, len(off))
    print(f"MDT: cobertura {100*np.isfinite(Z).mean():.1f}% "
          f"em {n*len(off)} pontos")

    # ---- estaca 0 na ponta esquerda; canal centrado no eixo
    sta = off + meia
    lb, rb = meia - larg / 2.0, meia + larg / 2.0
    dentro = (sta > lb) & (sta < rb)

    # ---- tabelas
    import csv
    p1 = os.path.join(SAIDA, "canal_secoes.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["rs", "estaca", "x", "y", "z", "origem", "z_mdt_lamina"])
        for i in range(n):
            rsi = (n - 1 - i) * dx
            for k in range(len(off)):
                P = pts[i * len(off) + k]
                z = Z[i, k]
                if dentro[k]:
                    w.writerow([f"{rsi:.1f}", f"{sta[k]:.2f}",
                                f"{P[0]:.2f}", f"{P[1]:.2f}", "",
                                "A LEVANTAR",
                                "" if not np.isfinite(z) else f"{z:.2f}"])
                else:
                    w.writerow([f"{rsi:.1f}", f"{sta[k]:.2f}",
                                f"{P[0]:.2f}", f"{P[1]:.2f}",
                                "" if not np.isfinite(z) else f"{z:.2f}",
                                "MDT SIG-SC 1 m", ""])
    print(f"\nperfis            -> {p1}")

    # ---- o pedido de levantamento: eixo e margens de cada secao
    p2 = os.path.join(SAIDA, "canal_batimetria_a_levantar.csv")
    with open(p2, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["rs", "ponto", "estaca", "x", "y",
                    "z_lamina_mdt", "z_leito_A_LEVANTAR"])
        for i in range(n):
            rsi = (n - 1 - i) * dx
            c = centros[i]
            nv = eixos_n[i]
            for rot, e in (("margem_esq", lb), ("eixo", meia),
                           ("margem_dir", rb)):
                P = c - (e - meia) * nv
                k = int(np.argmin(np.abs(sta - e)))
                z = Z[i, k]
                w.writerow([f"{rsi:.1f}", rot, f"{e:.2f}",
                            f"{P[0]:.2f}", f"{P[1]:.2f}",
                            "" if not np.isfinite(z) else f"{z:.2f}", ""])
    print(f"pedido de batimetria -> {p2}   ({3*n} pontos)")

    # ---- cutlines e eixo, para abrir no RAS Mapper / QGIS
    feats = []
    for i in range(n):
        c = centros[i]
        nv = eixos_n[i]
        P0 = c + meia * nv
        P1 = c - meia * nv
        feats.append({"type": "Feature",
                      "properties": {"rs": (n - 1 - i) * dx,
                                     "largura_canal": larg},
                      "geometry": mapping(LineString([P0, P1]))})
    p3 = os.path.join(SAIDA, "canal_cutlines.geojson")
    json.dump({"type": "FeatureCollection",
               "crs": {"type": "name",
                       "properties": {"name": "urn:ogc:def:crs:EPSG::31982"}},
               "features": feats}, open(p3, "w"), indent=1)
    p4 = os.path.join(SAIDA, "canal_eixo.geojson")
    json.dump({"type": "FeatureCollection",
               "crs": {"type": "name",
                       "properties": {"name": "urn:ogc:def:crs:EPSG::31982"}},
               "features": [{"type": "Feature",
                             "properties": {"nome": "Canal Retificado",
                                            "comprimento_m": canal.length},
                             "geometry": mapping(canal)}]},
              open(p4, "w"), indent=1)
    print(f"cutlines          -> {p3}")
    print(f"eixo              -> {p4}")

    # ---- resumo
    fora = Z[:, ~dentro]
    lam = Z[:, dentro]
    print("\nRESUMO")
    print(f"   secoes                 : {n}   RS de {(n-1)*dx:.0f} a 0")
    print(f"   pontos por secao       : {len(off)}   "
          f"({int(dentro.sum())} deles dentro do canal, sem cota)")
    print(f"   planicie (do MDT)      : mediana {np.nanmedian(fora):.2f} m   "
          f"p10 {np.nanpercentile(fora,10):.2f}   "
          f"p90 {np.nanpercentile(fora,90):.2f}")
    print(f"   lamina no canal (MDT)  : mediana {np.nanmedian(lam):.2f} m   "
          "<- TETO do leito, nao o leito")
    ka = int(np.argmin(np.abs(rs - rs_de(A))))
    kb = int(np.argmin(np.abs(rs - rs_de(B))))
    za = float(S[ka]["z"].min())
    zb = float(S[kb]["z"].min())
    print(f"   leito medido nas pontas: {za:+.2f} m (montante)  "
          f"{zb:+.2f} m (jusante)   queda {za-zb:.2f} m")
    print(f"   -> o levantamento precisa cobrir {canal.length/1000:.2f} km "
          f"entre essas duas cotas")
    return SAIDA


if __name__ == "__main__":
    main(sys.argv[1:])
