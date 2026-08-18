# -*- coding: utf-8 -*-
"""
Confere o que esta GRAVADO no .g01 contra o terreno, na mesma linha.

Buraco que este modulo fecha. O resto do hecras_qc audita a CUTLINE contra o
DEM: se ela esta no lugar certo, se cruza o canal, se e perpendicular. Mas as
cotas que vao para o HEC-RAS nao sao as do DEM -- passaram por calha imposta,
perfil condicionado, pilot channel e parede vertical. Nenhuma dessas etapas era
verificada.

O custo disso foi concreto: um pilot channel com o limite mal escrito rebaixou
as 1.232 secoes do modelo a uma cota unica, e a auditoria continuou dizendo
"0 saltos de area, 0 secoes rasas" -- uma bacia chata com parede vertical passa
nos dois testes. Quem viu foi o usuario, olhando UMA figura.

O que se mede aqui, secao a secao:

  ESCAVACAO      quanto o modelo esta abaixo do terreno, e onde. Escavar a
                 calha e legitimo; escavar a planicie inteira nao.
  ATERRO         quanto o modelo esta acima do terreno. Fora das paredes das
                 pontas, isso e terreno inventado.
  ACHATAMENTO    fracao dos pontos numa mesma cota. Terreno real nao tem 93%
                 da secao na mesma cota; calha imposta larga demais tem.
  TERRENO PRESERVADO fracao dos pontos em que o modelo AINDA e o terreno
                 (diferenca abaixo de meio metro). Diz se a forma do vale
                 sobreviveu ao condicionamento. Comparar DESNIVEL nao serve:
                 escavar rebaixa o minimo e o desnivel do modelo fica MAIOR
                 que o do terreno, entao a metrica dava 1,10 de mediana e nao
                 media o que promete.

Uso:
    python -m hecras_qc.modelo modelo.g01 terreno.tif [EPSG:31982]
"""
import sys

import numpy as np

from .dem import DEM
from .ras_geometry import ler_g01

ESCAVACAO_ALERTA = 8.0     # m; abaixo do terreno, fora da calha
ATERRO_ALERTA = 2.0        # m; acima do terreno, fora das pontas
ACHATAMENTO_ALERTA = 0.40  # fracao da secao numa cota so
PRESERVADO_MIN = 0.50      # fracao da secao que deve continuar sendo terreno
PONTAS = 2                 # pontos de cada ponta que podem ser parede


def comparar(secao_g01, dem):
    """Uma secao do .g01 contra o terreno na mesma linha."""
    sta = np.asarray(secao_g01["sta"], float)
    zm = np.asarray(secao_g01["z"], float)
    g = secao_g01["geometry"]
    if len(sta) < 5 or g is None:
        return None
    # amostra o DEM nas MESMAS estacas da tabela do modelo, para a comparacao
    # ser ponto a ponto e nao entre malhas diferentes
    L = float(g.length)
    if L <= 0:
        return None
    frac = np.clip(sta / (sta[-1] or 1.0), 0.0, 1.0)
    pts = [g.interpolate(float(f) * L) for f in frac]
    zt = dem.cota([p.x for p in pts], [p.y for p in pts])

    ok = np.isfinite(zt) & np.isfinite(zm)
    if ok.sum() < 5:
        return None
    d = zm - zt                                   # positivo = modelo acima
    miolo = np.ones(len(d), bool)
    miolo[:PONTAS] = miolo[-PONTAS:] = False      # as pontas podem ser parede
    m = ok & miolo

    vals, cnt = np.unique(np.round(zm, 2), return_counts=True)
    achat = float(cnt.max() / len(zm))
    rel_t = float(np.nanmax(zt[ok]) - np.nanmin(zt[ok]))
    preservado = float(np.mean(np.abs(d[m]) < 0.5)) if m.any() else 0.0

    return {
        "river": secao_g01["river"], "reach": secao_g01["reach"],
        "rs": secao_g01["rs"],
        "escavacao_max": float(-np.nanmin(d[m])) if m.any() else 0.0,
        "escavacao_media": float(-np.nanmean(np.clip(d[m], None, 0.0)))
                           if m.any() else 0.0,
        "aterro_max": float(np.nanmax(d[m])) if m.any() else 0.0,
        "achatamento": achat,
        "relevo_terreno": rel_t,
        "preservado": preservado,
        "cotas_distintas": int(len(vals)),
        "largura": float(sta[-1] - sta[0]),
    }


def avaliar(r):
    """Motivos pelos quais esta secao do modelo nao representa o terreno."""
    m = []
    if r["achatamento"] > ACHATAMENTO_ALERTA:
        m.append(f"{100*r['achatamento']:.0f}% da secao numa cota so")
    if r["preservado"] < PRESERVADO_MIN and r["relevo_terreno"] > 5.0:
        m.append(f"so {100*r['preservado']:.0f}% da secao ainda e terreno")
    if r["escavacao_max"] > ESCAVACAO_ALERTA:
        m.append(f"escavado {r['escavacao_max']:.1f} m abaixo do terreno")
    if r["aterro_max"] > ATERRO_ALERTA:
        m.append(f"aterrado {r['aterro_max']:.1f} m acima do terreno")
    return m


def rodar(g01, tif, crs="EPSG:31982"):
    import geopandas as gpd
    dem = DEM(tif)
    _, secoes = ler_g01(g01, com_perfil=True)
    if not secoes:
        print("nenhuma secao lida do .g01")
        return 1
    gdf = gpd.GeoDataFrame(secoes, crs=crs)
    if dem.crs is not None:
        gdf = gdf.to_crs(dem.crs_metrico)
    for s, geom in zip(secoes, gdf.geometry):
        s["geometry"] = geom

    linhas = [comparar(s, dem) for s in secoes]
    linhas = [r for r in linhas if r]
    print(f"{g01}: {len(linhas)} secoes conferidas contra o terreno\n")

    for nome, k in (("achatamento (fracao numa cota so)", "achatamento"),
                    ("terreno preservado (fracao)", "preservado"),
                    ("escavacao maxima (m)", "escavacao_max"),
                    ("aterro maximo no miolo (m)", "aterro_max")):
        v = np.array([r[k] for r in linhas], float)
        print(f"  {nome:<36} mediana {np.median(v):7.2f}   "
              f"p90 {np.percentile(v, 90):7.2f}   max {v.max():7.2f}")

    ruins = [(r, avaliar(r)) for r in linhas]
    ruins = [(r, m) for r, m in ruins if m]
    print(f"\n  secoes que NAO representam o terreno: {len(ruins)} de {len(linhas)}")
    if ruins:
        print(f"\n{'rio':<16}{'rch':<5}{'RS':>10}  motivo")
        for r, m in sorted(ruins,
                           key=lambda t: -t[0]["achatamento"])[:25]:
            print(f"{r['river']:<16}{r['reach']:<5}{r['rs']:>10.1f}  "
                  f"{'; '.join(m)[:70]}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(rodar(sys.argv[1], sys.argv[2],
                   sys.argv[3] if len(sys.argv) > 3 else "EPSG:31982"))
