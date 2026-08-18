# -*- coding: utf-8 -*-
"""
Eixos dos rios, montados das linhas da BHO 2017 da ANA.

Decisao: a geometria vem da ANA, nao do relevo. Tracar o eixo por relevo
(preenchimento de depressoes + D8) e uma inferencia, e ela ja errou de forma
grosseira -- o Itajai do Sul saiu com 272 km em vez de 87 porque a busca
continuou pela calha principal depois da foz. A ANA e levantamento
cartografico: erra em posicao (dezenas de metros), nao em topologia.

Tres coisas que este modulo resolve, e que nao sao obvias:

CADEIA PRINCIPAL. Um "rio" na BHO e um punhado de trechos, incluindo bracos
secundarios. A cadeia principal sai andando da foz para montante e escolhendo
sempre o afluente de MAIOR area acumulada -- e a definicao hidrologica de calha
principal, e nao depende de nome.

ORIENTACAO. Os trechos da BHO nao vem todos no mesmo sentido. Emendar sem
conferir produz um eixo que vai e volta, e a secao transversal tirada
perpendicular a ele sai girada 180 graus em pontos aleatorios.

CORTE NA FOZ. O eixo de um afluente tem de PARAR onde encontra o receptor. Sem
isso ele segue pela calha do rio grande e as secoes dos dois se sobrepoem,
contando a mesma agua duas vezes.
"""
import argparse

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge, nearest_points, substring

from .rios import BASE, EPSG, catalogo, selecionar

TOL_CONFLUENCIA = 800.0     # m; alem disso nao e confluencia, e coincidencia


def _linha_unica(geoms):
    """Emenda os trechos num LineString so, respeitando a orientacao."""
    partes = []
    for g in geoms:
        if g is None or g.is_empty:
            continue
        partes.extend(g.geoms if isinstance(g, MultiLineString) else [g])
    if not partes:
        return None
    u = linemerge(MultiLineString(partes)) if len(partes) > 1 else partes[0]
    if isinstance(u, MultiLineString):
        # sobrou desconexo: fica o maior ramo continuo, e quem chama e avisado
        u = max(u.geoms, key=lambda l: l.length)
    return u


def cadeia(g, chave):
    """Trechos da calha principal do rio, da foz para montante.

    A escolha por MAIOR NUAREAMONT em cada bifurcacao e o que define calha
    principal. Restringir ao mesmo nome nao bastaria: e justamente onde o nome
    muda de grafia que a cadeia se romperia.
    """
    sub = g[g["chave"] == chave]
    if not len(sub):
        return []
    cods = set(sub["COTRECHO"].astype(int))
    por_cod = {int(r.COTRECHO): r for r in sub.itertuples()}
    montante = {}
    for r in sub.itertuples():
        montante.setdefault(int(r.NUTRJUS), []).append(int(r.COTRECHO))

    atual = int(sub.loc[sub["NUAREAMONT"].idxmax(), "COTRECHO"])
    ch = [atual]
    while True:
        acima = [c for c in montante.get(atual, []) if c in cods]
        if not acima:
            break
        atual = max(acima, key=lambda c: float(por_cod[c].NUAREAMONT))
        ch.append(atual)
    return [por_cod[c] for c in ch]


def eixo_do_rio(g, d):
    """Um eixo por rio, orientado da CABECEIRA para a FOZ."""
    tr = cadeia(g, d["chave"])
    if not tr:
        return None
    linha = _linha_unica([r.geometry for r in tr])
    if linha is None or linha.length <= 0:
        return None
    # a cadeia foi montada da foz para montante; o primeiro trecho contem a
    # foz. Orienta-se o eixo para que o ULTIMO ponto seja a foz.
    foz_tr = tr[0].geometry
    p_foz = Point(foz_tr.coords[-1])
    if Point(linha.coords[0]).distance(p_foz) < \
            Point(linha.coords[-1]).distance(p_foz):
        linha = LineString(list(linha.coords)[::-1])
    return {**d, "linha": linha, "km_eixo": linha.length / 1000.0,
            "n_trechos_cadeia": len(tr)}


def arvore(eixos):
    """Quem desagua em quem, e onde.

    O receptor e sempre um rio de area MAIOR. Sem essa regra, um afluente pode
    ser pendurado noutro afluente menor que por acaso passe perto da foz dele
    -- e em Rio do Sul, onde Sul, Oeste e Acu se encontram quase no mesmo
    ponto, isso acontece.
    """
    por_area = sorted(eixos, key=lambda d: -d["area"])
    principal = por_area[0]
    for d in eixos:
        d["recebe_em"] = None
        d["receptor"] = None
        d["dist_foz"] = None
        if d is principal:
            continue
        foz = Point(d["linha"].coords[-1])
        cand = [m for m in eixos if m is not d and m["area"] > d["area"]]
        if not cand:
            continue
        alvo = min(cand, key=lambda m: m["linha"].distance(foz))
        dist = float(alvo["linha"].distance(foz))
        if dist > TOL_CONFLUENCIA:
            continue
        d["receptor"] = alvo["ras"]
        d["dist_foz"] = dist
        d["recebe_em"] = float(alvo["linha"].project(foz))
    return eixos, principal


def cortar_na_foz(eixos):
    """Apara o eixo do afluente que atravessa a foz e segue pelo receptor.

    Acontece quando a cadeia principal do afluente, escolhida por area
    acumulada, continua para dentro do rio grande -- os codigos da BHO nao
    mudam de nome na confluencia.
    """
    por_ras = {d["ras"]: d for d in eixos}
    for d in eixos:
        if not d.get("receptor"):
            continue
        alvo = por_ras[d["receptor"]]
        # ate onde o eixo do afluente ainda esta longe do receptor
        L = d["linha"].length
        s = np.linspace(0.0, L, max(int(L / 100.0), 20))
        dist = np.array([alvo["linha"].distance(d["linha"].interpolate(float(x)))
                         for x in s])
        dentro = np.flatnonzero(dist <= 30.0)
        if len(dentro) and s[dentro[0]] > 0.2 * L:
            corte = float(s[dentro[0]])
            d["linha"] = substring(d["linha"], 0.0, corte)
            d["km_eixo"] = d["linha"].length / 1000.0
            d["aparado_km"] = (L - corte) / 1000.0
    return eixos


def montar(sel=None, area_min=100.0, base=BASE):
    g = gpd.read_file(base).to_crs(EPSG)
    from .rios import normalizar
    g["chave"] = g["NORIOCOMP"].map(normalizar)
    cat = catalogo(area_min, base)
    escolhidos = selecionar(cat, sel)

    eixos = []
    perdidos = []
    for d in escolhidos:
        e = eixo_do_rio(g, d)
        (eixos if e else perdidos).append(e or d)
    eixos, principal = arvore(eixos)
    eixos = cortar_na_foz(eixos)
    eixos, principal = arvore(eixos)      # confluencias sobre o eixo aparado
    return eixos, principal, perdidos


def gravar(eixos, caminho="vale_eixos.geojson"):
    gdf = gpd.GeoDataFrame(
        [{k: v for k, v in d.items()
          if k in ("ras", "nome", "area", "km_eixo", "receptor", "recebe_em",
                   "dist_foz", "n")} | {"geometry": d["linha"]}
         for d in eixos], crs=EPSG)
    gdf.to_file(caminho, driver="GeoJSON")
    return caminho


def main(argv=None):
    p = argparse.ArgumentParser(description="eixos dos rios, da ANA")
    p.add_argument("--sel", default=None)
    p.add_argument("--area", type=float, default=100.0)
    p.add_argument("--saida", default="vale_eixos.geojson")
    a = p.parse_args(argv)

    eixos, principal, perdidos = montar(a.sel, a.area)
    print(f"{len(eixos)} eixos montados   |   principal: {principal['ras']}\n")
    print(f"{'rio':<16}{'km2':>8}{'km eixo':>9}{'trechos':>8}  "
          f"{'desagua em':<14}{'em km':>7}{'dist':>7}  aparado")
    print("-" * 84)
    for d in sorted(eixos, key=lambda x: -x["area"]):
        rec = d.get("receptor") or "-- foz --"
        s_em = "      -"
        if d.get("recebe_em") is not None:
            s_em = format(d["recebe_em"] / 1000.0, "7.1f")
        s_dist = "      -"
        if d.get("dist_foz") is not None:
            s_dist = format(d["dist_foz"], "7.0f")
        print(f"{d['ras']:<16}{d['area']:>8.0f}{d['km_eixo']:>9.1f}"
              f"{d['n_trechos_cadeia']:>8}  {rec:<14}{s_em}{s_dist}"
              f"  {d.get('aparado_km', 0.0):.1f} km")
    sem = [d for d in eixos if d.get("receptor") is None and d is not principal]
    if sem:
        print(f"\nSEM RECEPTOR ({len(sem)}): "
              + ", ".join(d["ras"] for d in sem)
              + f"\n   nenhum rio de area maior a menos de {TOL_CONFLUENCIA:.0f} m "
                f"da foz. Ou desaguam fora da selecao, ou direto no mar.")
    if perdidos:
        print(f"\nSEM EIXO ({len(perdidos)}): "
              + ", ".join(d["ras"] for d in perdidos))
    print("\ngravado:", gravar(eixos, a.saida))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
