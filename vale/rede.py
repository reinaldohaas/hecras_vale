# -*- coding: utf-8 -*-
"""
Rede: parte cada rio em trechos nas confluencias e monta as juncoes.

No HEC-RAS um rio e cortado em TRECHOS a cada confluencia, e cada confluencia
vira uma juncao com N trechos de montante e um de jusante. Tres regras que o
programa tem de respeitar, e que custaram rodadas:

REACH XY E DO TRECHO. Gravar o rio inteiro como geometria de cada trecho faz a
juncao virar teia de aranha no RAS Mapper e as bank lines saem erradas.

JUNC L&A E O CAMINHO ATRAVES DA JUNCAO -- da ultima secao de montante ate a
primeira de jusante. Gravar um valor fixo (500 m) onde a geometria da 75
desequilibra a continuidade exatamente na secao que o solver depois reporta.

JUNCAO COM UM ENTRANDO E UM SAINDO NAO EXISTE. O HEC-RAS recusa a geometria
antes de computar. Acontece quando o afluente entra tao perto da cabeceira do
receptor que o trecho de cima fica sem secoes -- e o caso do Iraputa, que entra
no Itajai_Norte a 5 km da cabeceira dele.
"""
import numpy as np
from shapely.ops import substring

MIN_SECOES_TRECHO = 3
NASCENTE_TOL = 50.0        # m; confluencia ate aqui e a NASCENTE do receptor


def _confluencias(eixos):
    """Por rio receptor, a lista de confluencias ordenada rio abaixo."""
    conf = {d["ras"]: [] for d in eixos}
    for d in eixos:
        if d.get("receptor") and d.get("recebe_em") is not None:
            conf[d["receptor"]].append({"de": d["ras"],
                                        "s": float(d["recebe_em"])})
    for k in conf:
        conf[k].sort(key=lambda c: c["s"])
    return conf


def montar(eixos, xs_por_rio, log=print):
    """Devolve (trechos, juncoes).

    trecho: rio, reach, linha, xs, a, b  (a e b sao estacas ao longo do eixo)
    juncao: nome, x, y, up [(rio,reach)], dn (rio,reach), dists
    """
    conf = _confluencias(eixos)
    por_ras = {d["ras"]: d for d in eixos}
    trechos, juncoes = [], []

    for d in eixos:
        ras = d["ras"]
        linha = d["linha"]
        L = linha.length
        xs = list(xs_por_rio.get(ras) or [])
        if not xs:
            continue
        # cortes validos: confluencia que deixe secoes dos dois lados
        cortes = []
        for c in conf.get(ras, []):
            rs_c = L - c["s"]
            acima = [x for x in xs if x["rs"] > rs_c]
            abaixo = [x for x in xs if x["rs"] <= rs_c]
            if len(acima) >= MIN_SECOES_TRECHO and \
                    len(abaixo) >= MIN_SECOES_TRECHO:
                cortes.append(c)
            else:
                log(f"      {ras}: confluencia de {c['de']} em "
                    f"{c['s']/1000:.1f} km nao divide o rio "
                    f"({len(acima)} secoes acima, {len(abaixo)} abaixo) -- "
                    f"os dois entram no MESMO trecho")
                c["sem_corte"] = True

        limites = [0.0] + [c["s"] for c in cortes] + [L]
        for i in range(len(limites) - 1):
            a, b = limites[i], limites[i + 1]
            sub = [x for x in xs if (L - b) < x["rs"] <= (L - a) + 1e-6]
            if len(sub) < 2:
                continue
            trechos.append({"rio": ras, "reach": f"R{i+1}",
                            "linha": substring(linha, a, b),
                            "xs": sorted(sub, key=lambda x: -x["rs"]),
                            "a": a, "b": b, "area": d["area"]})

    por_rio = {}
    for t in trechos:
        por_rio.setdefault(t["rio"], []).append(t)
    for v in por_rio.values():
        v.sort(key=lambda t: t["a"])

    # ------------------------------------------------ confluencia na NASCENTE
    # Um afluente pode entrar na quilometragem ZERO do receptor: o rio grande
    # NASCE da juncao. E o caso do Itajai-Acu, que comeca em Rio do Sul do
    # encontro do Sul com o Oeste.
    #
    # Nao ha o que dividir -- nao existe trecho do receptor acima do ponto --,
    # entao a regra geral (dividir o receptor e juntar as duas metades) nao se
    # aplica, e sem tratamento os dois afluentes ficam SEM JUNCAO: desligados
    # da rede, com a agua deles nao entrando em lugar nenhum e o receptor
    # recebendo um contorno de cabeceira que nao deveria ter.
    #
    # A juncao certa aqui tem os afluentes como montante e o PRIMEIRO trecho do
    # receptor como jusante -- dois entrando, um saindo, que e o que o HEC-RAS
    # exige. Varios afluentes na mesma nascente entram na MESMA juncao.
    for ras, cs in conf.items():
        if ras not in por_rio:
            continue
        nasc = [c for c in cs
                if c["s"] <= NASCENTE_TOL and c["de"] in por_rio]
        if not nasc:
            continue
        dn = por_rio[ras][0]
        ups = [por_rio[c["de"]][-1] for c in nasc]
        if len(ups) < 2 and len(ups) + 1 < 2:
            continue
        p_j = por_ras[ras]["linha"].interpolate(0.0)
        juncoes.append({
            "nome": f"Nasce_{ras}"[:16],
            "x": p_j.x, "y": p_j.y,
            "up": [(t["rio"], t["reach"]) for t in ups],
            "dn": (dn["rio"], dn["reach"]),
            # mesma correcao da juncao de foz: a distancia geometrica do fim do
            # eixo ate a juncao e ZERO, porque o eixo ja foi aparado ali. O que
            # o HEC-RAS quer e da ultima SECAO ate a juncao, que e a RS dela.
            "dists": [max(float(t["xs"][-1]["rs"]), 1.0) for t in ups]})
        for c in nasc:
            c["na_nascente"] = True
        log(f"      {ras} nasce da juncao de "
            f"{', '.join(c['de'] for c in nasc)}")

    # ------------------------------------------------------------ juncoes
    for ras, cs in conf.items():
        if ras not in por_rio:
            continue
        for c in cs:
            if (c.get("sem_corte") or c.get("na_nascente")
                    or c["de"] not in por_rio):
                continue
            recept = por_rio[ras]
            # trecho de jusante: o que COMECA nesta confluencia
            dn = next((t for t in recept if abs(t["a"] - c["s"]) < 1.0), None)
            up_rec = next((t for t in recept if abs(t["b"] - c["s"]) < 1.0), None)
            if dn is None:
                continue
            afl = por_rio[c["de"]][-1]           # ultimo trecho do afluente
            ups = [t for t in (up_rec, afl) if t is not None]
            if len(ups) < 2:
                # um entrando e um saindo: o HEC-RAS recusa. Melhor nao gravar
                # a juncao do que gravar uma invalida.
                log(f"      juncao de {c['de']} em {ras} descartada: "
                    f"so {len(ups)} trecho(s) de montante")
                continue
            p = por_ras[ras]["linha"].interpolate(c["s"])
            dists = []
            for t in ups:
                # caminho ATRAVES da juncao: da ultima secao do trecho de
                # montante ate a primeira do de jusante
                # DO AFLUENTE, A RS DA ULTIMA SECAO -- e nao a distancia
                # geometrica do fim do eixo ate a juncao. O eixo do afluente ja
                # foi APARADO para terminar exatamente na confluencia
                # (eixos.cortar_na_foz), entao essa distancia e zero e o
                # max(...,1.0) a transformava em 1 METRO. Todas as dez juncoes
                # gravavam "Junc L&A=1.00" para o afluente.
                #
                # Um trecho computacional de 1 m encostado em trechos de 25 a
                # 150 m, com dt=5 s, da Courant da ordem de 10, e o solver
                # implicito responde batendo o teto de iteracoes: o modelo
                # integrado instabilizava aos 15 SEGUNDOS em Benedito R2,
                # Norte R1/R2, Cedros e Acu R3 ao mesmo tempo, com erros de
                # nivel de apenas 0,10 a 0,45 m. Rio isolado nao tem juncao e
                # nunca exercitava este numero -- por isso todos completavam
                # sozinhos e morriam juntos.
                #
                # A RS e medida da foz (rs = L - s), entao a RS da ultima secao
                # E a distancia dela ate a confluencia.
                dists.append(max(abs(t["xs"][-1]["rs"] - dn["xs"][0]["rs"]),
                                 1.0) if t is up_rec else
                             max(float(t["xs"][-1]["rs"]), 1.0))
            juncoes.append({"nome": f"Foz_{c['de']}"[:16],
                            "x": p.x, "y": p.y,
                            "up": [(t["rio"], t["reach"]) for t in ups],
                            "dn": (dn["rio"], dn["reach"]),
                            "dists": dists})

    log(f"      {len(trechos)} trechos, {len(juncoes)} juncoes")
    return trechos, juncoes


def contornos_necessarios(trechos, juncoes):
    """Trechos que precisam de contorno de montante.

    Criterio do HEC-RAS: todo trecho que NAO e o 'dn' de alguma juncao precisa
    de condicao de contorno a montante. Heuristicas do tipo "nasceu de juncao"
    conflitam com o recuo das secoes e produzem o erro
    "R1 needs a upstream boundary condition".
    """
    com_dn = {(j["dn"][0], j["dn"][1]) for j in juncoes}
    return [t for t in trechos if (t["rio"], t["reach"]) not in com_dn]
