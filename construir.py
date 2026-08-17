# -*- coding: utf-8 -*-
"""
Constroi o modelo da bacia do Itajai, do relevo ao projeto HEC-RAS.

Reescrita limpa. O gerador anterior tinha 1.376 linhas com quinze correcoes
empilhadas que interagiam entre si -- o resultado oscilava entre 1 e 48 dos 192
passos conforme se mexia num parametro, sem direcao. Aqui cada etapa e um
modulo pequeno, verificavel sozinho:

    itajai/config.py      nome do projeto, EPSG, projecao
    itajai/terreno.py     relevo Copernicus: UTM, amostragem, terreno .hdf
    itajai/tracado.py     eixo dos rios pelo relevo (priority-flood + D8)
    itajai/topologia.py   quem desagua em quem, e area de drenagem (ANA)
    itajai/secao.py       onde e como cortar, e onde esta a calha
    itajai/perfil.py      condicionamento longitudinal e Manning
    itajai/escrita.py     .g01/.u01/.p01/.prj/.rasmap, com validacao

Uso:  python construir.py                 projeto Tajai, cheia sintetica
      python construir.py --terreno       refaz tambem o terreno do RAS Mapper
"""
import datetime
import os
import sys
import time

import numpy as np
from shapely.geometry import Point
from shapely.ops import substring

from itajai import config, terreno, tracado, topologia, secao, perfil, escrita

# cheia de referencia (2008: 5.700 m3/s em Itajai), rateada por area
Q_REF_FOZ = 5700.0
N_HORAS = 97
INICIO = datetime.datetime(2026, 8, 1)
MARE_MEDIA, MARE_AMP, MARE_PERIODO = 0.30, 0.50, 12.42
NOME_JUNCAO = {0: "Rio_do_Sul", 1: "Ibirama", 2: "Indaial", 3: "Itajai"}
FRACAO_CABECEIRA = 0.05   # da area propria; o resto entra como vazao lateral


def hidrograma(pico, base=None, n=N_HORAS, tp=26, te=46):
    """Cheia sintetica: subida rapida, recessao longa."""
    base = base if base is not None else max(pico * 0.15, 20.0)
    t = np.arange(n, dtype=float)
    # np.where avalia os DOIS ramos, entao a recessao recebe (t-tp) negativo e
    # a potencia 1,5 devolve NaN. O clip nao muda o resultado -- so evita que o
    # ramo descartado seja calculado sobre valor invalido.
    subida = base + (pico - base) * (t / tp) ** 2
    recessao = base + (pico - base) * np.exp(-(np.maximum(t - tp, 0.0) / te) ** 1.5)
    return np.round(np.where(t <= tp, subida, recessao), 1)


def mare(n=N_HORAS):
    t = np.arange(n, dtype=float)
    return np.round(MARE_MEDIA + MARE_AMP * np.sin(2 * np.pi * t / MARE_PERIODO), 3)


def series_evento(evento, barragens):
    """Hidrogramas reais do evento, por sub-bacia. None = cheia sintetica."""
    if not evento:
        return None, N_HORAS
    from itajai import hidrologia
    q, n = hidrologia.hidrogramas(evento, barragens=barragens)
    print(f"    evento {evento}, barragens "
          f"{'ATIVAS' if barragens else 'ABERTAS (sem obras)'}: {n} h")
    for k, v in q.items():
        print(f"      {k:<10} pico {v.max():8.1f} m3/s na hora {int(v.argmax()):3d}")
    return q, n


def fatia(q_ev, rede, receptor, k, area_km2, a_ref):
    """Serie do evento para uma area, tirada da fonte que cobre o rio k.

    Os afluentes de 2a ordem nao tem serie propria -- a area deles ja esta
    dentro da do rio que os recebe. Cada um leva a fatia proporcional a area, e
    o receptor fica com o resto: o volume do evento nao muda ao detalhar a rede.
    """
    if q_ev is None:
        return None
    alvo, base = k, rede[k]["area"]
    while alvo is not None and alvo not in q_ev:
        alvo = receptor.get(alvo)
        if alvo is not None:
            base = rede[alvo]["area"]
    if alvo is None:
        alvo, base = "acu_incr", a_ref
    if alvo not in q_ev:
        return None
    return q_ev[alvo] * (area_km2 / max(base, 1.0))


def main():
    t0 = time.time()
    projeto = config.PROJETO
    evento = next((a for a in sys.argv[1:] if a.isdigit()), None)
    barragens = "--sem-barragens" not in sys.argv
    if evento:
        projeto = f"{config.PROJETO}_{evento}"
        if not barragens:
            projeto += "_sb"
    print("=" * 70)
    print(f"{projeto}  |  construcao a partir do relevo")
    print("=" * 70)
    print()
    q_ev, n_horas = series_evento(evento, barragens)
    a_ref = None

    # ------------------------------------------------------------ topologia
    print("\n[1] topologia (ANA BHO 2017)")
    rede = topologia.carregar()
    receptor, _ = topologia.arvore(rede)
    area_total = rede[topologia.PRINCIPAL]["area"]
    for k, v in sorted(rede.items(), key=lambda x: -x[1]["area"]):
        alvo = receptor.get(k)
        print(f"    {v['nome']:<14} {v['area']:8.1f} km2"
              + (f"  -> {rede[alvo]['nome']}" if alvo else "   (calha principal)"))

    # --------------------------------------------------------------- eixos
    print("\n[2] eixo dos rios pelo relevo (priority-flood + D8)")
    terreno.preparar_utm()
    tr = tracado.Tracador()
    for k, v in sorted(rede.items(), key=lambda x: -x[1]["area"]):
        ln = tr.eixo(v["cabeceira"], v["foz"])
        if ln is None:
            raise SystemExit(f"nao consegui tracar {v['nome']}")
        v["linha"] = ln
        d = np.array([ln.distance(Point(p))
                      for p in list(v["linha_ana"].coords)[::20]])
        print(f"    {v['nome']:<14} {ln.length/1000:7.1f} km "
              f"(ANA {v['linha_ana'].length/1000:6.1f})   "
              f"afastamento mediano {np.median(d):5.0f} m")

    # -------------------------------------------------------------- secoes
    print("\n[3] secoes do relevo")
    am = terreno.Amostrador()
    # a arvore veio da ANA; a ESTACA da confluencia sai do eixo em uso
    conf = {m: [] for m in rede}
    for k, m in receptor.items():
        L = rede[m]["linha"].length
        s = rede[m]["linha"].project(Point(rede[k]["linha"].coords[-1]))
        # Uma confluencia na CABECEIRA do receptor nao deixa trecho a montante,
        # e a juncao fica com uma entrada e uma saida -- o HEC-RAS recusa:
        #   "Junction ... has only one reach flowing in and one flowing out.
        #    Junctions are for flow confluences and splits."
        # e aborta a simulacao na validacao, antes de calcular. Acontece com o
        # Taio (entra no Oeste) e o Iraputa (no Norte), porque a cadeia da ANA
        # para o receptor comeca exatamente ali. Afastar a confluencia de um
        # espacamento cria o trecho de montante que falta; o deslocamento e da
        # ordem do proprio espacamento das secoes.
        s_bruto = s
        s = float(np.clip(s, secao.ESPACAMENTO, max(L - secao.ESPACAMENTO,
                                                    secao.ESPACAMENTO)))
        # s_bruto guarda a estaca ANTES do afastamento. E ela que diz se o rio
        # NASCE da juncao -- o teste feito sobre a estaca ja afastada nunca
        # encontra nada abaixo de 1 m, e Acu, Oeste e Norte passavam a receber
        # hidrograma de cabeceira ALEM da vazao que chega pela juncao. O leitor
        # do ras-commander flagrou: 13 contornos, um por rio mais a mare,
        # quando deveriam ser 10.
        conf[m].append({"k": k, "s": s, "s_bruto": s_bruto,
                        "pt": rede[m]["linha"].interpolate(s)})
    for m in conf:
        conf[m].sort(key=lambda c: c["s"])

    ordem = sorted(rede, key=lambda k: -rede[k]["area"])
    for k in ordem:                       # o maior primeiro: os afluentes
        v = rede[k]                       # ancoram no leito de quem os recebe
        a_cab = v["area"] - sum(rede[c["k"]]["area"] for c in conf[k])
        v["xs"] = perfil.condicionar(
            secao.cortar_trecho(v["linha"], am, v["area"],
                                area_cabeceira=max(a_cab, v["area"] * 0.05)),
            v["nome"])
        extra = ""
        if k in receptor:
            m = receptor[k]
            s = next(c["s"] for c in conf[m] if c["k"] == k)
            rs_alvo = rede[m]["linha"].length - s
            leito = {d["rs"]: perfil.cota_talvegue(d) for d in rede[m]["xs"]}
            alvo = leito[min(leito, key=lambda r: abs(r - rs_alvo))]
            desl = perfil.ancorar(v["xs"], alvo)
            extra = f"   ancorado em {alvo:7.1f} m ({desl:+.1f})"
        # ESCAVA POR ULTIMO, com o perfil ja definido. Escavar no corte e
        # reajustar a cada passo do condicionamento reaplicava o trapezio sobre
        # um perfil ja alterado: a calha degenerava num pico de um ponto so
        # (uma fenda de 3,4 m num rio de calha de 106 m) e a secao deixava de
        # conduzir. Era a causa comum dos saltos de area de ate 18x entre
        # vizinhas, e do erro de balanco que o solver acusava.
        for d in v["xs"]:
            secao.escavar(d)
        # CONFERE o talvegue final, depois de escavar. O condicionamento
        # trabalha sobre z_alvo; se a geometria resultante ainda tiver degrau,
        # e porque algo entre uma coisa e outra o reintroduziu. Medir aqui e
        # barato e evita descobrir isso no log do solver.
        zt = np.array([d["z"].min() for d in v["xs"]])
        rs_v = np.array([d["rs"] for d in v["xs"]])
        subida = np.diff(zt) > 0.5           # leito subindo rio abaixo
        if subida.any():
            print(f"        ! {v['nome']}: {int(subida.sum())} degraus de "
                  f"subida no leito (max {np.diff(zt).max():.1f} m)")
        n_j = perfil.manning(v["xs"])
        print(f"    {v['nome']:<14} {len(v['xs']):4d} secoes"
              f"   Jarrett em {n_j:3d}{extra}")

    # ------------------------------------------------------------- trechos
    print("\n[4] trechos e juncoes")
    trechos, por_rio = [], {}
    for k in ordem:
        v = rede[k]
        L = v["linha"].length
        cortes = sorted({0.0} | {c["s"] for c in conf[k] if c["s"] > 1.0} | {L})
        por_rio[k] = []
        for i in range(len(cortes) - 1):
            a, b = cortes[i], cortes[i + 1]
            if b - a < secao.ESPACAMENTO:
                continue
            rs_hi, rs_lo = L - a, L - b
            sel = ([d for d in v["xs"] if rs_lo < d["rs"] <= rs_hi]
                   if por_rio[k] else
                   [d for d in v["xs"] if rs_lo <= d["rs"] <= rs_hi])
            if len(sel) < 2:
                continue
            # SEGMENTO do eixo, nao o rio inteiro. Dando a mesma linha a todos
            # os trechos, os quatro do Acu saem com os mesmos 174 km de
            # Reach XY: o RAS Mapper desenha um sobre o outro e liga os pontos
            # de juncao entre si, formando uma teia de linhas tracejadas. As
            # bank lines, que se apoiam no tracado do trecho, herdam o erro.
            t = {"rio": v["nome"], "reach": f"R{len(por_rio[k])+1}", "k": k,
                 "linha": substring(v["linha"], a, b), "xs": sel,
                 "a": a, "b": b}
            trechos.append(t)
            por_rio[k].append(t)
        print(f"    {v['nome']:<14} {len(por_rio[k])} trecho(s)")

    juncoes, n_princ = [], 0
    for m in ordem:
        for s in sorted({c["s"] for c in conf[m]}):
            entra = [c for c in conf[m] if abs(c["s"] - s) < 1.0]
            ups = [por_rio[c["k"]][-1] for c in entra if por_rio.get(c["k"])]
            up_m = [t for t in por_rio[m] if abs(t["b"] - s) < 1.0]
            dn_m = [t for t in por_rio[m] if abs(t["a"] - s) < 1.0]
            if not dn_m or not ups:
                continue
            nome = (NOME_JUNCAO.get(n_princ, f"J{n_princ+1}")
                    if m == topologia.PRINCIPAL
                    else f"Foz_{rede[entra[0]['k']]['nome']}"[:16])
            if m == topologia.PRINCIPAL:
                n_princ += 1
            rs_j = rede[m]["linha"].length - s
            dn_d = max(rs_j - dn_m[0]["xs"][0]["rs"], 0.0)
            dists = [max((t["xs"][-1]["rs"] - rs_j) if t["k"] == m
                         else t["xs"][-1]["rs"], 0.0) + dn_d for t in ups + up_m]
            juncoes.append({"nome": nome, "x": entra[0]["pt"].x,
                            "y": entra[0]["pt"].y,
                            "up": [(t["rio"], t["reach"]) for t in ups + up_m],
                            "dists": [max(d, 1.0) for d in dists],
                            "dn": (dn_m[0]["rio"], dn_m[0]["reach"])})
            print(f"      {nome:<16} {[u[0] for u in juncoes[-1]['up']]}"
                  f" -> {dn_m[0]['rio']}/{dn_m[0]['reach']}")

    # ------------------------------------------------------------ contornos
    print("\n[5] contornos")

    def area_ate(k, s):
        """Area de drenagem que ja chegou a estaca s do rio k.

        Propria do rio mais os afluentes que entraram a montante -- e,
        recursivamente, o que veio nos afluentes deles.
        """
        propria = rede[k]["area"] - sum(rede[c["k"]]["area"] for c in conf[k])
        return propria + sum(rede[c["k"]]["area"] for c in conf[k]
                             if c["s"] <= s + 1.0)

    # CRITERIO: precisa de contorno todo trecho que nao e o 'dn' de alguma
    # juncao. E o que o HEC-RAS exige, literalmente:
    #   "River: Itajai_Oeste  Reach: R1  needs a upstream boundary condition."
    # A heuristica anterior -- "o rio nasce de juncao, logo nao tem contorno" --
    # entrava em conflito com o afastamento da confluencia: mover a foz do Taio
    # 1 km rio abaixo CRIA um trecho de 1 km acima dela, que fica sem fonte.
    # Marcar o rio inteiro como alimentado pela juncao deixava esse trecho orfao.
    alimentados = {(j["dn"][0], j["dn"][1]) for j in juncoes}
    a_ref = area_total - sum(rede[c["k"]]["area"] for c in conf[topologia.PRINCIPAL]
                             if q_ev and c["k"] in q_ev)
    cabeceiras = []
    for k in ordem:
        if not por_rio[k]:
            continue
        t = por_rio[k][0]
        if (t["rio"], t["reach"]) in alimentados:
            print(f"    {rede[k]['nome']:<14} alimentado por juncao")
            continue
        # area que chega ao INICIO deste trecho: a propria do rio ate ali,
        # mais os afluentes ja incorporados. Num trecho de cabeceira curto
        # (o vao acima de uma confluencia afastada) isso e pequeno, e deve ser.
        # a cabeceira leva so a fracao que NAO vai para a lateral
        propria = (rede[k]["area"] - sum(rede[c["k"]]["area"] for c in conf[k]))                   * FRACAO_CABECEIRA
        propria += sum(rede[c["k"]]["area"] for c in conf[k]
                       if c["s"] <= t["a"] + 1.0)
        t["q_pico"] = Q_REF_FOZ * max(propria, 1.0) / area_total
        s_ev = fatia(q_ev, rede, receptor, k, propria, a_ref)
        t["serie"] = s_ev if s_ev is not None else hidrograma(t["q_pico"])
        t["q_pico"] = float(np.max(t["serie"]))
        t["q_base"] = float(t["serie"][0])
        cabeceiras.append(t)
        print(f"    {t['rio']:<14} Q pico {t['q_pico']:7.1f} m3/s")

    # vazao inicial de TODO trecho, pela area que ja chegou nele. Sem isto os
    # trechos de jusante partem indefinidos e o balanco estoura no passo 1.
    for t in trechos:
        if "q_base" in t:
            continue
        a = area_ate(t["k"], t["a"])
        t["q_pico"] = Q_REF_FOZ * max(a, 1.0) / area_total
        t["q_base"] = max(t["q_pico"] * 0.15, 20.0)
    print(f"    vazao inicial em {len(trechos)} trechos "
          f"({min(t['q_base'] for t in trechos):.0f} a "
          f"{max(t['q_base'] for t in trechos):.0f} m3/s)")

    # --- area propria de cada rio, distribuida ao longo dele
    laterais = []
    for k in ordem:
        # 95% da area propria vai para a lateral; os 5% restantes ficam na
        # cabeceira (ver FRACAO_CABECEIRA). Sem essa divisao a mesma area entra
        # duas vezes -- no contorno de montante E na lateral -- e o total
        # injetado passava de 5.700 para 8.468 m3/s.
        propria = (rede[k]["area"] - sum(rede[c["k"]]["area"] for c in conf[k]))                   * (1.0 - FRACAO_CABECEIRA)
        if propria < 1.0 or not por_rio[k]:
            continue
        L = sum(t["b"] - t["a"] for t in por_rio[k]) or 1.0
        for t in por_rio[k]:
            a = propria * (t["b"] - t["a"]) / L
            rss = sorted((d["rs"] for d in t["xs"]), reverse=True)
            if len(rss) < 4:          # a faixa nao pode tocar as extremas
                continue
            pico = Q_REF_FOZ * a / area_total
            s_ev = fatia(q_ev, rede, receptor, k, a, a_ref)
            laterais.append({"rio": t["rio"], "reach": t["reach"],
                             "rs_hi": rss[1], "rs_lo": rss[-2],
                             "serie": s_ev if s_ev is not None
                                      else hidrograma(pico)})
    q_lat = sum(float(np.max(l["serie"])) for l in laterais)
    q_cab = sum(t["q_pico"] for t in cabeceiras)
    print(f"    lateral em {len(laterais)} trechos: Q pico {q_lat:7.1f} m3/s")
    print(f"    TOTAL injetado: {q_cab + q_lat:7.1f} de {Q_REF_FOZ:.0f} m3/s")

    # -------------------------------------------------------------- escrita
    print("\n[6] escrita")
    saida = por_rio[topologia.PRINCIPAL][-1]
    print("   ", escrita.geometria(projeto, trechos, juncoes,
                                   f"{projeto} - eixo do relevo Copernicus"))
    print("   ", escrita.fluxo(projeto, trechos, cabeceiras, saida, mare(n_horas),
                                n_horas, laterais))
    print("   ", escrita.plano(projeto, INICIO, n_horas))
    if "--terreno" in sys.argv:
        hdf = terreno.preparar_hdf(config.WKT)
    else:
        p = os.path.join(terreno.PASTA, "Terreno.hdf")
        hdf = p if os.path.exists(p) else None
    print("   ", escrita.rasmap(projeto, hdf))

    v = escrita.validar(projeto)
    if v:
        print(f"\n[7] conferido pelo ras-commander: {v['secoes']} secoes, "
              f"{v['pontos']} pontos na primeira, margens {v['margens']}, "
              f"{v['zonas_manning']} zonas de Manning")
    print(f"\nconcluido em {time.time()-t0:.0f} s   ->  abra {projeto}.prj")


if __name__ == "__main__":
    main()
