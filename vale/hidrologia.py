# -*- coding: utf-8 -*-
"""
Vazoes: hidrogramas de cabeceira, contribuicao lateral, mare e a condicao
inicial acumulada pela rede.

Tres erros que este arquivo evita, todos medidos num modelo que rodava:

AREA CONTADA DUAS VEZES. A area propria do rio entrava no contorno de montante
E na contribuicao lateral: o total injetado passava de 5.700 para 8.468 m3/s.
Aqui a area propria e repartida -- FRACAO_CABECEIRA no contorno, o resto
distribuido ao longo do rio.

AREA QUE NAO ENTRAVA. So os contornos de cabeceira eram gravados, e a area
incremental de cada rio ficava de fora: entravam 2.788 de 5.700 m3/s. O modelo
partia cheio e DRENAVA (204.777 -> 179.238 hm3) em vez de encher.

CONDICAO INICIAL CHUTADA. Havia duas regras -- cabeceira usava o primeiro valor
da propria serie, os demais trechos um max(0,15*pico, 20) inventado -- e
nenhuma somava a vazao lateral, que em t=0 e o termo dominante. O Itajai_Sul
partia com 3 m3/s recebendo 50, e o Itajai_Acu com 466 quando a bacia inteira
injetava 333. O HEC-RAS monta o remanso inicial com esses numeros e no passo 1
recebe outros: o sistema leva um choque. Acumular pela rede garante
continuidade em toda juncao no instante inicial, que e o que o remanso
pressupoe.
"""
import datetime

import numpy as np

Q_REF_FOZ = 5700.0          # m3/s; pico de 1983 na foz, referencia do sintetico
INICIO = datetime.datetime(2026, 8, 1, 0, 0)

# (cota de espera em m, vazao maxima do vertedouro em m3/s)
BARRAGENS = {"sul": (110.0, 1200.0), "oeste": (110.0, 1500.0),
             "norte": (357.0, 3000.0)}
# Em 1983 a Oeste estava EM CONSTRUCAO; so Sul e Norte operavam.
BARRAGENS_POR_EVENTO = {"1983": {"sul", "norte"},
                        "2008": {"sul", "oeste", "norte"},
                        "2011": {"sul", "oeste", "norte"},
                        "2023": {"sul", "oeste", "norte"}}


def hidrograma(pico, horas=192, t_pico=72.0, base_frac=0.02, m=4.0):
    """Hidrograma sintetico suave, com base e um pico bem definido."""
    t = np.arange(horas, dtype=float)
    x = t / max(t_pico, 1.0)
    f = (x ** m) * np.exp(m * (1.0 - x))
    f = f / (f.max() or 1.0)
    return np.round(pico * (base_frac + (1.0 - base_frac) * f), 2)


def mare(horas=192, amplitude=0.6, media=0.3, periodo=12.42):
    """Mare na foz, como Stage Hydrograph."""
    t = np.arange(horas, dtype=float)
    return np.round(media + amplitude * np.sin(2 * np.pi * t / periodo), 3)


def puls_barragem(q_in, volume_hm3, q_vert_max, dt_h=1.0):
    """Amortecimento por reservatorio (Puls simplificado).

    Com o reservatorio cheio a saida e o MAIOR entre a vazao afluente e a
    saida anterior, limitada pelo vertedouro -- nunca mais que a afluente.
    Escrito como 'q + 2*excedente', o modelo fazia a barragem AMPLIFICAR o
    pico: a Norte saia com 3.000 m3/s recebendo 2.649.
    """
    q_in = np.asarray(q_in, float)
    vol = 0.0
    vmax = volume_hm3 * 1e6
    q_out = np.zeros_like(q_in)
    q = float(q_in[0])
    for t in range(len(q_in) - 1):
        entra = 0.5 * (q_in[t] + q_in[t + 1]) * dt_h * 3600.0
        vol += entra - q * dt_h * 3600.0
        if vol < 0.0:
            vol = 0.0
        if vol >= vmax:
            q = min(max(q, float(q_in[t + 1])), q_vert_max)
            vol = vmax
        else:
            q = min(q_vert_max, max(0.1 * float(q_in[t + 1]), q))
        q_out[t + 1] = q
    q_out[0] = q_in[0]
    return np.round(q_out, 2)


def series(op, eixos, xs_por_rio, arvore, log=print):
    """Monta cabeceiras, laterais e a vazao inicial de cada trecho.

    xs_por_rio: {ras -> lista de secoes}, ja condicionadas.
    arvore: {ras -> receptor_ras or None}, e {ras -> chainage da confluencia}.
    """
    area_total = max(d["area"] for d in eixos)
    por_ras = {d["ras"]: d for d in eixos}
    filhos = {}
    for d in eixos:
        if d.get("receptor"):
            filhos.setdefault(d["receptor"], []).append(d["ras"])

    cab, lat = [], []
    for d in eixos:
        ras = d["ras"]
        xs = xs_por_rio.get(ras) or []
        if len(xs) < 4:
            continue
        propria = d["area"] - sum(por_ras[f]["area"] for f in filhos.get(ras, []))
        propria = max(propria, 1.0)
        pico_total = Q_REF_FOZ * propria / area_total

        q_cab = hidrograma(pico_total * op.fracao_cabeceira, op.horas)
        cab.append({"rio": ras, "reach": "R1", "serie": q_cab,
                    "q_pico": float(q_cab.max()), "xs": xs})

        rss = sorted((x["rs"] for x in xs), reverse=True)
        lat.append({"rio": ras, "reach": "R1",
                    "rs_hi": rss[1], "rs_lo": rss[-2],
                    "serie": hidrograma(pico_total * (1.0 - op.fracao_cabeceira),
                                        op.horas)})

    if op.barragens and op.evento:
        alvo = BARRAGENS_POR_EVENTO.get(op.evento, set(BARRAGENS))
        for c in cab:
            chave = c["rio"].lower().replace("itajai_", "")
            if chave in alvo and chave in BARRAGENS:
                _, qmax = BARRAGENS[chave]
                antes = float(c["serie"].max())
                c["serie"] = puls_barragem(c["serie"], 100.0, qmax)
                log(f"      barragem {chave}: pico {antes:.0f} -> "
                    f"{c['serie'].max():.0f} m3/s")

    return cab, lat


def inicial_por_trecho(trechos, juncoes, cabs, lats, log=print):
    """Vazao em t=0 na secao de MONTANTE de cada trecho.

    POR TRECHO, E NAO POR RIO. A versao anterior (inicial()) devolvia UM valor
    por rio -- a vazao acumulada na FOZ dele -- e o passo de escrita gravava
    esse mesmo valor em todo trecho, inclusive no "Initial RS" da secao mais de
    MONTANTE. No Iraputa isso declarava 4 m3/s na cabeceira, onde o contorno
    entrega 0,21: o modelo era montado cheio e esvaziado no primeiro passo por
    um contorno 19 vezes menor, e o solver batia as 40 iteracoes com 2 cm de
    lamina. Acontecia em todos os doze rios, de 19 a 38 vezes.

    Os dois numeros estavam certos, cada um no seu lugar: 0,21 m3/s e o que
    entra na cabeceira (a cabeceira leva 5% da area do rio, o resto entra como
    vazao lateral ao longo do trecho) e 4 m3/s e o que sai na foz. O erro era
    escrever o segundo na posicao do primeiro.

    Percorre de montante para jusante somando o que chega: contorno de
    cabeceira, o trecho anterior do mesmo rio, os afluentes que desaguam na
    cabeceira do trecho, e a lateral que corre dentro dele.
    """
    q_cab = {(c["rio"], c["reach"]): float(c["serie"][0]) for c in cabs}
    q_lat = {(l["rio"], l["reach"]): float(l["serie"][0]) for l in lats}
    entra_de = {}
    for j in juncoes:
        entra_de.setdefault(tuple(j["dn"]), []).extend(
            tuple(u) for u in j.get("up", []))

    entrada, saida = {}, {}
    # (area, a): afluente antes do receptor, e dentro do rio de montante p/
    # jusante -- o receptor tem sempre area maior que quem desagua nele
    for t in sorted(trechos, key=lambda t: (t["area"], t["a"])):
        k = (t["rio"], t["reach"])
        q = q_cab.get(k, 0.0)
        # A juncao do HEC-RAS lista TODOS os trechos de montante, inclusive o
        # do proprio rio. Somar a juncao E o trecho anterior conta o mesmo
        # escoamento duas vezes, e o erro COMPOE trecho a trecho: no Acu, com
        # cinco trechos, a vazao inicial da foz saia 2.861 m3/s em vez de 114
        # -- 25 vezes. Onde ha juncao ela ja e a soma completa; o trecho
        # anterior so entra quando nao ha juncao alimentando este.
        ups = entra_de.get(k)
        if ups:
            q += sum(saida.get(u, 0.0) for u in ups)
        else:
            anteriores = [s for s in trechos
                          if s["rio"] == t["rio"] and s["b"] <= t["a"] + 1e-6]
            if anteriores:
                ant = max(anteriores, key=lambda s: s["b"])
                q += saida.get((ant["rio"], ant["reach"]), 0.0)
        entrada[k] = q
        saida[k] = q + q_lat.get(k, 0.0)

    v = list(entrada.values())
    if v:
        log(f"      vazao inicial na cabeceira de cada trecho: "
            f"{min(v):.2f} a {max(v):.2f} m3/s")
    return entrada


def inicial(eixos, cab, lat, log=print):
    """Vazao inicial de cada rio na FOZ dele. NAO SERVE PARA "Initial RS".

    Mantida porque o numero e util para conferencia -- e o que o rio entrega no
    fim --, mas nao e o que vai no .u01: la o valor tem de ser a vazao na secao
    onde a linha e escrita, que e a de MONTANTE do trecho. Usar este valor la
    foi o que pos 4 m3/s na cabeceira do Iraputa contra 0,21 do contorno.
    Use inicial_por_trecho().

    'eixos' vem ordenado por area; percorre-se do MENOR para o maior, porque o
    receptor tem sempre area maior que o afluente e precisa dos afluentes
    prontos antes.
    """
    q0_lat = {l["rio"]: float(l["serie"][0]) for l in lat}
    q0_cab = {c["rio"]: float(c["serie"][0]) for c in cab}
    saida = {}
    for d in sorted(eixos, key=lambda x: x["area"]):
        ras = d["ras"]
        entra = q0_cab.get(ras, 0.0) + q0_lat.get(ras, 0.0)
        entra += sum(saida.get(m["ras"], 0.0) for m in eixos
                     if m.get("receptor") == ras)
        saida[ras] = entra
    foz = max(eixos, key=lambda d: d["area"])["ras"]
    log(f"      vazao inicial acumulada: {min(saida.values()):.1f} a "
        f"{max(saida.values()):.1f} m3/s (foz {saida[foz]:.1f})")
    return saida
