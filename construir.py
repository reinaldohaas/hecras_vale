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

from itajai import config, terreno, tracado, topologia, secao, perfil, escrita

# cheia de referencia (2008: 5.700 m3/s em Itajai), rateada por area
Q_REF_FOZ = 5700.0
N_HORAS = 97
INICIO = datetime.datetime(2026, 8, 1)
MARE_MEDIA, MARE_AMP, MARE_PERIODO = 0.30, 0.50, 12.42
NOME_JUNCAO = {0: "Rio_do_Sul", 1: "Ibirama", 2: "Indaial", 3: "Itajai"}


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


def main():
    t0 = time.time()
    projeto = config.PROJETO
    print("=" * 70)
    print(f"{projeto}  |  construcao a partir do relevo")
    print("=" * 70)

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
        s = rede[m]["linha"].project(Point(rede[k]["linha"].coords[-1]))
        conf[m].append({"k": k, "s": s, "pt": rede[m]["linha"].interpolate(s)})
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
            t = {"rio": v["nome"], "reach": f"R{len(por_rio[k])+1}", "k": k,
                 "linha": v["linha"], "xs": sel, "a": a, "b": b}
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
    cabeceiras = []
    for k in ordem:
        # rio que NASCE de juncao nao pode ter contorno proprio: seria vazao
        # contada duas vezes
        if any(abs(c["s"]) < 1.0 for c in conf[k]) or not por_rio[k]:
            continue
        t = por_rio[k][0]
        propria = rede[k]["area"] - sum(rede[c["k"]]["area"] for c in conf[k])
        t["q_pico"] = Q_REF_FOZ * max(propria, 1.0) / area_total
        t["serie"] = hidrograma(t["q_pico"])
        t["q_base"] = float(t["serie"][0])
        cabeceiras.append(t)
        print(f"    {t['rio']:<14} Q pico {t['q_pico']:7.1f} m3/s")

    # -------------------------------------------------------------- escrita
    print("\n[6] escrita")
    saida = por_rio[topologia.PRINCIPAL][-1]
    print("   ", escrita.geometria(projeto, trechos, juncoes,
                                   f"{projeto} - eixo do relevo Copernicus"))
    print("   ", escrita.fluxo(projeto, cabeceiras, saida, mare(), N_HORAS))
    print("   ", escrita.plano(projeto, INICIO, N_HORAS))
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
