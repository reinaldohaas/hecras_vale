# -*- coding: utf-8 -*-
"""
Condicionamento do perfil longitudinal.

O solver precisa de um talvegue que desca de forma monotonica e sem degrau. O
DEM bruto nao da isso: tem contrapendentes, pocos e a lamina achatada sobre a
agua.

A regra que custou caro aprender: CONDICIONAR O PERFIL NAO PODE MOVER A
PLANICIE. A versao anterior deslocava a secao inteira (z = z - delta) para
impor o perfil. Como delta chega a 5,8 m no baixo Itajai, margens e planicie
desciam junto -- a margem em Itajai ficava a -2,75 m e a de Ilhota a -0,60 m,
abaixo do nivel do mar. Isso nao vinha do relevo, vinha do deslocamento. Aqui
so a CALHA se move, com o mesmo trapezio da escavacao.

Tambem: o talvegue de referencia e sempre o da calha junto ao eixo (i_thal),
nunca z.min(). Usar o minimo global fazia o perfil e a conducao descreverem
canais diferentes quando a secao cruzava outro leito.
"""
import numpy as np

DECL_MINIMA = 1e-4       # m/m; abaixo disso o trecho vira lago e o solver oscila
DECL_MAXIMA = 0.008      # m/m; acima disso o escoamento fica transcritico


def cota_talvegue(d):
    """Cota da calha junto ao eixo. Sobrevive aos ajustes, porque e por indice."""
    i = d.get("i_thal")
    return float(d["z"][i]) if i is not None else float(np.nanmin(d["z"]))


def mover_calha(d, delta):
    """Move so a calha em 'delta'; o terreno em volta fica onde o DEM o pos."""
    if abs(delta) < 1e-6:
        return
    sta = d["sta"]
    larg = d.get("larg_canal", 150.0)
    centro = sta[d["i_thal"]]
    dist = np.abs(sta - centro)
    meia = larg / 2.0
    talude = max(larg * 0.25, 30.0)
    frac = np.clip(1.0 - (dist - meia) / talude, 0.0, 1.0)
    frac[dist <= meia] = 1.0
    d["z"] = d["z"] + delta * frac


def condicionar(xs, rotulo=""):
    """Talvegue monotonico, com declividade entre DECL_MINIMA e DECL_MAXIMA.

    1. apara a cabeceira enquanto a declividade passar do maximo (torrente de
       montanha: o Benedito chega a 9,4% no DEM bruto);
    2. impoe decrescimento rio abaixo (tira contrapendente);
    3. limita a declividade ancorando a JUSANTE, para preservar a cota da foz
       ou da confluencia, que e o que amarra o trecho a rede.
    """
    if len(xs) < 3:
        return xs
    corte = 0
    while corte < len(xs) - 2:
        dx = xs[corte]["rs"] - xs[corte + 1]["rs"]
        dz = cota_talvegue(xs[corte]) - cota_talvegue(xs[corte + 1])
        if dx > 0 and abs(dz) / dx > DECL_MAXIMA:
            corte += 1
        else:
            break
    if corte:
        print(f"        {rotulo}: aparadas {corte} secoes de cabeceira "
              f"({(xs[0]['rs']-xs[corte]['rs'])/1000:.1f} km acima de "
              f"{100*DECL_MAXIMA:.1f}%)")
        xs = xs[corte:]

    for i in range(1, len(xs)):
        dx = xs[i - 1]["rs"] - xs[i]["rs"]
        teto = cota_talvegue(xs[i - 1]) - DECL_MINIMA * dx
        atual = cota_talvegue(xs[i])
        if atual > teto:
            mover_calha(xs[i], teto - atual)

    for i in range(len(xs) - 2, -1, -1):
        dx = xs[i]["rs"] - xs[i + 1]["rs"]
        lim = cota_talvegue(xs[i + 1]) + DECL_MAXIMA * dx
        atual = cota_talvegue(xs[i])
        if atual > lim:
            mover_calha(xs[i], lim - atual)
    return xs


def ancorar(xs, cota_alvo, degrau=0.5):
    """Desloca o trecho inteiro para casar a foz com o receptor na juncao.

    Aqui o deslocamento uniforme e legitimo: e mudanca de referencia de um
    afluente inteiro, entao o relevo relativo dentro dele nao se altera.
    """
    desl = (cota_alvo + degrau) - cota_talvegue(xs[-1])
    for d in xs:
        d["z"] = d["z"] + desl
    return desl


def manning(xs, n_canal=0.035, razao_planicie=1.8, n_max=0.10):
    """n de Manning por secao, com Jarrett (1984) nas gargantas.

    Corredeira em rocha nao tem a rugosidade de rio de planicie. Usar 0,035 no
    trecho do Salto Pilao (8 m/km) da Froude ~0,9: transcritico, e o solver
    diverge no pico. Jarrett, valido para 0,002 <= S <= 0,052:

        n = 0,39 S^0,38 R^-0,16

    que em S = 0,008 e R ~ 3 m da n ~ 0,052 -- o triplo do valor de planicie.
    Nao e ajuste para estabilizar: e a rugosidade que a literatura mede nesses
    trechos, e ela sozinha derruba o Froude para ~0,5, porque a lamina engrossa
    e a velocidade cai.
    """
    n_ajust = 0
    for i, d in enumerate(xs):
        d["n"] = n_canal
        viz = xs[min(i + 1, len(xs) - 1)]
        dx = d["rs"] - viz["rs"]
        if dx <= 0:
            continue
        S = abs(cota_talvegue(d) - cota_talvegue(viz)) / dx
        if S < 0.002:
            continue
        R = max(d.get("prof_canal", 3.0), 0.5)
        n_j = min(max(0.39 * S ** 0.38 * R ** -0.16, n_canal), n_max)
        if n_j > n_canal + 1e-4:
            d["n"] = round(float(n_j), 4)
            n_ajust += 1
    for d in xs:
        d["n_planicie"] = round(d["n"] * razao_planicie, 3)
    return n_ajust
