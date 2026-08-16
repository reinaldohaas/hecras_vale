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
DECL_MAXIMA = 0.008      # m/m; teto de PARTIDA, valido para rio de planicie
DECL_TETO = 0.05         # m/m; teto absoluto (limite de validade de Jarrett)
ESCAVACAO_MAXIMA = 12.0  # m; o quanto o talvegue pode se afastar do terreno


def teto_declividade(xs):
    """Teto de declividade do PROPRIO rio, tirado do terreno.

    Um teto unico para rio de planicie e rio de serra e a origem do pior
    defeito desta reescrita. O Rio Benedito cai 515,5 -> 53,0 m em 43,7 km:
    1,06% de media, 1,33% de mediana local, 13,6% de maxima. Forcando 0,8%,
    ancorado na foz e subindo, o alvo se afasta do terreno em 128,6 m de
    mediana e 279,2 m no pior ponto -- o modelo cavava um canion ficticio sob
    a serra. Dai vinham o erro de volume astronomico, a falha no assentamento
    (nao ha como estabelecer lamina num canal de 280 m) e a falha mudando de
    rio a cada correcao: cada rio de serra tinha o seu canion.

    O teto sai do percentil 90 da declividade real, com folga. Rio de planicie
    continua em 0,008; rio de serra ganha o que o relevo dele pede, e a
    rugosidade de Jarrett e o refino de espacamento ja acompanham.
    """
    S = np.array([d.get("S_terreno", 0.0) for d in xs], float)
    S = S[np.isfinite(S)]
    if not len(S):
        return DECL_MAXIMA
    return float(np.clip(np.percentile(S, 90) * 1.2, DECL_MAXIMA, DECL_TETO))


def cota_talvegue(d):
    """Cota ALVO do talvegue -- a que o perfil longitudinal definiu.

    Enquanto o condicionamento roda, a secao ainda e o terreno cru: quem carrega
    a decisao e z_alvo, um escalar. A geometria so muda no fim, quando
    secao.escavar() aplica o trapezio uma unica vez. Antes disso cada passo
    reescavava sobre o resultado do anterior e a calha degenerava num pico.
    """
    if "z_alvo" in d:
        return float(d["z_alvo"])
    i = d.get("i_thal")
    base = float(d["z"][i]) if i is not None else float(np.nanmin(d["z"]))
    return base - d.get("prof_canal", 0.0)


def mover_calha(d, delta):
    """Ajusta a cota ALVO do talvegue. Nao toca na geometria.

    Antes isto reaplicava o trapezio da escavacao a cada chamada. Como o
    condicionamento chama tres vezes (monotonia, limite de declividade,
    ancoragem) e cada uma incidia sobre o perfil ja alterado, os deslocamentos
    se acumulavam em pontos diferentes e a calha virava um pico de um ponto so.
    """
    if abs(delta) < 1e-6:
        return
    d["z_alvo"] = cota_talvegue(d) + delta


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
    dmax = teto_declividade(xs)
    if dmax > DECL_MAXIMA * 1.05:
        print(f"        {rotulo}: rio de serra, teto de declividade "
              f"{100*dmax:.2f}% (padrao {100*DECL_MAXIMA:.1f}%)")
    # o alvo parte do terreno menos a profundidade da calha
    for d in xs:
        d.setdefault("z_alvo", float(d["z"][d["i_thal"]]) - d.get("prof_canal", 0.0))
    corte = 0
    while corte < len(xs) - 2:
        dx = xs[corte]["rs"] - xs[corte + 1]["rs"]
        dz = cota_talvegue(xs[corte]) - cota_talvegue(xs[corte + 1])
        if dx > 0 and abs(dz) / dx > dmax:
            corte += 1
        else:
            break
    if corte:
        print(f"        {rotulo}: aparadas {corte} secoes de cabeceira "
              f"({(xs[0]['rs']-xs[corte]['rs'])/1000:.1f} km acima de "
              f"{100*dmax:.1f}%)")
        xs = xs[corte:]

    for i in range(1, len(xs)):
        dx = xs[i - 1]["rs"] - xs[i]["rs"]
        teto = cota_talvegue(xs[i - 1]) - DECL_MINIMA * dx
        atual = cota_talvegue(xs[i])
        if atual > teto:
            mover_calha(xs[i], teto - atual)

    # Limite de declividade e piso de escavacao no MESMO laco. Aplicar o piso
    # depois, como passe separado, sobrescreve z_alvo sem que as secoes
    # vizinhas saibam, e o perfil ganha degraus que ele proprio acabou de
    # remover -- o Itajai-Mirim saiu com 17,5% de degrau local sendo um rio de
    # teto 0,8%. Dentro do laco, cada secao ja enxerga a vizinha corrigida.
    fundo = 0
    for i in range(len(xs) - 2, -1, -1):
        dx = xs[i]["rs"] - xs[i + 1]["rs"]
        lim = cota_talvegue(xs[i + 1]) + dmax * dx
        piso = float(xs[i]["z"][xs[i]["i_thal"]]) - ESCAVACAO_MAXIMA
        atual = cota_talvegue(xs[i])
        novo = min(atual, lim)              # nao mais ingreme que o teto
        if novo < piso:                     # nem mais fundo que o terreno
            novo = piso
            fundo += 1
        if abs(novo - atual) > 1e-9:
            xs[i]["z_alvo"] = novo
    if fundo:
        print(f"        {rotulo}: {fundo} secoes no piso de "
              f"{ESCAVACAO_MAXIMA:.0f} m de escavacao")
    return xs


def ancorar(xs, cota_alvo, degrau=0.5):
    """Casa a foz do afluente com o leito do receptor, na juncao.

    Eu tinha escrito aqui que o deslocamento uniforme da secao inteira era
    legitimo, "porque o relevo relativo dentro do afluente nao se altera".
    Preserva o relativo e afunda o ABSOLUTO -- e num rio que termina no mar
    isso poe a planicie abaixo de zero. No Itajai-Mirim o deslocamento era de
    -6,0 m e as secoes dos ultimos 12 km saiam assim:

        RS 12.137   leito -8,85   topo -0,96
        RS  8.137   leito -9,25   topo -2,22

    O ponto MAIS ALTO da secao a -2 m: nao havia cota acima do mar em 1.474 m
    de largura. Qualquer lamina extrapola a tabela de conducao, e o HEC-RAS
    falhava ja no assentamento (200 passos de warm-up), antes do primeiro passo
    de tempo -- ele nao conseguia sequer estabelecer a lamina inicial.

    Duas mudancas:
      - so a CALHA se move, como no condicionamento. A planicie fica onde o
        DEM a pos.
      - o ajuste e ATENUADO ao longo do trecho: inteiro na foz, nulo na
        cabeceira. O desencontro esta na confluencia; nao ha razao para
        propaga-lo 114 km rio acima.
    """
    for d in xs:                       # garante z_alvo em todas
        d.setdefault("z_alvo", cota_talvegue(d))
    desl = (cota_alvo + degrau) - cota_talvegue(xs[-1])
    if abs(desl) < 1e-6 or len(xs) < 2:
        return desl
    rs = np.array([d["rs"] for d in xs], float)
    faixa = max(rs.max() - rs.min(), 1.0)
    for d in xs:
        mover_calha(d, desl * (rs.max() - d["rs"]) / faixa)
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
        # a declividade do TERRENO, nao a do perfil ja condicionado. Tres
        # afluentes de serra (dos Cedros, Taio, Benedito) saem do
        # condicionamento com a declividade cravada em DECL_MAXIMA: derivar n
        # dali da o mesmo valor para todos e subestima a rugosidade justamente
        # onde o escoamento e mais rapido. O terreno guarda a queda real.
        S = d.get("S_terreno")
        if S is None:
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
