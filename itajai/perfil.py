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
ESCAVACAO_MAXIMA = 12.0   # m; o quanto o talvegue pode se afastar do terreno
CORTE_MAX_FRACAO = 0.35   # ate quanto do trecho a cabeceira pode ser aparada
JANELA_CORTE = 3          # secoes na janela que mede a queda da cabeceira
SECOES_ACIMA_JUNCAO = 4   # secoes que sempre ficam acima da 1a confluencia


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


def condicionar(xs, rotulo="", rs_limite=None):
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
    # Aparar a cabeceira de torrente. O criterio antigo parava no PRIMEIRO par
    # abaixo do teto, entao um patamar isolado tres secoes abaixo do topo
    # interrompia o corte e deixava o resto da torrente dentro do modelo -- que
    # e de onde vinham os degraus de 10,75% no Benedito e 10,37% no dos Cedros,
    # e os erros de 20 a 25 m que o solver reportava no Taio e no alto Acu.
    #
    # Agora procura-se a secao mais a JUSANTE que ainda esta acima do teto,
    # medindo a queda numa janela de algumas secoes para que um par plano no
    # meio da descida nao mascare o conjunto, e corta-se tudo acima dela.
    #
    # Cortar nao perde agua: a area de drenagem do trecho aparado continua
    # entrando como vazao lateral, distribuida nas secoes que sobram. O que se
    # descarta e a pretensao de representar como rio 1D um trecho que desce 10%
    # -- ali o escoamento nao e gradualmente variado, e o modelo so tem a
    # perder tentando.
    # Nunca aparar ate engolir a confluencia mais de montante. Aparando por
    # fracao apenas, a cabeceira do Itajai_Norte acima da foz do Iraputa
    # desapareceu inteira, e a juncao ficou com UM trecho entrando e um saindo
    # -- o HEC-RAS recusa a geometria antes de computar ("Junctions are for
    # flow confluences and splits"), e a rodada morre sem sequer gerar o HDF.
    rs_v = np.array([x["rs"] for x in xs], float)
    zt_v = np.array([cota_talvegue(x) for x in xs], float)
    jan = JANELA_CORTE
    corte = 0
    lim_corte = max(int(CORTE_MAX_FRACAO * len(xs)), 0)
    if rs_limite is not None:
        acima = int(np.searchsorted(-rs_v, -float(rs_limite)))
        lim_corte = min(lim_corte, max(acima - SECOES_ACIMA_JUNCAO, 0))
    for i in range(min(lim_corte, len(xs) - jan - 1)):
        dx = rs_v[i] - rs_v[i + jan]
        if dx > 0 and (zt_v[i] - zt_v[i + jan]) / dx > dmax:
            corte = i + 1
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
    alisar_perfil(xs, dmax, rotulo=rotulo)
    return xs


def alisar_perfil(xs, dmax, n_iter=400, rotulo=""):
    """Distribui a queda ao longo do trecho, em vez de empilhar degraus.

    Os dois limites de condicionamento sao rigidos: um laco impoe a
    declividade MINIMA indo rio abaixo, outro a MAXIMA voltando. Entre os dois
    nao sobra valor intermediario, e o perfil sai alternando exatamente os
    extremos. No Itajai-Mirim, com secoes de 150 m:

        RS 110274,7   dz = -2,40      (= dmax, 1,6%)
        RS 110124,7   dz = -0,02      (= dmin, 0,01%)
        RS 109974,7   dz = -0,02
        RS 109824,7   dz = -2,02
        RS 109674,7   dz = -2,41

    Isso e uma escada de pocos e quedas: dezenas de degraus seguidos, cada um
    um salto transcritico. O solver nao permanente nao converge nisso -- as
    iteracoes batiam o teto de 40 em TODA linha desde o aquecimento, com as
    vazoes ainda de base, e o aborto vinha justamente aqui
    (Itajai_Mirim R1 RS 104738,3).

    A queda total do trecho e dada pelo terreno e nao muda; o que muda e como
    ela se reparte. Alisar e reimpor os limites, alternadamente, espalha a
    mesma queda de forma continua e so encosta nos limites onde e inevitavel.

    Os dois limites viram acumulacoes de minimo, com a distancia rio abaixo
    como variavel: para a declividade minima, z + dmin*t tem de ser nao
    crescente; para a maxima, z + dmax*t tem de ser nao decrescente.
    """
    if len(xs) < 5:
        return
    rs = np.array([x["rs"] for x in xs], float)
    t = rs[0] - rs                     # distancia rio abaixo, crescente
    z = np.array([cota_talvegue(x) for x in xs], float)
    terreno = np.array([float(x["z"][x["i_thal"]]) for x in xs], float)
    piso = terreno - ESCAVACAO_MAXIMA
    z_foz = z[-1]                      # amarra o trecho a rede: nao se mexe

    for _ in range(n_iter):
        m = z.copy()
        m[1:-1] = 0.25 * z[:-2] + 0.5 * z[1:-1] + 0.25 * z[2:]
        z = np.clip(m, piso, terreno)
        z[-1] = z_foz
        v = np.minimum.accumulate(z + DECL_MINIMA * t)      # declividade minima
        z = v - DECL_MINIMA * t
        u = np.minimum.accumulate((z + dmax * t)[::-1])[::-1]   # maxima
        z = u - dmax * t
        z = np.maximum(z, piso)
        z[-1] = z_foz
    # ------ e agora no dominio da DECLIVIDADE, nao da cota.
    # Limitar a declividade entre um minimo e um maximo nao impede o JOELHO: o
    # perfil pode respeitar os dois limites e mesmo assim passar de 0,07% para
    # 0,70% num vao de 150 m. Foi o que sobrou no Itajai-Acu na garganta do
    # Salto Pilao -- remanso quase plano ate RS 92.115, depois 0,70% em uma
    # secao, plano de novo, e 0,86% sustentado de RS 86.352 a 84.852. Fator de
    # doze, com as transicoes numa secao so.
    #
    # Hidraulicamente isso e o escoamento saindo de remanso profundo e lento
    # para quase-critico e voltando, dezenas de vezes seguidas, que e onde um
    # solver 1D nao permanente nao converge. O log do HEC-RAS reclama
    # exatamente ai (Itajai_Acu R3 82.796 a 86.883 e R2 91.117 a 105.042).
    #
    # A queda total do trecho nao muda: suaviza-se a declividade com media
    # movel ponderada pelo vao e reintegra-se a partir da FOZ, que e a cota
    # amarrada a rede. O desnivel se redistribui, a garganta continua sendo
    # garganta, mas a transicao para ela passa a ocupar varias secoes.
    dt = np.diff(t)
    ok = dt > 0
    if ok.sum() >= 5:
        decl = np.zeros_like(dt)
        decl[ok] = (z[:-1] - z[1:])[ok] / dt[ok]
        jan = 7
        k = np.ones(jan) / jan
        peso = np.convolve(np.pad(dt, jan // 2, mode="edge"), k, "valid")
        suave = np.convolve(np.pad(decl * dt, jan // 2, mode="edge"), k, "valid")
        decl = np.divide(suave, peso, out=decl.copy(), where=peso > 0)
        decl = np.clip(decl, DECL_MINIMA, dmax)
        zn = np.empty_like(z)
        zn[-1] = z_foz
        for i in range(len(z) - 2, -1, -1):
            zn[i] = zn[i + 1] + decl[i] * dt[i]
        z = np.minimum(zn, terreno)          # nunca acima do terreno
        z = np.maximum(z, piso)
        z[-1] = z_foz
        z = np.maximum(z, piso)
        z[-1] = z_foz

    # MONOTONIA POR ULTIMO, e sem reaplicar o piso depois dela.
    # O piso de escavacao (12 m abaixo do terreno) e uma restricao boa, mas ela
    # LEVANTA o leito onde o DEM tem um alto local -- e o Copernicus e modelo
    # de SUPERFICIE, entao ponte, mata densa e soleira aparecem como alto. Ao
    # ser aplicado depois da monotonia, o piso desfazia o que ela garantira:
    # no Itajai_Sul R1 RS 41.340 o leito ficou em 375,09 m entre vizinhas em
    # 366,18 e 365,36 -- um corcovo de 9 m em 300 m, isto e, uma barragem
    # dentro do modelo, com o solver tendo de empurrar agua ladeira acima. Era
    # a maior fonte de erro do run (4,83 m em RS 41.190).
    #
    # Entre respeitar o piso e nao ter contrapendente, nao ter contrapendente
    # ganha: escavar 15 m num ponto e uma aproximacao; um degrau ao contrario
    # e um erro de fisica.
    v = np.minimum.accumulate(z + DECL_MINIMA * t)
    z = v - DECL_MINIMA * t
    z[-1] = z_foz
    fundo_extra = int(np.sum(z < piso - 0.01))
    if fundo_extra:
        print(f"        {rotulo}: {fundo_extra} secoes abaixo do piso de "
              f"escavacao, para nao criar contrapendente")

    for x, zi in zip(xs, z):
        x["z_alvo"] = float(zi)


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
