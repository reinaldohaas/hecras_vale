# -*- coding: utf-8 -*-
"""
Onde cortar as secoes transversais, e com que largura.

Este modulo produz secoes de TERRENO puro. A calha e escavada depois
(vale/calha.py), uma unica vez, com o perfil longitudinal ja resolvido
(vale/perfil.py). A ordem importa: escavando aqui e reajustando a cada passo do
condicionamento, o trapezio e reaplicado sobre um perfil ja alterado e o que
sobra e um pico isolado -- uma fenda de 3,4 m num rio de calha de 106 m, que
nao conduz nada.

Quatro decisoes, cada uma paga com uma rodada de depuracao:

ESPACAMENTO PELA DECLIVIDADE. O Acu cai 195 m em 13 km na garganta do Salto
Pilao. A 1 km de espacamento sao 8 m de queda entre secoes vizinhas e o solver
falha no primeiro passo; o criterio dx <~ 0,15*D/S da ~75 m ali.

LIMITE DE CURVATURA. Numa curva de raio R, duas perpendiculares vizinhas se
encontram a R do eixo, do lado CONCAVO. Passar disso cruza as cutlines -- 24%
dos pares na primeira tentativa -- e a interpolacao da lamina entre elas deixa
de ter significado.

JANELA DESLOCADA, e nao so aparada. O limite de curvatura corta um lado so, e o
alargamento estica o outro atras de encosta: a calha acaba na borda do corte.
Houve secao com Bank Sta 511/617 num corte que terminava em 620 -- tres metros
de planicie de um lado. Deslocar a janela preserva a largura total; aparar so o
lado largo e o recurso quando a curvatura nao deixa deslocar.

ESPACAMENTO UNIFORME DOS PONTOS. Repartir N pontos ao meio entre dois lados de
larguras diferentes da 1,4 m de espacamento de um lado e 5,0 m do outro. A
conducao fica enviesada para o lado denso, e o proprio talvegue e procurado num
indice que nao corresponde ao eixo.
"""
import numpy as np


def tirar_picos(z, janela=5, limite=3.0):
    """Filtro de mediana no talvegue lido do terreno.

    Monotonia entre pares vizinhos nao pega pico isolado: um valor 19 m acima
    dos dois vizinhos passa nos dois testes de par. A mediana movel pega.
    """
    z = np.asarray(z, float)
    if len(z) < janela:
        return z.copy()
    k = janela // 2
    pad = np.pad(z, k, mode="edge")
    med = np.array([np.median(pad[i:i + janela]) for i in range(len(z))])
    fora = np.abs(z - med) > limite
    saida = z.copy()
    saida[fora] = med[fora]
    return saida


def direcao(linha, s, suaviza):
    """Direcao local do eixo, medida numa janela de +-suaviza metros.

    Entre dois vertices consecutivos a direcao oscila com o serrilhado da
    digitalizacao da ANA, e as perpendiculares saem tortas -- o RAS acusa
    "edge lines have self intersections".
    """
    a = linha.interpolate(max(0.0, s - suaviza))
    b = linha.interpolate(min(linha.length, s + suaviza))
    tx, ty = b.x - a.x, b.y - a.y
    n = float(np.hypot(tx, ty)) or 1.0
    return tx / n, ty / n


def estacas(linha, amostrador, op, area_foz=None):
    """Posicoes de corte, adensadas onde o leito e ingreme."""
    L = linha.length
    passo = min(op.espacamento / 4.0, 250.0)
    d = np.arange(0.0, L + passo, passo)
    P = [linha.interpolate(float(x)) for x in d]
    zb = amostrador.cota([p.x for p in P], [p.y for p in P])
    ok = np.isfinite(zb)
    if ok.sum() < 3:
        return list(np.arange(0.0, L, op.espacamento))
    zb = tirar_picos(np.interp(d, d[ok], zb[ok]))
    S = np.abs(np.gradient(zb, d))
    S = np.convolve(S, np.ones(5) / 5.0, "same")
    f = np.clip((S - op.decl_plano) / (op.decl_ingreme - op.decl_plano), 0.0, 1.0)
    dx = op.espacamento + (op.espacamento_min - op.espacamento) * f

    # TETO DE SAMUELS (1989), dx <= k*D/S. A interpolacao acima satura em
    # `espacamento_min` e para de responder: de 0,6% para cima ela devolve
    # 150 m tanto para 0,6% quanto para 5%, e o criterio pede 37 m e 4,5 m.
    # Era o mesmo espacamento para uma encosta oito vezes mais ingreme.
    if getattr(op, "samuels", False):
        # D E A PROFUNDIDADE DE CALHA CHEIA, e ela nao e a mesma no Itajai-Acu
        # e num ribeirao de cabeceira. Com D=1,5 m fixo o criterio ficava 5,3
        # vezes severo demais no Acu (D real ~8 m) e o modelo saiu com uma
        # secao a cada 100 m ou menos em rio de planicie -- 11.684 secoes nos
        # 12 rios, contra ~5.000 com o D certo. D fixo tambem nao respondia a
        # rio nenhum: o mesmo numero para 14.871 km2 e para 240 km2.
        #
        # Leopold, h = kh*A^eh, com a MESMA area de drenagem que o corte usa
        # (rampa da cabeceira ate a foz). Sem area, cai no valor fixo.
        D = np.full_like(S, float(op.samuels_D))
        if area_foz and getattr(op, "samuels_leopold", True):
            a_cab = max(float(area_foz) * op.fracao_cabeceira, 1.0)
            a = a_cab + (float(area_foz) - a_cab) * (d / max(L, 1.0))
            D = op.canal_kh * np.maximum(a, 1.0) ** op.canal_eh
        lim = op.samuels_k * D / np.maximum(S, 1e-6)
        dx = np.maximum(np.minimum(dx, lim), op.espacamento_piso)
    # nao encosta nos extremos: secao em cima da juncao conflita com o
    # comprimento declarado em Junc L&A e trava o solver
    recuo = op.espacamento_min * 0.5
    piso = float(getattr(op, "espacamento_piso", 25.0))
    ss, s = [recuo], recuo
    while s < L - recuo:
        s += _passo(s, d, dx, piso)
        if s < L - recuo:
            ss.append(s)
    # trecho final minusculo e pior que trecho final ausente: dois cortes a
    # poucos metros um do outro dao comprimento de trecho quase zero, e o
    # solver divide por ele.
    if ss and (L - recuo) - ss[-1] < 0.5 * piso:
        ss.pop()
    ss.append(L - recuo)
    return ss


def _minimo_no_intervalo(a, b, d, dx):
    """Menor dx exigido em QUALQUER ponto de [a, b], e nao so nas pontas."""
    i0 = int(np.searchsorted(d, a, "right"))
    i1 = int(np.searchsorted(d, b, "left"))
    v = [float(np.interp(a, d, dx)), float(np.interp(b, d, dx))]
    if i1 > i0:
        v.append(float(dx[i0:i1].min()))
    return min(v)


def _passo(s, d, dx, piso):
    """Maior passo que respeita dx em TODO o intervalo que ele cobre.

    O laco antigo fazia `s += interp(s, d, dx)`: media o espacamento exigido no
    ponto de partida e depois andava as cegas ate o fim do passo. Onde a
    declividade cresce DENTRO do passo, o corte seguinte cai muito alem do que
    o criterio permitia -- e como o criterio so e reavaliado no proximo ponto
    de partida, o trecho grande ja estava criado.

    Nao e detalhe: no Itajai do Norte deixou 420 dos 2.069 trechos acima do
    limite E acima do piso, o pior com 79 m onde o criterio pedia 4 m. Ou seja,
    a regra estava escrita e nao valia justamente nas transicoes -- exatamente
    onde ela importa, porque e ali que a declividade muda depressa.

    Ponto fixo em vez de tentativa unica: encurtar o passo pode trazer para
    dentro dele um trecho ainda mais ingreme. Converge porque `h` so diminui e
    tem o piso por baixo.
    """
    h = float(np.interp(s, d, dx))
    for _ in range(12):
        if h <= piso + 1e-9:
            return piso
        m = _minimo_no_intervalo(s, s + h, d, dx)
        if m >= h - 1e-6:
            break
        h = max(m, piso)
    return h


def desnivel(d):
    """Desnivel do TERRENO dentro da secao, em metros.

    E o que decide se a secao alcanca as encostas do vale. Fracao plana nao
    serve para isso: uma secao pode ter 40% de pontos numa cota so e ainda
    assim ter 60 m de parede nas pontas.
    """
    z = np.asarray(d["z"], float)
    z = z[np.isfinite(z)]
    return float(z.max() - z.min()) if z.size else 0.0


def meias_para_cheia(d, op):
    """Meia-largura de cada lado que a CHEIA precisa, e nao a que a area sugere.

    O criterio antigo era so o porte do rio -- 180*sqrt(A/100), com piso de
    500 m de meia-largura --, e o piso e que mandava: 129 das 148 secoes do
    Benedito estavam nele. Medido, a cheia de pico molhava 13% da largura na
    mediana e 6% na cabeceira: 62 m de agua dentro de uma secao de 966 m.

    Os 87% restantes nao ficam so sobrando, eles ESTRAGAM. A secao sobe as duas
    encostas do vale e vira bacia fechada -- 40 m de profundidade mediana,
    206 m no pior caso --, e numa bacia a conducao depende brutalmente da cota:
    entre duas vizinhas com fundos diferentes ela variava por um fator de
    2.809, o momento nao fechava e a vazao invertia de sinal no primeiro passo.

    Aqui a seccao vai ate onde o terreno passa da cota de cheia mais uma folga,
    com margem, limitada abaixo por meia_largura_min. Onde o vale e largo e a
    cheia espalha de verdade (a foz do Benedito usa 57% da secao) nao encolhe.
    """
    from .calha import altura_para_vazao, vazao_projeto
    sta = np.asarray(d["sta"], float)
    z = np.asarray(d["z"], float)
    i = int(d["i_thal"])
    h = altura_para_vazao(sta, z, 0.05, d.get("S_terreno"),
                          vazao_projeto(d["area_km2"]))
    alvo = float(z[i]) + h + op.folga_secao
    # ate onde a agua chega ANDANDO a partir do talvegue: a primeira subida
    # acima da cota fecha o lado. Procurar o minimo global pegaria uma
    # depressao do outro lado do morro.
    e = i
    while e > 0 and z[e - 1] <= alvo:
        e -= 1
    r = i
    while r < len(z) - 1 and z[r + 1] <= alvo:
        r += 1
    me = max((sta[i] - sta[e]) * op.margem_secao, op.meia_largura_min)
    md = max((sta[r] - sta[i]) * op.margem_secao, op.meia_largura_min)
    return float(me), float(md)


def largura_base(area_km2):
    """Meia-largura conforme o porte do rio.

    Alargar isto nao e de graca: com secao larga demais o corte atravessa
    meandros do proprio rio e o escoamento e contado duas vezes (subir o
    coeficiente de 180 para 440 derrubou a simulacao de 30 passos para 2).
    """
    return float(np.clip(180.0 * np.sqrt(max(area_km2, 1.0) / 100.0),
                         500.0, 2500.0))


def limites_curvatura(linha, ss, meia, op):
    """Ate onde cada semi-secao pode ir sem cruzar as vizinhas."""
    n = len(ss)
    dirs = [direcao(linha, s, op.janela_direcao) for s in ss]
    ang = np.unwrap([np.arctan2(t[1], t[0]) for t in dirs])
    esq = (np.full(n, float(meia)) if np.isscalar(meia)
           else np.asarray(meia, float).copy())
    dir_ = esq.copy()
    for i in range(n):
        for j in (i - 1, i + 1):
            if j < 0 or j >= n:
                continue
            ds = abs(ss[j] - ss[i])
            dth = abs(ang[j] - ang[i])
            if ds <= 0 or dth < 1e-6:
                continue
            R = op.folga_curva * ds / dth
            if (ang[j] - ang[i]) * (1 if j > i else -1) > 0:
                esq[i] = min(esq[i], R)
            else:
                dir_[i] = min(dir_[i], R)

    # Alisar ao longo do trecho. Sem isto o limite apara a secao i e nao a i+1,
    # e a largura sai em dente de serra (1000, 941, 746, 737, 1000...). O que
    # derruba o solver nao e a largura, e o SALTO entre vizinhas: area e
    # conducao mudam de degrau a cada secao.
    def alisar(v, jan=5):
        if len(v) < jan:
            return v
        m = np.array([v[max(0, i - jan // 2):i + jan // 2 + 1].min()
                      for i in range(len(v))])
        return np.convolve(np.pad(m, 1, mode="edge"), np.ones(3) / 3.0, "valid")

    e = np.maximum(alisar(esq), op.minimo_lado)
    d = np.maximum(alisar(dir_), op.minimo_lado)
    return e, d


def equilibrar(hw_e, hw_d, cap_e, cap_d, op):
    """Desloca a janela do corte para a calha nao ficar no canto.

    O que falta de um lado sai do outro, ate o limite de curvatura -- e la que
    as cutlines se cruzam. Se ainda sobrar desequilibrio, apara-se o lado
    largo, porque estreitar nunca cria cruzamento.
    """
    minlado = (hw_e + hw_d) / (1.0 + op.razao_lados)
    for falta, pode, lado in ((minlado - hw_d, cap_d - hw_d, "d"),
                              (minlado - hw_e, cap_e - hw_e, "e")):
        outro = hw_e if lado == "d" else hw_d
        mv = np.clip(np.minimum(np.minimum(falta, pode),
                                outro - op.minimo_lado), 0.0, None)
        if lado == "d":
            hw_d, hw_e = hw_d + mv, hw_e - mv
        else:
            hw_e, hw_d = hw_e + mv, hw_d - mv
    return (np.minimum(hw_e, hw_d * op.razao_lados),
            np.minimum(hw_d, hw_e * op.razao_lados))


def cortar(linha, s, amostrador, area_km2, he, hd, op):
    """Uma secao perpendicular ao eixo, so de terreno."""
    tx, ty = direcao(linha, s, op.janela_direcao)
    rx, ry = ty, -tx                                # normal, para a direita
    p = linha.interpolate(s)
    he, hd = float(he), float(hd)

    # NUMERO DE PONTOS PELO TAMANHO DA SECAO, e nao fixo. Com 280 pontos numa
    # secao de 175 m -- que e o que o recorte pela cota de cheia passou a
    # produzir -- o espacamento cai para 0,62 m, mais fino que o pixel de 30 m
    # do terreno: sao pontos inventados por interpolacao. E foi assim que o
    # HEC-RAS recusou a geometria: a estaca da margem calhou exatamente sobre
    # uma amostra, o construtor inseriu o ponto de novo e sobrou uma duplicata
    # ("Station and elevation data contains duplicate points").
    n_pts = int(np.clip(round((he + hd) / op.espacamento_pontos),
                        op.n_pontos_min, op.n_pontos))
    n_e = int(round(n_pts * he / max(he + hd, 1e-6)))
    n_e = min(max(n_e, 2), n_pts - 2)
    off = np.concatenate([np.linspace(-he, 0.0, n_e, endpoint=False),
                          np.linspace(0.0, hd, n_pts - n_e)])
    i_eixo = n_e
    z = amostrador.cota(p.x + off * rx, p.y + off * ry)
    if not np.isfinite(z).any():
        return None
    if not np.isfinite(z).all():
        ok = np.isfinite(z)
        z = np.interp(np.arange(len(z)), np.flatnonzero(ok), z[ok])

    sta = off + he
    # talvegue PROXIMO AO EIXO, nao o minimo global: o corte pode atravessar
    # outro canal mais fundo (o leito antigo do Mirim, um meandro do proprio
    # rio) e a calha iria parar la, com as bank lines em estrela.
    janela = max(op.canal_kw * max(area_km2, 1.0) ** op.canal_ew, 150.0)
    # E A JANELA NUNCA ALCANCA A BORDA. O piso de 150 m foi escrito pensando em
    # secao larga, mas o recorte pela cota de cheia deixa secoes de 120 m de
    # largura TOTAL -- e ai a janela cobre a secao inteira, a restricao ao eixo
    # deixa de existir e vence o minimo global, que pode estar na ponta.
    # Medido nos 12 rios: 3% das secoes com o talvegue na BORDA, e e onde a
    # calha e cavada -- o canal ia parar fora do rio. No Trombudo RS 39925 o
    # i_thal era 0, o primeiro ponto do corte.
    janela = min(janela, 0.35 * float(sta[-1] - sta[0]))
    m = np.abs(sta - sta[i_eixo]) <= janela
    idx = np.flatnonzero(m)
    i0 = int(idx[np.argmin(z[idx])]) if len(idx) else i_eixo
    # UMA calha por secao: o que estiver mais fundo que o talvegue do eixo e
    # outro canal que o corte cruzou. Num 1D isso nao e um segundo caminho de
    # escoamento, e um poco dentro da mesma secao.
    z = np.maximum(z, z[i0])

    return {"sta": sta, "z": z, "i_thal": i0,
            "cut": (p.x - he * rx, p.y - he * ry, p.x + hd * rx, p.y + hd * ry),
            "area_km2": float(area_km2), "s": float(s),
            "z_terreno": float(z[i0])}


def ajustar_ao_talvegue(linha, amostrador, op, log=print):
    """Puxa o eixo para o fundo do vale, sem mexer nas pontas.

    O eixo vem da BHO 2017 da ANA, que e um tracado ESQUEMATICO: no Mirim ele
    fica a 16 m do talvegue lido do terreno na mediana, 42 m no p90 e 296 m no
    pior caso -- contra uma meia-calha de 26 m. Em 28% das secoes o rio real
    esta inteiro fora da calha declarada, e em 12% as duas margens caem do
    MESMO lado do ponto do eixo.

    A consequencia nao aparece nas secoes, aparece nas linhas que o HEC-RAS
    deriva delas: a bank line liga os pontos de margem, segue o rio de verdade
    e por isso cruza o eixo toda vez que o rio troca de lado -- 320 vezes no
    Mirim. Com o eixo em cima do talvegue, secao, calha e bank line passam a
    concordar por construcao, em vez de por acaso.

    AS PONTAS NAO SE MOVEM. Elas sao a conexao com a juncao, e o snapping nas
    confluencias e exato (0,0 m); deslocar a foz de um afluente desliga o
    afluente da rede. A correcao e afunilada ate zero nos `eixo_taper` metros
    de cada extremidade.

    O deslocamento tambem e limitado e alisado: sem isso o eixo pula para um
    meandro abandonado ou para uma vala de drenagem vizinha, que no MDS sao
    tao fundos quanto o rio, e sai em dente de serra.
    """
    L = float(linha.length)
    passo = float(getattr(op, "eixo_passo", 50.0))
    jan = float(getattr(op, "eixo_janela", 80.0))
    res = float(getattr(op, "eixo_res", 5.0))
    dmax = float(getattr(op, "eixo_desloc_max", 120.0))
    taper = float(getattr(op, "eixo_taper", 300.0))
    if L < 4 * taper or passo <= 0:
        return linha, None

    # AS ESTACAS TEM DE INCLUIR O FIM. Com `arange(0, L, passo)` a ultima cai
    # antes de L e a polilinha nova TERMINA ali: a foz andava 45 m mesmo com o
    # afunilamento correto, porque o ponto final simplesmente nao existia.
    ss = np.append(np.arange(0.0, L, passo), L)
    offs = np.arange(-jan, jan + res * 0.5, res)
    px, py, nx, ny = [], [], [], []
    for s in ss:
        tx, ty = direcao(linha, float(s), op.janela_direcao)
        p = linha.interpolate(float(s))
        px.append(p.x); py.append(p.y)
        nx.append(ty); ny.append(-tx)          # normal, para a direita
    px = np.array(px); py = np.array(py)
    nx = np.array(nx); ny = np.array(ny)

    X = px[:, None] + offs[None, :] * nx[:, None]
    Y = py[:, None] + offs[None, :] * ny[:, None]
    Z = np.asarray(amostrador.cota(X.ravel(), Y.ravel()), float).reshape(X.shape)

    # PENALIDADE POR DISTANCIA, para nao trocar de rio. Um canal vizinho 5 cm
    # mais fundo puxaria o eixo 80 m se o criterio fosse so o minimo. A
    # penalidade e pequena perto do centro e cresce nas bordas da janela.
    pen = float(getattr(op, "eixo_penalidade", 0.02)) * np.abs(offs)
    custo = Z + pen[None, :]
    custo[~np.isfinite(Z)] = np.inf
    bom = np.isfinite(custo).any(axis=1)
    d = np.zeros(len(ss))
    d[bom] = offs[np.argmin(custo[bom], axis=1)]

    # alisar e limitar
    janela_s = max(3, int(round(float(getattr(op, "eixo_alisar", 250.0)) / passo)))
    if janela_s % 2 == 0:
        janela_s += 1
    d = np.convolve(np.pad(d, janela_s // 2, mode="edge"),
                    np.ones(janela_s) / janela_s, "valid")
    d = np.clip(d, -dmax, dmax)
    # afunilar nas pontas: a juncao nao pode se mover
    f = np.clip(np.minimum(ss, L - ss) / max(taper, 1e-9), 0.0, 1.0)
    d = d * f

    from shapely.geometry import LineString
    novo = LineString(np.c_[px + d * nx, py + d * ny])
    return novo, d


def cortar_rio(eixo, amostrador, op, log=print, prog=None):
    """Todas as secoes de um rio, da cabeceira para a foz."""
    linha = eixo["linha"]
    if getattr(op, "eixo_talvegue", False):
        novo, d = ajustar_ao_talvegue(linha, amostrador, op, log)
        if d is not None:
            a = np.abs(d)
            log(f"   {eixo['ras']}: eixo puxado para o talvegue "
                f"(deslocamento mediano {np.median(a):.1f} m, p90 "
                f"{np.percentile(a, 90):.1f} m, maximo {a.max():.1f} m; "
                f"comprimento {linha.length/1000:.1f} -> {novo.length/1000:.1f} km)")
            linha = novo
            eixo["linha"] = novo          # densificar() usa o MESMO eixo
    L = linha.length
    ss = estacas(linha, amostrador, op, eixo.get("area"))

    # declividade do terreno NAO condicionado, guardada por secao: e ela que
    # alimenta o Manning de Jarrett. Derivar n do perfil ja condicionado usa a
    # declividade do CLAMP e subestima a rugosidade justamente nos afluentes de
    # serra, que sao os que saturam o clamp.
    P = [linha.interpolate(float(s)) for s in ss]
    zt = amostrador.cota([p.x for p in P], [p.y for p in P])
    ok = np.isfinite(zt)
    if ok.sum() >= 2:
        zt = tirar_picos(np.interp(np.arange(len(zt)), np.flatnonzero(ok), zt[ok]))
        S_terr = np.convolve(np.abs(np.gradient(zt, np.asarray(ss, float))),
                             np.ones(3) / 3.0, "same")
    else:
        S_terr = np.zeros(len(ss))

    a_foz = eixo["area"]
    a_cab = max(a_foz * op.fracao_cabeceira, 1.0)
    areas = [a_cab + (a_foz - a_cab) * (s / max(L, 1.0)) for s in ss]
    hw = [largura_base(a) for a in areas]

    hw_e, hw_d = limites_curvatura(linha, ss, hw, op)
    cap_e, cap_d = limites_curvatura(linha, ss, np.maximum(hw_e, hw_d), op)
    hw_e, hw_d = equilibrar(hw_e, hw_d, cap_e, cap_d, op)

    xs = []
    for i, s in enumerate(ss):
        r = cortar(linha, s, amostrador, areas[i], hw_e[i], hw_d[i], op)
        if prog is not None:
            prog.passo(extra=f"{eixo['ras']} RS {L - s:,.0f}")
        if r is None:
            continue
        r["rs"] = round(L - s, 2)          # RS decresce para jusante
        r["S_terreno"] = float(S_terr[i])
        r["rio"] = eixo["ras"]
        r["_i"] = i
        xs.append(r)
    xs.sort(key=lambda d: -d["rs"])

    # ------------------------------------------------------ secao sem desnivel
    # Uma secao inteiramente dentro do fundo plano do vale nao contem cheia
    # nenhuma: no Itajai do Oeste RS 83.278 o terreno tinha 0,0 m de desnivel
    # em 575 m. Nao e a calha que a achatou -- e o vale que e plano ali --, mas
    # uma secao sem parede nao serve para conduzir, e a largura estava presa em
    # 288 m pelo limite de curvatura do meandro, nao pelo porte do rio.
    #
    # SO NAS PLANAS, e nunca em todas. Medido: ampliar 2x TODAS as secoes leva o
    # Acu de 0 para 213 cutlines cruzadas e o Oeste de 8 para 160 -- e o RAS ja
    # avisa de auto-interseccao com a largura atual. Ampliando so as planas o
    # custo e de UM cruzamento por rio.
    # enumerate, e nao xs.index(d): a secao e um dict com arrays numpy dentro,
    # entao list.index compara com == e o numpy devolve um array de booleanos
    # -- "truth value of an array with more than one element is ambiguous".
    largos = 0
    for k, d in enumerate(xs):
        if d.get("_i") is None or desnivel(d) > op.desnivel_minimo:
            continue
        i = d["_i"]
        r = cortar(linha, ss[i], amostrador, areas[i],
                   hw_e[i] * op.fator_alargar, hw_d[i] * op.fator_alargar, op)
        if r is None or desnivel(r) <= desnivel(d):
            continue                        # nao melhorou: fica como estava
        r.update({"rs": d["rs"], "S_terreno": d["S_terreno"],
                  "rio": d["rio"], "_i": i})
        xs[k] = r
        largos += 1
    if largos:
        log(f"   {eixo['ras']}: {largos} secoes sem desnivel alargadas "
            f"{op.fator_alargar:.0f}x (vale plano, largura presa pelo meandro)")

    # ------------------------------------------------- recorte pela cheia
    # SEGUNDA PASSADA, e nao um criterio novo na primeira: a cota de cheia so
    # se conhece depois de ter a secao. Corta-se generoso, mede-se, recorta-se.
    if op.recortar_secao:
        antes = np.array([float(d["sta"][-1]) for d in xs])
        n_rec = 0
        for k, d in enumerate(xs):
            if d.get("_i") is None:
                continue
            i = d["_i"]
            me, md = meias_para_cheia(d, op)
            # so encolhe, e so se valer a pena -- recortar 10% nao paga o corte
            if me > hw_e[i] * 0.9 and md > hw_d[i] * 0.9:
                continue
            r = cortar(linha, ss[i], amostrador, areas[i],
                       min(me, hw_e[i]), min(md, hw_d[i]), op)
            if r is None:
                continue
            r.update({"rs": d["rs"], "S_terreno": d["S_terreno"],
                      "rio": d["rio"], "_i": i})
            xs[k] = r
            n_rec += 1
        if n_rec:
            dep = np.array([float(d["sta"][-1]) for d in xs])
            log(f"   {eixo['ras']}: {n_rec} secoes recortadas pela cota de "
                f"cheia (largura mediana {np.median(antes):.0f} -> "
                f"{np.median(dep):.0f} m)")
    for d in xs:
        d.pop("_i", None)
    n_cort = len(xs)
    xs = densificar(xs, op, eixo.get("area"), log, eixo.get("linha"))
    log(f"   {eixo['ras']:<16} {n_cort:>4} secoes cortadas do terreno"
        + (f" + {len(xs)-n_cort} interpoladas = {len(xs)}"
           if len(xs) > n_cort else "")
        + f"   largura {min(d['sta'][-1] for d in xs):.0f} a "
          f"{max(d['sta'][-1] for d in xs):.0f} m")
    return xs


def _suavizar_alvo(alvo, op):
    """Limita o SALTO de espacamento entre vaos vizinhos.

    O criterio de Samuels e `dx <= k*D/S`: aperta onde o rio e ingreme e NAO
    DIZ NADA onde e plano -- com S tendendo a zero o limite tende ao infinito.
    A planicie herda o espacamento maximo enquanto a serra do mesmo rio fica no
    piso, e a razao entre vaos vizinhos chega a 14x.

    Medido no Cedros: 15 vaos acima de 400 m, um de 944 m, todos com
    declividade de 1e-4 a 3e-4, contra mediana de 25 m no resto do rio. E foi
    no ultimo quilometro plano, com vaos de 642 e 722 m, que o par
    Benedito+Cedros instabilizou -- com 2 a 3 CENTIMETROS de lamina numa secao
    de 500 m de largura, onde a conducao fica mal definida. Salto brusco de
    espacamento e fonte conhecida de instabilidade em 1D nao permanente,
    independente da declividade.

    Iterativo porque baixar um vao aperta o vizinho seguinte.
    """
    razao = float(getattr(op, "razao_dx", 2.0))
    if len(alvo) < 3 or razao <= 1.0:
        return alvo
    v = alvo.copy()
    for _ in range(60):
        antes = v.copy()
        viz = np.minimum(np.r_[v[1:], v[-1]], np.r_[v[0], v[:-1]])
        v = np.minimum(v, razao * viz)
        if np.allclose(v, antes):
            break
    return np.maximum(v, float(getattr(op, "espacamento_piso", 25.0)))



def densificar(xs, op, area_foz=None, log=print, eixo=None):
    """Insere secoes INTERPOLADAS onde o criterio numerico pede mais.

    Duas coisas diferentes estavam sendo confundidas numa so:

      GEOMETRIA -- de quantos em quantos metros o terreno precisa ser AMOSTRADO
      para que o vale esteja representado. Num rio de planicie isso e cada
      poucas centenas de metros; amostrar de 50 em 50 nao acrescenta terreno
      nenhum, so repete o mesmo vale.

      RESOLUCAO NUMERICA -- de quantos em quantos metros o solver precisa de um
      no para resolver Saint-Venant (Samuels). Isso pode ser dezenas de metros.

    Cortar do terreno na densidade NUMERICA e o que produzia 1.553 secoes no
    Mirim, e a declividade que exigia isso vinha do Copernicus: modelo de
    SUPERFICIE, com o dossel dentro. Medido no proprio Mirim -- o terreno
    "sobe" rio abaixo em 34% dos trechos, e o p90 da declividade lida e 3,6
    vezes a declividade media real do rio. O criterio respondia ao dossel.

    O jeito do HEC-RAS e o outro: poucas secoes reais e as intermediarias
    INTERPOLADAS. Quem interpola aqui e a biblioteca --
    `GeomCrossSection.interpolate_station_elevation`, com posicao lateral
    normalizada, estacas de margem interpoladas e o limite de 500 pontos.
    Nao ha geometria inventada por nos: a intermediaria e combinacao das duas
    vizinhas medidas.
    """
    if not getattr(op, "interpolar", True) or len(xs) < 2:
        return xs
    try:
        import pandas as pd

        from ras_commander.geom import GeomCrossSection as G
    except ImportError:
        log("      ras-commander indisponivel; sem interpolacao de secoes")
        return xs

    kh, eh = op.canal_kh, op.canal_eh
    piso = float(getattr(op, "espacamento_piso", 25.0))
    # ESPACAMENTO ALVO DE CADA PAR, calculado antes de inserir nada. O teto de
    # vizinhanca tem de agir SOBRE ELE, e nao sobre os vaos do corte: no corte
    # todos os vaos sao grandes e parecidos (>= 150 m), a razao entre vizinhos
    # e ~1 e o limite nao morde. O desnivel nasce DEPOIS -- o trecho ingreme e
    # subdividido a 25 m e o plano fica em 1.000 m, e ai a razao vira 14x.
    alvos = []
    for a_, b_ in zip(xs, xs[1:]):
        dxr_ = float(a_["rs"]) - float(b_["rs"])
        if dxr_ <= 0:
            alvos.append(float(op.espacamento))
            continue
        dz_ = abs(float(a_.get("z_terreno", 0.0))
                  - float(b_.get("z_terreno", 0.0)))
        S_ = max(dz_ / dxr_, 1e-6)
        A_ = 0.5 * (float(a_.get("area_km2", 1.0))
                    + float(b_.get("area_km2", 1.0)))
        D_ = kh * max(A_, 1.0) ** eh if getattr(op, "samuels_leopold", True) \
            else float(op.samuels_D)
        alvos.append(max(min(op.samuels_k * D_ / S_, dxr_), piso))
    alvos = _suavizar_alvo(np.array(alvos, float), op)

    saida, n_novas = [], 0
    for i, (a, b) in enumerate(zip(xs, xs[1:])):  # xs vem de montante p/ jusante
        saida.append(a)
        dxr = float(a["rs"]) - float(b["rs"])
        if dxr <= 0:
            continue
        A = 0.5 * (float(a.get("area_km2", 1.0)) + float(b.get("area_km2", 1.0)))
        lim = float(alvos[i])
        n = int(np.ceil(dxr / lim)) - 1
        if n <= 0:
            continue
        n = min(n, int(getattr(op, "interp_max", 40)))
        up = pd.DataFrame({"Station": a["sta"], "Elevation": a["z"]})
        dn = pd.DataFrame({"Station": b["sta"], "Elevation": b["z"]})
        for k in range(1, n + 1):
            t = k / (n + 1.0)
            # AS MARGENS AINDA NAO EXISTEM AQUI. lb/rb sao definidas ao escavar
            # a calha, um passo depois; pedi-las incondicionalmente levantava
            # KeyError, e o `except` mudo abaixo devolvia zero secoes
            # interpoladas como se nada fosse preciso. Falha silenciosa que
            # imita sucesso e o defeito mais caro deste projeto.
            bl = br = None
            if "lb" in a and "lb" in b:
                bl = (1 - t) * float(a["lb"]) + t * float(b["lb"])
                br = (1 - t) * float(a["rb"]) + t * float(b["rb"])
            try:
                df = G.interpolate_station_elevation(
                    up, dn, ratio=t, bank_left=bl, bank_right=br,
                    max_points=int(op.n_pontos))
            except Exception as e:                           # noqa: BLE001
                log(f"      interpolacao falhou entre RS {a['rs']:.0f} e "
                    f"{b['rs']:.0f}: {type(e).__name__} {e}")
                break
            r = dict(a)
            r["sta"] = np.asarray(df["Station"], float)
            r["z"] = np.asarray(df["Elevation"], float)
            # TALVEGUE NA JANELA DO EIXO, como o cortar() faz -- e nao o minimo
            # global. Com argmin cru, um meandro do proprio rio cruzado pelo
            # corte, ou uma clareira lida como depressao pelo MDS, rouba o
            # talvegue: medido nos 12 rios, 3% das secoes ficaram com ele na
            # BORDA e 41% fora do terco central. E o talvegue e onde a calha e
            # cavada, entao o canal ia parar fora do rio.
            #
            # A posicao do eixo nao e guardada na secao, mas o talvegue das
            # duas vizinhas ja esta dentro da janela dele (o cortar() garante),
            # e a interpolacao aqui e por posicao lateral NORMALIZADA -- entao
            # interpolar a posicao relativa das duas da uma referencia valida.
            def _rel(x):
                s = np.asarray(x["sta"], float)
                larg = float(s[-1] - s[0]) or 1.0
                return (float(s[int(x.get("i_thal", 0))]) - float(s[0])) / larg
            s_ = r["sta"]
            larg_ = float(s_[-1] - s_[0]) or 1.0
            eixo_ = float(s_[0]) + ((1 - t) * _rel(a) + t * _rel(b)) * larg_
            jan_ = max(op.canal_kw * max(A, 1.0) ** op.canal_ew, 150.0)
            m_ = np.flatnonzero(np.abs(s_ - eixo_) <= jan_)
            r["i_thal"] = int(m_[np.argmin(r["z"][m_])]) if len(m_) \
                else int(np.argmin(np.abs(s_ - eixo_)))
            r["rs"] = round(float(a["rs"]) - t * dxr, 2)
            r["interpolada"] = True
            for c in ("lb", "rb", "area_km2", "z_terreno", "S_terreno",
                      "n", "n_planicie", "s"):
                if c in a and c in b:
                    try:
                        r[c] = (1 - t) * float(a[c]) + t * float(b[c])
                    except (TypeError, ValueError):
                        pass
            # A CUTLINE SAI DA FAIXA DE ESTACAS, e nao da interpolacao dos
            # quatro numeros das vizinhas. Interpolar as pontas linearmente
            # parece obvio e esta errado: a biblioteca interpola as estacas por
            # POSICAO LATERAL NORMALIZADA, e quando as vizinhas tem larguras
            # diferentes -- 288 m e 658 m no Cedros -- os dois calculos
            # divergem. A cutline fica com um comprimento e as estacas com
            # outro.
            #
            # Isso importa porque o HEC-RAS mapeia estaca -> posicao ao longo
            # da CUTLINE: se os comprimentos nao batem, a secao inteira e
            # esticada ou comprimida no mapa. As edge lines e as bank lines,
            # que saem dai, deixam de coincidir com a secao -- que e
            # exatamente o que se ve no RAS Mapper.
            #
            # Medido: 284 secoes cortadas do terreno com descasamento ZERO,
            # 1.030 interpoladas com 46% acima de 1 m e pior caso de 370 m.
            #
            # O centro e a direcao vem das vizinhas (interpolar ponto e vetor e
            # legitimo); o COMPRIMENTO vem das estacas desta secao.
            if "cut" in a and "cut" in b:
                ca = np.asarray(a["cut"], float)
                cb = np.asarray(b["cut"], float)
                mi = (1 - t) * np.array([(ca[0] + ca[2]) / 2,
                                         (ca[1] + ca[3]) / 2]) \
                    + t * np.array([(cb[0] + cb[2]) / 2, (cb[1] + cb[3]) / 2])
                # O CENTRO TEM DE ESTAR NO EIXO, e nao na corda entre os
                # centros vizinhos. A secao cortada nasce em
                # `linha.interpolate(s)` -- um ponto DO eixo -- e por isso o
                # eixo cruza exatamente no meio dela (mediana medida: -0,00 m).
                # A interpolada saia da media dos dois centros, que num
                # meandro e a CORDA: passa por dentro da curva, longe do rio.
                #
                # Medido no Mirim: 160 secoes (11%) nao alcancavam o eixo, e
                # TODAS eram interpoladas -- nenhuma cortada. O RAS Mapper
                # reprova essas com "Some cross-sections do not cross a river
                # line - invalid geometry", nao consegue montar a superficie
                # de interpolacao delas ("Stored Interpolation Surface does
                # not contain XS'(s) at:") e o buraco se propaga: 236 pares
                # consecutivos ficaram sem superficie, de 1.417 esperadas.
                #
                # A DIRECAO continua vindo das vizinhas -- interpolar vetor e
                # legitimo, e o 382b95c ja tratou disso. So o centro muda.
                if eixo is not None and r.get("s") is not None:
                    q = eixo.interpolate(float(r["s"]))
                    mi = np.array([q.x, q.y])
                v = (1 - t) * np.array([ca[2] - ca[0], ca[3] - ca[1]]) \
                    + t * np.array([cb[2] - cb[0], cb[3] - cb[1]])
                # NAO chamar isto de `n`: `n` e o numero de secoes a inserir, e
                # `t = k / (n + 1)` e recalculado A CADA VOLTA do laco. Com o
                # nome repetido, da segunda secao em diante `t` passava a sair
                # do COMPRIMENTO DA CUTLINE (centenas de metros) em vez da
                # contagem, e as intermediarias desabavam para 0,07 a 3 m da
                # secao de montante -- enquanto o vao de 150 m seguia sem
                # refino nenhum. Era a origem das secoes coladas.
                norma = float(np.hypot(v[0], v[1])) or 1.0
                u = v / norma
                meia = 0.5 * float(r["sta"][-1] - r["sta"][0])
                r["cut"] = (float(mi[0] - meia * u[0]),
                            float(mi[1] - meia * u[1]),
                            float(mi[0] + meia * u[0]),
                            float(mi[1] + meia * u[1]))

                # NAO FORCAR AQUI O LADO DA ESTACA 0. Tentei: 35 secoes tinham
                # a estaca 0 na margem direita (33 delas interpoladas, porque
                # o vetor `v` e a SOMA dos vetores das vizinhas e num meandro
                # fechado a soma pode apontar para o outro lado), e inverter a
                # secao inteira parecia o conserto obvio.
                #
                # Piorou, nas duas medidas, com as mesmas opcoes:
                #     Validate Geometry           452 -> 612
                #     auto-intersecoes das edges  170 -> 214
                #
                # O motivo: a edge line liga PONTA A PONTA entre vizinhas. O
                # que importa nao e cada secao estar do lado "certo" em
                # absoluto, e sim concordar com a vizinha. Inverter uma secao
                # isolada troca a ponta que se conecta dos dois lados dela e
                # cria DOIS lacos onde havia um desalinhamento. Otimizar a
                # secao sozinha piora a linha que passa por ela.
            saida.append(r)
            n_novas += 1
    saida.append(xs[-1])
    saida.sort(key=lambda d: -d["rs"])
    if getattr(op, "curva_pos", True) and eixo is not None:
        saida = limitar_por_curvatura(saida, op, eixo, log)
    return saida


def limitar_por_curvatura(xs, op, eixo, log=print):
    """Aperta a largura das secoes onde o meandro fecha -- DEPOIS de densificar.

    O limitador que ja existia (`R = folga_curva * ds / dth`, dentro do
    cortar()) age em duas coisas que nao sao a geometria final: nas secoes
    CORTADAS, e nas larguras de ANTES do recorte pela cota de cheia. Entre ele
    e o arquivo gravado entram a densificacao -- que insere secoes entre as
    cortadas, com espacamento outro -- e o recorte, que muda toda a largura.

    O resultado e que a condicao de nao se cruzarem nunca era imposta sobre as
    secoes que o HEC-RAS realmente recebe. Nos meandros fechados do Mirim o
    raio de curvatura fica menor que a meia-largura (mediana 66 m), as
    vizinhas se cruzam pelo lado de dentro da curva, e a edge line que liga as
    pontas da o laco:

        "The generated edge lines have self intersections, the interpolation
         surface may not generate correctly because of this."

    Aqui o mesmo criterio e aplicado ao conjunto final, com o angulo lido da
    CUTLINE de cada secao (e nao do eixo) e a distancia medida entre elas.
    Apara so o que passa; nunca alarga.
    """
    if len(xs) < 3:
        return xs
    xs = sorted(xs, key=lambda d: -float(d["rs"]))       # montante -> jusante
    n = len(xs)
    s = np.array([float(d["s"]) if d.get("s") is not None else np.nan
                  for d in xs])
    if not np.all(np.isfinite(s)):
        return xs
    ang = np.unwrap([np.arctan2(d["cut"][3] - d["cut"][1],
                                d["cut"][2] - d["cut"][0]) for d in xs])

    # ONDE O EIXO CRUZA CADA SECAO, em estaca. Nao e o meio: a cortada nasce
    # com meias-larguras diferentes dos dois lados.
    c = np.empty(n)
    for i, d in enumerate(xs):
        A = np.array(d["cut"][:2], float)
        B = np.array(d["cut"][2:], float)
        st = np.asarray(d["sta"], float)
        from shapely.geometry import LineString
        g = LineString([A, B]).intersection(eixo)
        if g.is_empty:
            c[i] = 0.5 * (st[0] + st[-1])
            continue
        p = g if g.geom_type == "Point" else list(g.geoms)[0]
        ab = B - A
        t = float(np.dot(np.asarray(p.coords[0]) - A, ab)
                  / max(float(np.dot(ab, ab)), 1e-9))
        c[i] = st[0] + np.clip(t, 0.0, 1.0) * (st[-1] - st[0])

    esq = c - np.array([float(d["sta"][0]) for d in xs])
    dire = np.array([float(d["sta"][-1]) for d in xs]) - c
    le, ld = esq.copy(), dire.copy()
    for i in range(n):
        for j in (i - 1, i + 1):
            if j < 0 or j >= n:
                continue
            ds = abs(s[j] - s[i])
            dth = abs(ang[j] - ang[i])
            if ds <= 0 or dth < 1e-6:
                continue
            R = op.folga_curva * ds / dth
            if (ang[j] - ang[i]) * (1 if j > i else -1) > 0:
                le[i] = min(le[i], R)
            else:
                ld[i] = min(ld[i], R)

    # Alisar, pelo mesmo motivo do cortar(): o que derruba o solver nao e a
    # largura, e o SALTO de area e conducao entre vizinhas.
    def alisar(v, jan=5):
        if len(v) < jan:
            return v
        m = np.array([v[max(0, i - jan // 2):i + jan // 2 + 1].min()
                      for i in range(len(v))])
        return np.convolve(np.pad(m, 1, mode="edge"), np.ones(3) / 3.0, "valid")

    piso = float(getattr(op, "curva_piso", 30.0))
    le = np.maximum(alisar(le), piso)
    ld = np.maximum(alisar(ld), piso)

    saida, n_ap, antes, depois = [], 0, [], []
    for i, d in enumerate(xs):
        antes.append(float(d["sta"][-1] - d["sta"][0]))
        r = _aparar(d, c[i], le[i], ld[i])
        if r is not d:
            n_ap += 1
        depois.append(float(r["sta"][-1] - r["sta"][0]))
        saida.append(r)
    if n_ap:
        log(f"      {n_ap} secoes apertadas pelo raio do meandro "
            f"(largura mediana {np.median(antes):.0f} -> "
            f"{np.median(depois):.0f} m)")
    return saida


def _aparar(d, c, le, ld):
    """Corta a secao em [c-le, c+ld], com as pontas na cota interpolada.

    Estacas, cotas, cutline, margens e talvegue saem coerentes entre si -- e
    o comprimento da cutline continua sendo a faixa de estacas, que e o que o
    RAS Mapper exige ("The polyline length must match the last station minus
    the first station").
    """
    st = np.asarray(d["sta"], float)
    z = np.asarray(d["z"], float)
    a = max(float(st[0]), float(c - le))
    b = min(float(st[-1]), float(c + ld))
    # O CANAL NAO SE APARA. Cortando so pela curvatura, 5 a 10 secoes saiam
    # com a margem direita antes da esquerda -- a apara passava por dentro das
    # estacas de margem e o `clip` colapsava as duas no mesmo ponto. Aperta-se
    # a planicie; o que esta entre as margens fica.
    if d.get("lb") is not None and d.get("rb") is not None:
        a = min(a, float(d["lb"]))
        b = max(b, float(d["rb"]))
    if b - a < 5.0:                       # nao apara ate a secao sumir
        return d
    if a <= st[0] + 1e-6 and b >= st[-1] - 1e-6:
        return d                          # nada a fazer
    m = (st > a) & (st < b)
    ns = np.concatenate(([a], st[m], [b]))
    nz = np.concatenate(([float(np.interp(a, st, z))], z[m],
                         [float(np.interp(b, st, z))]))
    A = np.array(d["cut"][:2], float)
    B = np.array(d["cut"][2:], float)
    u = B - A
    u = u / max(float(np.hypot(u[0], u[1])), 1e-9)
    pa = A + (a - float(st[0])) * u
    pb = A + (b - float(st[0])) * u
    r = dict(d)
    r["sta"] = ns - ns[0]
    r["z"] = nz
    r["cut"] = (float(pa[0]), float(pa[1]), float(pb[0]), float(pb[1]))
    if d.get("lb") is not None and d.get("rb") is not None:
        r["lb"] = float(np.clip(float(d["lb"]) - a, 0.0, r["sta"][-1]))
        r["rb"] = float(np.clip(float(d["rb"]) - a, 0.0, r["sta"][-1]))
    # O TALVEGUE E O MESMO PONTO, se ele sobreviveu ao corte -- e nao o menor
    # da secao nova. Trocar por argmin joga o canal para onde a apara deixou a
    # cota mais baixa, que pode ser a propria ponta recem-criada.
    s_thal = float(st[int(d.get("i_thal", 0))])
    if a <= s_thal <= b:
        r["i_thal"] = int(np.argmin(np.abs(r["sta"] - (s_thal - a))))
    else:
        r["i_thal"] = int(np.argmin(nz))
    return r
