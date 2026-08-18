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


def estacas(linha, amostrador, op):
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
    # nao encosta nos extremos: secao em cima da juncao conflita com o
    # comprimento declarado em Junc L&A e trava o solver
    recuo = op.espacamento_min * 0.5
    ss, s = [recuo], recuo
    while s < L - recuo:
        s += float(np.interp(s, d, dx))
        if s < L - recuo:
            ss.append(s)
    ss.append(L - recuo)
    return ss


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
    dirs = [direcao(linha, s, op.espacamento_min * 1.7) for s in ss]
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
    tx, ty = direcao(linha, s, op.espacamento_min * 1.7)
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


def cortar_rio(eixo, amostrador, op, log=print, prog=None):
    """Todas as secoes de um rio, da cabeceira para a foz."""
    linha = eixo["linha"]
    L = linha.length
    ss = estacas(linha, amostrador, op)

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
    log(f"   {eixo['ras']:<16} {len(xs):>4} secoes   "
        f"largura {min(d['sta'][-1] for d in xs):.0f} a "
        f"{max(d['sta'][-1] for d in xs):.0f} m")
    return xs
