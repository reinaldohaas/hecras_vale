# -*- coding: utf-8 -*-
"""
Secoes transversais: onde cortar, como cortar, onde esta a calha.

Tres decisoes, cada uma paga com uma rodada de depuracao no modelo anterior:

ESPACAMENTO ADAPTATIVO. O Acu cai 195 m em 13 km na garganta do Salto Pilao
(confirmado em terreno de 1 m; nao e ruido de DEM). A 1 km de espacamento sao
8 m de queda ENTRE SECOES VIZINHAS, e o solver falha no primeiro passo. O
criterio usual para regime nao permanente e dx <~ 0,15*D/S, que com D~4 m e
S=0,008 da ~75 m. Entao o espacamento sai da declividade local.

CALHA NO EIXO. Localizar o talvegue pelo minimo GLOBAL da secao poe a calha no
lugar errado sempre que o corte atravessa outro canal mais fundo -- o leito
antigo do Mirim, um meandro do proprio rio. Com o eixo vindo do relevo isso e
raro, mas a busca fica restrita a uma janela em torno do eixo de qualquer jeito.

UMA CALHA POR SECAO. Depois de escavar, qualquer ponto AINDA mais fundo e outro
canal que o corte cruzou. Num 1D isso nao e um segundo caminho de escoamento: e
um poco mais fundo que o canal principal dentro da mesma secao, e a conducao
calculada em cima disso nao tem sentido.
"""
import numpy as np
from shapely.geometry import LineString

# --- corte
ESPACAMENTO = 1000.0     # m, no vale plano
ESPACAMENTO_MIN = 150.0  # m, na garganta
DECL_PLANO = 0.0010      # ate aqui, ESPACAMENTO
DECL_INGREME = 0.0060    # daqui pra cima, ESPACAMENTO_MIN
SUAVIZA_DIR = 250.0      # janela p/ a direcao local (evita cutlines cruzadas)
N_PONTOS = 280           # pontos por secao (limite do HEC-RAS: 450)

# --- calha, por geometria hidraulica (Leopold & Maddock)
#     ATENCAO: estes dois coeficientes NAO sao medidos. Sao relacao regional
#     aplicada a area de drenagem. Em Blumenau dao 6,8 m de escavacao abaixo da
#     lamina que o DEM mostra, o que poe o leito a -4 m. Se houver batimetria
#     real do Itajai-Acu, e aqui que ela entra.
CANAL_KH, CANAL_EH = 0.277, 0.35     # profundidade = KH * A^EH
CANAL_KW, CANAL_EW = 5.0, 0.40       # largura     = KW * A^EW
ALTURA_MARGEM = 3.0                  # folga acima do topo da calha


def largura_secao(area_km2):
    """Meia-largura da secao, conforme o porte do rio.

    Alargar isto nao e trivial: com secao larga demais o corte atravessa
    meandros do proprio rio e o escoamento e contado duas vezes (subir o
    coeficiente de 180 para 440 derrubou a simulacao de 30 para 2 passos).
    So funciona junto com area de escoamento inefetivo.
    """
    return float(np.clip(180.0 * np.sqrt(max(area_km2, 1.0) / 100.0), 500.0, 2500.0))


def canal(area_km2):
    """(profundidade, largura) da calha para a area de drenagem dada."""
    a = max(area_km2, 1.0)
    return CANAL_KH * a ** CANAL_EH, CANAL_KW * a ** CANAL_EW


def estacas(linha, amostrador):
    """Posicoes de corte, adensadas onde o leito e ingreme."""
    L = linha.length
    passo = min(ESPACAMENTO / 4.0, 250.0)
    d = np.arange(0.0, L + passo, passo)
    P = [linha.interpolate(float(x)) for x in d]
    zb = amostrador.talvegue([p.x for p in P], [p.y for p in P])
    ok = np.isfinite(zb)
    if ok.sum() < 3:
        return list(np.arange(0.0, L, ESPACAMENTO))
    zb = np.interp(d, d[ok], zb[ok])
    S = np.abs(np.gradient(zb, d))
    S = np.convolve(S, np.ones(5) / 5.0, "same")          # tira ruido do DEM
    f = np.clip((S - DECL_PLANO) / (DECL_INGREME - DECL_PLANO), 0.0, 1.0)
    dx = ESPACAMENTO + (ESPACAMENTO_MIN - ESPACAMENTO) * f
    # nao encosta nos extremos: secao em cima da juncao conflita com o
    # comprimento declarado em Junc L&A e trava o solver
    recuo = ESPACAMENTO_MIN * 0.5
    ss, s = [recuo], recuo
    while s < L - recuo:
        s += float(np.interp(s, d, dx))
        if s < L - recuo:
            ss.append(s)
    ss.append(L - recuo)
    return ss


def indice_eixo(sta, z, janela):
    """Indice do talvegue PROXIMO AO EIXO -- nao o minimo global."""
    i = len(sta) // 2
    m = np.abs(np.asarray(sta) - sta[i]) <= janela
    idx = np.flatnonzero(m & np.isfinite(z))
    return int(idx[np.nanargmin(np.asarray(z)[idx])]) if len(idx) else i


def _escavar(sta, z, prof, larg, i0):
    """Trapezio de largura 'larg' e profundidade 'prof', centrado em i0."""
    d = np.abs(sta - sta[i0])
    meia = larg / 2.0
    talude = max(larg * 0.25, 30.0)
    frac = np.clip(1.0 - (d - meia) / talude, 0.0, 1.0)
    frac[d <= meia] = 1.0
    return z - prof * frac, frac


def margens(sta, z, i0, prof_canal):
    """Estacas das margens: do talvegue ate o topo da calha + folga.

    Medir a folga a partir do TALVEGUE poe a margem dentro do canal escavado, e
    o modelo passa a achar que tudo extravasa. A margem real e o topo da calha.
    """
    lim = z[i0] + prof_canal + ALTURA_MARGEM
    e = i0
    while e > 0 and z[e - 1] < lim:
        e -= 1
    d = i0
    while d < len(z) - 1 and z[d + 1] < lim:
        d += 1
    e = min(max(e, 1), len(sta) - 3)
    d = max(min(d, len(sta) - 2), e + 1)
    return round(float(sta[e]), 2), round(float(sta[d]), 2)


FOLGA_CURVA = 0.70       # fracao do raio de curvatura ate onde a secao pode ir


def direcao(linha, s):
    """Direcao local do eixo, suavizada em +-SUAVIZA_DIR.

    Com janela de +-1 m as cutlines se cruzam nas curvas e o RAS avisa
    "edge lines have self intersections".
    """
    a = linha.interpolate(max(0.0, s - SUAVIZA_DIR))
    b = linha.interpolate(min(linha.length, s + SUAVIZA_DIR))
    tx, ty = b.x - a.x, b.y - a.y
    n = np.hypot(tx, ty) or 1.0
    return tx / n, ty / n


def limites_por_curvatura(linha, estacas, meia_largura):
    """Ate onde cada semi-secao pode ir sem cruzar as vizinhas.

    Numa curva de raio R, duas perpendiculares vizinhas convergem e se
    encontram a R do eixo, do lado CONCAVO. Passar disso e o que produz
    cutlines cruzadas -- 24% dos pares vizinhos, no primeiro corte desta
    reescrita. E dai que sai a mancha continua e sem sentido no Depth do RAS
    Mapper: ele interpola a superficie d'agua ENTRE as cutlines, e onde elas se
    cruzam a interpolacao nao tem significado.

    R sai da variacao de angulo entre estacas consecutivas: R = ds / |dtheta|.
    Aparar so o lado concavo preserva a largura nos trechos retos, que e onde
    a planicie precisa dela.
    """
    n = len(estacas)
    dirs = [direcao(linha, s) for s in estacas]
    ang = np.unwrap([np.arctan2(t[1], t[0]) for t in dirs])
    esq = np.full(n, float(meia_largura) if np.isscalar(meia_largura)
                  else 0.0)
    dir_ = esq.copy()
    if not np.isscalar(meia_largura):
        esq = np.asarray(meia_largura, float).copy()
        dir_ = esq.copy()
    for i in range(n):
        for j in (i - 1, i + 1):
            if j < 0 or j >= n:
                continue
            ds = abs(estacas[j] - estacas[i])
            dth = abs(ang[j] - ang[i])
            if ds <= 0 or dth < 1e-6:
                continue
            R = FOLGA_CURVA * ds / dth
            # dtheta > 0 = curva a esquerda: as perpendiculares convergem la
            if (ang[j] - ang[i]) * (1 if j > i else -1) > 0:
                esq[i] = min(esq[i], R)
            else:
                dir_[i] = min(dir_[i], R)
    # Suaviza ao longo do trecho. Sem isto o limite apara a secao i e nao a
    # i+1, e a largura sai em dente de serra -- no Rio do Testo dava
    # 1000, 941, 746, 737, 1000, 846, ... O que importa nao e a largura em si,
    # e o SALTO entre vizinhas: area e conducao mudam de degrau e o solver ve
    # uma contracao seguida de expansao a cada secao. Um minimo movel garante
    # que a secao nunca seja mais larga que as vizinhas precisam, e a media
    # movel tira o degrau que sobra.
    def alisar(v, jan=5):
        n_ = len(v)
        if n_ < jan:
            return v
        m = np.array([v[max(0, i - jan // 2):i + jan // 2 + 1].min()
                      for i in range(n_)])
        k = np.ones(3) / 3.0
        return np.convolve(np.pad(m, 1, mode="edge"), k, "valid")

    minimo = 120.0
    return (np.maximum(alisar(esq), minimo), np.maximum(alisar(dir_), minimo))


def cortar(linha, s, amostrador, meia_largura, area_km2, hw_esq=None,
           hw_dir=None):
    """Uma secao perpendicular ao eixo na posicao s."""
    tx, ty = direcao(linha, s)
    rx, ry = ty, -tx                               # normal a direita
    p = linha.interpolate(s)
    he = float(hw_esq if hw_esq is not None else meia_largura)
    hd = float(hw_dir if hw_dir is not None else meia_largura)
    off = np.concatenate([np.linspace(-he, 0, N_PONTOS // 2, endpoint=False),
                          np.linspace(0, hd, N_PONTOS - N_PONTOS // 2)])
    z = amostrador.cota(p.x + off * rx, p.y + off * ry)
    if np.isnan(z).all():
        return None
    if np.isnan(z).any():
        ok = ~np.isnan(z)
        z = np.interp(np.arange(len(z)), np.flatnonzero(ok), z[ok])
    sta = off + he

    prof, larg = canal(area_km2)
    i0 = indice_eixo(sta, z, max(larg, 150.0))
    z, _ = _escavar(sta, z, prof, larg, i0)
    # UMA calha: sobe o que ficou mais fundo que ela (outro canal cruzado)
    i0 = indice_eixo(sta, z, max(larg, 150.0))
    z = np.maximum(z, z[i0])
    lb, rb = margens(sta, z, i0, prof)
    cut = (p.x - he * rx, p.y - he * ry, p.x + hd * rx, p.y + hd * ry)
    return {"sta": sta, "z": z, "i_thal": i0, "lb": lb, "rb": rb,
            "cut": cut, "area_km2": area_km2, "prof_canal": prof,
            "larg_canal": larg}


def cortar_trecho(linha, amostrador, area_foz, rs0=0.0, area_cabeceira=None):
    """Todas as secoes de um rio, de montante para jusante.

    A largura cresce rio abaixo com a area de drenagem: usar a largura da foz
    na cabeceira desperdica os pontos na encosta e deixa o canal com 1-3 deles.
    """
    L = linha.length
    a0 = area_cabeceira if area_cabeceira is not None else area_foz * 0.05
    ss = estacas(linha, amostrador)
    # declividade do terreno NAO condicionado, guardada por secao. E ela que
    # deve alimentar o Manning de Jarrett: derivar n do perfil ja condicionado
    # usa a declividade do CLAMP (0,008 cravado), nao a do relevo, e portanto
    # subestima a rugosidade justamente nos afluentes de serra -- que sao os
    # que saturam o clamp.
    P = [linha.interpolate(float(s)) for s in ss]
    zt = amostrador.talvegue([p.x for p in P], [p.y for p in P])
    ok = np.isfinite(zt)
    if ok.sum() >= 2:
        zt = np.interp(np.arange(len(zt)), np.flatnonzero(ok), zt[ok])
        s_arr = np.asarray(ss, float)
        S_terr = np.abs(np.gradient(zt, s_arr))
        S_terr = np.convolve(S_terr, np.ones(3) / 3.0, "same")
    else:
        S_terr = np.zeros(len(ss))

    areas = [a0 + (area_foz - a0) * (s / max(L, 1.0)) for s in ss]
    hw = [largura_secao(a) for a in areas]
    # apara o lado concavo de cada secao ANTES de cortar: e o que impede as
    # cutlines de se cruzarem nas curvas
    hw_e, hw_d = limites_por_curvatura(linha, ss, hw)
    xs = []
    for i, s in enumerate(ss):
        r = cortar(linha, s, amostrador, hw[i], areas[i],
                   hw_esq=hw_e[i], hw_dir=hw_d[i])
        if r is None:
            continue
        r["rs"] = round(rs0 + (L - s), 2)
        r["S_terreno"] = float(S_terr[i])
        xs.append(r)
    xs.sort(key=lambda d: -d["rs"])
    fin, visto = [], set()
    for d in xs:                                   # o RAS exige RS unico
        if d["rs"] in visto:
            continue
        visto.add(d["rs"])
        fin.append(d)
    return fin
