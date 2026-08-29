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
CALHA_SINTETICA = False  # escavar batimetria que o DEM nao mostra?
                         # O Copernicus achata a lamina d'agua na cota da
                         # SUPERFICIE, entao o leito real esta abaixo do que
                         # ele mostra e a calha era escavada para compensar.
                         # Mas com isso o modelo parte SECO e o assentamento
                         # tem de encher 382 hm3 de canal inventado -- e onde
                         # ele falha. Desligado, o leito e a lamina do DEM:
                         # aproximacao grosseira, porem real, e o modelo parte
                         # praticamente cheio.
CANAL_KH, CANAL_EH = 0.277, 0.35     # profundidade = KH * A^EH
CANAL_KW, CANAL_EW = 5.0, 0.40       # largura     = KW * A^EW
ALTURA_MARGEM = 3.0                  # folga acima do topo da calha

# --- pilot channel (o mesmo recurso do HEC-RAS: Geometric Data > Pilot
#     Channels, com largura, n e inverts interpolados). Um entalhe estreito e
#     raso no talvegue.
#     Por que e necessario aqui: com CALHA_SINTETICA desligada o fundo da secao
#     e CHATO, com ~100 m de largura na cota do talvegue do DEM. Em lamina
#     baixa -- que e o aquecimento, exatamente onde o modelo falha, entre
#     00:07 e 00:25 de simulacao -- a area molhada e o raio hidraulico ficam
#     mal definidos: alguns centimetros de agua espalhados por 100 m. A
#     conducao calculada sobre isso oscila, e o solver nao converge.
#     O entalhe da um caminho de escoamento continuo e bem condicionado desde
#     a primeira iteracao, sem alterar a capacidade de cheia (25 m x 1,5 m sao
#     37 m2 num rio que conduz milhares de m3/s).
PILOT_ATIVO = True
PILOT_LARGURA = 25.0                 # m
PILOT_PROF = 1.5                     # m abaixo do fundo imposto


def largura_secao(area_km2):
    """Meia-largura da secao, conforme o porte do rio.

    Alargar isto nao e trivial: com secao larga demais o corte atravessa
    meandros do proprio rio e o escoamento e contado duas vezes (subir o
    coeficiente de 180 para 440 derrubou a simulacao de 30 para 2 passos).
    So funciona junto com area de escoamento inefetivo.
    """
    return float(np.clip(180.0 * np.sqrt(max(area_km2, 1.0) / 100.0), 500.0, 2500.0))


def canal(area_km2):
    """(profundidade, largura) da calha para a area de drenagem dada.

    Com CALHA_SINTETICA desligado a profundidade e zero: o leito passa a ser a
    lamina que o DEM mostra. A largura continua valendo, porque define a zona
    de canal do Manning e a busca das margens.
    """
    a = max(area_km2, 1.0)
    prof = CANAL_KH * a ** CANAL_EH if CALHA_SINTETICA else 0.0
    return prof, CANAL_KW * a ** CANAL_EW


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
    zb = tirar_picos(np.interp(d, d[ok], zb[ok]))
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


def indice_eixo(sta, z, janela, i_eixo=None):
    """Indice do talvegue PROXIMO AO EIXO -- nao o minimo global."""
    i = len(sta) // 2 if i_eixo is None else int(i_eixo)
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
RAZAO_LADOS = 2.5        # largura maxima de um lado em relacao ao outro
MINIMO_LADO = 120.0      # m, meia-largura minima de cada lado


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

    minimo = MINIMO_LADO
    e = np.maximum(alisar(esq), minimo)
    d = np.maximum(alisar(dir_), minimo)
    # Centrar a calha no corte. O limite de curvatura apara SO o lado concavo,
    # entao numa curva fechada a secao sai com 90 m de um lado e 700 do outro:
    # a calha encostada na borda, sem planicie de um lado (foi o Bank Sta
    # 2,58/92,99 num corte de 800 m no Mirim). Estreitar o lado LARGO nunca
    # cria cruzamento -- cruzamento vem de largura a mais --, entao limitar a
    # razao entre os lados centra a calha sem afrouxar o criterio de curvatura.
    e, d = np.minimum(e, d * RAZAO_LADOS), np.minimum(d, e * RAZAO_LADOS)
    return np.maximum(e, minimo), np.maximum(d, minimo)


def cortar(linha, s, amostrador, meia_largura, area_km2, hw_esq=None,
           hw_dir=None):
    """Uma secao perpendicular ao eixo na posicao s."""
    tx, ty = direcao(linha, s)
    rx, ry = ty, -tx                               # normal a direita
    p = linha.interpolate(s)
    he = float(hw_esq if hw_esq is not None else meia_largura)
    hd = float(hw_dir if hw_dir is not None else meia_largura)
    # ESPACAMENTO UNIFORME. Dividir os pontos ao meio entre os dois lados so
    # funciona se eles tiverem a mesma largura -- e desde o limitador de
    # curvatura nao tem: com he=200 m e hd=700 m saia 1,4 m de espacamento a
    # esquerda e 5,0 m a direita. Conducao calculada sobre pontos desigualmente
    # espacados fica enviesada para o lado denso, e o proprio talvegue e
    # procurado num indice que nao corresponde ao eixo. Aqui os pontos sao
    # repartidos na PROPORCAO das larguras, e o eixo continua sendo um ponto da
    # tabela.
    n_e = int(round(N_PONTOS * he / max(he + hd, 1e-6)))
    n_e = min(max(n_e, 2), N_PONTOS - 2)
    off = np.concatenate([np.linspace(-he, 0, n_e, endpoint=False),
                          np.linspace(0, hd, N_PONTOS - n_e)])
    i_eixo = n_e
    z = amostrador.cota(p.x + off * rx, p.y + off * ry)
    if np.isnan(z).all():
        return None
    if np.isnan(z).any():
        ok = ~np.isnan(z)
        z = np.interp(np.arange(len(z)), np.flatnonzero(ok), z[ok])
    sta = off + he

    # A secao sai do TERRENO. A calha e escavada depois, uma unica vez, com o
    # perfil longitudinal ja definido -- ver escavar(). Escavar aqui e depois
    # reajustar a cada passo do condicionamento reaplicava o trapezio sobre um
    # perfil ja modificado, e o que sobrava era um pico isolado: uma fenda de
    # 3,4 m num rio de calha de 106 m, que nao conduz nada.
    prof, larg = canal(area_km2)
    i0 = indice_eixo(sta, z, max(larg, 150.0), i_eixo)
    # UMA calha por secao: sobe o que estiver mais fundo que o talvegue do eixo
    # (o leito antigo, um meandro que o corte cruzou)
    z = np.maximum(z, z[i0])
    cut = (p.x - he * rx, p.y - he * ry, p.x + hd * rx, p.y + hd * ry)
    return {"sta": sta, "z": z, "i_thal": i0, "cut": cut,
            "area_km2": area_km2, "prof_canal": prof, "larg_canal": larg,
            "z_terreno": float(z[i0])}


def escavar(d):
    """Escava a calha UMA vez, com o fundo na cota que o perfil definiu.

    Chamada depois do condicionamento e da ancoragem, quando d["z_alvo"] ja
    tem a cota final do talvegue. O trapezio e aplicado sobre o TERRENO, entao
    e sempre integro -- largura cheia, taludes suaves -- em vez de resultar da
    soma de varios ajustes parciais.
    """
    sta, z = d["sta"], np.array(d["z"], float)
    larg = d["larg_canal"]
    i0 = d["i_thal"]
    alvo = d.get("z_alvo", z[i0] - d["prof_canal"])
    dist = np.abs(sta - sta[i0])
    meia = larg / 2.0
    talude = max(larg * 0.25, 30.0)

    # A calha e IMPOSTA, nao subtraida. Subtrair uma profundidade constante
    # preserva a forma do terreno: vale em V da calha em V, vale em U da calha
    # em U, e a area molhada salta entre secoes vizinhas -- eram 121 pares com
    # mais de 3x de diferenca a 1 m de lamina, contra 31 a 8 m, ou seja o
    # defeito estava na CALHA e nao na planicie.
    # Imposta, ela e um trapezio de fundo plano em z_alvo e largura larg, com
    # talude subindo ate encontrar o terreno. Como larg e z_alvo variam
    # suavemente ao longo do rio (area de drenagem e perfil condicionado), a
    # conducao de estiagem passa a ser continua por construcao.
    # O minimo garante que a calha so CORTA o terreno, nunca o preenche.
    subida = np.clip((dist - meia) / talude, 0.0, 1.0)
    z_canal = alvo + subida * np.maximum(z - alvo, 0.0)
    z = np.minimum(z, z_canal)

    # PILOT CHANNEL. O invert acompanha z_alvo, que ja vem do perfil alisado,
    # entao o entalhe e continuo rio abaixo por construcao -- que e o que o
    # HEC-RAS obtem interpolando os inverts de montante e jusante na ferramenta
    # dele.
    if PILOT_ATIVO and PILOT_PROF > 0:
        meia_p = PILOT_LARGURA / 2.0
        talude_p = max(PILOT_LARGURA * 0.6, 10.0)
        # O limite tem de VOLTAR AO TERRENO ao subir, como o trapezio
        # principal. Escrito como
        #     min(z, (alvo - PILOT_PROF) + sobe * PILOT_PROF)
        # o limite satura em 'alvo' longe do canal, e min(z, alvo) rebaixa a
        # secao INTEIRA ate o fundo: as 1.232 secoes do modelo viraram bacias
        # chatas com um entalhe no meio e parede vertical nas pontas. No
        # Itajai_Acu R4 RS 34.956 o terreno real vai de 0,50 a 51,54 m, com a
        # encosta subindo a 50 m, e a secao gravada tinha QUATRO cotas em
        # 2.908 m de largura.
        base_p = alvo - PILOT_PROF
        sobe = np.clip((dist - meia_p) / talude_p, 0.0, 1.0)
        z = np.minimum(z, base_p + sobe * np.maximum(z - base_p, 0.0))

    # SECAO RASA: onde o vale e plano de verdade, nem 10x a largura acha
    # terreno alto -- no Taio havia secoes de 1.000 m com 2,25 m de desnivel
    # TOTAL, e qualquer lamina extrapolava a tabela de conducao. O log do
    # solver as lista nominalmente ("Extrapolated above Cross Section Table").
    #
    # Nesse caso a saida e fechar a secao com parede vertical nas pontas, que
    # e o recurso padrao do proprio HEC-RAS ("glass wall"): a agua fica contida
    # e a tabela cobre a faixa toda. So entra em jogo se a lamina chegar la;
    # onde o terreno ja e alto, nada muda. Diferente da tentativa anterior,
    # aqui e aplicado APENAS as secoes que precisam, nao a todas.
    precisa = altura_para_vazao(sta, z, d.get("n"), d.get("S_terreno"),
                                vazao_projeto(d.get("area_km2", 10.0)))
    minima = float(np.clip(FOLGA_ALTURA * precisa,
                           ALTURA_MINIMA_SECAO, ALTURA_MAX_SECAO))
    util = float(z.max() - z[i0])
    if util < minima:
        alvo_topo = z[i0] + minima
        z[0] = max(z[0], alvo_topo)
        z[-1] = max(z[-1], alvo_topo)
        d["parede"] = round(minima - util, 2)
    d["h_precisa"] = round(precisa, 2)

    d["z"] = z
    # A profundidade que a margem enxerga e a ESCAVACAO -- terreno menos o
    # fundo --, nao a diferenca para o ponto vizinho. Com o fundo agora plano o
    # vizinho tambem esta em z_alvo, isso dava zero, e as margens eram
    # procuradas DENTRO da propria calha: o Rio Benedito saia com profundidade
    # de calha 0,00 m e area zero em todas as secoes.
    d["lb"], d["rb"] = margens(sta, z, i0,
                               max(d.get("z_terreno", z[i0]) - alvo, 0.5))
    return d


ALTURA_ALVO = 15.0    # m de desnivel que a secao deve ter acima do talvegue
ALTURA_MINIMA_SECAO = 12.0   # piso absoluto, quando a vazao pede menos
ALTURA_MAX_SECAO = 30.0      # teto: parede mais alta que isto nao e fisica
FOLGA_ALTURA = 1.4           # borda livre sobre a lamina de projeto
# Q = K * A^EXP, ancorado nos ~5.700 m3/s de 1983 na foz com 15.000 km2. Serve
# so para DIMENSIONAR a secao; a vazao que roda vem da hidrologia.
Q_K, Q_EXP = 2.6, 0.80
FATOR_MAX = 10.0      # ate quantas vezes a largura base pode ser esticada


def tirar_picos(z, janela=5, limite=3.0):
    """Remove picos isolados do talvegue amostrado no DEM.

    O condicionamento impoe decrescimento PAR A PAR, e um pico para cima
    seguido de descida satisfaz esse teste em cada par. No Itajai-Mirim
    sobrava um degrau de 19 m entre secoes a 150 m uma da outra:

        RS 91807   leito 126,99
        RS 91657   leito 146,01   <- sobe 19 m
        RS 91507   leito 123,48   <- desce 22 m

    e era exatamente ali que o solver reportava erro maximo em 19 das 63
    linhas do log. Aqui o valor e trocado pela mediana da vizinhanca sempre
    que se afasta dela mais que 'limite' metros -- filtro de mediana classico,
    que preserva o degrau REAL (uma queda sustentada) e remove o isolado.
    """
    z = np.asarray(z, float).copy()
    n = len(z)
    if n < janela:
        return z
    m = janela // 2
    med = np.array([np.median(z[max(0, i - m):min(n, i + m + 1)])
                    for i in range(n)])
    fora = np.abs(z - med) > limite
    z[fora] = med[fora]
    return z


def vazao_projeto(area_km2):
    """Vazao de pico de referencia para uma area de drenagem."""
    return Q_K * max(float(area_km2), 1.0) ** Q_EXP


def altura_para_vazao(sta, z, n, S, Q, h_max=ALTURA_MAX_SECAO):
    """Menor lamina sobre o talvegue que conduz Q, por Manning na PROPRIA secao.

    O criterio anterior era um numero fixo (12 m de parede, 15 m de alvo) igual
    para o Itajai-Acu e para um afluente de 30 km2. Fixo demais para o afluente
    e de menos para o rio grande: o solver listou nominalmente as secoes que
    passaram do topo da tabela ("Extrapolated above Cross Section Table"), e
    todas estavam travadas em 12,00 m -- o piso, nao o terreno.
    """
    z0 = float(np.min(z))
    S = max(float(S or 0.0), 1e-4)
    n = max(float(n or 0.035), 0.02)
    for h in np.arange(1.0, h_max + 0.01, 0.5):
        prof = np.clip(z0 + h - z, 0.0, None)
        A = float(np.trapezoid(prof, sta)) if hasattr(np, "trapezoid")             else float(np.trapz(prof, sta))
        if A <= 0.0:
            continue
        molh = (prof[:-1] + prof[1:]) > 0.0
        P = float(np.sum(np.hypot(np.diff(sta), np.diff(z))[molh]))
        if P <= 0.0:
            continue
        if (A / n) * (A / P) ** (2.0 / 3.0) * S ** 0.5 >= Q:
            return float(h)
    return float(h_max)


def alvo_por_area(area_km2, S):
    """Altura alvo do ALARGAMENTO, antes de existir secao para medir.

    Canal largo: R ~ h, entao h = (Q n / (w sqrt(S)))^(3/5), com w tomado como
    tres vezes a largura de calha da relacao hidraulica.
    """
    Q = vazao_projeto(area_km2)
    w = max(3.0 * CANAL_KW * max(float(area_km2), 1.0) ** CANAL_EW, 20.0)
    S = max(float(S or 0.0), 1e-4)
    h = (Q * 0.045 / (w * S ** 0.5)) ** 0.6
    return float(np.clip(FOLGA_ALTURA * h, ALTURA_ALVO, ALTURA_MAX_SECAO))


def alargar_ate_conter(linha, s, amostrador, he, hd, area_km2,
                       alvo=ALTURA_ALVO):
    """Estica a secao ate encontrar terreno 'alvo' metros acima do talvegue.

    A largura vinha so da area de drenagem. Onde o vale e raso isso nao alcanca
    terreno alto: no Taio, no Trombudo, no Benedito e no Mirim a secao tinha
    menos de 8 m de desnivel util, e na cheia de 1983 -- 2,4x maior que a
    sintetica -- o HEC-RAS extrapolava a tabela de conducao justamente nessas
    secoes ("Extrapolated above Cross Section Table"), com erro chegando a
    21 m. Sao os mesmos rios que a auditoria ja marcava por altura util.

    Estica em passos, so ate conseguir a altura ou atingir FATOR_MAX. Nao e
    parede vertical: e ir buscar o terreno que existe mais longe.
    """
    tx, ty = direcao(linha, s)
    rx, ry = ty, -tx
    p = linha.interpolate(s)
    for f in (1.0, 1.5, 2.0, 3.0, 5.0, 7.0, FATOR_MAX):
        e, d = he * f, hd * f
        off = np.linspace(-e, d, 60)
        z = amostrador.cota(p.x + off * rx, p.y + off * ry)
        if not np.isfinite(z).any():
            break
        if float(np.nanmax(z) - np.nanmin(z)) >= alvo:
            return he * f, hd * f
    return he * FATOR_MAX, hd * FATOR_MAX


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
        zt = tirar_picos(np.interp(np.arange(len(zt)), np.flatnonzero(ok), zt[ok]))
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
    # alarga onde o vale e raso, ANTES do corte definitivo
    for i, s in enumerate(ss):
        hw_e[i], hw_d[i] = alargar_ate_conter(
            linha, s, amostrador, hw_e[i], hw_d[i], areas[i],
            alvo=alvo_por_area(areas[i], S_terr[i]))
    # e volta a limitar por curvatura, agora com as larguras novas
    hw_e2, hw_d2 = limites_por_curvatura(linha, ss, np.maximum(hw_e, hw_d))
    hw_e = np.minimum(hw_e, hw_e2)
    hw_d = np.minimum(hw_d, hw_d2)
    # DESLOCA a janela do corte para a calha nao ficar no canto. Limitar a
    # razao antes do alargamento nao adianta: e o alargamento que joga a calha
    # para a borda, porque ele estica o lado onde houver encosta ate achar
    # altura. No Itajai_Mirim RS 104738,3 -- onde o solver abortou -- a calha
    # saiu em 511,54/616,95 num corte que termina em 620: TRES metros de
    # planicie a direita do rio, contra 511 do outro lado.
    #
    # Deslocar e melhor que aparar, porque preserva a largura total: o que
    # falta de um lado sai do outro. O que nao se pode e passar do limite de
    # curvatura, que e onde as cutlines se cruzam. Entao desloca-se ate esse
    # limite e, se ainda sobrar desequilibrio, ai sim apara-se o lado largo
    # (estreitar nunca cria cruzamento).
    minlado = (hw_e + hw_d) / (1.0 + RAZAO_LADOS)
    for falta, pode, dar in ((minlado - hw_d, hw_d2 - hw_d, "d"),
                             (minlado - hw_e, hw_e2 - hw_e, "e")):
        outro = hw_e if dar == "d" else hw_d
        mv = np.clip(np.minimum(np.minimum(falta, pode), outro - MINIMO_LADO),
                     0.0, None)
        if dar == "d":
            hw_d, hw_e = hw_d + mv, hw_e - mv
        else:
            hw_e, hw_d = hw_e + mv, hw_d - mv
    hw_e, hw_d = (np.minimum(hw_e, hw_d * RAZAO_LADOS),
                  np.minimum(hw_d, hw_e * RAZAO_LADOS))
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
