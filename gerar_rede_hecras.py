# -*- coding: utf-8 -*-
"""
Gerador da REDE 1D real da Bacia do Itajai para o HEC-RAS 7.0.1 (Fase 1).

Monta a rede dendritica a partir da TOPOLOGIA OFICIAL da ANA (campo NUTRJUS
da BHO 2017) e do relevo real (DEM Copernicus 30 m):

    Itajai do Sul   ─┐
                     ├─[Rio do Sul   km   0,0]──► Acu R1
    Itajai do Oeste ─┘
    Itajai do Norte ──[Ibirama      km  39,3]──► Acu R2
    Rio Benedito    ──[Indaial      km  93,0]──► Acu R3
    Itajai-Mirim    ──[Itajai       km 180,8]──► Acu R4 ──► foz

Diferencas em relacao a gerar_geometria_hecras.py (1 trecho unico):
  - a rede vem de NUTRJUS, nao de linemerge/max(length) — os afluentes
    deixam de ser descartados;
  - o Itajai-Acu passa a ser a calha principal ate a foz;
  - o Itajai-Mirim entra como AFLUENTE, o que elimina a duplicacao de
    secoes no canal de Itajai (antes Mirim e Acu tinham secoes no mesmo
    ponto com leitos divergentes de 1,5 m);
  - as vazoes de contorno sao rateadas por area de drenagem (NUAREAMONT).

Formatos validados contra os projetos-exemplo oficiais do HEC-RAS:
  series posicionais em colunas de 8 caracteres (10/linha), Boundary
  Location em 6 campos com padding, Use Restart + Initial Flow Loc,
  juncao por NOME (Junct Name / Up River,Reach / Dn River,Reach / Junc L&A).

Uso:   python gerar_rede_hecras.py
Depois: python run_hecras.py Itajai_Rede
"""
import datetime
import os
import unicodedata
import numpy as np
import geopandas as gpd
import rasterio
from shapely.geometry import LineString, Point, shape
from shapely.ops import substring
from pyproj import Transformer

# ------------------------------------------------------------------ PARAMETROS
PROJECT   = "Itajai_Rede"      # recebe sufixo _<EVENTO> se houver evento
GEOJSON   = "rios_itajai.geojson"
DEM       = "dem_itajai.tif"
# MDT do SIG-SC a 1 m (terreno, nao superficie). Cobre 100% dos rios
# modelados no mesmo CRS. Contra o Copernicus GLO-30, nas encostas a
# diferenca chega a 15 m -- que e a COPA DA MATA sendo tratada como terreno.
USAR_SIGSC = True
# Corta a secao onde ela reencontraria o rio. Uma secao 1D tem de cruzar o
# canal UMA vez; em trecho meandrante 17% das cutlines cruzavam 2 ou 3 vezes
# e a mesma agua era contada vezes repetidas.
CORTAR_NO_REENCONTRO = True
UTM_EPSG  = 31982              # SIRGAS 2000 / UTM 22S

SPACING   = 1000.0             # espacamento entre secoes (m) no trecho PLANO
SPACING_MIN = 150.0            # espacamento nas gargantas. O Acu cai 195 m em
                               # 13 km entre Lontras e Ibirama (Salto Pilao) --
                               # confirmado no terreno SIG-SC de 1 m, nao e
                               # ruido de DEM. A 1 km de espacamento sao 8 m de
                               # queda ENTRE SECOES VIZINHAS: o solver unsteady
                               # falha no primeiro passo ("Solution Solver
                               # Failed" em Itajai_Acu R1 150611.9). O criterio
                               # usual para regime nao permanente e
                               # dx <~ 0,15*D/S; com D~4 m e S=0,008 da ~75 m.
DECL_PLANO  = 0.0010           # ate aqui, espacamento SPACING
DECL_INGREME = 0.0060          # daqui pra cima, espacamento SPACING_MIN
HALFWIDTH = 2500.0             # meia-largura MAXIMA da secao (m). A secao TEM
                               # de conter a cheia: com 700 m o HEC-RAS avisava
                               # "Extrapolated above Cross Section Table" nos
                               # ultimos 27 km do Acu (agua ate 2,45 m acima do
                               # topo da secao), assumia paredes verticais e
                               # estrangulava a conducao no baixo vale.
NPTS      = 280                # pontos por secao (limite HEC-RAS: 450).
                               # Denso o bastante p/ resolver o canal
                               # dentro de uma secao larga o suficiente
                               # para conter a cheia.
SMOOTH    = 250.0              # janela p/ direcao local (evita cutlines cruzadas)
MIN_SLOPE = 1e-4               # declividade minima imposta ao talvegue
MAX_SLOPE = 0.008              # declividade maxima (0,8%): acima disso o
                               # escoamento fica supercritico e o solver 1D
                               # nao converge. As bacias reais modeladas neste
                               # vale usam 0,15%-1,0% (reaches_bacia_completa)
MIN_AREA  = 200.0              # area de drenagem minima (km2) p/ iniciar um
                               # rio: abaixo disso e torrente de montanha
# n de Manning do CANAL por zona, lido de
# itajai_flood_model/data/reaches_bacia_completa.csv. As zonas sao delimitadas
# por uma cidade de referencia, cuja estaca e calculada projetando as
# coordenadas sobre o eixo do rio (as distancias do CSV nao batem com a
# geometria da ANA: o Acu tem 153 km la e 187,6 km aqui).
#   rio -> (cidade, lat, lon, n acima da cidade, n abaixo)
MANNING = {
    "Itajai_Acu":   ("Blumenau", -26.9180, -49.0660, 0.035, 0.030),
    "Itajai_Mirim": ("Brusque",  -27.0980, -48.9120, 0.038, 0.028),
}
N_CANAL_PADRAO   = 0.035
RAZAO_PLANICIE   = 1.8         # n da planicie = n do canal x isto
# --- escavacao da calha (bathymetry) -------------------------------------
# O DEM e de SUPERFICIE: sobre o rio ele mede a lamina d'agua, nao o leito.
CAVAR_CANAL = True
CANAL_KH  = 0.277              # h = KH * A^EH   -> ~8,0 m na foz (14.871 km2)
CANAL_EH  = 0.35               #                    ~1,8 m com 200 km2
CANAL_KW  = 5.0                # w = KW * A^EW   -> ~233 m na foz
CANAL_EW  = 0.40               #                    ~35 m com 200 km2
BANK_H    = 3.0                # altura acima do talvegue p/ definir a margem
# Contorno de jusante. O trecho final e ESTUARIO sob mare, nao canal com
# declividade: com Friction Slope numa foz quase plana a profundidade normal
# fica alta, gera remanso e o baixo vale estoca agua demais.
MARE          = True
MARE_MEDIA    = 0.30           # nivel medio do mar (m)
MARE_AMPLITUDE= 0.50           # semi-amplitude (m) -> variacao de ~1,0 m
MARE_PERIODO  = 12.42          # componente M2 (h)
DS_SLOPE  = 0.0005             # declividade (usada so se MARE = False)
# EVENTO: usa a chuva REAL observada (via hidrologia_evento) no lugar do
# hidrograma triangular sintetico. None = sintetico.
EVENTO      = None       # "2008"/"2011"/... usa chuva real (ver nota abaixo)
BARRAGENS   = True             # False = cenario "sem obras" (comportas abertas)
DATA_INICIO = None             # datetime do inicio da simulacao. Vem da data
                               # REAL do evento (primeiro registro do arquivo
                               # de chuva observada); None = cenario sintetico.
NHORAS    = 97                 # ordinatas horarias (h = 0..96). 48 h nao
                               # bastam: o pico leva mais de 48 h para
                               # percorrer os 187 km ate a foz, e a
                               # simulacao truncava antes de ele chegar.

# rios modelados: chave -> (padrao de nome na ANA, nome HEC-RAS)
RIOS = {
    "acu":      ("Itajaí-açu",                "Itajai_Acu"),
    "sul":      ("Itajaí do Sul",             "Itajai_Sul"),
    "oeste":    ("Itajaí do Oeste",           "Itajai_Oeste"),
    "norte":    ("Itajaí do Norte|Hercílio",  "Itajai_Norte"),
    "benedito": ("Benedito",                  "Rio_Benedito"),
    "mirim":    ("Itajaí-mirim",              "Itajai_Mirim"),
    # Afluentes de 2a ordem. Nenhum deles desagua no Acu, exceto Luis Alves e
    # do Testo: o Trombudo, o Taio e o das Pombas entram no Oeste; o Krauel e o
    # Iraputa no Norte; o dos Cedros no Benedito. Sao os oito rios acima de
    # 240 km2 que faltavam na rede -- juntos, 4.000 km2 que ate aqui entravam
    # como vazao incremental distribuida, sem geometria propria (portanto sem
    # mancha de inundacao em Luis Alves, Taio, Trombudo Central ou Mirim Doce).
    "luisalves": ("Luís Alves",   "Rio_Luis_Alves"),
    "trombudo":  ("Trombudo",     "Rio_Trombudo"),
    "taio":      ("Taió",         "Rio_Taio"),
    "pombas":    ("das Pombas",   "Rio_das_Pombas"),
    "cedros":    ("dos Cedros",   "Rio_dos_Cedros"),
    "krauel":    ("Krauel",       "Rio_Krauel"),
    "iraputa":   ("Iraputã",      "Rio_Iraputa"),
    "testo":     ("do Testo",     "Rio_do_Testo"),
}
# Canais retificados que SUBSTITUEM o curso natural. A base da ANA traz o
# leito antigo meandrante; onde houve retificacao, o rio real corre pelo
# canal. No Itajai-Mirim isso troca 19,23 km de meandros por 7,55 km de
# canal -- muda declividade e tempo de viagem de forma relevante.
CANAIS = {"mirim": os.path.join("dados_estruturas", "canal_itajai_mirim.geojson")}
MAIN = "acu"
# ESCOPO: quais afluentes entram na rede. Reduzir o escopo diminui o numero
# de juncoes (cada uma e um ponto potencial de instabilidade) e permite
# validar a cadeia inteira ate a mancha antes de voltar os demais rios.
#   completo -> ["sul","oeste","norte","benedito","mirim"]
#   reduzido -> ["mirim"]   (Acu + Mirim, 1 juncao)
ESCOPO = ["sul", "oeste", "norte", "benedito", "mirim"]
# Os oito afluentes de 2a ordem (Luis Alves, Trombudo, Taio, das Pombas, dos
# Cedros, Krauel, Iraputa, do Testo) ficam FORA por ora: a area deles continua
# entrando como vazao incremental distribuida, que e o que ja acontecia antes.
# Basta acrescenta-los a lista acima para voltarem -- o codigo da arvore ja
# resolve sozinho em que rio cada um desagua.
# Area de drenagem que entra pela CABECEIRA do Acu. No escopo reduzido os
# afluentes ausentes sao somados aqui para que as vazoes a jusante (Blumenau,
# Itajai) fiquem na ordem de grandeza correta.
LATERAIS = []
# Afluentes injetados como VAZAO LATERAL na estaca da confluencia, sem virar
# trecho nem juncao. Motivo: uma varredura mostrou que 1 juncao roda limpa
# (0 falhas, erro de volume 0,02-0,18%) mas 2 ou mais divergem -- inclusive
# combinando afluentes que funcionam isolados. A vazao lateral entrega a
# agua no ponto certo sem adicionar juncao.
INCREMENTAL = True             # distribui a area de drenagem do PROPRIO Acu
                               # (a que nao pertence a nenhum afluente nomeado)
                               # como Uniform Lateral Inflow ao longo da calha.
                               # Sem isso faltam ~3.300 km2 de contribuicao e a
                               # vazao em Blumenau sai ~25% baixa.
AREA_CABECEIRA_ACU = 0.0    # Sul + Oeste (km2), que formam o Acu
                               # na confluencia de Rio do Sul
# nomes das juncoes por ordem de km ao longo do Acu
NOME_JUNCAO = {0: "Rio_do_Sul", 1: "Ibirama", 2: "Indaial", 3: "Itajai"}

# cheia de referencia (2008: 5.700 m3/s em Itajai) rateada por area
Q_REF_FOZ  = 5700.0
EDIT_TIME  = "Node Last Edited Time= Aug/15/2026 00:00:00"


# ------------------------------------------------------------------- UTILIDADES
def p16(s):
    return f"{str(s)[:16]:<16}"


def f8(v):
    return f"{v:8.2f}"


def serie8(vals):
    """Serie posicional: colunas de 8 caracteres, 10 valores por linha."""
    return "\n".join("".join(f8(v) for v in vals[i:i + 10])
                     for i in range(0, len(vals), 10))


WARMUP = 8                     # horas de vazao constante antes da cheia:
                               # sem aquecimento o solver diverge nos
                               # primeiros minutos (partida a frio nas juncoes)


def mare(n=None):
    """Nivel do mar na foz: onda semidiurna M2."""
    n = n or NHORAS          # NAO usar default: NHORAS muda com o evento
    return [MARE_MEDIA + MARE_AMPLITUDE * np.sin(2 * np.pi * h / MARE_PERIODO)
            for h in range(n)]


def hidrograma(pico, base=None, n=None, tp=26, te=46):
    n = n or NHORAS          # idem: o default congelaria o valor antigo
    base = base if base is not None else max(pico * 0.15, 20.0)
    v = []
    for h in range(n):
        if h <= WARMUP:
            q = base
        elif h <= tp:
            q = base + (pico - base) * ((h - WARMUP) / (tp - WARMUP))
        elif h <= te:
            q = pico - (pico - base) * ((h - tp) / (te - tp))
        else:
            q = base
        v.append(q)
    return v


class Dem:
    """Amostrador do DEM. Le a banda inteira em memoria: ds.sample() derruba
    o processo em rasterio 1.5.x (crash nativo, sem traceback)."""

    def __init__(self, path):
        self.ds = rasterio.open(path)
        self.arr = self.ds.read(1)
        self.rows, self.cols = self.arr.shape
        self.nodata = self.ds.nodata
        self.inv = ~self.ds.transform
        self.tr = Transformer.from_crs(UTM_EPSG, self.ds.crs.to_epsg(),
                                       always_xy=True)

    def sample(self, xs, ys):
        lon, lat = self.tr.transform(np.asarray(xs), np.asarray(ys))
        a, b, c, d, e, f = (self.inv.a, self.inv.b, self.inv.c,
                            self.inv.d, self.inv.e, self.inv.f)
        col = np.floor(a * lon + b * lat + c).astype(int)
        row = np.floor(d * lon + e * lat + f).astype(int)
        ok = ((row >= 0) & (row < self.rows) & (col >= 0) & (col < self.cols))
        out = np.full(lon.shape, np.nan)
        out[ok] = self.arr[row[ok], col[ok]]
        if self.nodata is not None:
            out[out == self.nodata] = np.nan
        out[out < -500] = np.nan
        return out


# ------------------------------------------------------------------- TOPOLOGIA
def aplicar_canal(linha, caminho):
    """Substitui o trecho do eixo entre as pontas do canal pelo proprio canal."""
    import json as _json
    from shapely.ops import linemerge as _lm
    if not os.path.exists(caminho):
        print(f"      ! canal nao encontrado: {caminho}")
        return linha
    feats = _json.load(open(caminho, encoding="utf-8"))["features"]
    cu = gpd.GeoSeries([shape(f["geometry"]) for f in feats],
                       crs=4326).to_crs(UTM_EPSG)
    canal = _lm(list(cu))
    if canal.geom_type != "LineString":
        canal = max(canal.geoms, key=lambda g: g.length)
    a, b = Point(canal.coords[0]), Point(canal.coords[-1])
    sa, sb = linha.project(a), linha.project(b)
    if sa > sb:                       # canal digitalizado ao contrario
        canal = LineString(list(canal.coords)[::-1])
        sa, sb = sb, sa
    if sb - sa < canal.length * 0.5:  # nao encurta: nao e retificacao
        return linha
    montante = list(substring(linha, 0, sa).coords)
    jusante = list(substring(linha, sb, linha.length).coords)
    nova = LineString(montante + list(canal.coords) + jusante)
    print(f"      canal retificado: {(sb-sa)/1000:.2f} km de meandros -> "
          f"{canal.length/1000:.2f} km de canal "
          f"(eixo {linha.length/1000:.1f} -> {nova.length/1000:.1f} km)")
    return nova


def montar_rede():
    g = gpd.read_file(GEOJSON).to_crs(UTM_EPSG)
    g["NORIOCOMP"] = g["NORIOCOMP"].astype(str)
    by = {int(r.COTRECHO): r for r in g.itertuples()}
    pred = {}
    for r in g.itertuples():
        pred.setdefault(int(r.NUTRJUS), []).append(int(r.COTRECHO))

    def cadeia(pat):
        """Calha principal do rio: da foz para montante, sempre pelo ramo de
        maior area de drenagem."""
        sub = g[g["NORIOCOMP"].str.contains(pat, case=False, na=False)]
        out = int(sub.loc[sub["NUAREAMONT"].idxmax(), "COTRECHO"])
        ch = [out]
        while True:
            # str() obrigatorio: nem todo trecho da BHO tem nome, e o valor
            # vem como NaN (float). Sem isso a busca quebra com AttributeError
            # ao subir por qualquer rio cujo caminho passe por trecho sem nome.
            ups = [c for c in pred.get(ch[-1], [])
                   if c in by and any(a.lower() in str(by[c].NORIOCOMP).lower()
                                      for a in pat.split("|"))]
            if not ups:
                break
            melhor = max(ups, key=lambda c: float(by[c].NUAREAMONT or 0))
            if float(by[melhor].NUAREAMONT or 0) < MIN_AREA:
                break          # torrente de cabeceira: nao modela
            ch.append(melhor)
        return ch[::-1]                       # cabeceira -> foz

    def eixo(ch):
        pts = []
        for c in ch:
            gg = by[c].geometry
            for l in ([gg] if gg.geom_type == "LineString" else list(gg.geoms)):
                cc = list(l.coords)
                if pts and (Point(pts[-1]).distance(Point(cc[0])) >
                            Point(pts[-1]).distance(Point(cc[-1]))):
                    cc = cc[::-1]
                pts += cc if not pts else cc[1:]
        return LineString(pts)

    rede = {}
    for k, (pat, nome) in RIOS.items():
        if k != MAIN and k not in ESCOPO and k not in LATERAIS:
            continue
        ch = cadeia(pat)
        ln = eixo(ch)
        if k in CANAIS:
            ln = aplicar_canal(ln, CANAIS[k])
        rede[k] = {"nome": nome, "linha": ln,
                   "area": float(by[ch[-1]].NUAREAMONT)}
    return rede


# --------------------------------------------------------------------- SECOES
DS_ILHA = 400.0               # ate esta diferenca de estaca, o outro braco
                               # ainda e "o mesmo lugar do rio" -> ilha


def _ate_reencontro(p, rx, ry, hw, eixo, s, linha, folga=40.0):
    """Ate onde a semi-secao pode ir sem reencontrar o rio.

    Nem todo reencontro e igual, e tratar os dois do mesmo jeito custava caro:

    CURVA (meandro) -- o raio alcanca o MESMO rio muitas estacas adiante ou
    atras. Aquela agua ja e contabilizada na secao de la; incluir aqui conta o
    mesmo escoamento duas vezes e ainda cria cutlines cruzadas. Aqui a secao
    para: e a barreira.

    ILHA / BRACO SECUNDARIO -- o outro braco esta na MESMA altura do rio
    (diferenca de estaca pequena) e corre no MESMO sentido. Ele conduz de
    verdade: e por ele que parte da cheia passa, como no leito antigo do
    Itajai-Mirim. Parar ali amputa a secao justamente onde a agua vai.
    Entao a secao ATRAVESSA e segue.

    A distincao e por estaca ao longo do eixo e pelo produto escalar das
    direcoes locais -- nao por distancia em linha reta, que confunde as duas.
    """
    raio = LineString([(p.x, p.y), (p.x + hw * rx, p.y + hw * ry)])
    it = raio.intersection(eixo)
    if it.is_empty:
        return hw
    pts = [it] if it.geom_type == "Point" else list(getattr(it, "geoms", []))
    pts = [q for q in pts if getattr(q, "geom_type", "") == "Point"]

    def direcao(est):
        a = linha.interpolate(max(0.0, est - SMOOTH))
        b = linha.interpolate(min(linha.length, est + SMOOTH))
        v = np.array([b.x - a.x, b.y - a.y])
        return v / (np.linalg.norm(v) or 1.0)

    u_aqui = direcao(s)
    corte = []
    for q in pts:
        d = float(np.hypot(q.x - p.x, q.y - p.y))
        if d <= folga:                     # o proprio cruzamento no eixo
            continue
        s_la = linha.project(q)
        mesmo_lugar = abs(s_la - s) <= DS_ILHA
        mesmo_sentido = float(np.dot(u_aqui, direcao(s_la))) > 0.5
        if mesmo_lugar and mesmo_sentido:
            continue                       # ilha/braco: conduz, pode atravessar
        corte.append(d)                    # meandro: barreira
    if not corte:
        return hw
    return max(min(corte) - folga, 120.0)


def cortar(linha, s, dem, hw=HALFWIDTH, eixo=None):
    """Secao perpendicular ao eixo na posicao s, amostrada no DEM.
    A direcao usa uma janela de +-SMOOTH m: com +-1 m as cutlines se cruzam
    nas curvas (o RAS avisa 'edge lines have self intersections')."""
    a = linha.interpolate(max(0.0, s - SMOOTH))
    b = linha.interpolate(min(linha.length, s + SMOOTH))
    tx, ty = b.x - a.x, b.y - a.y
    n = np.hypot(tx, ty) or 1.0
    rx, ry = ty / n, -tx / n                  # normal a direita
    p = linha.interpolate(s)
    hw_e = hw_d = hw
    if CORTAR_NO_REENCONTRO and eixo is not None:
        hw_e = _ate_reencontro(p, -rx, -ry, hw, eixo, s, linha)
        hw_d = _ate_reencontro(p,  rx,  ry, hw, eixo, s, linha)
    off = np.concatenate([np.linspace(-hw_e, 0, NPTS // 2, endpoint=False),
                          np.linspace(0, hw_d, NPTS - NPTS // 2)])
    z = dem.sample(p.x + off * rx, p.y + off * ry)
    if np.isnan(z).all():
        return None
    if np.isnan(z).any():                     # tapa buracos por interpolacao
        ok = ~np.isnan(z)
        z = np.interp(np.arange(len(z)), np.flatnonzero(ok), z[ok])
    sta = off + hw_e                          # 0 .. hw_e+hw_d
    cut = (p.x - hw_e * rx, p.y - hw_e * ry,
           p.x + hw_d * rx, p.y + hw_d * ry)
    return np.round(sta, 2), z, cut


def canal_geometria(area_km2):
    """Profundidade e largura da calha por geometria hidraulica.

    O DEM Copernicus e modelo de SUPERFICIE: sobre o rio ele registra a
    lamina d'agua no instante da aquisicao, nao o leito. Sem escavar a calha,
    toda a secao de escoamento abaixo dessa lamina fica faltando -- a area
    molhada e subestimada e a cota calculada nao se apoia em nada.

    Usa-se a forma classica de Leopold & Maddock, com a area de drenagem como
    substituto da vazao de margens plenas:
        h = CANAL_KH * A^CANAL_EH     (m)
        w = CANAL_KW * A^CANAL_EW     (m)
    Os coeficientes foram fixados para reproduzir a ordem de grandeza
    conhecida do Itajai: ~8 m de profundidade e ~230 m de largura na foz
    (canal do porto de Itajai, dragado), caindo a ~1,8 m e ~35 m nos
    afluentes de 200 km2.
    """
    a = max(float(area_km2), 1.0)
    h = CANAL_KH * a ** CANAL_EH
    w = CANAL_KW * a ** CANAL_EW
    return float(h), float(w)


def indice_eixo(sta, z, janela):
    """Indice do talvegue PROXIMO AO EIXO do rio, nao o minimo global.

    O eixo da cutline esta no offset 0, que por construcao e sta[len//2].
    Usar o minimo GLOBAL da secao poe a calha no lugar errado sempre que o
    corte atravessa outro canal mais fundo. E o que acontecia no Itajai-Mirim:
    a secao cruza o LEITO ANTIGO, que no DEM e mais baixo que o canal
    retificado, entao o "talvegue" caia la. O RAS Mapper mostrou o efeito --
    as bank lines desenhadas sobre o leito antigo, e no Acu cruzando-se em
    estrela, porque a margem saltava de um canal para o outro entre secoes
    vizinhas. Junto com elas iam a zona de Manning do canal e a propria
    escavacao, ou seja, a calha condutora ficava fora do caminho da agua.

    A busca fica restrita a uma janela em volta do eixo: pega o talvegue real
    (o tracado da ANA nao passa exatamente no fundo) sem pular de canal.
    """
    i_eixo = len(sta) // 2
    m = np.abs(np.asarray(sta) - sta[i_eixo]) <= janela
    idx = np.flatnonzero(m & np.isfinite(z))
    if not len(idx):
        return i_eixo
    return int(idx[np.nanargmin(np.asarray(z)[idx])])


def zt(d):
    """Cota do talvegue da secao -- SEMPRE a do canal junto ao eixo.

    Nao pode ser z.min(): onde a secao atravessa o leito antigo do Mirim (ou
    um meandro do proprio rio), o minimo global fica noutro canal. Ate aqui o
    condicionamento do perfil longitudinal usava esse minimo enquanto a calha
    condutora, as margens e a escavacao ja apontavam para o eixo -- perfil e
    conducao descreviam canais diferentes, e a simulacao caiu de 30 para 2
    passos. O indice e guardado no corte e sobrevive aos deslocamentos que o
    condicionamento aplica a secao inteira.
    """
    i = d.get("i_thal")
    return float(d["z"][i]) if i is not None else float(np.nanmin(d["z"]))


def cavar_canal(sta, z, area_km2):
    """Rebaixa a calha no talvegue junto ao EIXO do rio.

    Escava um trapezio de largura w e profundidade h centrado ali, com
    taludes para nao criar degrau vertical.
    """
    if not CAVAR_CANAL:
        return z
    h, w = canal_geometria(area_km2)
    z = np.asarray(z, dtype=float).copy()
    i0 = indice_eixo(sta, z, max(w, 150.0))
    centro = sta[i0]
    d = np.abs(sta - centro)
    meia = w / 2.0
    talude = max(w * 0.25, 30.0)          # transicao suave ate o terreno
    # perfil de escavacao: 1 no fundo, decaindo linearmente ate 0 no talude
    frac = np.clip(1.0 - (d - meia) / talude, 0.0, 1.0)
    frac[d <= meia] = 1.0
    return z - h * frac


def margens(sta, z, h_canal=0.0, altura_margem=BANK_H, larg_canal=150.0):
    """Margens topograficas: do talvegue ate o TOPO DA CALHA + BANK_H.

    Com a calha escavada, medir BANK_H acima do talvegue coloca a margem
    DENTRO do canal -- o modelo passa a achar que tudo extravasa. A margem
    real e o topo da calha (o terreno original, h_canal acima do fundo
    escavado); BANK_H e a folga a partir dali.
    O valor DEVE coincidir com um sta da tabela, na mesma precisao (.2f).

    Parte do talvegue JUNTO AO EIXO (ver indice_eixo): partir do minimo global
    punha as margens em volta do leito antigo quando a secao o atravessava."""
    i = indice_eixo(sta, z, max(larg_canal, 150.0))
    lim = z[i] + h_canal + altura_margem
    li = i
    while li > 0 and z[li] < lim:
        li -= 1
    ri = i
    while ri < len(z) - 1 and z[ri] < lim:
        ri += 1
    li = min(max(li, 1), len(sta) - 3)
    ri = max(min(ri, len(sta) - 2), li + 1)
    return round(float(sta[li]), 2), round(float(sta[ri]), 2)


def largura(area_km2):
    """Meia-largura da secao proporcional ao porte do rio. Usar 700 m num
    afluente de montanha desperdicava quase todos os pontos na encosta e
    deixava o canal com 1-3 pontos.

    ALARGAR ISTO NAO E TRIVIAL. O leito antigo do Itajai-Mirim corre a 1.528 m
    do canal retificado (mediana medida entre os dois tracados) e fica fora dos
    737 m que esta escala da ao Mirim -- ou seja, o modelo ignora esse caminho
    da agua. Mas subir o coeficiente para 440 derrubou a simulacao de 30 para
    2 passos: a secao larga passa a atravessar meandros do proprio rio, e ai o
    escoamento e contado duas vezes. Alargar so vai funcionar junto com area de
    escoamento inefetivo, que precisa do formato #Ineffective do .g01.
    """
    return float(np.clip(180.0 * np.sqrt(max(area_km2, 1.0) / 100.0), 500.0, HALFWIDTH))


def estacoes(linha, dem):
    """Estacas ao longo do eixo, com espacamento adaptado a declividade.

    Faz uma varredura barata do talvegue (minimo numa janela de 100 m em
    volta do eixo, sem montar cutline) so para saber ONDE o rio e ingreme,
    e adensa as secoes ali. Sai de SPACING no vale plano para SPACING_MIN na
    garganta, interpolando linearmente entre DECL_PLANO e DECL_INGREME.
    """
    L = linha.length
    passo = min(SPACING / 4.0, 250.0)
    d = np.arange(0.0, L + passo, passo)
    P = [linha.interpolate(float(x)) for x in d]
    xs, ys = [], []
    for p in P:
        for a in (-50.0, 0.0, 50.0):
            for b in (-50.0, 0.0, 50.0):
                xs.append(p.x + a)
                ys.append(p.y + b)
    z = dem.sample(np.array(xs), np.array(ys)).reshape(len(P), 9)
    with np.errstate(all="ignore"):
        zb = np.nanmin(z, axis=1)
    ok = np.isfinite(zb)
    if ok.sum() < 3:
        return list(np.arange(0.0, L, SPACING))
    zb = np.interp(d, d[ok], zb[ok])
    S = np.abs(np.gradient(zb, d))
    S = np.convolve(S, np.ones(5) / 5.0, mode="same")     # tira o ruido do DEM
    f = np.clip((S - DECL_PLANO) / (DECL_INGREME - DECL_PLANO), 0.0, 1.0)
    dx = SPACING + (SPACING_MIN - SPACING) * f
    ss, s = [0.0], 0.0
    while s < L:
        s += float(np.interp(s, d, dx))
        if s < L:
            ss.append(s)
    return ss


def secoes(linha, dem, rs0, hw=HALFWIDTH, area=None):
    # hw pode ser um numero OU uma funcao da fracao percorrida (0=montante,
    # 1=jusante). A largura TEM de crescer rio abaixo: usar a largura da foz
    # na cabeceira do Acu (4.390 m para 120 m3/s) cria uma lamina de papel,
    # o trecho vira lago estagnado (vazoes negativas) e o solver instabiliza.
    """Corta as secoes de um trecho. rs0 = RS do extremo de JUSANTE."""
    L = linha.length
    ss = estacoes(linha, dem)
    # NAO cria secao colada em NENHUM dos extremos: uma secao em cima da
    # juncao conflita com o comprimento declarado em 'Junc L&A' e trava o
    # solver. O extremo de jusante ja era protegido; faltava o de MONTANTE --
    # o Acu nasce da juncao de Rio do Sul e tinha a primeira secao em RS
    # 187611,90 num rio de 187.611,9 m, ou seja, exatamente sobre ela.
    recuo = SPACING_MIN * 0.5
    ss = [s for s in ss if s >= recuo] or [recuo]
    if ss[0] > recuo * 1.5:
        ss.insert(0, recuo)
    if L - ss[-1] > SPACING_MIN * 0.6:
        ss.append(L - recuo)
    xs = []
    for s in ss:
        hw_s = hw(s / max(L, 1.0)) if callable(hw) else hw
        r = cortar(linha, s, dem, hw_s, linha)
        if r is None:
            continue
        sta, z, cut = r
        a_km2 = area(s / max(L, 1.0)) if callable(area) else (area or 1000.0)
        z = cavar_canal(sta, z, a_km2)     # escava a calha no talvegue
        _, w_c = canal_geometria(a_km2)
        # UM canal por secao. Depois de escavar no eixo, qualquer ponto ainda
        # MAIS FUNDO e outro canal que o corte atravessou -- o leito antigo do
        # Mirim, um meandro do proprio rio. Num modelo 1D isso nao e um segundo
        # caminho de escoamento: e um poco mais fundo que o canal principal
        # dentro da MESMA secao, e a conducao calculada em cima disso nao tem
        # sentido (a simulacao caiu de 30 para 2 passos assim que a escavacao
        # passou a ser no eixo). Subir esses pontos ate o fundo da calha deixa
        # uma calha so; o resto da secao nao e tocado.
        i_t = indice_eixo(sta, z, max(w_c, 150.0))
        z = np.maximum(z, z[i_t])
        z[i_t] = min(z[i_t], float(np.nanmin(z)))
        xs.append({"rs": round(rs0 + (L - s), 2), "sta": sta, "z": z,
                   "cut": cut, "area_km2": a_km2,
                   "i_thal": indice_eixo(sta, z, max(w_c, 150.0))})
    xs.sort(key=lambda d: -d["rs"])           # montante -> jusante
    # remove RS repetido (o RAS exige unicidade)
    fin, visto = [], set()
    for d in xs:
        if d["rs"] in visto:
            continue
        visto.add(d["rs"])
        fin.append(d)
    return fin


def ajustar_talvegue(d, delta):
    """Move o TALVEGUE em 'delta', deixando a planicie onde o DEM a pos.

    O condicionamento deslocava a secao INTEIRA (z = z - delta) para impor o
    perfil longitudinal. Como delta chega a 5,8 m no baixo Itajai, a planicie e
    as margens desciam junto: a margem em Itajai ficava a -2,75 m e a de Ilhota
    a -0,60 m, ou seja, abaixo do nivel do mar -- e isso nao veio do relevo,
    veio do deslocamento. Aqui so a calha se move, com o mesmo perfil
    trapezoidal da escavacao, e o terreno fora dela fica intacto.
    """
    if abs(delta) < 1e-6:
        return
    sta, z = d["sta"], d["z"]
    _, w = canal_geometria(d.get("area_km2", 1000.0))
    centro = sta[d.get("i_thal", len(sta) // 2)]
    dist = np.abs(sta - centro)
    meia = w / 2.0
    talude = max(w * 0.25, 30.0)
    frac = np.clip(1.0 - (dist - meia) / talude, 0.0, 1.0)
    frac[dist <= meia] = 1.0
    z = z + delta * frac
    if delta > 0:            # ao SUBIR a calha, nao pode passar do terreno
        z = np.minimum(z, d["z"] + np.maximum(delta, 0.0) * frac)
    d["z"] = z


def condicionar(xs, rotulo=""):
    """Prepara o talvegue para o solver unsteady, em tres passos:

    1. apara a CABECEIRA enquanto a declividade local passar de MAX_SLOPE
       (torrentes de montanha: o Benedito chega a 9,4% no DEM bruto);
    2. impoe decrescimento monotonico rio abaixo (tira as contrapendentes);
    3. limita a declividade a MAX_SLOPE ancorando no extremo de JUSANTE
       (rebaixa o lado de montante), preservando a cota da foz/confluencia.
    """
    if len(xs) < 3:
        return xs
    # -- 1. apara cabeceira ingreme
    corte = 0
    while corte < len(xs) - 2:
        dx = xs[corte]["rs"] - xs[corte + 1]["rs"]
        dz = zt(xs[corte]) - zt(xs[corte + 1])
        if dx > 0 and abs(dz) / dx > MAX_SLOPE:
            corte += 1
        else:
            break
    if corte:
        print(f"        {rotulo}: aparadas {corte} secoes de cabeceira "
              f"({(xs[0]['rs']-xs[corte]['rs'])/1000:.1f} km, decl > "
              f"{100*MAX_SLOPE:.1f}%)")
        xs = xs[corte:]
    # -- 2. monotonico rio abaixo
    for i in range(1, len(xs)):
        dx = xs[i - 1]["rs"] - xs[i]["rs"]
        teto = zt(xs[i - 1]) - MIN_SLOPE * dx
        atual = zt(xs[i])
        if atual > teto:
            ajustar_talvegue(xs[i], teto - atual)
    # -- 3. limita declividade ancorando a JUSANTE
    for i in range(len(xs) - 2, -1, -1):
        dx = xs[i]["rs"] - xs[i + 1]["rs"]
        lim = zt(xs[i + 1]) + MAX_SLOPE * dx
        atual = zt(xs[i])
        if atual > lim:
            ajustar_talvegue(xs[i], lim - atual)
    return xs


# ---------------------------------------------------------------- ESCRITA RAS
def bl(river, reach, rs):
    return (f"Boundary Location={p16(river)},{p16(reach)},{str(rs)[:8]:<8}"
            f",        ,                ,                ")


def escrever(trechos, juncoes):
    # 'Viewing Rectangle' e a extensao geografica da geometria, e estava
    # gravado como o placeholder "0 , 1 , 1 , 0" -- um quadrado de 1 m na
    # origem. E dele que o HEC-RAS deriva o atributo Extents do .g01.hdf, e e
    # o Extents que o RAS Mapper usa para saber ONDE desenhar. Com [0,1,0,1] a
    # janela do RAS Mapper abria vazia mesmo com a geometria marcada: ele
    # procurava a bacia inteira dentro de um metro quadrado na origem.
    # A ordem dos quatro campos e xmin, xmax, ymax, ymin.
    xs, ys = [], []
    for t in trechos:
        for x, y in t["linha"].coords:
            xs.append(x); ys.append(y)
        for d in t["xs"]:
            cu = d["cut"]
            xs += [cu[0], cu[2]]
            ys += [cu[1], cu[3]]
    folga = 0.02 * max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
    x0, x1 = min(xs) - folga, max(xs) + folga
    y0, y1 = min(ys) - folga, max(ys) + folga
    g = [f"Geom Title={PROJECT} - rede real ANA + relevo DEM",
         "Program Version=7.01",
         f"Viewing Rectangle= {x0:.6f} , {x1:.6f} , {y1:.6f} , {y0:.6f} "]
    wkt = ('PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",'
           'DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",'
           'SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],'
           'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
           'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],'
           'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
           'PARAMETER["false_northing",10000000],UNIT["metre",1]]')
    g.append(f"Spatial Reference System={wkt}")
    g.append("")

    for j in juncoes:
        g.append(f"Junct Name={p16(j['nome'])}")
        g.append(f"Junct Desc={j['desc']}, 0 , 0 , 0 ,0")
        g.append(f"Junct X Y & Text X Y={j['x']:.2f},{j['y']:.2f},"
                 f"{j['x'] + 800:.2f},{j['y'] + 800:.2f}")
        for r, rc in j["up"]:
            g.append(f"Up River,Reach={p16(r)},{p16(rc)}")
        g.append(f"Dn River,Reach={p16(j['dn'][0])},{p16(j['dn'][1])}")
        for d in j.get("dists") or [500.0] * len(j["up"]):
            g.append(f"Junc L&A={d:.2f},0")
        g.append("")

    for t in trechos:
        g.append(f"River Reach={p16(t['rio'])},{p16(t['reach'])}")
        # Eixo em RESOLUCAO CHEIA. Decimar para 400 pontos cortava as curvas
        # e o rio aparecia poligonal no RAS Mapper e no app -- as cutlines
        # sempre usaram a linha completa, mas o traçado desenhado nao.
        c = list(t["linha"].coords)
        g.append(f"Reach XY= {len(c)} ")
        for i in range(0, len(c), 2):
            par = c[i:i + 2]
            g.append("".join(f"{x:16.4f}{y:16.4f}" for x, y in par))
        g.append("")
        xs = t["xs"]
        for i, d in enumerate(xs):
            dx = (round(d["rs"] - xs[i + 1]["rs"], 2)
                  if i < len(xs) - 1 else 0.0)
            g.append(f"Type RM Length L Ch R = 1 ,{d['rs']:.2f},"
                     f"{dx:.2f},{dx:.2f},{dx:.2f}")
            cu = d["cut"]
            g.append("XS GIS Cut Line=2")
            g.append(f"{cu[0]:16.4f}{cu[1]:16.4f}{cu[2]:16.4f}{cu[3]:16.4f}")
            g.append(EDIT_TIME)
            sta, z = d["sta"], d["z"]
            g.append(f"#Sta/Elev= {len(sta)} ")
            par = [v for p in zip(sta, z) for v in p]
            g += [ "".join(f8(v) for v in par[i:i + 10])
                   for i in range(0, len(par), 10) ]
            h_c, w_c = canal_geometria(d.get('area_km2', 1000.0))
            lb, rb = margens(sta, z, h_c, larg_canal=w_c)
            n_ch = d.get("n", N_CANAL_PADRAO)
            n_fp = round(n_ch * RAZAO_PLANICIE, 3)
            g.append("#Mann= 3 ,-1,0")
            g.append("".join(f"{v:>8}" for v in
                             [f"{sta[0]:.2f}", f"{n_fp:.3f}", "0",
                              f"{lb:.2f}", f"{n_ch:.3f}", "0",
                              f"{rb:.2f}", f"{n_fp:.3f}", "0"]))
            g.append(f"Bank Sta={lb:.2f},{rb:.2f}")
            g.append("XS Rating Curve= 0 ,0")
            g.append("Exp/Cntr=0.3,0.1")
            g.append("")

    open(f"{PROJECT}.g01", "w", encoding="ascii", errors="replace").write(
        "\n".join(g) + "\n")
    n_xs = sum(len(t["xs"]) for t in trechos)
    print(f"  [OK] {PROJECT}.g01  ({len(trechos)} trechos, {len(juncoes)} juncoes, {n_xs} secoes)")


def escrever_fluxo(trechos, cabeceiras, saida, laterais=(), uniformes=()):
    u = [f"Flow Title=Cheia_Rede_Real", "Program Version=7.01", "Use Restart= 0 "]
    for t in trechos:
        u.append(f"Initial Flow Loc={p16(t['rio'])},{p16(t['reach'])},"
                 f"{t['xs'][0]['rs']:<8.0f},{t['q_base']:.0f}")
    for t in cabeceiras:
        u.append(bl(t["rio"], t["reach"], f"{t['xs'][0]['rs']:.2f}"))
        u.append("Interval=1HOUR")
        u.append(f"Flow Hydrograph= {NHORAS} ")
        u.append(serie8(t.get("serie") if t.get("serie") is not None
                        else hidrograma(t["q_pico"])))
    for lt in laterais:
        u.append(bl(lt["rio"], lt["reach"], f"{lt['rs']:.2f}"))
        u.append("Interval=1HOUR")
        u.append(f"Lateral Inflow Hydrograph= {NHORAS} ")
        u.append(serie8(lt.get("serie") if lt.get("serie") is not None
                        else hidrograma(lt["q_pico"])))
    for un in uniformes:
        # Uniform Lateral Inflow: o campo 4 do Boundary Location recebe o RS
        # de JUSANTE do intervalo (formato do exemplo oficial UngagedAreaInflows)
        rs_hi = f"{un['rs_hi']:.2f}"[:8].ljust(8)
        rs_lo = f"{un['rs_lo']:.2f}"[:8].ljust(8)
        u.append(f"Boundary Location={p16(un['rio'])},{p16(un['reach'])},"
                 f"{rs_hi},{rs_lo},                ,                ")
        u.append("Interval=1HOUR")
        u.append(f"Uniform Lateral Inflow Hydrograph= {NHORAS} ")
        u.append(serie8(un.get("serie") if un.get("serie") is not None
                        else hidrograma(un["q_pico"])))
    u.append(bl(saida["rio"], saida["reach"], f"{saida['xs'][-1]['rs']:.2f}"))
    if MARE:
        u.append("Interval=1HOUR")
        u.append(f"Stage Hydrograph= {NHORAS} ")
        u.append(serie8(mare(NHORAS)))
    else:
        u.append(f"Friction Slope={DS_SLOPE}")
    open(f"{PROJECT}.u01", "w", encoding="ascii", errors="replace").write(
        "\n".join(u) + "\n")
    print(f"  [OK] {PROJECT}.u01  ({len(cabeceiras)} hidrogramas de cabeceira, "
          f"{len(laterais)} vazoes laterais)")


MES_RAS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]


def data_ras(dt):
    """datetime -> '05JUL1983,0000', o formato de data do HEC-RAS."""
    return f"{dt.day:02d}{MES_RAS[dt.month-1]}{dt.year},{dt.hour:02d}00"


def escrever_plano_prj():
    # A janela tem de terminar EXATAMENTE na ultima ordinata da serie
    # (h = NHORAS-1). Terminar depois disso aborta o calculo com
    # "Time series data ends before the end of the simulation".
    #
    # A data vem do EVENTO. Estava fixa em 01AUG2026, entao a cheia de julho de
    # 1983 aparecia rotulada como agosto de 2026 no RAS, no RAS Mapper, no HDF e
    # em tudo que le dali -- inclusive na interface web.
    total_h = NHORAS - 1
    ini = DATA_INICIO or datetime.datetime(2026, 8, 1)
    fim = ini + datetime.timedelta(hours=total_h)
    open(f"{PROJECT}.p01", "w", encoding="ascii").write("\n".join([
        "Plan Title=Rede_Real_Bacia_Itajai", "Program Version=7.01",
        "Short Identifier=REDE", "Geom File=g01", "Flow File=u01",
        f"Simulation Date={data_ras(ini)},{data_ras(fim)}",
        "Mixed Flow Regime",
        # --- opcoes do solver 1D, que ate aqui ficaram TODAS no padrao.
        #
        # 'Mixed Flow Regime' e o LPI sao coisas separadas: o primeiro so
        # permite regime misto, o segundo amortece os termos de inercia
        # quando o Froude se aproxima de 1. O amortecimento vem DESLIGADO de
        # fabrica (Froude Reduction=False), e e exatamente do que a garganta
        # do Salto Pilao precisa -- o diagnostico acusa 51 secoes com Fr > 0,9.
        "UNET Froude Reduction=True",
        "UNET Froude Limit= 0.8 ",
        "UNET Froude Power= 4 ",
        # O log mostrava '20' no fim de cada linha de iteracao: o solver batia
        # o teto de iteracoes com erro de 1,06 m e desistia sem convergir.
        "UNET MxIter= 40 ",
        "UNET Max Iter WO Improvement= 20 ",
        "UNET Theta= 1 ",
        "UNET Theta Warmup= 1 ",
        "UNET ZTol= 0.01 ",
        "UNET ZSATol= 0.01 ",
        "UNET DZMax Abort= 30 ",
        # Passos so para assentar a condicao inicial, antes de comecar a
        # contar o tempo. E o que tira o transiente de partida que secava os
        # trechos planos (257 secoes chegavam a menos de 5 cm de lamina).
        "UNET MaxInSteps= 200 ",
        "UNET DtIC= 0 ",
        "Flow Smoothing Iterations=10",
        "Unsteady Friction Slope Method= 2 ",
        "UNET 1D Methodology=Finite Difference",
        # 15SEC era necessario para a rede instavel de varias juncoes; com
        # uma juncao so, 1MIN converge e roda ~4x mais rapido (relevante para
        # eventos reais de 168-312 h).
        "Computation Interval=1MIN", "Output Interval=1HOUR",
        "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
        "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
        # Mapeamento de planicie ligado. Ele so PRODUZ os rasters se houver um
        # terreno importado no RAS Mapper (Project > New Terrain), porque a
        # profundidade e a cota d'agua menos o terreno. Sem terreno o HEC-RAS
        # simplesmente nao escreve os mapas, sem erro.
        "Run RASMapper=-1"]) + "\n")
    open(f"{PROJECT}.prj", "w", encoding="ascii").write("\n".join([
        f"Proj Title={PROJECT}", "Current Plan=p01",
        "Default Exp/Contr=0.3,0.1", "SI Units", "Geom File=g01",
        "Unsteady File=u01", "Plan File=p01", "Y Axis Title=Elevation",
        "X Axis Title(PR)=Distance", "X Axis Title(CS)=Station",
        f"RASMap Filename={PROJECT}.rasmap"]) + "\n")

    # O terreno NAO e referenciado aqui de proposito: o RAS Mapper espera
    # terreno no formato HDF dele (gerado pelo import de terreno), e apontar
    # para o GeoTIFF cru produz uma cascata de HDF5-DIAG "file signature not
    # found". Para ver o relevo, importe o DEM pelo RAS Mapper uma vez.
    # Sem estes dois arquivos o RAS termina com
    #   "Error executing: StoreAllMaps / RasMapFilename '' does not exist."
    # Alem de silenciar o erro, fazem o projeto abrir georreferenciado no
    # RAS Mapper, com o DEM ja carregado como terreno.
    wkt = ('PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",'
           'DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",'
           'SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],'
           'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
           'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],'
           'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
           'PARAMETER["false_northing",10000000],UNIT["metre",1]]')
    with open(f"{PROJECT}.projection", "w", encoding="utf-8") as f:
        f.write(wkt)
    # O .rasmap so declarava a projecao, com <Terrains/> e <Results/> vazios e
    # NENHUM bloco <Geometries>. Sem ele o RAS Mapper nao sabe que camada
    # desenhar. A estrutura abaixo e a que o proprio RAS Mapper grava, com os
    # caminhos no formato ".\arquivo" que ele usa.
    #
    # O terreno so e referenciado se ja tiver sido IMPORTADO pelo RAS Mapper
    # (Project > New Terrain), que o converte para o .hdf dele. Apontar para o
    # GeoTIFF cru produz uma cascata de "HDF5-DIAG: file signature not found".
    # Enquanto nao existir, a entrada fica vazia e o resto do .rasmap funciona.
    # O .hdf tem prioridade: e o formato proprio, criado em Project > New
    # Terrain, e o unico que serve para o HEC-RAS calcular profundidade. O .tif
    # ja reprojetado (preparar_terreno.py) entra como alternativa para pelo
    # menos VER o relevo, e e tambem o arquivo que se aponta no import.
    terr = None
    for cand in ("Terrain/Terreno.hdf", "Terrain/Terrain.hdf",
                 f"Terrain/{PROJECT}.hdf", "Terrain/Terreno_Copernicus.tif"):
        if os.path.exists(cand):
            terr = cand.replace("/", "\\")
            break
    terrenos = ('  <Terrains Checked="True" Expanded="True">\n'
                f'    <Layer Name="{os.path.splitext(os.path.basename(terr))[0]}" '
                f'Type="TerrainLayer" Checked="True" '
                f'Filename=".\\{terr}" />\n'
                '  </Terrains>\n') if terr else '  <Terrains Checked="True" />\n'

    # limites da planicie de inundacao, se ja gerados (gerar_planicie.py)
    shp = f"{PROJECT}_planicie.shp"
    mapas = ('  <MapLayers Checked="True" Expanded="True">\n'
             f'    <Layer Name="Planicie de inundacao" '
             f'Type="PolygonFeatureLayer" Checked="True" '
             f'Filename=".\\{shp}" />\n'
             '  </MapLayers>\n') if os.path.exists(shp) else \
            '  <MapLayers Checked="True" />\n'
    with open(f"{PROJECT}.rasmap", "w", encoding="utf-8") as f:
        f.write(
            '<?xml version="1.0" encoding="utf-8"?>\n<RASMapper>\n'
            '  <Version>2.00</Version>\n'
            f'  <RASProjectionFilename Filename=".\\{PROJECT}.projection" />\n'
            # A geometria NAO e declarada aqui. O RAS Mapper ja a descobre pelo
            # projeto, com a arvore completa (Rivers, Junctions, Bank Lines,
            # Cross Sections...) e usando o Geom Title como nome. Declarar
            # tambem, com outro nome, fazia aparecerem DUAS entradas apontando
            # para o mesmo .g01.hdf -- uma com sub-camadas e a minha vazia, o
            # que so confundia. O que faltava para ela desenhar era o
            # Viewing Rectangle, nao esta declaracao.
            '  <Geometries Checked="True" Expanded="True" />\n'
            '  <Results Checked="True" Expanded="True">\n'
            f'    <Layer Name="REDE" Type="RASResults" Checked="True" '
            f'Expanded="True" Filename=".\\{PROJECT}.p01.hdf">\n'
            f'      <Layer Name="Event Conditions" Type="RASEventConditions" '
            f'Filename=".\\{PROJECT}.p01.hdf" />\n'
            # Mapas armazenados: e o que faz o HEC-RAS gerar a MANCHA. Sem
            # estas tres camadas declaradas, 'Run RASMapper=-1' nao tem o que
            # produzir. ProfileName="Max" pede o envelope maximo da cheia --
            # a mancha de inundacao do evento inteiro.
            '      <Layer Name="Depth" Type="RASResultsMap" Checked="True">\n'
            '        <MapParameters MapType="depth" '
            'ProfileIndex="2147483647" ProfileName="Max" />\n'
            '      </Layer>\n'
            '      <Layer Name="WSE" Type="RASResultsMap" Checked="True">\n'
            '        <MapParameters MapType="elevation" '
            'ProfileIndex="2147483647" ProfileName="Max" />\n'
            '      </Layer>\n'
            '      <Layer Name="Velocity" Type="RASResultsMap" Checked="True">\n'
            '        <MapParameters MapType="velocity" '
            'ProfileIndex="2147483647" ProfileName="Max" />\n'
            '      </Layer>\n'
            '    </Layer>\n'
            '  </Results>\n'
            + terrenos + mapas +
            '</RASMapper>\n')
    print(f"  [OK] {PROJECT}.p01 / {PROJECT}.prj / {PROJECT}.rasmap")


# ------------------------------------------------------------------- PIPELINE
def main():
    print("=" * 68)
    print("REDE 1D REAL DA BACIA DO ITAJAI — topologia ANA + relevo DEM")
    print("=" * 68)

    global NHORAS, PROJECT, DATA_INICIO
    q_ev = None
    if EVENTO:
        from hidrologia_evento import hidrogramas
        q_ev, NHORAS = hidrogramas(EVENTO, barragens=BARRAGENS)
        PROJECT = f"Itajai_Rede_{EVENTO}"   # projeto separado do sintetico
        # data real do evento, do proprio arquivo de chuva observada
        for cam in (os.path.join("itajai_flood_model", "data",
                                 "rainfall_events", f"chuva_real_{EVENTO}.csv"),):
            if os.path.exists(cam):
                with open(cam, encoding="utf-8") as fh:
                    fh.readline()                     # cabecalho
                    t0 = fh.readline().split(",")[0].strip()
                try:
                    DATA_INICIO = datetime.datetime.fromisoformat(t0)
                    print(f"      inicio da simulacao: {data_ras(DATA_INICIO)} "
                          f"(data real do evento)")
                except ValueError:
                    pass
        print(f"[0/4] Evento {EVENTO}: chuva real observada, "
              f"{NHORAS} h, barragens "
              f"{'ATIVAS' if BARRAGENS else 'ABERTAS (sem obras)'}")
        for k, v in q_ev.items():
            print(f"      {k:<10} Q pico {v.max():8.1f} m3/s")

    rede = montar_rede()
    acu = rede[MAIN]["linha"]
    area_total = rede[MAIN]["area"]
    print(f"\n[1/4] Topologia — calha principal: {acu.length/1000:.1f} km, "
          f"{area_total:.0f} km2")

    # --- ARVORE da rede: em qual rio cada afluente desagua
    # Antes todo afluente era projetado no eixo do Acu, o que so vale para quem
    # entra no proprio Acu. Dos afluentes maiores da bacia, a maioria NAO
    # entra: o Trombudo, o Taio e o das Pombas desaguam no Oeste, o Krauel e o
    # Iraputa no Norte, o dos Cedros no Benedito. Projetar esses no Acu punha a
    # confluencia a dezenas de km do lugar certo. Agora a rede e uma arvore e o
    # Acu e so a raiz.
    ativos = [k for k in rede if k == MAIN or k in ESCOPO]
    receptor, filhos = {}, {k: [] for k in ativos}
    for k in ativos:
        if k == MAIN:
            continue
        foz = Point(list(rede[k]["linha"].coords)[-1])
        # o receptor e sempre um rio MAIOR: assim um afluente nunca e
        # pendurado noutro afluente menor que por acaso passe perto da foz
        cand = [m for m in ativos
                if m != k and rede[m]["area"] > rede[k]["area"]]
        if not cand:
            continue
        alvo = min(cand, key=lambda m: rede[m]["linha"].distance(foz))
        d = rede[alvo]["linha"].distance(foz)
        if d > 500.0:
            print(f"      ! {rede[k]['nome']}: foz a {d:.0f} m do rio mais "
                  f"proximo, fora da rede")
            continue
        receptor[k] = alvo
        filhos[alvo].append({"k": k, "s": rede[alvo]["linha"].project(foz),
                             "pt": foz})
    for m in filhos:
        filhos[m].sort(key=lambda d: d["s"])
    for m in sorted(filhos, key=lambda x: -rede[x]["area"]):
        for c in filhos[m]:
            Lm = rede[m]["linha"].length
            print(f"      {rede[c['k']]['nome']:<16} -> {rede[m]['nome']:<14} "
                  f"a {(Lm-c['s'])/1000:6.1f} km da foz dele")
    conf = filhos[MAIN]        # confluencias no proprio Acu

    dem = None
    if USAR_SIGSC:
        try:
            from dem_sigsc import DemSIGSC, DemHibrido
            # Reserva obrigatoria: os tiles do SIG-SC nao cobrem tudo. Sem
            # ela o Itajai do Sul (5% coberto) fica com 5 secoes em 87 km.
            dem = DemHibrido(DemSIGSC(), Dem(DEM))
            print(f"      relevo: {dem.nome}")
        except Exception as e:
            print(f"      ! SIG-SC indisponivel ({e}); usando {DEM}")
    if dem is None:
        dem = Dem(DEM)
        print(f"      relevo: Copernicus GLO-30 ({DEM})")
    print(f"\n[2/4] Cortando secoes do DEM (espacamento {SPACING:.0f} m, "
          f"largura {2*HALFWIDTH:.0f} m)...")

    # --- corta e condiciona CADA rio como um perfil continuo (evita degrau de
    #     leito nas juncoes internas) e so depois divide em trechos, nas
    #     confluencias dos SEUS proprios afluentes.
    def hw_area(frac):
        return AREA_CABECEIRA_ACU + (area_total - AREA_CABECEIRA_ACU) * frac

    def hw_acu(frac):
        return largura(hw_area(frac))

    DROP = 0.5                                # desnivel na juncao (m)
    trechos, acu_reaches = [], []
    reaches_de, bed_de = {}, {}

    def bed_em(k, rs_alvo):
        """cota do leito do rio k na estaca mais proxima de rs_alvo"""
        b = bed_de[k]
        return b[min(b, key=lambda r: abs(r - rs_alvo))]

    # Do MAIOR para o menor: quando um afluente e processado o leito do rio que
    # o recebe ja existe, e ele pode ser ancorado na cota certa da confluencia.
    for k in sorted(ativos, key=lambda x: -rede[x]["area"]):
        ln = rede[k]["linha"]
        L = ln.length
        if k == MAIN:
            xs = condicionar(secoes(ln, dem, 0.0, hw_acu, hw_area),
                             rede[k]["nome"])
            anc = ""
        else:
            a_k = rede[k]["area"]
            xs = condicionar(secoes(ln, dem, 0.0, largura(a_k), a_k),
                             rede[k]["nome"])
            m = receptor[k]
            s_conf = next(c["s"] for c in filhos[m] if c["k"] == k)
            alvo = bed_em(m, rede[m]["linha"].length - s_conf) + DROP
            desl = alvo - zt(xs[-1])
            for d in xs:                      # desloca o trecho inteiro
                d["z"] = d["z"] + desl
            anc = f"  ancorado em {alvo:.1f} m ({desl:+.1f} m)"
        bed_de[k] = {d["rs"]: zt(d) for d in xs}

        cortes = sorted({0.0} | {c["s"] for c in filhos[k] if c["s"] > 1.0}
                        | {L})
        reaches_de[k] = []
        for i in range(len(cortes) - 1):
            a, b = cortes[i], cortes[i + 1]
            if b - a < SPACING:
                continue
            rs_hi, rs_lo = L - a, L - b
            sel = ([d for d in xs if rs_lo < d["rs"] <= rs_hi] if reaches_de[k]
                   else [d for d in xs if rs_lo <= d["rs"] <= rs_hi])
            if len(sel) < 2:
                continue
            t = {"rio": rede[k]["nome"], "reach": f"R{len(reaches_de[k])+1}",
                 "linha": substring(ln, a, b), "xs": sel, "a": a, "b": b,
                 "rio_k": k}
            trechos.append(t)
            reaches_de[k].append(t)
            if k == MAIN:
                acu_reaches.append(t)
        # 'k' marca APENAS o trecho de jusante de um afluente -- e o que se
        # liga a juncao e o que recebe o hidrograma de cabeceira adiante.
        if k != MAIN and reaches_de[k]:
            c = next(c for c in filhos[receptor[k]] if c["k"] == k)
            reaches_de[k][-1].update({"k": k, "conf_s": c["s"], "pt": c["pt"]})
        print(f"      {rede[k]['nome']:<16}{L/1000:7.1f} km, {len(xs):4d} secoes,"
              f" {len(reaches_de[k])} trecho(s){anc}")

    # --- juncoes: em CADA rio, toda estaca que recebe afluente vira uma juncao
    #     ligando (afluentes que chegam ali + trecho de montante) -> trecho de
    #     jusante do proprio rio.
    print("\n[3/4] Juncoes")
    juncoes = []
    n_acu = 0
    for m in sorted(filhos, key=lambda x: -rede[x]["area"]):
        for s in sorted({c["s"] for c in filhos[m]}):
            entra = [c for c in filhos[m] if abs(c["s"] - s) < 1.0]
            ups = [reaches_de[c["k"]][-1] for c in entra
                   if reaches_de.get(c["k"])]
            up_m = [t for t in reaches_de[m] if abs(t["b"] - s) < 1.0]
            dn_m = [t for t in reaches_de[m] if abs(t["a"] - s) < 1.0]
            if not dn_m or not ups:
                continue
            # As juncoes do Acu tem nome de cidade; as dos afluentes levam o
            # nome do proprio afluente, senao NOME_JUNCAO passaria a rotular
            # confluencias erradas assim que a rede deixou de ser so o Acu.
            if m == MAIN:
                nome = NOME_JUNCAO.get(n_acu, f"J{n_acu+1}")
                n_acu += 1
            else:
                nome = f"Foz_{rede[entra[0]['k']]['nome']}"[:16].strip()
            # 'Junc L&A' e o caminho que a agua percorre ATRAVES da juncao:
            # da ultima secao do trecho de montante ate a primeira do de
            # jusante. Estava gravado como 500 m fixo para todas, enquanto a
            # geometria real da 75 m no Norte e 150 m no vao do Acu. Declarar
            # comprimento errado desequilibra a continuidade exatamente na
            # primeira secao abaixo da juncao -- que e onde o solver falhava
            # ("Solution Solver Failed" em Itajai_Acu R2 148108.3).
            rs_j = rede[m]["linha"].length - s        # estaca da juncao em m
            dn_dist = max(rs_j - dn_m[0]["xs"][0]["rs"], 0.0)
            dists = []
            for t in ups + up_m:
                if t["rio_k"] == m:                   # mesmo rio, a montante
                    d_up = max(t["xs"][-1]["rs"] - rs_j, 0.0)
                else:                                 # afluente: foz em RS 0
                    d_up = max(t["xs"][-1]["rs"], 0.0)
                dists.append(round(max(d_up + dn_dist, 1.0), 2))
            j = {"s": s, "nome": nome, "desc": "Confluencia",
                 "x": entra[0]["pt"].x, "y": entra[0]["pt"].y,
                 "up": [(t["rio"], t["reach"]) for t in ups + up_m],
                 "dists": dists,
                 "dn": (dn_m[0]["rio"], dn_m[0]["reach"])}
            juncoes.append(j)
            print(f"      {j['nome']:<18} {[u[0] for u in j['up']]} -> "
                  f"{j['dn'][0]}/{j['dn'][1]}")

    # --- vazoes rateadas por area de drenagem
    # cabeceira do Acu vira contorno de vazao no escopo reduzido
    acu_head = acu_reaches[0]
    acu_head["q_pico"] = Q_REF_FOZ * AREA_CABECEIRA_ACU / area_total
    acu_head["q_base"] = max(acu_head["q_pico"] * 0.15, 20.0)
    if q_ev:
        # a cabeceira do Acu e a confluencia Sul + Oeste (Rio do Sul)
        acu_head["serie"] = q_ev["sul"] + q_ev["oeste"]
        acu_head["q_base"] = float(acu_head["serie"][0])

    # Um rio que nasce de uma juncao NAO pode ter contorno proprio na cabeceira
    # -- seria vazao contada duas vezes. Vale para o Acu (juncao Sul+Oeste) e
    # agora tambem para qualquer afluente que receba outro na estaca zero.
    def nasce_de_juncao(k):
        return any(abs(c["s"]) < 1.0 for c in filhos.get(k, []))

    acu_nasce_de_juncao = nasce_de_juncao(MAIN)
    cabeceiras = [] if acu_nasce_de_juncao else [acu_head]
    if acu_nasce_de_juncao:
        print("      cabeceira do Acu = juncao Sul+Oeste (sem contorno proprio)")

    # area de drenagem que ja chegou ao inicio de um trecho, em QUALQUER rio:
    # a area propria do rio (o que nao pertence a nenhum afluente nomeado dele)
    # mais os afluentes que ja entraram a montante. Antes so o Acu acumulava;
    # com a arvore o Oeste, o Norte e o Benedito tambem tem afluentes proprios.
    def incremental_de(k):
        """Area do rio k que NAO pertence a nenhum afluente nomeado dele."""
        return max(rede[k]["area"]
                   - sum(rede[c["k"]]["area"] for c in filhos.get(k, []))
                   - (AREA_CABECEIRA_ACU if k == MAIN else 0.0), 0.0)

    def area_ate(k, s):
        propria = AREA_CABECEIRA_ACU if k == MAIN else incremental_de(k)
        return propria + sum(rede[c["k"]]["area"] for c in filhos.get(k, [])
                             if c["s"] <= s + 1.0)

    for t in trechos:
        a = area_ate(t["rio_k"], t["a"])
        if t is not acu_head:
            t["q_pico"] = Q_REF_FOZ * a / area_total
            t["q_base"] = max(t["q_pico"] * 0.15, 20.0)
            if q_ev and t["rio_k"] == MAIN:
                # Com evento, a vazao inicial dos trechos do Acu TEM de vir da
                # mesma fonte da cabeceira. Misturar serie real na cabeceira
                # com valor sintetico aqui cria um degrau de vazao na juncao
                # (R1 com 119 e R2 com 386 m3/s) e o solver nao converge.
                q0 = float(q_ev["sul"][0] + q_ev["oeste"][0])
                for c in conf:                    # afluentes ja incorporados
                    if c["s"] <= t["a"] + 1.0 and c["k"] in q_ev:
                        q0 += float(q_ev[c["k"]][0])
                for k_lat in LATERAIS:            # laterais ja incorporadas
                    if k_lat in rede and k_lat in q_ev:
                        p_lat = Point(list(rede[k_lat]["linha"].coords)[-1])
                        if acu.project(p_lat) <= t["a"] + 1.0:
                            q0 += float(q_ev[k_lat][0])
                if "acu_incr" in q_ev:            # incremental ja acumulado
                    q0 += float(q_ev["acu_incr"][0]) * (t["a"] / max(acu.length, 1.0))
                t["q_base"] = q0
    # O contorno de vazao entra na CABECEIRA de cada rio, isto e, no primeiro
    # trecho -- nao no ultimo. Enquanto cada afluente tinha um unico trecho os
    # dois eram o mesmo; agora que Oeste, Norte e Benedito se dividem em varios,
    # pendurar o hidrograma no trecho de jusante injetaria a agua ja na foz.
    # De onde vem o hidrograma de cada rio novo. O evento so traz serie para
    # Sul, Oeste, Norte, Benedito, Mirim e o incremental da bacia -- os oito
    # afluentes recem-adicionados nao tem serie propria. Mas a area deles JA
    # ESTA dentro da serie do rio que os recebe: o Trombudo esta dentro do
    # Oeste, o Krauel dentro do Norte, o Luis Alves dentro do incremental. Cada
    # um leva entao a fatia da serie do seu receptor proporcional a area, e o
    # receptor fica com o resto. Assim o volume do evento nao muda ao detalhar
    # a rede, e ninguem recebe agua sintetica misturada com serie real.
    a_ref = area_total - sum(rede[c["k"]]["area"] for c in filhos[MAIN]
                             if q_ev and c["k"] in q_ev)

    def fonte_serie(k):
        """(serie, area coberta) da fonte de hidrograma que cobre o rio k."""
        if q_ev and k in q_ev:
            return q_ev[k], rede[k]["area"]
        m = receptor.get(k)
        if m is not None:
            return fonte_serie(m)
        if q_ev and "acu_incr" in q_ev:
            return q_ev["acu_incr"], a_ref
        return None, 1.0

    def serie_de_area(k, area_km2):
        s, a = fonte_serie(k)
        return None if s is None else s * (area_km2 / max(a, 1.0))

    for k in sorted(ativos, key=lambda x: -rede[x]["area"]):
        if k == MAIN or not reaches_de.get(k):
            continue
        if nasce_de_juncao(k):
            print(f"      cabeceira do {rede[k]['nome']} = juncao "
                  f"(area propria entra como vazao lateral)")
            continue
        t = reaches_de[k][0]
        # a cabeceira leva a area PROPRIA do rio; o que pertence aos afluentes
        # dele entra nas respectivas confluencias, nao aqui
        own = incremental_de(k)
        t["q_pico"] = Q_REF_FOZ * own / area_total
        s = serie_de_area(k, own)
        if s is not None:
            t["serie"] = s
            t["q_base"] = float(s[0])
        cabeceiras.append(t)
    print("      vazao inicial acumulada na calha:")
    for t in acu_reaches:
        print(f"        {t['rio']}/{t['reach']}: Q inicial = {t['q_base']:6.1f} m3/s")
    print("\n[4/4] Contornos (rateio da cheia de 2008 por area de drenagem)")
    for t in cabeceiras:
        print(f"      {t['rio']:<14} Q pico = {t['q_pico']:7.1f} m3/s")

    # --- afluentes injetados como vazao lateral (sem juncao)
    laterais = []
    for k in LATERAIS:
        if k not in rede:
            continue
        p_ = Point(list(rede[k]["linha"].coords)[-1])
        s_ = acu.project(p_)
        rs_conf = acu.length - s_
        alvo = None
        for t in acu_reaches:
            rss = [d["rs"] for d in t["xs"]]
            if min(rss) <= rs_conf <= max(rss):
                alvo = t
                break
        if alvo is None:
            print(f"      ! {rede[k]['nome']}: confluencia fora dos trechos, ignorado")
            continue
        rs_sec = min((d["rs"] for d in alvo["xs"]), key=lambda r: abs(r - rs_conf))
        pico = Q_REF_FOZ * rede[k]["area"] / area_total
        laterais.append({"rio": alvo["rio"], "reach": alvo["reach"], "rs": rs_sec,
                         "nome": rede[k]["nome"], "q_pico": pico,
                         "serie": q_ev.get(k) if q_ev else None})
        print(f"      lateral: {rede[k]['nome']:<14} entra em {alvo['reach']} "
              f"RS {rs_sec/1000:6.1f} km, Q pico {pico:7.1f} m3/s")

    # --- area de drenagem PROPRIA de cada rio, distribuida ao longo da calha
    #
    # Duas correcoes que a arvore tornou necessarias:
    #
    # 1. O contab antigo somava rede[k]["area"] de TODO o ESCOPO. Como a area
    #    do Trombudo ja esta dentro da do Oeste, os netos eram contados duas
    #    vezes: a soma passava de 15.500 km2 numa bacia de 14.871 e o
    #    incremental do Acu zerava. So os afluentes DIRETOS do Acu contam.
    #
    # 2. Cada afluente tem incremental proprio, nao so o Acu. O Oeste tem
    #    1.593 km2 que nao pertencem ao Taio, ao Trombudo nem ao das Pombas;
    #    como ele passou a nascer da juncao do Taio, sem isto receberia apenas
    #    a agua do Taio e esses 1.593 km2 sumiriam do balanco.
    uniformes = []
    if INCREMENTAL:
        for k in sorted(ativos, key=lambda x: -rede[x]["area"]):
            rr = reaches_de.get(k) or []
            incr = incremental_de(k)
            # So quem NASCE DE JUNCAO recebe a area propria distribuida. Quem
            # tem contorno de cabeceira ja a recebeu la; repetir aqui dobrava a
            # vazao (o Sul aparecia com 776 m3/s na cabeceira e 776 de lateral).
            if not rr or incr < 1.0 or not nasce_de_juncao(k):
                continue
            L_tot = sum(t["linha"].length for t in rr)
            print(f"      incremental {rede[k]['nome']:<16}{incr:7.0f} km2 "
                  f"({100*incr/area_total:4.1f}% da bacia) em {L_tot/1000:6.1f} km")
            for t in rr:
                a = incr * t["linha"].length / L_tot
                pico = Q_REF_FOZ * a / area_total
                # O intervalo NAO pode tocar as secoes extremas do trecho:
                # "Uniform lateral inflows cannot start on the upstream cross
                #  section of a reach" (idem para a de jusante). Recua uma secao.
                ordenadas = sorted((d["rs"] for d in t["xs"]), reverse=True)
                if len(ordenadas) < 4:
                    continue                  # trecho curto demais p/ uniforme
                serie = serie_de_area(k, a)
                uniformes.append({"rio": t["rio"], "reach": t["reach"],
                                  "rs_hi": ordenadas[1], "rs_lo": ordenadas[-2],
                                  "q_pico": pico, "serie": serie})
                print(f"        {t['reach']}: {a:7.0f} km2  ->  Q pico "
                      f"{pico:6.1f} m3/s  (RS {ordenadas[1]/1000:.1f} a "
                      f"{ordenadas[-2]/1000:.1f} km)")

    # --- n de Manning por secao, a partir da zona de cada rio
    for t in trechos:
        cfg = MANNING.get(t["rio"])
        if cfg:
            cidade, la_, lo_, n_cima, n_baixo = cfg
            pc = gpd.GeoSeries([Point(lo_, la_)], crs=4326).to_crs(UTM_EPSG).iloc[0]
            s_c = t["linha"].project(pc)
            # RS da cidade no sistema do trecho: RS decresce rio abaixo
            rs_topo = max(d["rs"] for d in t["xs"])
            rs_cidade = rs_topo - s_c
            for d in t["xs"]:
                d["n"] = n_cima if d["rs"] > rs_cidade else n_baixo
            dentro = 0 <= s_c <= t["linha"].length
            print(f"      Manning {t['rio']}/{t['reach']}: {cidade} em RS "
                  f"{rs_cidade/1000:.1f} km -> n {n_cima} acima / {n_baixo} abaixo"
                  f"{'' if dentro else '  (cidade fora do trecho)'}")
        else:
            for d in t["xs"]:
                d["n"] = N_CANAL_PADRAO

    # --- n de Manning nas GARGANTAS, por Jarrett (1984)
    # Corredeira em rocha nao tem a rugosidade de rio de planicie. Usar
    # n = 0,035 no trecho do Salto Pilao (8 m/km) da Froude ~0,9: o
    # escoamento fica transcritico e o solver diverge no pico da cheia.
    #
    #     Jarrett (1984), para 0,002 <= S <= 0,052:   n = 0,39 S^0,38 R^-0,16
    #
    # que em S = 0,008 e R ~ 3 m da n ~ 0,052 -- o triplo do valor de
    # planicie. Isso nao e um ajuste para estabilizar: e a rugosidade que a
    # literatura mede nesses trechos, e ela por si so derruba o Froude para
    # ~0,5, porque a lamina engrossa e a velocidade cai.
    n_ajust = 0
    for t in trechos:
        xs = t["xs"]
        for i, d in enumerate(xs):
            viz = xs[min(i + 1, len(xs) - 1)]
            dx = d["rs"] - viz["rs"]
            if dx <= 0:
                continue
            S = abs(zt(d) - zt(viz)) / dx
            if S < 0.002:                       # fora da faixa de Jarrett
                continue
            R = max(canal_geometria(d.get("area_km2", 1000.0))[0], 0.5)
            n_j = 0.39 * S ** 0.38 * R ** -0.16
            n_j = float(min(max(n_j, d["n"]), 0.10))   # nunca abaixa o n
            if n_j > d["n"] + 1e-4:
                d["n"] = round(n_j, 4)
                n_ajust += 1
    if n_ajust:
        print(f"      Manning de garganta (Jarrett 1984): {n_ajust} secoes "
              f"com S >= 0,002 -> n ate 0,10")

    escrever(trechos, juncoes)
    escrever_fluxo(trechos, cabeceiras, acu_reaches[-1], laterais, uniformes)
    escrever_plano_prj()
    print(f"\nPronto.  Rode:  python run_hecras.py {PROJECT}")


if __name__ == "__main__":
    main()
