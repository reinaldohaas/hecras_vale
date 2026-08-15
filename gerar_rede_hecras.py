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
import os
import unicodedata
import numpy as np
import geopandas as gpd
import rasterio
from shapely.geometry import LineString, Point
from shapely.ops import substring
from pyproj import Transformer

# ------------------------------------------------------------------ PARAMETROS
PROJECT   = "Itajai_Rede"
GEOJSON   = "rios_itajai.geojson"
DEM       = "dem_itajai.tif"
UTM_EPSG  = 31982              # SIRGAS 2000 / UTM 22S

SPACING   = 1000.0             # espacamento entre secoes (m)
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
}
MAIN = "acu"
# ESCOPO: quais afluentes entram na rede. Reduzir o escopo diminui o numero
# de juncoes (cada uma e um ponto potencial de instabilidade) e permite
# validar a cadeia inteira ate a mancha antes de voltar os demais rios.
#   completo -> ["sul","oeste","norte","benedito","mirim"]
#   reduzido -> ["mirim"]   (Acu + Mirim, 1 juncao)
ESCOPO = ["mirim"]
# Area de drenagem que entra pela CABECEIRA do Acu. No escopo reduzido os
# afluentes ausentes sao somados aqui para que as vazoes a jusante (Blumenau,
# Itajai) fiquem na ordem de grandeza correta.
LATERAIS = ["norte", "benedito"]
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
AREA_CABECEIRA_ACU = 5033.0    # Sul + Oeste (km2), que formam o Acu
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
            ups = [c for c in pred.get(ch[-1], [])
                   if c in by and any(a.lower() in by[c].NORIOCOMP.lower()
                                      for a in pat.split("|"))]
            if not ups:
                break
            melhor = max(ups, key=lambda c: by[c].NUAREAMONT)
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
        rede[k] = {"nome": nome, "linha": eixo(ch),
                   "area": float(by[ch[-1]].NUAREAMONT)}
    return rede


# --------------------------------------------------------------------- SECOES
def cortar(linha, s, dem, hw=HALFWIDTH):
    """Secao perpendicular ao eixo na posicao s, amostrada no DEM.
    A direcao usa uma janela de +-SMOOTH m: com +-1 m as cutlines se cruzam
    nas curvas (o RAS avisa 'edge lines have self intersections')."""
    a = linha.interpolate(max(0.0, s - SMOOTH))
    b = linha.interpolate(min(linha.length, s + SMOOTH))
    tx, ty = b.x - a.x, b.y - a.y
    n = np.hypot(tx, ty) or 1.0
    rx, ry = ty / n, -tx / n                  # normal a direita
    p = linha.interpolate(s)
    off = np.linspace(-hw, hw, NPTS)
    z = dem.sample(p.x + off * rx, p.y + off * ry)
    if np.isnan(z).all():
        return None
    if np.isnan(z).any():                     # tapa buracos por interpolacao
        ok = ~np.isnan(z)
        z = np.interp(np.arange(len(z)), np.flatnonzero(ok), z[ok])
    sta = off + hw                            # 0 .. 2*hw
    cut = (p.x - hw * rx, p.y - hw * ry,
           p.x + hw * rx, p.y + hw * ry)
    return np.round(sta, 2), z, cut


def margens(sta, z):
    """Margens topograficas: do talvegue ate BANK_H m de cada lado.
    O valor DEVE coincidir com um sta da tabela, na mesma precisao (.2f)."""
    i = int(np.nanargmin(z))
    lim = z[i] + BANK_H
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
    deixava o canal com 1-3 pontos."""
    return float(np.clip(180.0 * np.sqrt(max(area_km2, 1.0) / 100.0), 500.0, HALFWIDTH))


def secoes(linha, dem, rs0, hw=HALFWIDTH):
    """Corta as secoes de um trecho. rs0 = RS do extremo de JUSANTE."""
    L = linha.length
    ss = list(np.arange(0.0, L, SPACING))
    # NAO cria secao colada no extremo: uma secao a 1 m da juncao conflita
    # com o comprimento declarado em 'Junc L&A' e trava o solver.
    if L - ss[-1] > SPACING * 0.6:
        ss.append(L - SPACING * 0.5)
    xs = []
    for s in ss:
        r = cortar(linha, s, dem, hw)
        if r is None:
            continue
        sta, z, cut = r
        xs.append({"rs": round(rs0 + (L - s), 2), "sta": sta, "z": z,
                   "cut": cut})
    xs.sort(key=lambda d: -d["rs"])           # montante -> jusante
    # remove RS repetido (o RAS exige unicidade)
    fin, visto = [], set()
    for d in xs:
        if d["rs"] in visto:
            continue
        visto.add(d["rs"])
        fin.append(d)
    return fin


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
        dz = xs[corte]["z"].min() - xs[corte + 1]["z"].min()
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
        teto = xs[i - 1]["z"].min() - MIN_SLOPE * dx
        atual = xs[i]["z"].min()
        if atual > teto:
            xs[i]["z"] = xs[i]["z"] - (atual - teto)
    # -- 3. limita declividade ancorando a JUSANTE
    for i in range(len(xs) - 2, -1, -1):
        dx = xs[i]["rs"] - xs[i + 1]["rs"]
        lim = xs[i + 1]["z"].min() + MAX_SLOPE * dx
        atual = xs[i]["z"].min()
        if atual > lim:
            xs[i]["z"] = xs[i]["z"] - (atual - lim)
    return xs


# ---------------------------------------------------------------- ESCRITA RAS
def bl(river, reach, rs):
    return (f"Boundary Location={p16(river)},{p16(reach)},{str(rs)[:8]:<8}"
            f",        ,                ,                ")


def escrever(trechos, juncoes):
    g = [f"Geom Title={PROJECT} - rede real ANA + relevo DEM",
         "Program Version=7.01",
         "Viewing Rectangle= 0 , 1 , 1 , 0 "]
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
        for _ in j["up"]:
            g.append("Junc L&A=500,0")
        g.append("")

    for t in trechos:
        g.append(f"River Reach={p16(t['rio'])},{p16(t['reach'])}")
        c = list(t["linha"].coords)
        if len(c) > 400:
            c = [c[i] for i in np.linspace(0, len(c) - 1, 400).astype(int)]
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
            lb, rb = margens(sta, z)
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


def escrever_plano_prj():
    # A janela tem de terminar EXATAMENTE na ultima ordinata da serie
    # (h = NHORAS-1). Terminar depois disso aborta o calculo com
    # "Time series data ends before the end of the simulation".
    total_h = NHORAS - 1
    DIA_FIM = 1 + total_h // 24
    HORA_FIM = f"{(total_h % 24) * 100:04d}"
    open(f"{PROJECT}.p01", "w", encoding="ascii").write("\n".join([
        "Plan Title=Rede_Real_Bacia_Itajai", "Program Version=7.01",
        "Short Identifier=REDE", "Geom File=g01", "Flow File=u01",
        f"Simulation Date=01AUG2026,0000,{DIA_FIM:02d}AUG2026,{HORA_FIM}",
        "Mixed Flow Regime",
        # 15SEC era necessario para a rede instavel de varias juncoes; com
        # uma juncao so, 1MIN converge e roda ~4x mais rapido (relevante para
        # eventos reais de 168-312 h).
        "Computation Interval=1MIN", "Output Interval=1HOUR",
        "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
        "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
        "Run RASMapper=0"]) + "\n")
    open(f"{PROJECT}.prj", "w", encoding="ascii").write("\n".join([
        f"Proj Title={PROJECT}", "Current Plan=p01",
        "Default Exp/Contr=0.3,0.1", "SI Units", "Geom File=g01",
        "Unsteady File=u01", "Plan File=p01", "Y Axis Title=Elevation",
        "X Axis Title(PR)=Distance", "X Axis Title(CS)=Station",
        f"RASMap Filename={PROJECT}.rasmap"]) + "\n")

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
    with open(f"{PROJECT}.rasmap", "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n<RASMapper>\n'
                '  <Version>2.00</Version>\n'
                f'  <RASProjectionFilename Filename="{PROJECT}.projection" />\n'
                '  <Terrains>\n'
                f'    <Layer Name="Terreno" Filename="{DEM}" />\n'
                '  </Terrains>\n'
                '  <Results />\n</RASMapper>\n')
    print(f"  [OK] {PROJECT}.p01 / {PROJECT}.prj / {PROJECT}.rasmap")


# ------------------------------------------------------------------- PIPELINE
def main():
    print("=" * 68)
    print("REDE 1D REAL DA BACIA DO ITAJAI — topologia ANA + relevo DEM")
    print("=" * 68)

    global NHORAS
    q_ev = None
    if EVENTO:
        from hidrologia_evento import hidrogramas
        q_ev, NHORAS = hidrogramas(EVENTO, barragens=BARRAGENS)
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

    # posicao de cada afluente ao longo do Acu
    conf = []
    for k, v in rede.items():
        if k == MAIN or k not in ESCOPO:
            continue
        p = Point(list(v["linha"].coords)[-1])
        conf.append({"k": k, "s": acu.project(p), "pt": p})
    conf.sort(key=lambda d: d["s"])
    for c in conf:
        print(f"      {rede[c['k']]['nome']:<14} entra em {c['s']/1000:6.1f} km "
              f"({(acu.length-c['s'])/1000:5.1f} km da foz)")

    dem = Dem(DEM)
    print(f"\n[2/4] Cortando secoes do DEM (espacamento {SPACING:.0f} m, "
          f"largura {2*HALFWIDTH:.0f} m)...")

    # --- Acu: corta e condiciona como UM perfil continuo (evita degrau de
    #     leito nas juncoes internas), so depois divide em trechos
    acu_xs = condicionar(secoes(acu, dem, 0.0, largura(area_total)), rede[MAIN]["nome"])
    acu_bed = {d["rs"]: d["z"].min() for d in acu_xs}

    def bed_em(rs_alvo):
        """cota do leito do Acu na estaca mais proxima de rs_alvo"""
        k = min(acu_bed, key=lambda r: abs(r - rs_alvo))
        return acu_bed[k]

    cortes = sorted({0.0} | {c["s"] for c in conf if c["s"] > 1.0} |
                    {acu.length})
    trechos, acu_reaches = [], []
    for i in range(len(cortes) - 1):
        a, b = cortes[i], cortes[i + 1]
        if b - a < SPACING:
            continue
        seg = substring(acu, a, b)
        rs_hi, rs_lo = acu.length - a, acu.length - b
        xs = [d for d in acu_xs if rs_lo < d["rs"] <= rs_hi] if acu_reaches \
             else [d for d in acu_xs if rs_lo <= d["rs"] <= rs_hi]
        if len(xs) < 2:
            continue
        t = {"rio": rede[MAIN]["nome"], "reach": f"R{len(acu_reaches)+1}",
             "linha": seg, "xs": xs, "a": a, "b": b}
        trechos.append(t); acu_reaches.append(t)
        print(f"      {t['rio']}/{t['reach']}: {seg.length/1000:6.1f} km, "
              f"{len(xs):3d} secoes  (km {a/1000:.1f}–{b/1000:.1f})")

    # --- afluentes: condiciona e ANCORA o leito na cota do Acu na confluencia
    DROP = 0.5                                # desnivel na juncao (m)
    for c in conf:
        k = c["k"]; ln = rede[k]["linha"]
        xs = condicionar(secoes(ln, dem, 0.0, largura(rede[k]["area"])), rede[k]["nome"])
        alvo = bed_em(acu.length - c["s"]) + DROP
        desl = alvo - xs[-1]["z"].min()
        for d in xs:                          # desloca o trecho inteiro
            d["z"] = d["z"] + desl
        t = {"rio": rede[k]["nome"], "reach": "R1", "linha": ln, "xs": xs,
             "k": k, "conf_s": c["s"], "pt": c["pt"]}
        trechos.append(t)
        print(f"      {t['rio']}/R1: {ln.length/1000:6.1f} km, {len(xs):3d} secoes"
              f"  (ancorado em {alvo:.1f} m, deslocado {desl:+.1f} m)")

    # --- juncoes: cada confluencia liga (afluentes + Acu montante) -> Acu jusante
    print("\n[3/4] Juncoes")
    juncoes = []
    for i, c in enumerate(conf):
        ups = [t for t in trechos if t.get("k") == c["k"]]
        # trecho do Acu imediatamente a montante desta confluencia
        up_acu = [t for t in acu_reaches if abs(t["b"] - c["s"]) < 1.0]
        dn_acu = [t for t in acu_reaches if abs(t["a"] - c["s"]) < 1.0]
        if not dn_acu:
            continue
        # afluentes que entram exatamente no mesmo ponto (ex.: Sul + Oeste)
        ups += [t for t in trechos
                if t.get("k") and t["k"] != c["k"]
                and abs(t.get("conf_s", -1) - c["s"]) < 1.0
                and not any(u is t for u in ups)]
        if any(j["s"] == c["s"] for j in juncoes):
            continue
        j = {"s": c["s"], "nome": NOME_JUNCAO.get(len(juncoes), f"J{len(juncoes)+1}"),
             "desc": "Confluencia", "x": c["pt"].x, "y": c["pt"].y,
             "up": [(t["rio"], t["reach"]) for t in ups + up_acu],
             "dists": [ (t["xs"][-1]["rs"] if t.get("k")
                         else t["xs"][-1]["rs"] - (acu.length - c["s"]))
                        for t in ups + up_acu ],
             "dn": (dn_acu[0]["rio"], dn_acu[0]["reach"])}
        juncoes.append(j)
        print(f"      {j['nome']:<12} {[u[0] for u in j['up']]} -> {j['dn'][0]}/{j['dn'][1]}")

    # --- vazoes rateadas por area de drenagem
    # cabeceira do Acu vira contorno de vazao no escopo reduzido
    acu_head = acu_reaches[0]
    acu_head["q_pico"] = Q_REF_FOZ * AREA_CABECEIRA_ACU / area_total
    acu_head["q_base"] = max(acu_head["q_pico"] * 0.15, 20.0)
    if q_ev:
        # a cabeceira do Acu e a confluencia Sul + Oeste (Rio do Sul)
        acu_head["serie"] = q_ev["sul"] + q_ev["oeste"]
        acu_head["q_base"] = float(acu_head["serie"][0])

    cabeceiras = [acu_head]
    for t in trechos:
        if t.get("k"):                         # afluente: area propria
            a = rede[t["k"]]["area"]
        else:                                  # trecho do Acu: ACUMULA os
            # afluentes que ja entraram a montante deste trecho. Sem isso o
            # Acu comeca praticamente seco e a juncao instabiliza no warm-up.
            a = AREA_CABECEIRA_ACU + sum(rede[c["k"]]["area"]
                                         for c in conf if c["s"] <= t["a"] + 1.0)
        if t is not acu_head:
            t["q_pico"] = Q_REF_FOZ * a / area_total
            t["q_base"] = max(t["q_pico"] * 0.15, 20.0)
        if t.get("k"):
            if q_ev and t["k"] in q_ev:
                t["serie"] = q_ev[t["k"]]
                t["q_base"] = float(t["serie"][0])
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

    # --- area de drenagem do proprio Acu, distribuida ao longo da calha
    uniformes = []
    if INCREMENTAL:
        contab = AREA_CABECEIRA_ACU + sum(rede[k]["area"] for k in
                                          list(ESCOPO) + list(LATERAIS)
                                          if k in rede)
        incr = max(area_total - contab, 0.0)
        L_tot = sum(t["linha"].length for t in acu_reaches)
        print(f"      incremental do Acu: {incr:.0f} km2 "
              f"({100*incr/area_total:.1f}% da bacia), distribuida em "
              f"{L_tot/1000:.1f} km")
        for t in acu_reaches:
            a = incr * t["linha"].length / L_tot
            pico = Q_REF_FOZ * a / area_total
            # O intervalo NAO pode tocar as secoes extremas do trecho:
            # "Uniform lateral inflows cannot start on the upstream cross
            #  section of a reach" (idem para a de jusante). Recua uma secao.
            ordenadas = sorted((d["rs"] for d in t["xs"]), reverse=True)
            if len(ordenadas) < 4:
                continue                      # trecho curto demais p/ uniforme
            rs_hi = ordenadas[1]
            rs_lo = ordenadas[-2]
            frac = t["linha"].length / L_tot
            uniformes.append({"rio": t["rio"], "reach": t["reach"],
                              "rs_hi": rs_hi, "rs_lo": rs_lo, "q_pico": pico,
                              "serie": (q_ev["acu_incr"] * frac)
                                       if (q_ev and "acu_incr" in q_ev) else None})
            print(f"        {t['reach']}: {a:7.0f} km2  ->  Q pico {pico:6.1f} m3/s"
                  f"  (RS {rs_hi/1000:.1f} a {rs_lo/1000:.1f} km)")

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

    escrever(trechos, juncoes)
    escrever_fluxo(trechos, cabeceiras, acu_reaches[-1], laterais, uniformes)
    escrever_plano_prj()
    print(f"\nPronto.  Rode:  python run_hecras.py {PROJECT}")


if __name__ == "__main__":
    main()
