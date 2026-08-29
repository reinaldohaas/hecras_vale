# -*- coding: utf-8 -*-
"""QC das cross sections EXISTENTES de um projeto HEC-RAS, contra o DEM.

    python scripts/qc_secoes.py modelo/so_mirim.prj
    python scripts/qc_secoes.py modelo/so_mirim.prj --figuras 12
    python scripts/qc_secoes.py modelo/so_mirim.prj --rs 96865        # uma so

DIAGNOSTICO. NAO ESCREVE NA GEOMETRIA. Le o .prj, descobre o .gNN pelo
"Geom File=", le as secoes como estao, amostra o DEM ao longo de cada cutline
NAS MESMAS ESTACAS do HEC-RAS, e compara.

Premissas, explicitas:

  A BATIMETRIA E DADO, NAO ERRO. Entre as estacas de margem o perfil do
  HEC-RAS deve mesmo ficar ABAIXO do DEM -- o DEM de superficie enxerga a
  lamina d'agua, nao o fundo. Canal abaixo do DEM e registrado como
  'batimetria', nunca como defeito.

  O OVERBANK E QUE SE CONFERE. Fora das margens o perfil deveria seguir o
  terreno; divergencia ali e inconsistencia de verdade.

  TALVEGUE PERTO DA PONTA E DIAGNOSTICO, NAO ACAO. Relata-se; nao se desloca
  nem se estende secao nenhuma por causa disso.

Por River Station devolve: perfil HEC-RAS, perfil DEM e a diferenca; talvegue
do HEC-RAS e do DEM com o deslocamento horizontal entre eles; largura; cota
das margens; declividade local; e a comparacao com as duas vizinhas. Detecta
spikes, inversoes, descontinuidades e deslocamento horizontal.
"""
import os
import sys

import numpy as np

TOL_OVERBANK = 1.0      # m; divergencia mediana no overbank que ja conta
TOL_CANAL = 0.3         # m; abaixo disto o canal esta colado no DEM
PONTA = 0.05            # fracao das estacas que conta como extremidade
SPIKE = 3.0             # m acima da mediana movel para o ponto ser spike
DESCONT = 2.0           # m/m; degrau entre pontos vizinhos
DESLOC = 0.25           # fracao da largura: deslocamento horizontal que conta


# ------------------------------------------------------------------ leitura
def geometria_do_projeto(prj):
    """<projeto>.prj -> caminho do .gNN, pelo 'Geom File='."""
    ext = None
    for ln in open(prj, encoding="latin-1", errors="replace"):
        if ln.startswith("Geom File="):
            ext = ln.split("=", 1)[1].strip()
            break
    if not ext:
        raise SystemExit(f"{prj}: nao achei 'Geom File='")
    g = os.path.splitext(prj)[0] + "." + ext
    if not os.path.exists(g):
        raise SystemExit(f"{prj}: 'Geom File={ext}' aponta para {g}, inexistente")
    return g


def _bloco(linhas, i, largura):
    n = int(linhas[i].split("=")[1])
    v, i = [], i + 1
    while len(v) < 2 * n and i < len(linhas):
        l = linhas[i]
        if not l.strip() or l[:1].isalpha() or l[:1] == "#":
            break
        v += [float(l[j:j + largura]) for j in range(0, len(l), largura)
              if l[j:j + largura].strip()]
        i += 1
    return v[:2 * n]


def ler_secoes(g01):
    linhas = open(g01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    rio = reach = None
    secoes, i = [], 0
    while i < len(linhas):
        l = linhas[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio = p[0].strip()
            reach = p[1].strip() if len(p) > 1 else ""
        elif l.startswith("Type RM Length L Ch R"):
            p = [x.strip() for x in l.split("=", 1)[1].split(",")]
            # estruturas (Type 5 = barragem inline) tem comprimentos
            # vazios: ",,," -- tolerar
            d = {"rio": rio, "reach": reach, "tipo": p[0], "rs": float(p[1]),
                 "len_ch": float(p[3]) if len(p) > 3 and p[3] else np.nan,
                 "lb": None, "rb": None}
            j = i + 1
            while j < len(linhas) and not linhas[j].startswith("Type RM Length"):
                s = linhas[j]
                if s.startswith("River Reach="):
                    # cabecalho do PROXIMO reach: devolve ao laco de fora,
                    # senao o rio/reach nunca atualiza e toda secao sai
                    # atribuida ao primeiro rio do arquivo
                    break
                if s.startswith("#Sta/Elev"):
                    v = _bloco(linhas, j, 8)
                    d["sta"] = np.array(v[0::2]); d["z"] = np.array(v[1::2])
                elif s.startswith("XS GIS Cut Line"):
                    d["cut"] = np.array(_bloco(linhas, j, 16)).reshape(-1, 2)
                elif s.startswith("Bank Sta="):
                    b = [float(x) for x in s.split("=", 1)[1].split(",")]
                    d["lb"], d["rb"] = b[0], b[1]
                j += 1
            if "sta" in d:
                secoes.append(d)
            i = j
            continue
        i += 1
    return secoes


class Dem:
    def __init__(self, caminho):
        import rasterio
        self.src = rasterio.open(caminho)
        self.nod = self.src.nodata
        self.banda = self.src.read(1)
        self.inv = ~self.src.transform

    def cota(self, xs, ys):
        xs = np.atleast_1d(np.asarray(xs, float))
        ys = np.atleast_1d(np.asarray(ys, float))
        c, r = self.inv * (xs, ys)
        c = np.floor(c).astype(int); r = np.floor(r).astype(int)
        ok = ((r >= 0) & (r < self.banda.shape[0])
              & (c >= 0) & (c < self.banda.shape[1]))
        z = np.full(xs.shape, np.nan)
        z[ok] = self.banda[r[ok], c[ok]]
        if self.nod is not None:
            z[z == self.nod] = np.nan
        return z


def amostrar_dem(d, dem):
    """DEM ao longo da cutline, NAS MESMAS ESTACAS do HEC-RAS."""
    from shapely.geometry import LineString
    st = d["sta"]
    if "cut" not in d or len(d["cut"]) < 2:
        return np.full(len(st), np.nan)
    ln = LineString(d["cut"])
    f = np.clip((st - st[0]) / max(st[-1] - st[0], 1e-9), 0.0, 1.0)
    P = [ln.interpolate(float(x), normalized=True) for x in f]
    return dem.cota([p.x for p in P], [p.y for p in P])


# ---------------------------------------------------------------- analise
def _mediana_movel(z, jan=7):
    n = len(z)
    if n < jan:
        return z.copy()
    m = np.empty(n)
    h = jan // 2
    for i in range(n):
        m[i] = np.median(z[max(0, i - h):min(n, i + h + 1)])
    return m


def diagnosticar(d, dem, ant=None, prox=None):
    st, z = d["sta"], d["z"]
    zd = amostrar_dem(d, dem)
    dif = z - zd
    lb, rb = d.get("lb"), d.get("rb")
    larg = float(st[-1] - st[0])

    canal = ((st >= lb) & (st <= rb)) if (lb is not None and rb is not None
                                          and rb > lb) else np.zeros(len(st), bool)
    esq = st < lb if lb is not None else np.zeros(len(st), bool)
    dire = st > rb if rb is not None else np.zeros(len(st), bool)

    def med(m):
        v = dif[m & np.isfinite(dif)]
        return float(np.median(v)) if v.size else np.nan

    def pior(m):
        v = dif[m & np.isfinite(dif)]
        return float(v[np.argmax(np.abs(v))]) if v.size else np.nan

    # ---- talvegues, e o deslocamento horizontal entre eles
    k = int(np.argmin(z))
    sta_ras, z_ras = float(st[k]), float(z[k])
    if np.isfinite(zd).any():
        kd = int(np.nanargmin(zd))
        sta_dem, z_dem = float(st[kd]), float(zd[kd])
    else:
        kd, sta_dem, z_dem = -1, np.nan, np.nan
    desloc = sta_dem - sta_ras

    # ---- spikes: ponto muito fora da mediana movel do proprio perfil
    mm = _mediana_movel(z)
    res = z - mm
    spikes = int(np.sum(np.abs(res) > SPIKE))
    spike_max = float(np.abs(res).max()) if len(res) else 0.0

    # ---- descontinuidade: degrau entre pontos vizinhos
    ds = np.diff(st); dz = np.diff(z)
    incl = np.abs(dz) / np.maximum(ds, 1e-6)
    descont = int(np.sum(incl > DESCONT))
    descont_max = float(incl.max()) if incl.size else 0.0

    r = {
        "rio": d["rio"], "reach": d["reach"], "rs": d["rs"], "n": len(st),
        "largura": larg, "len_ch": d.get("len_ch", np.nan),
        "lb": lb, "rb": rb,
        "larg_canal": (rb - lb) if (lb is not None and rb is not None) else np.nan,
        "z_lb": float(np.interp(lb, st, z)) if lb is not None else np.nan,
        "z_rb": float(np.interp(rb, st, z)) if rb is not None else np.nan,
        "sem_dem": float(np.mean(~np.isfinite(zd))),
        "d_canal": med(canal), "d_esq": med(esq), "d_dir": med(dire),
        "pior_esq": pior(esq), "pior_dir": pior(dire),
        "sta_thal_ras": sta_ras, "z_thal_ras": z_ras,
        "sta_thal_dem": sta_dem, "z_thal_dem": z_dem,
        "desloc_horiz": desloc,
        "desloc_frac": desloc / larg if larg > 0 else np.nan,
        "spikes": spikes, "spike_max": spike_max,
        "descont": descont, "descont_max": descont_max,
        "i_thal": k,
    }

    # ---- vizinhas: declividade local e saltos
    def dz_ate(o):
        if o is None:
            return np.nan, np.nan
        zo = float(np.min(o["z"]))
        dxr = abs(float(d["rs"]) - float(o["rs"]))
        return zo, dxr

    z_ant, dx_ant = dz_ate(ant)
    z_prox, dx_prox = dz_ate(prox)
    r["z_thal_ant"], r["z_thal_prox"] = z_ant, z_prox
    # declividade local (montante -> jusante), com o RS decrescendo
    if np.isfinite(z_ant) and dx_ant > 0:
        r["decl_montante"] = (z_ant - z_ras) / dx_ant
    else:
        r["decl_montante"] = np.nan
    if np.isfinite(z_prox) and dx_prox > 0:
        r["decl_jusante"] = (z_ras - z_prox) / dx_prox
    else:
        r["decl_jusante"] = np.nan
    r["larg_ant"] = float(ant["sta"][-1] - ant["sta"][0]) if ant is not None else np.nan
    r["larg_prox"] = float(prox["sta"][-1] - prox["sta"][0]) if prox is not None else np.nan
    razao = [x for x in (r["larg_ant"], r["larg_prox"]) if np.isfinite(x) and x > 0]
    r["salto_largura"] = (max(max(larg / x, x / larg) for x in razao)
                          if razao and larg > 0 else np.nan)

    ob = [v for v in (r["d_esq"], r["d_dir"]) if np.isfinite(v)]
    r["overbank_max"] = max(abs(v) for v in ob) if ob else np.nan

    # ---------------------------------------------------------- marcas
    m = []
    if r["sem_dem"] > 0.05:
        m.append("sem-dem")
    if np.isfinite(r["d_canal"]):
        if r["d_canal"] < -TOL_CANAL:
            m.append("batimetria")
        elif r["d_canal"] > TOL_CANAL:
            m.append("canal-acima-do-dem")
    if np.isfinite(r["overbank_max"]) and r["overbank_max"] > TOL_OVERBANK:
        m.append("overbank-diverge")
    if spikes:
        m.append("spike")
    if descont:
        m.append("descontinuidade")
    if np.isfinite(r["desloc_frac"]) and abs(r["desloc_frac"]) > DESLOC:
        m.append("desloc-horizontal")
    if k < PONTA * len(st) or k > (1 - PONTA) * len(st):
        m.append("talvegue-na-ponta")
    if np.isfinite(r["decl_jusante"]) and r["decl_jusante"] < 0:
        m.append("inversao")                       # leito sobe rio abaixo
    if np.any(np.diff(st) < 0):
        m.append("estacas-fora-de-ordem")
    if np.any(np.diff(st) == 0):
        m.append("estacas-repetidas")
    if lb is not None and rb is not None and rb <= lb:
        m.append("margens-invertidas")
    if np.isfinite(r["salto_largura"]) and r["salto_largura"] > 3.0:
        m.append("salto-largura")
    r["marcas"] = m
    r["_dif"] = dif
    r["_zd"] = zd
    return r


CAMPOS = ["rio", "reach", "rs", "n", "largura", "len_ch", "lb", "rb",
          "larg_canal", "z_lb", "z_rb", "sem_dem",
          "sta_thal_ras", "z_thal_ras", "sta_thal_dem", "z_thal_dem",
          "desloc_horiz", "desloc_frac",
          "d_canal", "d_esq", "d_dir", "pior_esq", "pior_dir", "overbank_max",
          "spikes", "spike_max", "descont", "descont_max",
          "z_thal_ant", "z_thal_prox", "decl_montante", "decl_jusante",
          "larg_ant", "larg_prox", "salto_largura"]


def figura(d, r, destino):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    st, z, zd = d["sta"], d["z"], r["_zd"]
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.4), height_ratios=[2.2, 1],
                           sharex=True)
    a = ax[0]
    a.plot(st, z, "-", lw=1.6, color="#1f4e79", label="HEC-RAS (station-elevation)")
    a.plot(st, zd, "-", lw=1.2, color="#8a6d3b", label="DEM ao longo da cutline")
    if r["lb"] is not None:
        for s_, rot in ((r["lb"], "margem esq"), (r["rb"], "margem dir")):
            a.axvline(s_, color="#1a7a4c", ls="--", lw=1)
        a.axvspan(r["lb"], r["rb"], color="#1a7a4c", alpha=.07)
    a.plot(r["sta_thal_ras"], r["z_thal_ras"], "v", ms=9, color="#1f4e79",
           label=f"talvegue HEC-RAS ({r['sta_thal_ras']:.0f} m)")
    if np.isfinite(r["sta_thal_dem"]):
        a.plot(r["sta_thal_dem"], r["z_thal_dem"], "v", ms=9, color="#8a6d3b",
               label=f"talvegue DEM ({r['sta_thal_dem']:.0f} m)")
    a.set_ylabel("cota (m)")
    a.legend(fontsize=8, loc="upper center", ncol=2)
    a.set_title(f"{r['rio']} {r['reach']}  RS {r['rs']:.2f}   "
                f"largura {r['largura']:.0f} m   canal {r['larg_canal']:.0f} m   "
                f"| {', '.join(r['marcas']) or 'sem marcas'}",
                fontsize=10, loc="left")
    b = ax[1]
    b.plot(st, r["_dif"], "-", lw=1.2, color="#c0392b")
    b.axhline(0, color="#333", lw=.8)
    b.set_ylabel("HEC-RAS - DEM (m)")
    b.set_xlabel("estaca (m)")
    b.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(destino)
    plt.close(fig)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    prj = argv[0]
    g01 = geometria_do_projeto(prj)
    raiz = os.path.dirname(prj) or "."
    nome = os.path.splitext(os.path.basename(prj))[0]
    import glob
    cand = glob.glob(os.path.join(raiz, "Terrain", f"{nome}_Terreno*.tif"))
    if not cand:
        raise SystemExit(f"nao achei o DEM em {raiz}/Terrain")
    tif = cand[0]

    print(f"projeto  : {prj}")
    print(f"geometria: {g01}   (pelo 'Geom File=')")
    print(f"DEM      : {tif}")
    secoes = ler_secoes(g01)
    print(f"secoes lidas: {len(secoes)}")
    dem = Dem(tif)
    res = []
    for i, d in enumerate(secoes):
        ant = secoes[i - 1] if i > 0 and secoes[i - 1]["reach"] == d["reach"] else None
        prox = (secoes[i + 1] if i + 1 < len(secoes)
                and secoes[i + 1]["reach"] == d["reach"] else None)
        res.append(diagnosticar(d, dem, ant, prox))

    from collections import Counter
    c = Counter(m for r in res for m in r["marcas"])
    limpas = sum(1 for r in res if not r["marcas"])
    print("\n" + "=" * 78)
    print("MARCAS (uma secao pode ter varias)")
    print("=" * 78)
    for k, v in c.most_common():
        print(f"   {k:<24} {v:5d}   {100*v/len(res):5.1f}%")
    print(f"   {'(sem marca nenhuma)':<24} {limpas:5d}   {100*limpas/len(res):5.1f}%")

    def col(k):
        v = np.array([r[k] for r in res], float)
        return v[np.isfinite(v)]

    print("\n" + "=" * 78)
    print("DISTRIBUICOES")
    print("=" * 78)
    for k, rot, un in (("d_canal", "HEC-RAS - DEM no canal", "m"),
                       ("d_esq", "HEC-RAS - DEM no overbank esq", "m"),
                       ("d_dir", "HEC-RAS - DEM no overbank dir", "m"),
                       ("desloc_horiz", "talvegue DEM - talvegue HEC", "m"),
                       ("largura", "largura da secao", "m"),
                       ("larg_canal", "largura do canal", "m"),
                       ("decl_jusante", "declividade local p/ jusante", "m/m"),
                       ("salto_largura", "salto de largura vs vizinha", "x")):
        v = col(k)
        if not v.size:
            continue
        print(f"   {rot:<32} p10 {np.percentile(v,10):+9.3f}  "
              f"mediana {np.median(v):+9.3f}  p90 {np.percentile(v,90):+9.3f}  "
              f"max|.| {np.abs(v).max():9.3f} {un}")

    print("\n" + "=" * 78)
    print("AS 15 COM MAIOR DIVERGENCIA NO OVERBANK")
    print("=" * 78)
    ordem = sorted([r for r in res if np.isfinite(r["overbank_max"])],
                   key=lambda r: -r["overbank_max"])
    print(f"   {'RS':>10} {'larg':>6} {'canal':>7} {'obEsq':>7} {'obDir':>7} "
          f"{'desloc':>7} {'spk':>4}  marcas")
    for r in ordem[:15]:
        print(f"   {r['rs']:>10.2f} {r['largura']:>6.0f} {r['d_canal']:>+7.2f} "
              f"{r['d_esq']:>+7.2f} {r['d_dir']:>+7.2f} "
              f"{r['desloc_horiz']:>+7.1f} {r['spikes']:>4d}  "
              f"{','.join(r['marcas'])}")

    import csv
    saida = os.path.join(raiz, f"qc_{nome}.csv")
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(CAMPOS + ["marcas"])
        for r in res:
            w.writerow([r[k] for k in CAMPOS] + ["|".join(r["marcas"])])
    print(f"\numa linha por River Station: {saida}")

    nfig = 0
    if "--figuras" in argv:
        nfig = int(argv[argv.index("--figuras") + 1])
    if "--rs" in argv:
        alvo = float(argv[argv.index("--rs") + 1])
        i = int(np.argmin([abs(r["rs"] - alvo) for r in res]))
        pasta = os.path.join(raiz, "figuras"); os.makedirs(pasta, exist_ok=True)
        p = os.path.join(pasta, f"qc_{nome}_RS{res[i]['rs']:.0f}.svg")
        figura(secoes[i], res[i], p)
        print(f"figura: {p}")
    elif nfig:
        pasta = os.path.join(raiz, "figuras"); os.makedirs(pasta, exist_ok=True)
        for r in ordem[:nfig]:
            i = res.index(r)
            p = os.path.join(pasta, f"qc_{nome}_RS{r['rs']:.0f}.svg")
            figura(secoes[i], r, p)
        print(f"{nfig} figuras em {pasta}")
    return res, secoes, dem


if __name__ == "__main__":
    main(sys.argv[1:])
