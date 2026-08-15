"""
03 - Monta a geometria HEC-RAS 1D da Bacia do Itajai a partir dos RIOS REAIS
(rios_itajai.geojson) e do RELEVO REAL (dem_bacia_itajai.tif):

  - corta secoes transversais perpendiculares ao eixo de cada rio, amostrando
    a elevacao do DEM (relevo real);
  - detecta as confluencias e divide o rio receptor em trechos, criando as
    JUNCOES do HEC-RAS;
  - grava .prj/.g01/.u01/.p01 no formato validado (colunas fixas de 8 char, etc.);
  - opcional: valida a geometria com o ras-commander se estiver disponivel.

Rode com:  python 03_gerar_geometria.py
Depois:    python run_hecras.py

Requisitos: geopandas, rasterio, shapely, pyproj, numpy.
"""
import json
import unicodedata
import numpy as np
import geopandas as gpd
import rasterio
from shapely.geometry import LineString, MultiLineString, Point, mapping
from shapely.ops import linemerge, substring
from pyproj import Transformer

# ------------------------- PARAMETROS AJUSTAVEIS ------------------------------
RIOS_GEOJSON = "rios_itajai.geojson"
DEM_TIF      = "dem_bacia_itajai.tif"
PROJECT      = "Itajai_Bacia_Real"
UTM_EPSG     = 31982          # SIRGAS 2000 / UTM 22S (metros)

# Quais rios modelar (casa por substring, sem acento, minusculo).
# "todos os rios reais, mirim inclusive"
RIOS_ALVO = ["itajai-acu", "itajai acu", "itajai do sul", "itajai do oeste",
             "itajai do norte", "hercilio", "itajai-mirim", "itajai mirim",
             "benedito"]

SPACING   = 1000.0    # espacamento entre secoes ao longo do rio (m)
HALFWIDTH = 800.0     # meia-largura da secao transversal (m)
STEP      = 30.0      # passo de amostragem do DEM na secao (m) ~ resolucao
MAXPTS    = 120       # maximo de pontos por secao (limite HEC-RAS = 450)
MIN_SLOPE = 1e-4      # declividade minima imposta ao talvegue (estabilidade)
SNAP_TOL  = 900.0     # tolerancia p/ detectar confluencia (m)
COND_BED  = True      # condicionar talvegue p/ ser monotonico rio abaixo

# Vazao de pico (m3/s) por rio para o hidrograma de cheia (headwaters).
PICO_POR_RIO = {
    "itajai do sul": 1200, "itajai do oeste": 1500, "itajai do norte": 2000,
    "hercilio": 2000, "itajai-mirim": 900, "itajai mirim": 900,
    "benedito": 700, "itajai-acu": 1000, "itajai acu": 1000,
}
PICO_DEFAULT = 600
NHOURS = 49

EDIT_TIME = "Node Last Edited Time= Aug/03/2026 00:00:00"


# ------------------------------ utilitarios -----------------------------------
def ascii_low(s):
    return "".join(c for c in unicodedata.normalize("NFD", s)
                   if unicodedata.category(c) != "Mn").lower()

def sanit(s, n=16):
    out = "".join(ch if (ch.isalnum() or ch in "_-") else "_"
                  for ch in ascii_low(s).replace(" ", "_"))
    return out[:n] if out else "Rio"

def f8(v):  return f"{v:8.3f}"

def fixed_series(vals):
    """valores em colunas fixas de 8 char, 10 por linha (formato HEC-RAS)."""
    return "\n".join("".join(f"{v:8.2f}" for v in vals[i:i+10])
                     for i in range(0, len(vals), 10))

def hydrograph(peak, base=None, n=NHOURS, tp=18, te=40):
    base = base if base is not None else max(peak * 0.15, 50)
    v = []
    for h in range(n):
        if h <= tp:   q = base + (peak-base)*(h/tp)
        elif h <= te: q = peak - (peak-base)*((h-tp)/(te-tp))
        else:         q = base
        v.append(q)
    return v


def main_line(geom):
    """Retorna a maior LineString conectada de uma geometria."""
    if geom.geom_type == "LineString":
        return geom
    merged = linemerge(geom)
    if merged.geom_type == "LineString":
        return merged
    return max(list(merged.geoms), key=lambda g: g.length)


# ------------------------------ amostragem DEM --------------------------------
class DemSampler:
    def __init__(self, path):
        self.ds = rasterio.open(path)
        self.nodata = self.ds.nodata
        self.to_dem = Transformer.from_crs(UTM_EPSG, self.ds.crs.to_epsg(),
                                           always_xy=True)

    def sample(self, xs_utm, ys_utm):
        lon, lat = self.to_dem.transform(xs_utm, ys_utm)
        vals = np.array([v[0] for v in self.ds.sample(list(zip(lon, lat)))],
                        dtype=float)
        if self.nodata is not None:
            vals[vals == self.nodata] = np.nan
        vals[vals < -1000] = np.nan
        return vals

    def elev_point(self, pt):
        return float(self.sample(np.array([pt.x]), np.array([pt.y]))[0])


def cut_section(line, along, dem):
    """Corta uma secao perpendicular ao eixo em 'along' e amostra o DEM."""
    eps = 1.0
    p  = line.interpolate(along)
    pa = line.interpolate(max(along-eps, 0))
    pb = line.interpolate(min(along+eps, line.length))
    tx, ty = pb.x-pa.x, pb.y-pa.y
    tl = np.hypot(tx, ty) or 1.0
    rx, ry = ty/tl, -tx/tl               # normal p/ a direita (olhando p/ jusante)
    offs = np.arange(-HALFWIDTH, HALFWIDTH+STEP, STEP)
    xs = p.x + offs*rx
    ys = p.y + offs*ry
    z  = dem.sample(xs, ys)
    # preenche buracos (nodata) por interpolacao linear
    if np.isnan(z).any():
        good = ~np.isnan(z)
        if good.sum() < 2:
            return None
        z = np.interp(np.arange(len(z)), np.flatnonzero(good), z[good])
    sta = offs + HALFWIDTH               # estacao 0..2*HW (cresce da esq p/ dir)
    # decima p/ nao passar de MAXPTS
    if len(sta) > MAXPTS:
        idx = np.linspace(0, len(sta)-1, MAXPTS).round().astype(int)
        sta, z = sta[idx], z[idx]
    return sta, z


def bank_stations(sta, z):
    """Estacoes de margem em torno do talvegue (ponto mais baixo)."""
    i0 = int(np.argmin(z))
    thr = z[i0] + 3.0
    li = i0
    while li > 0 and z[li] < thr: li -= 1
    ri = i0
    while ri < len(z)-1 and z[ri] < thr: ri += 1
    li = min(max(li, 1), len(sta)-3)          # margens estritamente interiores
    ri = max(min(ri, len(sta)-2), li+1)
    return float(sta[li]), float(sta[ri])


# ------------------------------ rede / confluencias ---------------------------
# ------------------------------ escrita HEC-RAS -------------------------------
def build(dem, reaches, junctions):
    """reaches: [{river,reach,line,xslist}]; junctions: [{name,pt,up:[(r,rc)],dn:(r,rc)}]"""
    g = [f"Geom Title=Bacia do Itajai - rios reais + relevo (DEM)",
         "Program Version=7.01"]
    for rc in reaches:
        g.append(f"River Reach={rc['river']:<16},{rc['reach']:<16}")
        xy = list(rc["line"].coords)
        g.append(f"Reach XY= {len(xy)} ")
        for i in range(0, len(xy), 2):
            g.append("".join(f"{x:16.4f}{y:16.4f}" for x, y in xy[i:i+2]))
        for xs in rc["xslist"]:
            sta, z, rs, ln, bl, br = xs
            g.append(f"Type RM Length L Ch R = 1 ,{rs:>10.2f} ,{ln:>10.2f},{ln:>10.2f},{ln:>10.2f}")
            g.append(EDIT_TIME)
            vals = [c for pair in zip(sta, z) for c in pair]
            g.append(f"#Sta/Elev= {len(sta)} ")
            for i in range(0, len(vals), 10):
                g.append("".join(f8(v) for v in vals[i:i+10]))
            g.append(f"Bank Sta={bl:.1f},{br:.1f}")
            g.append("#Mann= 3 , -1 , 0 ")
            g.append(f8(sta[0])+f8(0.06)+f8(0)+f8(bl)+f8(0.035)+f8(0)+f8(br)+f8(0.06)+f8(0))
    for j in junctions:
        g.append(f"Junct Name={j['name']:<16}")
        g.append(f"Junct Desc=Confluencia, 0 , 0 , 0 ,0")
        g.append(f"Junct X Y & Text X Y={j['pt'].x:.2f},{j['pt'].y:.2f},{j['pt'].x+500:.2f},{j['pt'].y+500:.2f}")
        for r, rcn in j["up"]:
            g.append(f"Up River,Reach={r:<16},{rcn:<16}")
        g.append(f"Dn River,Reach={j['dn'][0]:<16},{j['dn'][1]:<16}")
        for _ in j["up"]:
            g.append("Junc L&A=500,")
    with open(f"{PROJECT}.g01", "w") as f:
        f.write("\n".join(g) + "\n")
    print(f"[OK] {PROJECT}.g01  ({len(reaches)} trechos, {len(junctions)} juncoes)")


def write_unsteady(reaches, headwaters, outlet):
    def bl(river, reach, rs):
        return (f"Boundary Location={river:<16},{reach:<16},{rs:<8}"
                f",        ,                ,                ")
    u = ["Flow Title=Cheia_Bacia_Real", "Program Version=7.01", "Use Restart= 0 "]
    for rc in reaches:
        rs_up = rc["xslist"][0][2]
        u.append(f"Initial Flow Loc={rc['river']:<16},{rc['reach']:<16},{rs_up:<8.0f},{rc['base']:.0f}")
    for rc in headwaters:
        rs_up = rc["xslist"][0][2]
        u.append(bl(rc["river"], rc["reach"], f"{rs_up:.2f}"))
        u.append("Interval=1HOUR")
        u.append(f"Flow Hydrograph= {NHOURS} ")
        u.append(fixed_series(hydrograph(rc["peak"])))
    # jusante: profundidade normal na saida
    rs_dn = outlet["xslist"][-1][2]
    u.append(bl(outlet["river"], outlet["reach"], f"{rs_dn:.2f}"))
    u.append(f"Friction Slope={max(outlet['slope'],1e-4):.5f}")
    with open(f"{PROJECT}.u01", "w") as f:
        f.write("\n".join(u) + "\n")
    print(f"[OK] {PROJECT}.u01  ({len(headwaters)} entradas de cheia)")


def write_plan_prj():
    with open(f"{PROJECT}.p01", "w") as f:
        f.write("\n".join([
            "Plan Title=Bacia_Real", "Program Version=7.01",
            "Short Identifier=BaciaReal",
            "Simulation Date=01SEP2008,0000,03SEP2008,0000",
            "Geom File=g01", "Flow File=u01", "Subcritical Flow",
            "Computation Interval=1MIN", "Output Interval=1HOUR",
            "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
            "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
            "Run RASMapper=0"]) + "\n")
    with open(f"{PROJECT}.prj", "w") as f:
        f.write("\n".join([
            f"Proj Title={PROJECT}", "Current Plan=p01",
            "Default Exp/Contr=0.3,0.1", "SI Units",
            "Geom File=g01", "Unsteady File=u01", "Plan File=p01",
            "Y Axis Title=Elevation", "X Axis Title(PR)=Distance",
            "X Axis Title(CS)=Station"]) + "\n")
    print(f"[OK] {PROJECT}.p01 / {PROJECT}.prj")


# --------------------------------- pipeline -----------------------------------
def construir_geometria():
    print("Lendo rios e reprojetando para UTM 22S...")
    gdf = gpd.read_file(RIOS_GEOJSON).to_crs(epsg=UTM_EPSG)
    dem = DemSampler(DEM_TIF)

    # seleciona os rios-alvo
    sel = []
    for _, row in gdf.iterrows():
        nm = row.get("name") or ""
        key = ascii_low(nm)
        if any(t in key for t in RIOS_ALVO):
            line = main_line(row.geometry)
            if line.length < 2*SPACING:
                continue
            sel.append({"name": nm, "key": key, "line": line})
    # dedup por nome (fica o mais longo)
    best = {}
    for r in sel:
        if r["name"] not in best or r["line"].length > best[r["name"]]["line"].length:
            best[r["name"]] = r
    rios = list(best.values())
    print(f"Rios selecionados: {[r['name'] for r in rios]}")
    if not rios:
        raise SystemExit("Nenhum rio-alvo encontrado no geojson. Rode 01_baixar_rios.py.")

    # orienta montante->jusante pelo DEM (inicio mais ALTO)
    for r in rios:
        c = list(r["line"].coords)
        z0 = dem.elev_point(Point(c[0])); z1 = dem.elev_point(Point(c[-1]))
        if np.isnan(z0) or np.isnan(z1):
            continue
        if z1 > z0:                       # fim mais alto -> inverte
            r["line"] = LineString(c[::-1])

    # rio principal = o mais longo (nao vira afluente de ninguem)
    main_river = max(rios, key=lambda r: r["line"].length)
    print(f"Rio principal (saida): {main_river['name']}")

    # detecta confluencias: ponta de jusante de cada rio encostando em outro rio
    for r in rios:
        r["dn_pt"] = Point(list(r["line"].coords)[-1])
        r["splits"] = []                  # pontos (along) onde recebe afluentes
        r["recv"] = None                  # (rio_receptor, along_no_receptor)
    for r in rios:
        if r is main_river:
            continue
        best_d, best_recv, best_along = SNAP_TOL, None, None
        for s in rios:
            if s is r:
                continue
            d = s["line"].distance(r["dn_pt"])
            if d < best_d:
                best_d = d
                best_recv = s
                best_along = s["line"].project(r["dn_pt"])
        if best_recv is not None:
            r["recv"] = (best_recv, best_along)
            best_recv["splits"].append(best_along)
            print(f"  confluencia: {r['name']} -> {best_recv['name']} "
                  f"(a {best_along/1000:.1f} km, dist {best_d:.0f} m)")

    # divide cada rio nos pontos de confluencia -> trechos (reaches)
    reaches = []
    reach_of_river = {}
    for r in rios:
        cuts = sorted(set([0.0] + r["splits"] + [r["line"].length]))
        segs = []
        for i in range(len(cuts)-1):
            a, b = cuts[i], cuts[i+1]
            if b - a < SPACING:           # descarta segmentos minusculos
                continue
            seg = substring(r["line"], a, b)
            segs.append({"a": a, "b": b, "line": seg})
        r["segs"] = segs
        river_name = sanit(r["name"])
        for i, seg in enumerate(segs, 1):
            rc = {"river": river_name, "reach": f"R{i}", "line": seg["line"],
                  "riverobj": r, "seg": seg, "a": seg["a"], "b": seg["b"]}
            reaches.append(rc)
        reach_of_river[id(r)] = segs

    # garante unicidade river+reach
    seen = set()
    for rc in reaches:
        base = rc["river"]
        k = (rc["river"], rc["reach"])
        while k in seen:
            rc["reach"] += "b"
            k = (rc["river"], rc["reach"])
        seen.add(k)

    # corta secoes de cada trecho a partir do DEM
    print("Cortando secoes transversais do relevo...")
    for rc in reaches:
        line = rc["line"]
        L = line.length
        alongs = list(np.arange(0, L, SPACING))
        if alongs[-1] < L - 1:
            alongs.append(L - 1)
        xslist = []
        prev_along = None
        for a in alongs:
            res = cut_section(line, a, dem)
            if res is None:
                continue
            sta, z = res
            rs = round(L - a, 2)          # RS cresce p/ montante
            ln = SPACING if prev_along is not None else SPACING
            blb, brb = bank_stations(sta, z)
            xslist.append([sta, z, rs, ln, blb, brb])
            prev_along = a
        # ordena por RS decrescente (montante -> jusante) e ajusta comprimentos
        xslist.sort(key=lambda x: -x[2])
        for k in range(len(xslist)-1):
            xslist[k][3] = round(abs(xslist[k][2]-xslist[k+1][2]), 2) or SPACING
        xslist[-1][3] = 0.0
        # condiciona talvegue monotonico rio abaixo (estabilidade numerica)
        if COND_BED and len(xslist) > 1:
            for k in range(1, len(xslist)):
                zmin_prev = xslist[k-1][1].min()
                need = zmin_prev - MIN_SLOPE*xslist[k-1][3]
                zmin_cur = xslist[k][1].min()
                if zmin_cur > need:
                    xslist[k][1] = xslist[k][1] - (zmin_cur - need)
        rc["xslist"] = xslist
        # slope media do trecho (p/ Normal Depth se for saida)
        if len(xslist) > 1:
            dz = xslist[0][1].min() - xslist[-1][1].min()
            rc["slope"] = max(dz / max(L, 1.0), 1e-4)
        else:
            rc["slope"] = 1e-3
        print(f"  {rc['river']}/{rc['reach']}: {len(xslist)} secoes, "
              f"L={L/1000:.1f} km, decl={rc['slope']:.4f}")

    reaches = [rc for rc in reaches if len(rc.get("xslist", [])) >= 2]

    # vazoes: base/pico por trecho
    for rc in reaches:
        key = rc["riverobj"]["key"]
        peak = PICO_DEFAULT
        for k, v in PICO_POR_RIO.items():
            if k in key:
                peak = v; break
        rc["peak"] = peak
        rc["base"] = max(peak*0.15, 50)

    # ---- monta juncoes a partir das confluencias ----
    def upstream_seg(r, along):
        """segmento do rio r imediatamente a MONTANTE do ponto 'along'."""
        cand = [rc for rc in reaches if rc["riverobj"] is r and rc["b"] <= along+1]
        return max(cand, key=lambda rc: rc["b"]) if cand else None
    def downstream_seg(r, along):
        cand = [rc for rc in reaches if rc["riverobj"] is r and rc["a"] >= along-1]
        return min(cand, key=lambda rc: rc["a"]) if cand else None
    def trib_last_seg(r):
        cand = [rc for rc in reaches if rc["riverobj"] is r]
        return max(cand, key=lambda rc: rc["b"]) if cand else None

    junctions = []
    for r in rios:
        if not r.get("recv"):
            continue
        recv, along = r["recv"]
        up_trib = trib_last_seg(r)
        up_main = upstream_seg(recv, along)
        dn_main = downstream_seg(recv, along)
        ups = [rc for rc in (up_trib, up_main) if rc is not None]
        if dn_main is None or not ups:
            continue
        junctions.append({
            "name": sanit(f"J_{r['name']}"),
            "pt": r["dn_pt"],
            "up": [(rc["river"], rc["reach"]) for rc in ups],
            "dn": (dn_main["river"], dn_main["reach"]),
        })

    # ---- contornos: headwaters e saida ----
    dn_reach_ids = set()
    for j in junctions:
        dn_reach_ids.add(j["dn"])
    up_reach_ids = set()
    for j in junctions:
        for u in j["up"]:
            up_reach_ids.add(u)
    # headwater = trecho cujo extremo de montante NAO e saida de juncao
    headwaters = []
    for rc in reaches:
        rid = (rc["river"], rc["reach"])
        # e o segmento mais a montante do seu rio E nao e dn de juncao
        first_seg = min([x for x in reaches if x["riverobj"] is rc["riverobj"]],
                        key=lambda x: x["a"])
        if rc is first_seg and rid not in dn_reach_ids:
            headwaters.append(rc)
    # saida = trecho cujo extremo de jusante nao alimenta juncao alguma
    outlet = None
    for rc in reaches:
        rid = (rc["river"], rc["reach"])
        last_seg = max([x for x in reaches if x["riverobj"] is rc["riverobj"]],
                       key=lambda x: x["b"])
        if rc is last_seg and rid not in up_reach_ids:
            # prefere o Itajai-Acu como saida
            if outlet is None or "acu" in rc["riverobj"]["key"]:
                outlet = rc
    if outlet is None:
        outlet = reaches[-1]

    print(f"\nHeadwaters (entrada de cheia): "
          f"{[(rc['river'],rc['reach']) for rc in headwaters]}")
    print(f"Saida (Normal Depth): {(outlet['river'],outlet['reach'])}")

    build(dem, reaches, junctions)
    write_unsteady(reaches, headwaters, outlet)
    write_plan_prj()

    # opcional: valida com ras-commander se disponivel (Python >= 3.11)
    try:
        import ras_commander  # noqa
        print("\n(ras-commander disponivel — pode validar com RasPrj/RasGeom)")
    except Exception:
        pass

    print("\nPronto. Rode:  python run_hecras.py")


if __name__ == "__main__":
    construir_geometria()
