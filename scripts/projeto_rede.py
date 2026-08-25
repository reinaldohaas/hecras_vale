# -*- coding: utf-8 -*-
"""Monta o projeto HEC-RAS rodavel da REDE (varios reaches e juncoes).

    python scripts/projeto_rede.py modelo/itajai_rede/itajai_rede.g01

E o `projeto_rio_avulso.py` generalizado para a rede que o `construir_rede.py`
costura. A diferenca esta no contorno, que agora sai da TOPOLOGIA:

  CABECEIRAS   um reach cuja montante nao e alimentada por juncao (nao aparece
               como "Dn River" de nenhuma) recebe o hidrograma de vazao do
               legado, do rio de mesmo nome. Sao Sul, Oeste, Norte e Mirim.

  FOZ          o reach cuja jusante nao entra em juncao (nao aparece como "Up
               River") e a saida. Se o fundo da ultima secao alcanca a mare
               (< 2 m), leva a mare do legado; senao, profundidade normal.

  JUNCOES      nao levam contorno: a agua passa de um reach ao outro por elas.

O ACU NAO TEM HIDROGRAMA DE MONTANTE: ele nasce da juncao Sul+Oeste. Sua vazao
vem inteira das juncoes. A ENTRADA LATERAL do Acu (o deflivio da bacia entre
os afluentes, que o legado tem) fica de fora deste primeiro fechamento -- os
quatro afluentes ja trazem o grosso da cheia --, e entra como refinamento
depois que a rede fechar. Fica dito aqui para nao se tomar por esquecimento.
"""
import argparse
import os
import shutil
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from ras_io import escrever                       # noqa: E402
from projeto_rio_avulso import hidrograma, mare, col8, LIMITE_MARE  # noqa: E402
from construir_rede import ler_geom, ler_juncoes  # noqa: E402

HIDRO = "legado/Itajai_Rede_1983.u01"
# reach do modelo (river) -> nome no legado para buscar o hidrograma. Aqui sao
# iguais; o mapa existe para o dia em que divergirem (ex.: Canal_Retif).
FONTE_HIDRO = {}


def juncoes_da_geom(path):
    return ler_juncoes(path)


def montar(geom, hidro=HIDRO, terreno=None, q_override=None):
    """Monta prj/p01/u01 da rede. `q_override` = {(river,reach): Q} substitui
    o hidrograma de uma cabeceira (usado pelo roteamento das barragens)."""

    pasta = os.path.dirname(geom)
    nome = os.path.basename(geom).split(".")[0]
    ext = os.path.basename(geom).split(".")[-1]
    _, reaches = ler_geom(geom)
    juncts = juncoes_da_geom(geom)

    dn = {j["dn"] for j in juncts}                # reaches alimentados por juncao
    ups = {u for j in juncts for u in j["ups"]}   # reaches que entram em juncao
    R = {}
    for r in reaches:
        rs = [s[0] for s in r["secoes"]]
        R[(r["river"], r["reach"])] = {
            "rs_up": max(rs), "rs_dn": min(rs), "obj": r}
    cabeceiras = [k for k in R if k not in dn]
    fozes = [k for k in R if k not in ups]
    print(f"geometria : {geom}")
    print(f"reaches   : {len(R)}   juncoes: {len(juncts)}")
    print("cabeceiras (recebem hidrograma):",
          ", ".join(f"{k[0]},{k[1]}" for k in cabeceiras))
    print("foz(es)    :", ", ".join(f"{k[0]},{k[1]}" for k in fozes))

    # ---- vazao inicial de cada reach: acumula a montante pelas juncoes
    Q0 = {}
    Qhead = {}
    for k in cabeceiras:
        if q_override and k in q_override:
            Q = np.asarray(q_override[k], float)
        else:
            Q = hidrograma(hidro, FONTE_HIDRO.get(k[0], k[0]))
        if Q is None:
            raise SystemExit(f"sem hidrograma de {k[0]} em {hidro}")
        Qhead[k] = Q
        Q0[k] = float(Q[0])
    # propaga Q0 juncao a juncao, de montante para jusante
    for _ in range(len(juncts) + 1):
        for j in juncts:
            s = sum(Q0.get(u, 0.0) for u in j["ups"])
            if s > 0:
                Q0[j["dn"]] = max(Q0.get(j["dn"], 0.0), s)

    # ================= u01 =================
    u = [f"Flow Title={nome}", "Program Version=7.01", "Use Restart= 0 "]
    for k in sorted(R):
        rs_up = R[k]["rs_up"]
        u.append(f"Initial Flow Loc={k[0]:<16.16},{k[1]:<16.16},"
                 f"{rs_up:<8.0f},{Q0.get(k, 1.0):.0f}")
    u.append("")
    # cabeceiras: Flow Hydrograph
    for k in cabeceiras:
        Q = Qhead[k]
        rs_up = R[k]["rs_up"]
        u += [f"Boundary Location={k[0]:<16.16},{k[1]:<16.16},{rs_up:<8.2f},"
              f"{'':<8},{'':<16},{'':<16}",
              "Interval=1HOUR", f"Flow Hydrograph= {len(Q)} "]
        u += col8(Q)
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
              "Fixed Start Date/Time=01AUG2026,0000",
              "Is Critical Boundary=False", "Critical Boundary Flow=", ""]
    # foz: mare ou profundidade normal
    for k in fozes:
        rs_dn = R[k]["rs_dn"]
        z_foz = _z_min(R[k]["obj"]["secoes"][-1][1])
        u += [f"Boundary Location={k[0]:<16.16},{k[1]:<16.16},{rs_dn:<8.2f},"
              f"{'':<8},{'':<16},{'':<16}"]
        H = mare(hidro) if z_foz < LIMITE_MARE else None
        if H is not None:
            u += ["Interval=1HOUR", f"Stage Hydrograph= {len(H)} "]
            u += col8(H)
            u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
                  "Fixed Start Date/Time=01AUG2026,0000", ""]
            print(f"   foz {k[0]},{k[1]} RS {rs_dn:.0f} (fundo {z_foz:.2f} m): "
                  f"MARE {H.min():.2f} a {H.max():.2f} m")
        else:
            decl = _declive(R[k]["obj"]["secoes"])
            u += [f"Friction Slope={decl:.6f},0", ""]
            print(f"   foz {k[0]},{k[1]} RS {rs_dn:.0f} (fundo {z_foz:.2f} m): "
                  f"profundidade normal, decl {decl:.5f}")
    escrever(os.path.join(pasta, f"{nome}.u01"), "\n".join(u))

    # ================= prj / p01 =================
    prj = [f"Proj Title={nome}", "Current Plan=p01",
           "Default Exp/Contr=0.3,0.1", "SI Units", f"Geom File={ext}",
           "Unsteady File=u01", "Plan File=p01",
           "Y Axis Title=Elevation", "X Axis Title(PF)=Main Channel Distance",
           "X Axis Title(XS)=Station", "BEGIN DESCRIPTION:",
           "Rede Itajai (MDT + batimetria 1983), 5 rios e 3 juncoes",
           "END DESCRIPTION:", "DSS Start Date=", "DSS Start Time=",
           "DSS End Date=", "DSS End Time=", "DSS Export Filename=",
           "DSS Export Rating Curves= 0 ", "DSS Export Rating Curve Sorted= 0 ",
           "DSS Export Volume Flow Curves= 0 ", "DXF Filename=",
           "DXF OffsetX= 0 ", "DXF OffsetY= 0 ", "DXF ScaleX= 1 ",
           "DXF ScaleY= 10 ", "GIS Export Profiles= 0 "]
    escrever(os.path.join(pasta, f"{nome}.prj"), "\n".join(prj))

    p01 = [f"Plan Title={nome}", "Program Version=7.01",
           f"Short Identifier={nome[:16]:<16}",
           "Simulation Date=01AUG2026,0000,08AUG2026,2300",
           f"Geom File={ext}", "Flow File=u01", "Mixed Flow Regime",
           "UNET Froude Reduction=True", "UNET Froude Limit= 0.8 ",
           "UNET Froude Power= 4 ", "UNET ZTol= 0.02 ", "UNET ZSATol= 0.02 ",
           "UNET MxIter= 40 ", "UNET Theta Warmup= 1 ",
           # dt=5 s: com o espacamento adaptativo (ate 1000-1250 m no
           # estuario) o pico da cheia sobre o canal de mare foi
           # instavel a 15 s -- explosao as 65,25 h, na chegada do
           # pico (Q oscilando a +459.578 m3/s, WSE 135 m a 2 km da
           # foz, Out=0 em 192 h). A corrida boa de 01:23 aguentava
           # 15 s porque tinha secoes a 150 m ali. 5 s e tambem o
           # que o gerador do Antigravity usa, estavel.
           "Computation Interval=5SEC", "Output Interval=1HOUR",
           "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
           "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
           "Run RASMapper=-1", "UNET D1 Cores= 0 "]
    escrever(os.path.join(pasta, f"{nome}.p01"), "\n".join(p01))

    # ================= projecao + rasmap =================
    fonte = None
    for c in ("SIRGAS2000_UTM22S.prj", "modelo/SIRGAS2000_UTM22S.prj",
              "modelo/itajai_acu/SIRGAS2000_UTM22S.prj"):
        if os.path.exists(c):
            fonte = c
            break
    if fonte:
        shutil.copy2(fonte, os.path.join(pasta, "SIRGAS2000_UTM22S.prj"))
    terreno = a.terreno
    if terreno is None:
        for c in (os.path.join("modelo", "Terrain", "vale30_Terreno.hdf"),
                  os.path.join(os.path.dirname(pasta), "Terrain",
                               "vale30_Terreno.hdf")):
            if os.path.exists(c):
                terreno = c
                break
    P = np.vstack([r["xy"] for r in reaches if len(r["xy"])])
    bloco = ["  <Terrains />"]
    if terreno and os.path.exists(terreno):
        rel = os.path.relpath(terreno, pasta).replace("/", "\\")
        tn = os.path.basename(terreno)[:-4]
        bloco = ['  <Terrains Checked="True" Expanded="True">',
                 f'    <Layer Name="{tn}" Type="TerrainLayer" Checked="True" '
                 f'Filename="{rel}"', '      Expanded="False">',
                 '      <ResampleMethod>near</ResampleMethod>',
                 '      <Surface On="True" />', '    </Layer>',
                 '  </Terrains>']
    rasmap = ["<RASMapper>", "  <Version>2.0.0</Version>",
              '  <RASProjectionFilename Filename=".\\SIRGAS2000_UTM22S.prj" />',
              '  <Geometries Checked="True" Expanded="True">',
              f'    <Layer Name="{nome}" Type="RASGeometry" Checked="True" '
              f'Expanded="True" Filename=".\\{nome}.{ext}.hdf">',
              '      <Layer Type="RASRiver" Checked="True" />',
              '      <Layer Type="RASBankLines" Checked="True" />',
              '      <Layer Type="FlowPaths" Checked="True" />',
              '      <Layer Type="RiverStations" Checked="True" />',
              '      <Layer Type="RASXS" Checked="True" />',
              '      <Layer Type="RASErrors" Checked="True" />',
              "    </Layer>", "  </Geometries>"] + bloco + [
              "  <CurrentView>", f"    <MinX>{P[:,0].min():.2f}</MinX>",
              f"    <MaxX>{P[:,0].max():.2f}</MaxX>",
              f"    <MinY>{P[:,1].min():.2f}</MinY>",
              f"    <MaxY>{P[:,1].max():.2f}</MaxY>",
              "  </CurrentView>", "</RASMapper>"]
    escrever(os.path.join(pasta, f"{nome}.rasmap"), "\n".join(rasmap))
    t2 = open(os.path.join(pasta, f"{nome}.prj"), encoding="latin-1").read()
    if "RASMap Filename=" not in t2:
        escrever(os.path.join(pasta, f"{nome}.prj"),
                 t2.rstrip("\r\n") + f"\nRASMap Filename={nome}.rasmap")

    print("\nCONFERENCIA")
    import xml.etree.ElementTree as ET
    ET.parse(os.path.join(pasta, f"{nome}.rasmap"))
    print("   rasmap: XML valido")
    for e in ("prj", "p01", "u01", ext, "rasmap"):
        p = os.path.join(pasta, f"{nome}.{e}")
        ok = os.path.exists(p)
        print(f"   {'OK   ' if ok else 'FALTA'} {nome}.{e}"
              + (f"   {os.path.getsize(p)} bytes" if ok else ""))
    return pasta


def _z_min(bloco):
    for i, l in enumerate(bloco):
        if l.startswith("#Sta/Elev"):
            v, j = [], i + 1
            while j < len(bloco) and not bloco[j].lstrip()[:1].isalpha() \
                    and not bloco[j].startswith("#"):
                L = bloco[j]
                v += [float(L[c:c + 8]) for c in range(0, len(L), 8)
                      if L[c:c + 8].strip()]
                j += 1
            z = np.array(v[1::2])
            return float(z.min()) if len(z) else 1e9
    return 1e9


def _declive(secoes):
    z = np.array([_z_min(b) for _, b in secoes])
    # comprimento de canal da penultima coluna do header
    ch = []
    for _, b in secoes:
        try:
            ch.append(float(b[0].split(",")[3]))
        except Exception:
            ch.append(0.0)
    ch = np.array(ch)
    k = max(len(z) - 11, 0)
    d = (z[k] - z[-1]) / max(ch[k:-1].sum(), 1.0)
    return float(max(d, 1e-4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geom")
    ap.add_argument("--hidrograma", default=HIDRO)
    ap.add_argument("--terreno", default=None)
    a = ap.parse_args()
    return montar(a.geom, a.hidrograma, a.terreno)


if __name__ == "__main__":
    main()
