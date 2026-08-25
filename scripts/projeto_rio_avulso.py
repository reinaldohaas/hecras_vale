# -*- coding: utf-8 -*-
"""Monta um projeto HEC-RAS rodavel para um rio gerado do relevo.

    python scripts/projeto_rio_avulso.py modelo/benedito_mono/benedito_mono.g01 \
        --hidrograma legado/Itajai_Rede_1983.u01 --rio-fonte Rio_Benedito

Cria `.prj`, `.p01` e `.u01` ao lado da geometria. O rio gerado nao tem
contorno proprio, entao:

  MONTANTE   hidrograma de vazao COPIADO da rede legada, para o mesmo rio.
             E o unico dado de cheia que existe nesta bacia; gerar um
             sintetico so mudaria de lugar a invencao.

  JUSANTE    profundidade normal, com a declividade MEDIDA nos ultimos
             trechos da propria geometria. Rio avulso nao tem mare nem
             remanso conhecido -- quando ele entrar na rede, o contorno passa
             a ser a juncao, e este vai fora.

O `.p01` usa o mesmo passo e as mesmas tolerancias do resto do projeto
(15SEC, ZTol 0,02, 40 iteracoes), para que os numeros sejam comparaveis.
"""
import argparse
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes    # noqa: E402
from ras_io import escrever         # noqa: E402


def hidrograma(u01, rio, tipo=r"Flow Hydrograph"):
    """Vazao horaria daquele rio na rede legada.

    Procura primeiro o hidrograma de MONTANTE (`Flow Hydrograph`). O
    Itajai-Acu nao tem: na rede legada ele so recebe ENTRADA LATERAL ao longo
    dos quatro reaches, porque nasce da juncao dos outros. Para ele o
    chamador pede o tipo lateral, e o valor entra como vazao de montante --
    o que superestima a cabeceira e subestima o resto, e por isso fica
    impresso no relatorio.
    """
    t = open(u01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    ini = [i for i, l in enumerate(t) if l.startswith("Boundary Location=")]
    for a, b in zip(ini, ini[1:] + [len(t)]):
        p = t[a].split("=", 1)[1].split(",")
        if p[0].strip() != rio:
            continue
        for j in range(a, b):
            m = re.match(r"^%s=\s*(\d+)" % tipo, t[j])
            if not m:
                continue
            n = int(m.group(1))
            v, k = [], j + 1
            while k < b and len(v) < n:
                x = t[k]
                if not x.strip() or re.match(r"^[A-Za-z]", x):
                    break
                v += [float(x[c:c + 8]) for c in range(0, len(x), 8)
                      if x[c:c + 8].strip()]
                k += 1
            return np.array(v[:n])
    return None


def col8(v):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%8.2f" % x
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geom")
    ap.add_argument("--hidrograma", default="legado/Itajai_Rede_1983.u01")
    ap.add_argument("--rio-fonte", required=True)
    a = ap.parse_args()

    pasta = os.path.dirname(a.geom)
    nome = os.path.basename(a.geom).split(".")[0]
    S = ler_secoes(a.geom)
    S.sort(key=lambda d: -d["rs"])
    t = open(a.geom, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    rr = next(l for l in t if l.startswith("River Reach="))
    rio, rch = [x.strip() for x in rr.split("=", 1)[1].split(",")]

    Q = hidrograma(a.hidrograma, a.rio_fonte)
    origem_q = "Flow Hydrograph de montante"
    if Q is None:
        Q = hidrograma(a.hidrograma, a.rio_fonte,
                       r"Uniform Lateral Inflow Hydrograph")
        origem_q = ("Uniform Lateral Inflow -- este rio NAO tem hidrograma de "
                    "montante na rede legada")
    if Q is None:
        raise SystemExit(f"nao achei hidrograma de '{a.rio_fonte}' em "
                         f"{a.hidrograma}")
    print(f"geometria : {a.geom}")
    print(f"rio/reach : {rio} , {rch}   {len(S)} secoes")
    print(f"hidrograma: {len(Q)} horas   pico {Q.max():.1f} m3/s   "
          f"({a.rio_fonte}, {origem_q})")

    z = np.array([float(d["z"].min()) for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    k = max(len(S) - 11, 0)
    decl = (z[k] - z[-1]) / max(ch[k:-1].sum(), 1.0)
    decl = float(max(decl, 1e-4))
    print(f"jusante   : profundidade normal, declividade medida "
          f"{decl:.5f} nos ultimos {len(S)-k} trechos")

    prj = [f"Proj Title={nome}", "Current Plan=p01",
           "Default Exp/Contr=0.3,0.1", "SI Units", "Geom File=g01",
           "Unsteady File=u01", "Plan File=p01",
           "Y Axis Title=Elevation", "X Axis Title(PF)=Main Channel Distance",
           "X Axis Title(XS)=Station", "BEGIN DESCRIPTION:",
           f"{rio} gerado do relevo (MDT SIG-SC 1 m), sem perfil esculpido",
           "END DESCRIPTION:", "DSS Start Date=", "DSS Start Time=",
           "DSS End Date=", "DSS End Time=", "DSS Export Filename=",
           "DSS Export Rating Curves= 0 ", "DSS Export Rating Curve Sorted= 0 ",
           "DSS Export Volume Flow Curves= 0 ",
           "DXF Filename=", "DXF OffsetX= 0 ", "DXF OffsetY= 0 ",
           "DXF ScaleX= 1 ", "DXF ScaleY= 10 ", "GIS Export Profiles= 0 "]
    escrever(os.path.join(pasta, f"{nome}.prj"), "\n".join(prj))

    p01 = [f"Plan Title={nome}", "Program Version=7.01",
           f"Short Identifier={nome[:16]:<16}",
           "Simulation Date=01AUG2026,0000,08AUG2026,2300",
           "Geom File=g01", "Flow File=u01", "Mixed Flow Regime",
           "UNET Froude Reduction=True", "UNET Froude Limit= 0.8 ",
           "UNET Froude Power= 4 ", "UNET ZTol= 0.02 ", "UNET ZSATol= 0.02 ",
           "UNET MxIter= 40 ", "Computation Interval=15SEC",
           "Output Interval=1HOUR", "Instantaneous Interval=1HOUR",
           "Mapping Interval=1HOUR", "Run HTab=-1", "Run UNet=-1",
           "Run PostProcess=-1", "Run RASMapper=-1", "UNET D1 Cores= 0 "]
    escrever(os.path.join(pasta, f"{nome}.p01"), "\n".join(p01))

    rs0, rs1 = S[0]["rs"], S[-1]["rs"]
    u01 = [f"Flow Title={nome}", "Program Version=7.01", "Use Restart= 0 ",
           f"Initial RS={rio:<16.16},{rch:<16.16},{rs0:.2f},{Q[0]:.2f}",
           "",
           f"Boundary Location={rio:<16.16},{rch:<16.16},{rs0:<8.2f},"
           f"{'':<8},{'':<16},{'':<16}",
           "Interval=1HOUR", f"Flow Hydrograph= {len(Q)} "]
    u01 += col8(Q)
    u01 += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
            "Fixed Start Date/Time=01AUG2026,0000",
            "Is Critical Boundary=False", "Critical Boundary Flow=",
            "",
            f"Boundary Location={rio:<16.16},{rch:<16.16},{rs1:<8.2f},"
            f"{'':<8},{'':<16},{'':<16}",
            f"Friction Slope={decl:.6f},0"]
    escrever(os.path.join(pasta, f"{nome}.u01"), "\n".join(u01))

    # ---- PROJECAO E .rasmap. Sem eles o projeto "roda" mas nao tem sistema
    # de coordenadas: o RAS Mapper nao desenha, o CompleteGeometryCommand nao
    # monta o HDF da geometria (falha com "Geometry not found in
    # WriteAttributePreCheck") e nada pode ser validado sem antes simular.
    # O `.prj` de projecao e um arquivo ESRI, distinto do `.prj` do projeto.
    import shutil
    fonte = None
    for c in ("SIRGAS2000_UTM22S.prj", "modelo/SIRGAS2000_UTM22S.prj",
              "modelo/mirim_novo/SIRGAS2000_UTM22S.prj"):
        if os.path.exists(c):
            fonte = c
            break
    if fonte is None:
        raise SystemExit("nao achei SIRGAS2000_UTM22S.prj para a projecao")
    shutil.copy2(fonte, os.path.join(pasta, "SIRGAS2000_UTM22S.prj"))
    print(f"projecao  : {fonte} -> {nome}/SIRGAS2000_UTM22S.prj")

    P = np.vstack([np.asarray(d["cut"], float) for d in S])
    rasmap = [
        "<RASMapper>", "  <Version>2.0.0</Version>",
        '  <RASProjectionFilename Filename=".\\SIRGAS2000_UTM22S.prj" />',
        '  <Geometries Checked="True" Expanded="True">',
        f'    <Layer Name="{nome}" Type="RASGeometry" Checked="True" '
        f'Expanded="True" Filename=".\\{nome}.g01.hdf">',
        '      <Layer Type="RASRiver" Checked="True" />',
        '      <Layer Type="RASBankLines" Checked="True" />',
        '      <Layer Type="FlowPaths" Checked="True" />',
        '      <Layer Type="RiverStations" Checked="True" />',
        '      <Layer Type="RASXS" Checked="True" />',
        '      <Layer Type="RASErrors" Checked="True" />',
        "    </Layer>", "  </Geometries>", "  <Terrains />",
        "  <CurrentView>",
        f"    <MinX>{P[:,0].min():.2f}</MinX>",
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
    for e in ("prj", "p01", "u01", "g01", "rasmap"):
        p = os.path.join(pasta, f"{nome}.{e}")
        ok = os.path.exists(p)
        print(f"   {'OK   ' if ok else 'FALTA'} {nome}.{e}"
              + (f"   {os.path.getsize(p)} bytes" if ok else ""))
    return pasta


if __name__ == "__main__":
    main()
