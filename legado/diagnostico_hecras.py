# -*- coding: utf-8 -*-
"""
Diz ONDE e POR QUE a simulacao do HEC-RAS divergiu.

Rodando pela interface COM o HEC-RAS nao grava o .computeMsgs.txt, que e onde
ficam as mensagens boas ("Solution Solver Failed" com a secao exata). Este
modulo tira a mesma informacao do proprio .p01.hdf, que sempre e escrito:

  - o instante em que parou e o passo em que a lamina explodiu;
  - as secoes com maior salto de cota entre saidas consecutivas;
  - secoes com agua ACIMA do topo da secao (que e o "Extrapolated above Cross
    Section Table" do log, e costuma ser causa, nao sintoma);
  - secoes secas e vazoes negativas, que sinalizam trecho drenando;
  - Froude por secao no ultimo instante, para achar transicao critica.

Uso:  python diagnostico_hecras.py [PROJETO]
"""
import os
import sys

import numpy as np
import h5py

BASE = ("Results/Unsteady/Output/Output Blocks/Base Output/"
        "Unsteady Time Series/Cross Sections")


def carregar(projeto):
    with h5py.File(f"{projeto}.p01.hdf", "r") as f:
        sol = f["Results/Unsteady"].attrs.get("Solution")
        sol = sol.decode() if isinstance(sol, bytes) else str(sol)
        g = f[BASE]
        ws = g["Water Surface"][:]
        q = g["Flow"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = np.array([x["River"].decode().strip() for x in at])
        rch = np.array([x["Reach"].decode().strip() for x in at])
        rs = np.array([float(x["RS"].decode()) for x in at])
        se = f["Geometry/Cross Sections/Station Elevation Values"][:]
        info = f["Geometry/Cross Sections/Station Elevation Info"][:]
    z_min = np.array([se[a:a + n, 1].min() for a, n in info])
    z_max = np.array([se[a:a + n, 1].max() for a, n in info])
    larg = np.array([np.ptp(se[a:a + n, 0]) for a, n in info])
    return dict(sol=sol, ws=ws, q=q, riv=riv, rch=rch, rs=rs,
                z_min=z_min, z_max=z_max, larg=larg)


def rotulo(d, k):
    return f"{d['riv'][k]:<14}{d['rch'][k]:<4}RS {d['rs'][k]/1000:8.2f} km"


def relatorio(projeto):
    d = carregar(projeto)
    ws, q = d["ws"], d["q"]
    nt = ws.shape[0]
    ok = "Success" in d["sol"]
    print("=" * 72)
    print(f"{projeto}   Solution = {d['sol']}   {nt} saidas x {ws.shape[1]} secoes")
    print("=" * 72)
    if ok:
        print("  simulacao concluiu.")

    prof = ws - d["z_min"][None, :]

    # --- 1. salto de cota entre as duas ultimas saidas: onde explodiu
    if nt >= 2:
        salto = ws[-1] - ws[-2]
        j = np.argsort(-np.abs(salto))[:10]
        print(f"\n[1] MAIOR SALTO DE COTA na ultima saida (t={nt-1})")
        print(f"    {'secao':<34}{'salto m':>9}{'cota':>9}{'prof':>8}{'Q':>9}")
        for k in j:
            print(f"    {rotulo(d,k):<34}{salto[k]:9.2f}{ws[-1,k]:9.2f}"
                  f"{prof[-1,k]:8.2f}{q[-1,k]:9.0f}")

    # --- 2. agua acima do TOPO da secao (o 'Extrapolated' do log)
    acima = ws.max(axis=0) - d["z_max"]
    n_ac = int((acima > 0).sum())
    print(f"\n[2] AGUA ACIMA DO TOPO DA SECAO: {n_ac} de {len(acima)}")
    if n_ac:
        for k in np.argsort(-acima)[:8]:
            if acima[k] <= 0:
                break
            print(f"    {rotulo(d,k):<34}+{acima[k]:6.2f} m acima do topo"
                  f"   (largura {d['larg'][k]:6.0f} m)")

    # --- 3. secoes secas e vazao negativa
    seca = (prof < 0.05).any(axis=0)
    neg = (q < -1.0).any(axis=0)
    print(f"\n[3] SECOES QUE SECAM (prof < 5 cm): {int(seca.sum())}")
    for k in np.flatnonzero(seca)[:6]:
        t0 = int(np.argmin(prof[:, k]))
        print(f"    {rotulo(d,k):<34}min {prof[:,k].min():5.2f} m na saida {t0}")
    print(f"\n[4] SECOES COM VAZAO NEGATIVA (< -1 m3/s): {int(neg.sum())}")
    for k in np.flatnonzero(neg)[:6]:
        print(f"    {rotulo(d,k):<34}Q min {q[:,k].min():9.0f} m3/s")

    # --- 5. Froude aproximado no ultimo instante (prof media x largura)
    with np.errstate(all="ignore"):
        h = np.maximum(prof[-1], 0.01)
        b = np.maximum(d["larg"], 1.0)
        v = np.abs(q[-1]) / (b * h)
        fr = v / np.sqrt(9.80665 * h)
    j = np.argsort(-fr)[:8]
    print(f"\n[5] MAIOR FROUDE na ultima saida  "
          f"({int((fr>0.9).sum())} secoes com Fr > 0,9)")
    for k in j:
        print(f"    {rotulo(d,k):<34}Fr {fr[k]:6.2f}   V {v[k]:6.2f} m/s"
              f"   h {h[k]:6.2f} m")
    return d


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else "Itajai_Rede_1983"
    if not os.path.exists(f"{p}.p01.hdf"):
        raise SystemExit(f"{p}.p01.hdf nao encontrado")
    relatorio(p)
