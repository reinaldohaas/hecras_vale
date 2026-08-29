# -*- coding: utf-8 -*-
"""Mede largura de LAMINA e de CALHA PLENA dos rios no MDT 1 m do SIG-SC.

    python scripts/largura_do_sigsc.py taha_ai.g01 --cada 500 \
        --saida doc/larguras_sigsc

SO MEDE. Nao altera geometria nenhuma; sai um CSV por rio.

COMO

  Em cada transecto perpendicular ao eixo (passo `--cada` m, meia-largura
  de 250 m, amostra de 1 m):

    LAMINA       faixa continua em torno do talvegue com cota dentro de
                 20 cm do espelho local -- a agua do dia do voo. E o PISO
                 da largura real do canal.
    CALHA PLENA  a partir da lamina, sobe-se cada margem ate a QUEBRA DE
                 BARRANCO: primeiro ponto em que a inclinacao local cai
                 abaixo de 15% depois de ter subido ao menos 1 m acima do
                 espelho. E a largura de margens plenas -- a que os
                 `Bank Sta` do modelo deveriam declarar.

  Transectos sem espelho detectavel (mata fechada, sombra, eixo fora
  d'agua) saem marcados e nao entram nas medianas.

  CONFINAMENTO (auditoria de 26/08): lamina so vale se as DUAS margens
  sobem ao menos SUBIDA_MIN dentro de CONFINAMENTO m alem da borda
  d'agua. Sem isso o flat e arrozeira, banhado, reservatorio (Barragem
  Norte) ou estuario -- sai marcado "sem margem" e fora das medianas.

POR QUE

  O gerador declarou canais por formula de regime (w = 5*A^0.4): no
  Benedito deu 117-336 m onde o rio tem ~30-50. O SIG-SC de 1 m mede a
  largura de verdade, rio a rio, sem pedir nada a ninguem.
"""
import csv
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_geometria import ler_eixos     # noqa: E402

MDT = "taha_ai_novo/Terrain/taha_ai_corredor_1m_completo.tif"
MEIA = 250          # m para cada lado do eixo
TOL_LAMINA = 0.20   # m acima do espelho que ainda e agua
SUBIDA_MIN = 1.0    # m acima do espelho antes de aceitar quebra
DECL_QUEBRA = 0.15  # inclinacao local abaixo disto = topo de barranco
CONFINAMENTO = 30   # m alem da borda d'agua em que a margem TEM de subir

# TETO de plausibilidade (26/08, conhecimento de campo do Reinaldo): no
# alto vale ha arrozeiras, represas e canais LONGITUDINAIS paralelos ao
# rio que passam ate na regra do confinamento -- e dificilmente algum rio
# de la passa de 100 m. Lamina acima do teto sai censurada ("acima do
# teto"). Excecao: o baixo Acu estuarino e largo de verdade.
TETO_PADRAO = 100.0
TETO = {"Itajai_Acu": 400.0}


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def transecto(src, T, P, n):
    xs = P[0] + np.arange(-MEIA, MEIA + 1) * n[0]
    ys = P[1] + np.arange(-MEIA, MEIA + 1) * n[1]
    cc = ((xs - T.c) / T.a).astype(int)
    rr = ((ys - T.f) / T.e).astype(int)
    z = np.full(len(xs), np.nan)
    from rasterio.windows import Window
    ok = (rr >= 0) & (rr < src.height) & (cc >= 0) & (cc < src.width)
    for k in np.flatnonzero(ok):
        v = float(src.read(1, window=Window(cc[k], rr[k], 1, 1))[0, 0])
        z[k] = np.nan if v == src.nodata else v
    return z


def medir(z):
    """(larg_lamina, larg_calha, espelho) | None sem agua | "solto"
    se a agua nao e confinada por margens (arrozeira/banhado/acude)."""
    c = len(z) // 2
    jan = z[c - 25:c + 26]
    if np.all(np.isnan(jan)):
        return None
    esp = np.nanmin(jan)
    i0 = c - 25 + int(np.nanargmin(jan))
    agua = np.abs(z - esp) <= TOL_LAMINA
    a = i0
    while a > 0 and agua[a - 1]:
        a -= 1
    b = i0
    while b < len(z) - 1 and agua[b + 1]:
        b += 1
    lamina = b - a + 1

    def confinada(borda, passo):
        trecho = z[max(borda - CONFINAMENTO, 0):borda + 1] if passo < 0 \
            else z[borda:borda + CONFINAMENTO + 1]
        return np.nanmax(trecho) - esp >= SUBIDA_MIN \
            if not np.all(np.isnan(trecho)) else False

    if not (confinada(a, -1) and confinada(b, +1)):
        return "solto"

    def margem(ini, passo):
        j = ini
        while 0 < j < len(z) - 1:
            j += passo
            if np.isnan(z[j]):
                break
            if z[j] - esp >= SUBIDA_MIN:
                # inclinacao local numa janela de 5 m
                j0, j1 = max(j - 2, 0), min(j + 3, len(z))
                tramo = z[j0:j1]
                if np.all(~np.isnan(tramo)):
                    decl = abs(np.polyfit(np.arange(len(tramo)),
                                          tramo, 1)[0])
                    if decl < DECL_QUEBRA:
                        return j
        return j

    ja = margem(a, -1)
    jb = margem(b, +1)
    return lamina, jb - ja + 1, float(esp)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import rasterio
    g01 = argv[0]
    cada = _arg(argv, "--cada", 500.0, float)
    pasta = _arg(argv, "--saida", "doc/larguras_sigsc")
    os.makedirs(pasta, exist_ok=True)
    src = rasterio.open(MDT)
    T = src.transform

    E = ler_eixos(g01)
    por_rio = {}
    for (rio, reach), ls in E.items():
        por_rio.setdefault(rio, []).append((reach, ls))

    print(f"{'rio':16s} {'transectos':>10s} {'lamina med':>10s} "
          f"{'calha med':>10s} {'calha p90':>10s}")
    for rio, partes in sorted(por_rio.items()):
        partes.sort(key=lambda t: t[0])
        coords = []
        for _, ls in partes:
            coords += list(ls.coords)
        from shapely.geometry import LineString
        eixo = LineString(coords)
        linhas, lam, cal = [], [], []
        L = eixo.length
        for s in np.arange(200, L - 200, cada):
            P0 = np.asarray(eixo.interpolate(s).coords[0])
            P1 = np.asarray(eixo.interpolate(min(s + 30, L)).coords[0])
            t = P1 - P0
            t = t / max(np.hypot(*t), 1e-9)
            nvec = np.array([-t[1], t[0]])
            z = transecto(src, T, P0, nvec)
            m = medir(z)
            if m is None:
                linhas.append([f"{(L-s)/1000:.2f}", "", "", "", "sem agua"])
                continue
            if m == "solto":
                linhas.append([f"{(L-s)/1000:.2f}", "", "", "",
                               "sem margem"])
                continue
            la, ca, esp = m
            if la > TETO.get(rio, TETO_PADRAO):
                linhas.append([f"{(L-s)/1000:.2f}", "", "", f"{esp:.2f}",
                               "acima do teto"])
                continue
            linhas.append([f"{(L-s)/1000:.2f}", la, ca, f"{esp:.2f}", ""])
            lam.append(la)
            cal.append(ca)
        arq = os.path.join(pasta, f"{rio}.csv")
        with open(arq, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(["dist_foz_km", "largura_lamina_m",
                        "largura_calha_plena_m", "cota_espelho", "obs"])
            w.writerows(linhas)
        if lam:
            print(f"{rio:16s} {len(linhas):10d} {np.median(lam):10.0f} "
                  f"{np.median(cal):10.0f} {np.percentile(cal, 90):10.0f}")
        else:
            print(f"{rio:16s} {len(linhas):10d} {'-':>10s} {'-':>10s} "
                  f"{'-':>10s}")
    print(f"\nCSVs por rio em {pasta}/")


if __name__ == "__main__":
    main(sys.argv[1:])
