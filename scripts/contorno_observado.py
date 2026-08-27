# -*- coding: utf-8 -*-
"""Troca os contornos sinteticos pelo OBSERVADO de 1983 (ANA), por sistema.

    python scripts/contorno_observado.py taha_ai --series doc/ana_1983

O .u01 e o .p01 do projeto SAO EDITADOS, com backup `.antes_do_observado`.

COMO A AGUA E DISTRIBUIDA (sem contagem dupla)

  Cada SISTEMA de rio tem uma estacao de 1983 e uma razao de areas da
  propria ANA ate a foz do sistema:

    Norte+Iraputa   Ibirama (3330)          x 1.00  (estacao ~na foz)
    Oeste+Taio+
    Trombudo+Pombas Taio (1570)             x 1.47  (foz 2300, UHE SP Oeste)
    Sul             Ituporanga (1650)       x 1.21  (foz 1990, UHE SP Aurora)
    Mirim           Brusque (1240)          x 1.35  (foz ~1676, area da bacia)
    Benedito        Benedito Novo (717)     x 1.48  (foz 1600 MENOS os 536
                                                     do Cedros = 1064)
    Cedros          Arrozeira (536)         x 1.00  (estacao ~na foz)
    Testo           forma do Benedito Novo, pico de projeto mantido
    Acu (laterais)  forma de Blumenau, pico de projeto mantido

  O total do sistema (observado x razao) e repartido entre a cabeceira e
  as laterais do sistema NA PROPORCAO DOS PICOS SINTETICOS atuais -- a
  distribuicao espacial do projeto fica, o tempo e o volume passam a ser
  os observados. Series diarias entram como Interval=1DAY.

  A mare vira a M2 sintetica do gerador legado (media 0,30 m, semi-
  amplitude 0,50 m, periodo 12,42 h), horaria, na janela do evento.

  CONFERENCIA: imprime o fechamento em Blumenau -- soma de tudo que entra
  a montante vs os 5.274 m3/s observados (difere pelo amortecimento; e
  régua, nao meta).
"""
import csv
import datetime
import math
import os
import re
import shutil
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from ras_io import escrever            # noqa: E402


# 01/07 mesmo: aquecer desde 20/06 em vazao baixa foi TESTADO (27/08) e
# explodiu na hora 17 com erro de 5511% -- vazao baixa e o pior regime
# desta rede (lamina fina), nao um assentamento
INICIO = datetime.date(1983, 7, 1)
FIM = datetime.date(1983, 8, 5)

# sistema -> (rios do modelo no sistema, csv da estacao, razao de areas)
SISTEMAS = {
    "norte":    (["Itajai_Norte", "Rio_Iraputa"],
                 "83440000_Ibirama_vazao.csv", 1.00),
    "oeste":    (["Itajai_Oeste", "Rio_Taio", "Rio_Trombudo",
                  "Rio_das_Pombas"],
                 "83050000_Taio_vazao.csv", 2300.0 / 1570.0),
    "sul":      (["Itajai_Sul"],
                 "83250000_Ituporanga_vazao.csv", 1990.0 / 1650.0),
    "mirim":    (["Itajai_Mirim"],
                 "83900000_Brusque_vazao.csv", 1676.0 / 1240.0),
    "benedito": (["Rio_Benedito"],
                 "83660000_Benedito_Novo_vazao.csv", 1064.0 / 717.0),
    "cedros":   (["Rio_dos_Cedros"],
                 "83675000_Arrozeira_vazao.csv", 1.00),
}
# rios sem estacao em 1983: forma do doador, PICO DE PROJETO mantido
FORMA_DOADA = {
    "Rio_do_Testo": "83660000_Benedito_Novo_vazao.csv",
    "Itajai_Acu":   "83800002_Blumenau_vazao.csv",
}

# fracao do total do sistema que entra pela CABECEIRA do tronco principal
# = area de drenagem na estacao / area na foz (ANA). Repartir por picos
# sinteticos sub-alimentava as cabeceiras: Taio simulava 23 m3/s onde a
# regua marcou 922 (medido em 27/08, rodada de 212 h). O resto do sistema
# segue nos outros contornos na proporcao dos picos sinteticos.
FRACAO_TOPO = {
    "Itajai_Oeste":   1570.0 / 2300.0,   # Taio no topo do dominio
    "Itajai_Sul":      434.0 / 1990.0,   # Saltinho perto do topo
    "Itajai_Mirim":   1240.0 / 1676.0,   # topo amputado ~Brusque
    "Rio_Benedito":    717.0 / 1064.0,   # topo amputado ~Benedito Novo
    "Rio_dos_Cedros":            1.00,   # topo amputado ~Arrozeira (=foz)
}
MARE_MEDIA, MARE_AMP, MARE_T = 0.30, 0.50, 12.42


def serie_obs(pasta, arq):
    datas, vals = [], []
    for r in csv.reader(open(os.path.join(pasta, arq), encoding="utf-8"),
                        delimiter=";"):
        if r[0] == "data":
            continue
        d = datetime.date.fromisoformat(r[0])
        if INICIO <= d <= FIM:
            datas.append(d)
            vals.append(float(r[1]))
    dias = (FIM - INICIO).days + 1
    serie = np.full(dias, np.nan)
    for d, v in zip(datas, vals):
        serie[(d - INICIO).days] = v
    # falha pontual: interpola; borda: repete
    idx = np.arange(dias)
    ok = ~np.isnan(serie)
    if not ok.all():
        serie = np.interp(idx, idx[ok], serie[ok])
    return serie


def fmt_serie(vals):
    corpo, lin = [], ""
    for i, x in enumerate(vals):
        lin += "%8.2f" % x
        if (i + 1) % 10 == 0:
            corpo.append(lin)
            lin = ""
    if lin:
        corpo.append(lin)
    return corpo


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    base = argv[0].rstrip("/\\")
    pasta = argv[argv.index("--series") + 1] if "--series" in argv \
        else "doc/ana_1983"
    u01 = base + ".u01"
    p01 = base + ".p01"
    shutil.copy2(u01, u01 + ".antes_do_observado")
    shutil.copy2(p01, p01 + ".antes_do_observado")

    t = open(u01, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n")
    blocos = re.split(r"(?=^Boundary Location=)", t, flags=re.M)

    # picos sinteticos atuais por bloco
    info = []
    for k, b in enumerate(blocos):
        m = re.match(r"Boundary Location=([^,]+),", b)
        if not m:
            info.append(None)
            continue
        h = re.search(r"(Flow Hydrograph|Uniform Lateral Inflow Hydrograph"
                      r"|Stage Hydrograph)=\s*(\d+)", b)
        if not h:
            info.append(None)
            continue
        vals = []
        for l in b[h.end():].split("\n")[1:]:
            if not l.strip() or l[:1].isalpha():
                break
            vals += [float(l[i:i + 8]) for i in range(0, len(l), 8)
                     if l[i:i + 8].strip()]
        info.append({"rio": m.group(1).strip(), "tipo": h.group(1),
                     "pico": max(vals) if vals else 0.0})

    rio_sistema = {}
    for s, (rios, arq, razao) in SISTEMAS.items():
        for r in rios:
            rio_sistema[r] = s

    soma_sist = {}
    pico_topo = {}
    for i in info:
        if i and i["tipo"] != "Stage Hydrograph" \
                and i["rio"] in rio_sistema:
            s = rio_sistema[i["rio"]]
            soma_sist[s] = soma_sist.get(s, 0.0) + i["pico"]
            if i["tipo"] == "Flow Hydrograph" \
                    and i["rio"] == SISTEMAS[s][0][0]:
                pico_topo[s] = pico_topo.get(s, 0.0) + i["pico"]

    obs = {s: serie_obs(pasta, arq) * razao
           for s, (rios, arq, razao) in SISTEMAS.items()}
    formas = {r: serie_obs(pasta, a) for r, a in FORMA_DOADA.items()}
    dias = (FIM - INICIO).days + 1
    data_ini = INICIO.strftime("%d%b1983").upper()

    total_entrando = 0.0
    novo_blocos = []
    for k, b in enumerate(blocos):
        i = info[k]
        if i is None:
            novo_blocos.append(b)
            continue
        rio, tipo = i["rio"], i["tipo"]
        if tipo == "Stage Hydrograph":
            horas = dias * 24
            mare = [MARE_MEDIA + MARE_AMP *
                    math.sin(2 * math.pi * h / MARE_T) for h in range(horas)]
            corpo = ["Interval=1HOUR",
                     "Stage Hydrograph= %d " % horas] + fmt_serie(mare)
        else:
            if rio in rio_sistema:
                s = rio_sistema[rio]
                f = FRACAO_TOPO.get(SISTEMAS[s][0][0])
                topo = (tipo == "Flow Hydrograph"
                        and rio == SISTEMAS[s][0][0])
                if f is None:
                    fr = i["pico"] / max(soma_sist[s], 1e-9)
                elif topo:
                    fr = f
                else:
                    fr = (1 - f) * i["pico"] / max(
                        soma_sist[s] - pico_topo.get(s, 0.0), 1e-9)
                serie = obs[s] * fr
            elif rio in formas:
                f = formas[rio]
                serie = f / max(f.max(), 1e-9) * i["pico"]
            else:
                novo_blocos.append(b)
                continue
            total_entrando += float(serie.max())
            chave = ("Flow Hydrograph" if tipo == "Flow Hydrograph"
                     else "Uniform Lateral Inflow Hydrograph")
            corpo = ["Interval=1DAY",
                     "%s= %d " % (chave, dias)] + fmt_serie(serie)
        cab = b.split("\n")[0]
        # o RAS deste u01 ancora a serie no INICIO DO PLANO (o formato
        # com data fixa nao parseia aqui: "unable to resolve 01JUL1983")
        resto = ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=False",
                 "Fixed Start Date/Time=,"]
        if tipo == "Flow Hydrograph":
            resto += ["Flow Hydrograph Slope= 0.001 "]
        novo_blocos.append("\n".join([cab] + corpo + resto) + "\n\n")

    escrever(u01, "".join(novo_blocos))

    tp = open(p01, encoding="latin-1", errors="replace").read()
    tp = re.sub(r"Simulation Date=.*",
                "Simulation Date=%s,0000,31JUL1983,2300" % data_ini, tp, 1)
    escrever(p01, tp)

    print(f"contornos observados gravados em {u01}")
    print(f"janela: 01JUL1983 - 31JUL1983 (diario; mare M2 horaria)")
    print(f"\npicos por sistema (observado x razao de areas):")
    for s in SISTEMAS:
        print(f"   {s:9s}: {obs[s].max():7.1f} m3/s")
    print(f"\nsoma dos picos entrando no modelo: {total_entrando:.0f} m3/s")
    print("observado em Blumenau (11.803 km2) : 5274 m3/s  <- regua de "
          "fechamento (difere pelo amortecimento e defasagem)")


if __name__ == "__main__":
    main(sys.argv[1:])
