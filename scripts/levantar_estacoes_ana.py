# -*- coding: utf-8 -*-
"""Levanta as estacoes da ANA nos rios do vale: o que existe, onde, e quando.

    python scripts/levantar_estacoes_ana.py --saida doc/estacoes_ana_vale.csv

Consulta o inventario PUBLICO (HidroInventario, sem token) para cada rio do
modelo e tabula, por estacao fluviometrica:

    codigo, nome, rio, municipio, area de drenagem,
    escala (regua), registrador de nivel,
    DESCARGA LIQUIDA (medicoes de vazao -- e a materia-prima da CURVA-CHAVE:
    onde ha descarga liquida ha curva calibrada e largura/area medidas),
    sedimentos, qualidade d'agua, telemetria,
    periodos de operacao de escala e de descarga.

O perfil transversal e a propria curva-chave exigem o token do
HidroWebService (cadastro do usuario); este script mapeia ONDE pedir
quando o token sair.
"""
import csv
import os
import re
import sys
import urllib.request
import urllib.parse

URL = ("http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroInventario"
       "?codEstDE=&codEstATE=&tpEst=1&nmEst=&nmRio={rio}"
       "&codSubBacia=&codBacia=&nmMunicipio=&nmEstado="
       "&sgResp=&sgOper=&telemetrica=")

# FRAGMENTOS sem acento: o casamento do HidroInventario e sensivel a
# acento ("ITAJAI" nao acha "ITAJAÍ-AÇU"; "ITAJA" acha)
RIOS = ["ITAJA", "HERC", "BENEDITO", "CEDROS", "TROMBUDO", "TAIO",
        "POMBAS", "TESTO", "IRAPUT", "ALVES"]

CAMPOS = ["Codigo", "Nome", "RioNome", "nmMunicipio", "nmEstado",
          "AreaDrenagem", "TipoEstacaoEscala",
          "TipoEstacaoRegistradorNivel", "TipoEstacaoDescLiquida",
          "TipoEstacaoSedimentos", "TipoEstacaoQualAgua",
          "TipoEstacaoTelemetrica", "PeriodoEscalaInicio",
          "PeriodoEscalaFim", "PeriodoDescLiquidaInicio",
          "PeriodoDescLiquidaFim", "ResponsavelSigla", "OperadoraSigla",
          "Latitude", "Longitude"]


def consulta(rio):
    url = URL.format(rio=urllib.parse.quote(rio))
    with urllib.request.urlopen(url, timeout=90) as r:
        return r.read().decode("utf-8", errors="replace")


def blocos(xml):
    for b in re.findall(r"<Table[^>]*>(.*?)</Table>", xml, flags=re.S):
        d = {}
        for k in CAMPOS:
            m = re.search(r"<%s>([^<]*)</%s>" % (k, k), b)
            d[k] = (m.group(1).strip() if m else "")
        yield d


def main(argv):
    saida = "doc/estacoes_ana_vale.csv"
    if "--saida" in argv:
        saida = argv[argv.index("--saida") + 1]
    vistos = {}
    for rio in RIOS:
        try:
            xml = consulta(rio)
        except Exception as e:
            print(f"   {rio}: FALHOU ({e})")
            continue
        n = 0
        for d in blocos(xml):
            if d["nmEstado"].upper() != "SANTA CATARINA":
                continue
            vistos[d["Codigo"]] = d
            n += 1
        print(f"   {rio:16s}: {n} estacoes em SC")
    os.makedirs(os.path.dirname(saida) or ".", exist_ok=True)
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, CAMPOS, delimiter=";")
        w.writeheader()
        for d in sorted(vistos.values(), key=lambda d: d["Codigo"]):
            w.writerow(d)
    print(f"\n{len(vistos)} estacoes unicas -> {saida}")

    # resumo do que interessa a curva-chave
    com_dl = [d for d in vistos.values()
              if d["TipoEstacaoDescLiquida"] == "1"]
    print(f"\ncom DESCARGA LIQUIDA (base da curva-chave): {len(com_dl)}")
    for d in sorted(com_dl, key=lambda d: -float(d["AreaDrenagem"] or 0)):
        print(f"   {d['Codigo']} {d['Nome'][:28]:28s} {d['RioNome'][:22]:22s}"
              f" area {d['AreaDrenagem']:>7s} km2  descarga "
              f"{d['PeriodoDescLiquidaInicio'][:4]}-"
              f"{d['PeriodoDescLiquidaFim'][:4]}"
              f"  telem={'S' if d['TipoEstacaoTelemetrica']=='1' else 'n'}")


if __name__ == "__main__":
    main(sys.argv[1:])
