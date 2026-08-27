# -*- coding: utf-8 -*-
"""Baixa series diarias (vazao e cota) do servico PUBLICO da ANA.

    python scripts/baixar_series_ana.py --inicio 01/06/1983 --fim 31/08/1983 \
        --saida doc/ana_1983

Usa o HidroSerieHistorica (telemetriaws1, sem token). As series vem em
registros MENSAIS com colunas Vazao01..Vazao31 / Cota01..Cota31; aqui
viram um CSV longo (data;valor;consistencia) por estacao e tipo.

As estacoes padrao sao as 11 do vale com medicao de descarga ANTERIOR a
1983 (levantadas por levantar_estacoes_ana.py) -- as que enxergaram a
enchente de julho de 1983 com curva-chave propria.
"""
import csv
import os
import re
import sys
import urllib.request

URL = ("http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroSerieHistorica"
       "?codEstacao={cod}&dataInicio={ini}&dataFim={fim}"
       "&tipoDados={tipo}&nivelConsistencia=")

ESTACOES = {
    "83800002": "Blumenau",
    "83300200": "Rio_do_Sul_Novo",
    "83440000": "Ibirama",
    "83345000": "Barra_do_Prata",
    "83250000": "Ituporanga",
    "83105000": "Saltinho",
    "83050000": "Taio",
    "83900000": "Brusque",
    "83660000": "Benedito_Novo",
    "83675000": "Arrozeira",
    "83880000": "Luiz_Alves",
}
TIPOS = {"3": "vazao", "1": "cota"}


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def baixar(cod, tipo, ini, fim):
    url = URL.format(cod=cod, ini=ini.replace("/", "%2F"),
                     fim=fim.replace("/", "%2F"), tipo=tipo)
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("utf-8", errors="replace")


def registros(xml, rotulo):
    """(data_iso, valor) dos campos <RotuloDD> de cada bloco mensal."""
    out = []
    for b in re.findall(r"<SerieHistorica[^>]*>(.*?)</SerieHistorica>",
                        xml, flags=re.S):
        m = re.search(r"<DataHora>(\d{4})-(\d{2})-", b)
        if not m:
            continue
        ano, mes = m.group(1), m.group(2)
        for d in range(1, 32):
            v = re.search(r"<%s%02d>([^<]+)</%s%02d>"
                          % (rotulo, d, rotulo, d), b)
            if v and v.group(1).strip():
                out.append((f"{ano}-{mes}-{d:02d}",
                            v.group(1).strip().replace(",", ".")))
    # o servico repete o mes por nivel de consistencia; consistido (2)
    # vem depois e SOBRESCREVE o bruto na ordem do arquivo
    dedup = {}
    for data, v in out:
        dedup[data] = v
    return sorted(dedup.items())


def main(argv):
    ini = _arg(argv, "--inicio", "01/06/1983")
    fim = _arg(argv, "--fim", "31/08/1983")
    pasta = _arg(argv, "--saida", "doc/ana_1983")
    os.makedirs(pasta, exist_ok=True)

    picos = []
    for cod, nome in ESTACOES.items():
        linha = {"estacao": f"{cod} {nome}"}
        for tipo, rot in TIPOS.items():
            try:
                xml = baixar(cod, tipo, ini, fim)
            except Exception as e:
                print(f"   {nome} {rot}: FALHOU ({e})")
                continue
            rotulo = "Vazao" if rot == "vazao" else "Cota"
            regs = registros(xml, rotulo)
            arq = os.path.join(pasta, f"{cod}_{nome}_{rot}.csv")
            with open(arq, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f, delimiter=";")
                w.writerow(["data", rot])
                w.writerows(regs)
            if regs:
                vals = [(float(v), d) for d, v in regs]
                pico, dia = max(vals)
                linha[rot] = f"{pico:9.1f} em {dia}"
                linha[f"n_{rot}"] = len(regs)
            else:
                linha[rot] = "SEM DADO"
        picos.append(linha)
        print(f"   {nome:16s} vazao: {linha.get('vazao','-'):24s} "
              f"cota: {linha.get('cota','-')}")
    print(f"\nCSVs em {pasta}/")


if __name__ == "__main__":
    main(sys.argv[1:])
