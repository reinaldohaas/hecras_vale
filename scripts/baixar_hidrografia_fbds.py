# -*- coding: utf-8 -*-
"""Baixa a hidrografia vetorial da FBDS para os municipios do vale.

    python scripts/baixar_hidrografia_fbds.py --saida doc/fbds

FONTE  https://geo.fbds.org.br/SC/<MUNICIPIO>/HIDROGRAFIA/
       (mapeamento RapidEye 5 m, fev/2018, acesso publico)

Por municipio vem:
    RIOS_DUPLOS   poligono do rio onde ha margem dupla -- A LARGURA REAL
    MASSAS_DAGUA  lagos, acudes e REPRESAS (para excluir das larguras)

O codigo IBGE do arquivo (SC_4202404_...) e descoberto lendo o indice do
diretorio; nada e hardcoded. So baixa o que ainda nao existe no disco.
"""
import os
import re
import sys
import urllib.request

BASE = "https://geo.fbds.org.br/SC/{mun}/HIDROGRAFIA/"

MUNICIPIOS = [
    # baixo vale / Acu
    "ITAJAI", "NAVEGANTES", "ILHOTA", "LUIS_ALVES", "GASPAR", "BLUMENAU", "INDAIAL",
    "ASCURRA", "APIUNA",
    # Benedito / Cedros / Testo
    "TIMBO", "POMERODE", "RIO_DOS_CEDROS", "BENEDITO_NOVO",
    "DOUTOR_PEDRINHO",
    # Norte (Hercilio)
    "IBIRAMA", "JOSE_BOITEUX", "VITOR_MEIRELES", "DONA_EMMA",
    "PRESIDENTE_GETULIO", "WITMARSUM",
    # alto vale (Oeste / Sul / Trombudo / Taio / Pombas)
    "LONTRAS", "RIO_DO_SUL", "AGRONOMICA", "TROMBUDO_CENTRAL",
    "LAURENTINO", "RIO_DO_OESTE", "TAIO", "MIRIM_DOCE", "POUSO_REDONDO",
    "ITUPORANGA", "AURORA", "AGROLANDIA", "ATALANTA",
    # cabeceiras (Sul ate Alfredo Wagner/Bom Retiro; Oeste ate
    # Santa Cecilia/Rio do Campo; Norte ate Itaiopolis/Papanduva)
    "ALFREDO_WAGNER", "BOM_RETIRO", "CHAPADAO_DO_LAGEADO", "IMBUIA", "PETROLANDIA", "LEOBERTO_LEAL", "SANTA_CECILIA", "RIO_DO_CAMPO", "SALETE", "MONTE_CASTELO", "PAPANDUVA", "ITAIOPOLIS",
    # Mirim
    "VIDAL_RAMOS", "PRESIDENTE_NEREU", "BOTUVERA", "BRUSQUE", "GUABIRUBA",
]
CAMADAS = ["RIOS_DUPLOS", "MASSAS_DAGUA", "RIOS_SIMPLES"]
PARTES = [".shp", ".shx", ".dbf", ".prj", ".cpg"]


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def baixar(url, destino):
    with urllib.request.urlopen(url, timeout=180) as r, \
            open(destino, "wb") as f:
        f.write(r.read())


def main(argv):
    pasta = _arg(argv, "--saida", "doc/fbds")
    os.makedirs(pasta, exist_ok=True)
    ok, sem = [], []
    for mun in MUNICIPIOS:
        url_dir = BASE.format(mun=mun)
        try:
            with urllib.request.urlopen(url_dir, timeout=120) as r:
                indice = r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"   {mun:20s} FALHOU indice ({e})")
            sem.append(mun)
            continue
        m = re.search(r"SC_(\d+)_RIOS_DUPLOS\.shp", indice)
        if not m:
            print(f"   {mun:20s} sem RIOS_DUPLOS no indice")
            sem.append(mun)
            continue
        cod = m.group(1)
        dest_mun = os.path.join(pasta, mun)
        os.makedirs(dest_mun, exist_ok=True)
        n = 0
        for cam in CAMADAS:
            for ext in PARTES:
                arq = f"SC_{cod}_{cam}{ext}"
                destino = os.path.join(dest_mun, arq)
                if os.path.exists(destino) and os.path.getsize(destino) > 0:
                    continue
                if arq not in indice and ext != ".cpg":
                    continue
                try:
                    baixar(url_dir + arq, destino)
                    n += 1
                except Exception:
                    if ext in (".shp", ".dbf", ".shx"):
                        print(f"      {mun}: falta {arq}")
        ok.append(mun)
        print(f"   {mun:20s} cod {cod}  ({n} arquivos novos)")
    print(f"\n{len(ok)} municipios em {pasta}/; falharam: {sem or 'nenhum'}")


if __name__ == "__main__":
    main(sys.argv[1:])
