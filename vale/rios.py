# -*- coding: utf-8 -*-
"""
Catalogo de rios do Vale do Itajai, da BHO 2017 da ANA, e a selecao de quais
entram no modelo.

Da ANA vem o que ela tem de confiavel: QUEM desagua em QUEM e a AREA DE
DRENAGEM. A geometria (por onde o rio passa) vem do relevo -- misturar as duas
coisas quebra a rede, porque Sul, Oeste e Acu se encontram no MESMO ponto em
Rio do Sul e o criterio geometrico de "rio maior mais proximo" pendura o Sul no
Oeste em vez de no Acu.

NOME NORMALIZADO, e nao o nome cru. A base da ANA grava o mesmo rio com duas
grafias:

    "Rio Itajaí do Oeste"   3.007 km2   68 km   (jusante, COM acento)
    "Rio Itajai do Oeste"   1.103 km2   56 km   (montante, SEM acento)

Os dois se tocam -- distancia zero. Casar o nome por acento, como o gerador
anterior fazia, truncava o Itajai do Oeste em 68 km e deixava 56 km de
cabeceira fora do modelo. A area usada continuava certa (3.007 na foz), entao
a vazao entrava e nada acusava o erro: faltava so a geometria, silenciosamente.
O mesmo acontece com "Ribeirao Dollmann" e "Ribeirao Dollman".

Uso:
    python -m vale.rios                  # catalogo numerado
    python -m vale.rios --area 50        # outro limiar
    python -m vale.rios --sel 1,2,3      # confere uma selecao
    python -m vale.rios --sel todos
"""
import argparse
import re
import unicodedata

import geopandas as gpd

BASE = "rios_itajai.geojson"
EPSG = 31982
AREA_MINIMA = 100.0        # km2

# Os 12 rios do modelo atual, por nome normalizado. Servem de conjunto padrao
# ("atuais") para nao trocar o escopo sem que alguem tenha pedido.
ATUAIS = [
    "rio itajai-acu", "rio itajai do norte ou hercilio", "rio itajai do oeste",
    "rio itajai do sul", "rio itajai-mirim", "rio benedito", "rio dos cedros",
    "rio trombudo", "rio iraputa", "rio taio", "rio das pombas", "rio do testo",
]


def normalizar(nome):
    """Sem acento, minusculo, espacos colapsados. E a chave de agrupamento."""
    s = unicodedata.normalize("NFKD", str(nome))
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"\s+", " ", s).strip().lower()


def nome_ras(nome):
    """Nome aceito pelo HEC-RAS: ASCII, sem espaco, ate 16 caracteres.

    Os conectivos (do, da, dos, das, de) saem antes do corte, senao
    "Rio Itajai do Norte ou Hercilio" vira "Itajai_Do_Norte_" -- truncado no
    conectivo, com underscore solto no fim, e a 16 caracteres dois rios
    diferentes podem colidir no mesmo nome. Sem eles sai "Itajai_Norte", que
    e o nome usado no modelo atual.
    """
    s = normalizar(nome)
    s = re.sub(r"^(rio|ribeirao|corrego|arroio)\s+", "", s)
    s = re.sub(r"\s+ou\s+.*$", "", s)          # "X ou Y" -> X
    palavras = [w for w in re.split(r"[^a-z0-9]+", s)
                if w and w not in ("do", "da", "dos", "das", "de")]
    return "_".join(w.capitalize() for w in palavras)[:16].rstrip("_")


def catalogo(area_min=AREA_MINIMA, base=BASE):
    """Todos os rios nomeados com area de drenagem acima do limiar.

    Devolve lista de dicionarios, ordenada por area decrescente, com o numero
    que a selecao usa.
    """
    g = gpd.read_file(base).to_crs(EPSG)
    g["NORIOCOMP"] = g["NORIOCOMP"].astype(str)
    # filtra pela chave NORMALIZADA, nao pelo valor cru: a base traz "NaN"
    # com maiuscula, que passava pelo filtro e virava um "rio" de 141 km2
    g["chave"] = g["NORIOCOMP"].map(normalizar)
    g = g[~g["chave"].isin(["nan", "none", "null", ""])].copy()

    linhas = []
    for chave, sub in g.groupby("chave"):
        area = float(sub["NUAREAMONT"].max())
        if area < area_min:
            continue
        grafias = sorted({str(x) for x in sub["NORIOCOMP"]})
        linhas.append({
            "chave": chave,
            "nome": max(grafias, key=len),      # a grafia mais completa
            "grafias": grafias,
            "ras": nome_ras(chave),
            "area": area,
            "km": float(sub.length.sum() / 1000.0),
            "trechos": int(len(sub)),
            "atual": chave in ATUAIS,
        })
    linhas.sort(key=lambda d: -d["area"])
    for i, d in enumerate(linhas, 1):
        d["n"] = i
    return linhas


def _termo(s):
    """Chave de BUSCA: sem acento, sem o generico da frente, e sem separador.

    A tabela do catalogo mostra duas grafias do mesmo rio -- "Rio Itajai-mirim"
    e o nome do HEC-RAS "Itajai_Mirim" -- e a busca so casava com a primeira.
    Quem lia a tabela e copiava a coluna "hec-ras" recebia "nao encontrei".
    Aqui o hifen, o underscore e o espaco viram a mesma coisa, e o generico
    ("rio", "ribeirao"...) sai dos DOIS lados, entao as duas grafias casam.
    """
    s = normalizar(s)
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"^(rio|ribeirao|corrego|arroio)\s+", "", s)
    return re.sub(r"\s+", " ", s).strip()


def selecionar(cat, spec=None):
    """Resolve uma selecao contra o catalogo.

    spec aceita:
        None ou "atuais"  -> os 12 rios do modelo de hoje
        "todos"           -> tudo que passou do limiar
        "1,3,5" ou "1-6"  -> por numero do catalogo
        "acu,mirim"       -> por trecho do nome, na grafia da BHO ou na do
                             HEC-RAS ("Itajai_Mirim" e "Itajai-mirim" valem)
    Misturar formas e permitido: "1-6,luis alves".
    """
    if spec is None or str(spec).strip().lower() in ("", "atuais"):
        return [d for d in cat if d["atual"]]
    s = str(spec).strip().lower()
    if s == "todos":
        return list(cat)

    escolhidos, faltando, largos = [], [], []
    for parte in [p.strip() for p in s.split(",") if p.strip()]:
        if re.fullmatch(r"\d+", parte):
            achado = [d for d in cat if d["n"] == int(parte)]
        elif re.fullmatch(r"\d+\s*-\s*\d+", parte):
            a, b = (int(x) for x in re.split(r"-", parte))
            achado = [d for d in cat if a <= d["n"] <= min(a, b) or
                      min(a, b) <= d["n"] <= max(a, b)]
        else:
            alvo = _termo(parte)
            # TERMO VAZIO E PEDIDO DE TUDO, e nao de nada. "selecao=Rio" --
            # que e o que sobra quando o .bat quebra "selecao=Rio Itajai_Mirim"
            # no espaco -- casava com os 36 rios do catalogo, e a rodada de um
            # rio so virava o vale inteiro sem uma linha de aviso. Custou uma
            # execucao inteira para ser percebido. Quem quer tudo escreve
            # "todos"; qualquer outro caminho para o catalogo completo e erro.
            if not alvo:
                largos.append((parte, len(cat)))
                continue
            achado = [d for d in cat
                      if alvo in _termo(d["chave"]) or alvo in _termo(d["ras"])]
            if len(achado) > max(3, len(cat) // 4):
                largos.append((parte, len(achado)))
                continue
        if not achado:
            faltando.append(parte)
        for d in achado:
            if d not in escolhidos:
                escolhidos.append(d)
    if largos:
        det = "; ".join(f"'{p}' casa com {q} rios" for p, q in largos)
        raise ValueError(
            f"selecao larga demais: {det}. Seja especifico (o nome do HEC-RAS "
            f"da tabela serve: 'Itajai_Mirim') ou escreva 'todos' se e o vale "
            f"inteiro que voce quer. Veja 'python -m vale.rios'."
            + (" ATENCAO: valor com espaco tem de vir entre aspas -- o .bat "
               "corta no espaco." if any(" " in p for p, _ in largos) else ""))
    if faltando:
        raise ValueError(
            f"nao encontrei no catalogo: {', '.join(faltando)}. "
            f"Rode 'python -m vale.rios' para ver os numeros disponiveis.")
    escolhidos.sort(key=lambda d: -d["area"])
    return escolhidos


def imprimir(cat, escolhidos=None):
    marca = {d["n"] for d in (escolhidos or [])}
    print(f"{'n':>3}  {'rio':<38}{'km2':>8}{'km':>7}{'trechos':>8}  "
          f"{'hec-ras':<16} {'atual':<6}{'sel' if escolhidos else ''}")
    print("-" * 100)
    for d in cat:
        print(f"{d['n']:>3}  {d['nome'][:38]:<38}{d['area']:>8.0f}"
              f"{d['km']:>7.0f}{d['trechos']:>8}  {d['ras']:<16} "
              f"{'sim' if d['atual'] else '':<6}"
              f"{'X' if d['n'] in marca else ''}")
        if len(d["grafias"]) > 1:
            print(f"     grafias fundidas: {' | '.join(d['grafias'])}")
    if escolhidos is not None:
        a = sum(d["area"] for d in escolhidos)
        print(f"\nselecionados: {len(escolhidos)} rios   "
              f"(maior area {max(d['area'] for d in escolhidos):.0f} km2, "
              f"soma das areas proprias nao e aditiva -- ha aninhamento)")
        print("   " + ", ".join(d["ras"] for d in escolhidos))


def main(argv=None):
    p = argparse.ArgumentParser(description="catalogo de rios do Vale")
    p.add_argument("--area", type=float, default=AREA_MINIMA)
    p.add_argument("--sel", default=None,
                   help="'todos', 'atuais', numeros (1,3,5 / 1-6) ou nomes")
    p.add_argument("--base", default=BASE)
    a = p.parse_args(argv)
    cat = catalogo(a.area, a.base)
    imprimir(cat, selecionar(cat, a.sel))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
