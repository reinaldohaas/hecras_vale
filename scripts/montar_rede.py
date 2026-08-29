# -*- coding: utf-8 -*-
"""Costura a rede e APERTA sozinha ate o validador do HEC-RAS zerar.

    python scripts/montar_rede.py

Faz, em laco:

    1. construir_rede.py  --excluir <descartes>   costura os rios
    2. projeto_rede.py                            monta o projeto
    3. ler_erros_geometria.py                     valida SEM solver
    4. as secoes Fatais que cruzam o reach vizinho na confluencia entram na
       lista de descarte, e volta ao 1.

E o mesmo remedio do `construir_rio.py`, que aperta a taxa de largura sozinho:
aqui o que se remove sao as secoes que, na juncao, se sobrepoem a do outro rio.
Nao se inventa cota nem se mexe na secao -- so se tira a que o RAS recusa, e a
juncao conduz a agua dali. Cada descarte fica no CSV, para o relatorio dizer
quantas secoes de cada rio sairam e onde.
"""
import csv
import os
import re
import subprocess
import sys

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
PY = os.path.join(RAIZ, ".venv", "Scripts", "python.exe")
if not os.path.exists(PY):
    PY = sys.executable

GEOM = "modelo/itajai_rede/itajai_rede.g01"
DESCARTES = "doc/rede_descartes.csv"
MAX_ITER = 8


def roda(args, silenciar=True):
    p = subprocess.run([PY] + args, cwd=RAIZ, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    saida = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0:
        print(saida[-2000:])
        raise SystemExit(f"falhou: {' '.join(args)}")
    return saida


def fatais_xs(csv_erros):
    """(river, rs) das secoes com erro Fatal de cruzamento na costura."""
    out = []
    if not os.path.exists(csv_erros):
        return out
    for r in csv.DictReader(open(csv_erros, encoding="utf-8"), delimiter=";"):
        if r["nivel"].strip() != "Fatal":
            continue
        if "Cross Section" not in r["camada"] and "Cross Section" not in \
                r["onde"]:
            continue
        m = re.match(r"\s*([A-Za-z_]+)\s*,\s*\w+\s*\(([-\d.]+)\)", r["onde"])
        if m:
            out.append((m.group(1), round(float(m.group(2)), 2)))
    return out


def main():
    os.makedirs(os.path.dirname(DESCARTES), exist_ok=True)
    descartes = []       # [(river, rs)]
    if os.path.exists(DESCARTES):
        os.remove(DESCARTES)

    for it in range(1, MAX_ITER + 1):
        print(f"\n{'='*68}\nITERACAO {it}   ({len(descartes)} secao(oes) "
              f"descartada(s) ate aqui)\n{'='*68}")
        # 1. costura
        args = ["scripts/construir_rede.py", "--saida", GEOM]
        if descartes:
            args += ["--excluir", DESCARTES]
        roda(args)
        # 2. projeto
        roda(["scripts/projeto_rede.py", GEOM])
        # 3. valida (forca refazer o hdf)
        h = os.path.join(RAIZ, GEOM + ".hdf")
        if os.path.exists(h):
            os.remove(h)
        saida = roda(["scripts/ler_erros_geometria.py", GEOM])
        m = re.search(r"mensagens:\s*(\d+)", saida)
        nf = re.search(r"Fatal (\d+)", saida)
        n_msg = int(m.group(1)) if m else -1
        n_fat = int(nf.group(1)) if nf else 0
        print(f"   validador: {n_msg} mensagens, {n_fat} Fatal")
        # 4. novas fatais
        csv_erros = os.path.join(RAIZ, "modelo", "itajai_rede",
                                 "itajai_rede_erros.csv")
        novas = [x for x in fatais_xs(csv_erros) if x not in descartes]
        if n_fat == 0:
            print("\n   REDE LIMPA: 0 Fatal no validador do HEC-RAS")
            break
        if not novas:
            print(f"\n   parou com {n_fat} Fatal e nada novo a descartar "
                  "-- estes nao sao de secao cruzando reach; ver o CSV")
            break
        descartes += novas
        with open(DESCARTES, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(["river", "reach", "rs"])
            for riv, rs in descartes:
                w.writerow([riv, "", rs])
        print(f"   +{len(novas)} secao(oes) para descartar: "
              + ", ".join(f"{r}@{s:.0f}" for r, s in novas))

    # resumo dos descartes por rio
    if descartes:
        from collections import Counter
        c = Counter(r for r, _ in descartes)
        print("\ndescartes por rio (secoes na confluencia que o RAS recusou):")
        for riv, n in c.most_common():
            print(f"   {riv:14} {n}")
        print(f"lista completa -> {DESCARTES}")
    return GEOM


if __name__ == "__main__":
    main()
