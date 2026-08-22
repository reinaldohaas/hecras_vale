# -*- coding: utf-8 -*-
"""Le a lista do "Validate Geometry" do HEC-RAS, por tipo e por local.

    python scripts/ler_erros_geometria.py modelo/mirim_t33/mirim_t33.g01.hdf

Chama o MESMO motor da interface -- `RasMapperLib.RASGeometry.ValidateGeometry`,
que e o que esta por tras do botao -- e devolve as mensagens uma a uma, em vez
do numero solto que a janela mostra.

POR QUE ISTO PRECISOU EXISTIR

  O contador dizia 755 numa geometria e 1880 noutra, e nenhuma quantidade
  mensuravel no arquivo chegava perto: tentei discordancia entre estaca e
  cutline, cruzamento com o eixo, sobreposicao de vizinhas, HTab, auto-
  interseccao das edge lines e superficies de interpolacao faltando. Nenhuma
  batia. Adivinhar a composicao de um numero e trabalho perdido -- este script
  simplesmente pergunta.

COMO SE CHEGA LA

  `RASGeometry` aceita o caminho do `.gNN.hdf` no construtor, `LoadAll` monta
  as camadas e `ValidateGeometry` preenche a camada `Errors` (um `RASErrors`).
  Cada entrada e uma colecao com uma ou mais mensagens, e cada mensagem traz
  nivel (Info/Warning/Fatal) e texto. `GetFeatureNames()` diz sobre qual
  feicao -- e por ai se chega a River Station.

  Requer o pythonnet e o RasMapperLib da versao instalada; `vale.terreno.
  registrar_hecras` cuida do caminho, que na 7.0.1 a biblioteca nao acha
  sozinha.
"""
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

LIMITE = 20000       # trava contra laco infinito no GetErrors


def ler(hdf):
    from vale.terreno import registrar_hecras, HECRAS_DIR
    registrar_hecras(log=lambda *a: None)
    import clr
    sys.path.append(HECRAS_DIR)
    os.add_dll_directory(HECRAS_DIR)
    clr.AddReference("RasMapperLib")
    import System
    asm = [a for a in System.AppDomain.CurrentDomain.GetAssemblies()
           if a.GetName().Name == "RasMapperLib"][0]
    TG = asm.GetType("RasMapperLib.RASGeometry")
    G = System.Activator.CreateInstance(
        TG, System.Array[System.Object]([os.path.abspath(hdf)]))
    G.LoadAll(None)
    G.ValidateGeometry(True)
    E = G.Errors

    saida, i = [], 0
    while i < LIMITE:
        try:
            c = E.GetErrors(i)
        except Exception:
            break
        if c is None:
            break
        try:
            onde = str(c.GetFeatureNames())
        except Exception:
            onde = ""
        try:
            camada = str(c.GetLayerNames())
        except Exception:
            camada = ""
        try:
            for er in c.Errors:
                saida.append({"nivel": str(er.Level), "camada": camada,
                              "onde": onde, "mensagem": str(er.Message)})
        except Exception:
            saida.append({"nivel": "?", "camada": camada, "onde": onde,
                          "mensagem": str(c.GetErrorDescription())})
        i += 1
    return saida, i


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    hdf = argv[0]
    if not os.path.exists(hdf):
        raise SystemExit(f"nao achei {hdf}")
    msgs, n_col = ler(hdf)
    print(f"geometria: {hdf}")
    print(f"colecoes de erro: {n_col}   mensagens: {len(msgs)}\n")

    def norm(m):
        return re.sub(r"[-+]?\d+[.,]?\d*", "N", m)

    c = collections.Counter(norm(m["mensagem"]) for m in msgs)
    lv = collections.Counter(m["nivel"] for m in msgs)
    print("POR NIVEL: " + "   ".join(f"{k} {v}" for k, v in lv.most_common()))
    print("\nPOR TIPO:")
    for k, v in c.most_common():
        print(f"   {v:6d}  ({100*v/len(msgs):4.1f}%)  {k[:96]}")

    base = os.path.splitext(os.path.splitext(hdf)[0])[0]
    p = base + "_erros.csv"
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, ["nivel", "camada", "onde", "mensagem"],
                           delimiter=";")
        w.writeheader()
        w.writerows(msgs)
    print(f"\ntabela completa -> {p}")

    # onde doem: as River Stations mais citadas
    rs = collections.Counter()
    for m in msgs:
        for x in re.findall(r"\b\d{2,6}\.\d\b", m["onde"]):
            rs[x] += 1
    if rs:
        print("\nRiver Stations mais citadas:")
        for k, v in rs.most_common(10):
            print(f"   RS {k:>12}  {v} mensagem(ns)")
    return p


if __name__ == "__main__":
    main(sys.argv[1:])
