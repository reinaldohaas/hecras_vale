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
from ras_io import escrever   # noqa: E402

LIMITE = 20000       # trava contra laco infinito no GetErrors


def preparar_hdf(caminho):
    """Garante o `.gNN.hdf`, construindo-o SEM RODAR o solver se preciso.

    O `.hdf` da geometria costumava aparecer so depois de uma simulacao, e por
    isso a validacao vinha DEPOIS de rodar -- gastando dez minutos de solver
    para descobrir que a geometria tinha centenas de erros. Nao precisa:
    `RasMapperLib.Scripting.CompleteGeometryCommand` monta o HDF a partir do
    texto do `.gNN`, que e exatamente o que o RAS Mapper faz ao abrir.

    Aceita `.gNN` ou `.gNN.hdf`. Se o texto for mais novo que o HDF, refaz.
    """
    if caminho.lower().endswith(".hdf"):
        g = caminho[:-4]
        h = caminho
    else:
        g = caminho
        h = caminho + ".hdf"
    if not os.path.exists(g):
        if os.path.exists(h):
            return h
        raise SystemExit(f"nao achei {g}")
    if os.path.exists(h) and os.path.getmtime(h) >= os.path.getmtime(g):
        return h
    # `CompleteGeometryCommand` NAO serve para isto: ele COMPLETA um HDF que ja
    # exista e falha calado com "Geometry not found in WriteAttributePreCheck",
    # devolvendo `Execute ok` sem escrever arquivo nenhum. O caminho que
    # funciona e rodar o PREPROCESSADOR GEOMETRICO sozinho -- `Run HTab=-1`
    # com `Run UNet=0` --, que e barato: segundos, contra dez minutos de
    # solver.
    import re
    import shutil
    import tempfile
    from ras_commander import init_ras_project, RasCmdr
    from vale.terreno import HECRAS_DIR

    pasta = os.path.dirname(os.path.abspath(g)) or "."
    base = os.path.basename(g).split(".")[0]
    prj = os.path.join(pasta, base + ".prj")
    p01 = os.path.join(pasta, base + ".p01")
    if not (os.path.exists(prj) and os.path.exists(p01)):
        raise SystemExit(f"para montar o HDF preciso de {base}.prj e "
                         f"{base}.p01 ao lado da geometria")
    guarda = p01 + ".antes_do_htab"
    shutil.copy2(p01, guarda)
    try:
        t = open(p01, encoding="latin-1", errors="replace").read()
        for chave, val in (("Run HTab", "-1"), ("Run UNet", "0"),
                           ("Run PostProcess", "0"), ("Run RASMapper", "0")):
            if re.search(r"(?m)^%s=" % chave, t):
                t = re.sub(r"(?m)^%s=.*$" % chave, f"{chave}={val}", t)
            else:
                t = t.rstrip("\r\n") + f"\r\n{chave}={val}"
        escrever(p01, t)
        print(f"   montando {os.path.basename(h)} pelo preprocessador "
              "geometrico (sem solver)...")
        p = init_ras_project(prj, os.path.join(HECRAS_DIR, "Ras.exe"))
        RasCmdr.compute_plan("01", ras_object=p, force_rerun=True,
                             clear_geompre=True)
    finally:
        shutil.move(guarda, p01)
    if not os.path.exists(h):
        raise SystemExit(f"o preprocessador nao produziu {h}")
    return h


def ler(hdf):
    hdf = preparar_hdf(hdf)
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
