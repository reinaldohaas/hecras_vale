# -*- coding: utf-8 -*-
"""Roda os rios EM PARALELO e devolve uma tabela com o que interessa julgar.

    python scripts/rodar_rios.py                    # todos os de modelo/
    python scripts/rodar_rios.py itajai_acu itajai_mirim
    python scripts/rodar_rios.py --workers 3 --cores 4

POR QUE PARALELO

  Cada rio e um PROJETO separado, com pasta propria. Rodar um de cada vez
  deixa a maquina ociosa: o solver 1D do HEC-RAS nao escala bem em nucleos --
  passar de 2 para 8 nucleos num rio so nao divide o tempo por quatro --, mas
  seis rios em seis processos usam a maquina inteira. `ras-commander` da o
  `compute_parallel` para varios PLANOS de um mesmo projeto; aqui sao varios
  PROJETOS, entao a paralelizacao e um pool de threads sobre `compute_plan`,
  cada um com seu `RasPrj`.

  `--workers` e quantos rios ao mesmo tempo, `--cores` quantos nucleos cada
  um pede. Manter workers * cores <= nucleos da maquina.

O QUE A TABELA MOSTRA, E POR QUE ESSES NUMEROS

  erro de volume    e o veredito. Um modelo que perde metade da agua nao
                    esta "quase certo": esta errado. Vem do `.pNN.hdf`, ou
                    do texto do computeMsgs quando o HDF nao saiu.
  iteracoes max     40 em todo passo quer dizer que o solver nunca convergiu
                    e seguiu no limite -- resultado que nao vale nada, mesmo
                    sem mensagem de erro.
  tempo             para saber o que custa cada tentativa.

  NAO roda quem tem erro de geometria: isso ja foi decidido no
  `construir_rio.py`, e rodar geometria ruim so gasta o relogio.
"""
import argparse
import concurrent.futures as cf
import os
import re
import sys
import time

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)
sys.path.insert(0, RAIZ)


def um(pasta, cores):
    """Roda um projeto. Devolve um dicionario com o que julgar."""
    from ras_commander import init_ras_project, RasCmdr
    from vale.terreno import HECRAS_DIR
    nome = os.path.basename(pasta)
    prj = os.path.join(pasta, nome + ".prj")
    t0 = time.time()
    r = {"rio": nome, "ok": False, "erro": None}
    try:
        p = init_ras_project(prj, os.path.join(HECRAS_DIR, "Ras.exe"))
        RasCmdr.compute_plan("01", ras_object=p, num_cores=cores,
                             force_rerun=True)
        r["ok"] = True
    except Exception as e:                       # noqa: BLE001
        r["erro"] = str(e)[:200]
    r["min"] = (time.time() - t0) / 60.0
    r.update(medir(pasta, nome))
    return r


def medir(pasta, nome):
    """Erro de volume e iteracoes, do computeMsgs -- que sempre existe."""
    out = {"vol": None, "volpct": None, "iter": None, "falhou": None}
    msg = os.path.join(pasta, nome + ".p01.computeMsgs.txt")
    if not os.path.exists(msg):
        return out
    t = open(msg, encoding="latin-1", errors="replace").read()
    m = re.search(r"Volume Accounting Error in 1000 m\^3:\s*([-\d.]+)", t)
    if m:
        out["vol"] = float(m.group(1))
    m = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", t)
    if m:
        out["volpct"] = float(m.group(1))
    # a ultima coluna das linhas do solver e a contagem de iteracoes. Contar
    # so as linhas 0 e 1 -- que foi o meu erro antes -- da maximo 2 quando o
    # verdadeiro e 40.
    it = [int(x) for x in re.findall(r"\s(\d+)\s*$", t, re.M)]
    if it:
        out["iter"] = max(it)
    out["falhou"] = "Solution Solver Failed" in t
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rios", nargs="*")
    ap.add_argument("--modelo", default="modelo")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--cores", type=int, default=2)
    a = ap.parse_args()

    if a.rios:
        pastas = [os.path.join(a.modelo, r) for r in a.rios]
    else:
        pastas = sorted(
            os.path.join(a.modelo, d) for d in os.listdir(a.modelo)
            if os.path.isdir(os.path.join(a.modelo, d))
            and os.path.exists(os.path.join(a.modelo, d, d + ".prj")))
    pastas = [p for p in pastas if os.path.exists(p)]
    if not pastas:
        raise SystemExit(f"nenhum projeto em {a.modelo}")
    print(f"{len(pastas)} projeto(s), {a.workers} ao mesmo tempo, "
          f"{a.cores} nucleo(s) cada")
    for p in pastas:
        print("   " + p)

    t0 = time.time()
    res = []
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        fut = {ex.submit(um, p, a.cores): p for p in pastas}
        for f in cf.as_completed(fut):
            r = f.result()
            res.append(r)
            print(f"   [{len(res)}/{len(pastas)}] {r['rio']} "
                  f"{'ok' if r['ok'] else 'FALHOU'}  {r['min']:.1f} min")
    res.sort(key=lambda d: d["rio"])

    print(f"\n{'='*76}\nRESULTADO   ({(time.time()-t0)/60:.1f} min de relogio)"
          f"\n{'='*76}")
    print(f"{'rio':<16}{'min':>7}{'iter':>7}{'volume 1000m3':>16}{'erro %':>10}"
          f"   veredito")
    for r in res:
        v = "-" if r["vol"] is None else f"{r['vol']:.0f}"
        vp = "-" if r["volpct"] is None else f"{r['volpct']:.2f}"
        it = "-" if r["iter"] is None else str(r["iter"])
        if r["falhou"]:
            j = "SOLVER FALHOU"
        elif r["volpct"] is not None and abs(r["volpct"]) > 5:
            j = "roda, mas perde agua demais"
        elif not r["ok"]:
            j = "nao rodou: " + (r["erro"] or "?")[:40]
        else:
            j = "ok"
        print(f"{r['rio']:<16}{r['min']:>7.1f}{it:>7}{v:>16}{vp:>10}   {j}")
    return res


if __name__ == "__main__":
    main()
