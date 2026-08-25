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
  PROJETOS, entao a paralelizacao e um pool de processos sobre `compute_plan`,
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

# POR QUE PROCESSO, E NAO THREAD.
#
#   `ras_commander` guarda um objeto `ras` GLOBAL no modulo, e `init_ras_project`
#   mexe nele mesmo quando se passa `ras_object=`. Com varios rios em THREADS o
#   global e um so, compartilhado, e os projetos se atropelam: rodando os cinco
#   juntos, o init do `itajai_norte` foi ler `itajai_norte\itajai_acu.u01` -- o
#   arquivo de OUTRO rio -- e tres dos cinco terminaram SEM `.bco01`. Processo
#   separado tem seu proprio import de `ras_commander`, seu proprio global, e
#   nada vaza de um rio para o outro. O solver ja e um Ras.exe a parte; o pool
#   de processos so garante que o SETUP de cada rio tambem seja isolado.

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
    """Erro de volume e iteracoes, lidos do `.bco01`.

    O HEC-RAS NAO gera aqui o `.p01.computeMsgs.txt` que a versao antiga
    procurava -- por isso a tabela saia toda com `-`. O que interessa esta no
    `.bco01`, o log textual do solver:

      erro de volume    no bloco "Total Volume Accounting (for the entire
                        model)", na linha logo abaixo do cabecalho
                        "Error   Percent Error". Ha um bloco parecido so para
                        a "1D Flow Area"; queremos o do MODELO INTEIRO.
      iteracoes         cada passo do solver imprime uma ou mais linhas
                        `<i> <Reach> <valores...>`, onde `<i>` e a iteracao
                        (0,1,2,...) dentro daquele passo. O maximo desse `<i>`
                        em toda a corrida diz o pior caso de convergencia: 40
                        e o teto (nunca convergiu), 6 e folgado.
      falhou            se o bloco de volume nem existe, a corrida nao chegou
                        ao fim -- sinal mais confiavel que caçar uma frase.

    Serve a qualquer rio: o nome do reach entra como `\\w+`, nao fixo.
    """
    out = {"vol": None, "volpct": None, "iter": None, "falhou": None}
    bco = os.path.join(pasta, nome + ".bco01")
    if not os.path.exists(bco):
        out["falhou"] = True
        return out
    t = open(bco, encoding="latin-1", errors="replace").read()

    i = t.find("Total Volume Accounting")
    if i != -1:
        m = re.search(r"Error\s+Percent Error\s*\n\s*\*+\s+\*+\s*\n"
                      r"\s*([-\d.]+)\s+([-\d.]+)", t[i:i + 800])
        if m:
            out["vol"] = float(m.group(1))
            out["volpct"] = float(m.group(2))

    it = [int(x) for x in
          re.findall(r"^\s*(\d+)\s+\w+\s+[\d.]+\s+[\d.]+\s+[\d.]+", t, re.M)]
    if it:
        out["iter"] = max(it)

    out["falhou"] = out["vol"] is None or "Solution Solver Failed" in t
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
    with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
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
