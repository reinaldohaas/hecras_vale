# -*- coding: utf-8 -*-
"""
Roda um plano e devolve o log do solver.

Existe por dois motivos, os dois pagos com sessoes inteiras de depuracao:

PROJETO ERRADO. RasCmdr.compute_plan('01') resolve o plano dentro da PASTA, e a
pasta do repositorio tem tres projetos. Ele devolveu SUCCESS tendo computado
outro projeto, duas vezes. Aqui o projeto e copiado sozinho para uma pasta
isolada antes de rodar, entao '01' so pode ser o dele.

LOG CEGO. O Compute_CurrentPlan via COM nao escreve computeMsgs.txt -- so a GUI
escreve. Por muitas rodadas a unica forma de saber onde falhou foi o usuario
abrir o HEC-RAS e colar o log. As mensagens estao no .p01.hdf, e
HdfResultsPlan.get_compute_messages le de la.

Uso:  python rodar.py taha_ai_1983_sb [linhas_do_log]
"""
import os
import re
import shutil
import sys
import pathlib

RAS = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"
TMP = pathlib.Path(os.environ.get("TEMP", ".")) / "taha_ai_runs"


def isolar(projeto, raiz="."):
    """Copia so os arquivos DESTE projeto para uma pasta propria."""
    raiz = pathlib.Path(raiz)
    dest = TMP / projeto
    if dest.exists():
        shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir(parents=True, exist_ok=True)
    # NAO copiar resultados: o .p01.hdf e o .O01 da pasta do repositorio sao de
    # execucoes antigas, e get_compute_messages leria o log ERRADO. Aconteceu:
    # depois de corrigir a condicao inicial o resumo saiu identico ao da rodada
    # anterior, digito por digito, e a conclusao seria que a correcao nao fez
    # efeito -- quando o que houve foi ler o log de antes.
    # Inclui .u01.hdf e .g01.hdf: sao GERADOS pelo HEC-RAS a partir do texto, e
    # o solver le deles. Com um .u01.hdf antigo por perto a condicao inicial
    # nova do .u01 e ignorada -- o resumo saiu identico digito por digito
    # depois de trocar todas as vazoes iniciais.
    resultados = (".p01.hdf", ".u01.hdf", ".g01.hdf", ".O01", ".O02", ".r01",
                  ".x01", ".bco01", ".ic.o01", ".dss")
    for f in raiz.glob(f"{projeto}.*"):
        if any(f.name.lower().endswith(e.lower()) for e in resultados):
            continue
        shutil.copy2(f, dest / f.name)
    terreno = raiz / "Terrain"
    if terreno.is_dir() and not (dest / "Terrain").exists():
        shutil.copytree(terreno, dest / "Terrain")
    return dest / f"{projeto}.prj"


def rodar(projeto, raiz="."):
    from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan
    prj = isolar(projeto, raiz)
    p = init_ras_project(str(prj), RAS)
    r = RasCmdr.compute_plan("01", ras_object=p, force_rerun=True,
                             clear_geompre=True)
    hdf = prj.with_suffix(".p01.hdf")
    log = ""
    try:
        log = str(HdfResultsPlan.get_compute_messages(hdf))
    except Exception as e:                       # noqa: BLE001
        log = f"(sem mensagens: {e})"
    return r, log


def resumir(log):
    """O essencial: ate onde chegou, onde doeu, quanto de erro de volume."""
    fim = re.search(r"went unstable at:\s*(\S+\s+\S+)", log)
    onde = re.search(r"Minimum error exceeds allowable tolerance at\s+(\S+)\s+"
                     r"(\S+)\s*\n\s*\n(\S+)\s+(\S+)\s+([\d.]+)", log)
    vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", log)
    piores, ultimo = [], None
    for m in re.finditer(r"^(\d{2}\w{3}\d{4}\s+\d{2}:\d{2}:\d{2})\s+(\S+)\s+"
                         r"(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$",
                         log, re.M):
        ultimo = m.group(1)
        piores.append((float(m.group(6)), m.group(2), m.group(3),
                       m.group(4), m.group(1), int(m.group(7))))
    linhas = []
    linhas.append(f"instavel em: {fim.group(1) if fim else 'nao (completou)'}")
    linhas.append(f"ultimo passo com log: {ultimo}")
    linhas.append(f"erro de volume: {vol.group(1) + '%' if vol else 'n/d'}")
    if onde:
        linhas.append(f"ABORTOU EM: {onde.group(3)} {onde.group(4)} "
                      f"RS {onde.group(5)}  ({onde.group(1)} {onde.group(2)})")
    linhas.append("")
    linhas.append("maiores erros de nivel por passo:")
    for e, rio, rea, rs, t, it in sorted(piores, reverse=True)[:10]:
        linhas.append(f"   {e:8.2f} m  {rio:<16}{rea:<4}RS {rs:>10}  {t}  it={it}")
    ex = re.search(r"Extrapolated above Cross Section Table at:\s*\*+\s*\n(.*?)\n\n",
                   log, re.S)
    if ex:
        linhas.append("")
        linhas.append("extrapolou acima da tabela:")
        linhas += ["   " + l.strip() for l in ex.group(1).strip().splitlines()[:10]]
    return "\n".join(linhas)


if __name__ == "__main__":
    proj = sys.argv[1] if len(sys.argv) > 1 else "taha_ai_1983_sb"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    res, log = rodar(proj)
    print("=" * 70)
    print(f"{proj}: {res}")
    print("=" * 70)
    print(resumir(log))
    if n:
        print("\n" + "=" * 70 + f"\nultimas {n} linhas do log\n" + "=" * 70)
        print("\n".join(log.splitlines()[-n:]))
