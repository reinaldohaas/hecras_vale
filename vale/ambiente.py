# -*- coding: utf-8 -*-
"""
Confere o ambiente e diz EXATAMENTE o que falta e como instalar.

Existe porque o modo anterior de falhar era pessimo: o .bat tinha o caminho do
interpretador cravado, e quando o ambiente mudou ele disse apenas "nao
encontrei" -- sem dizer qual Python serviria, o que faltava nele, nem o que
rodar para consertar. Um programa que depende de dez bibliotecas tem de saber
dizer qual delas nao esta la.

Uso:
    python -m vale.ambiente          confere e mostra o que falta
    python -m vale.ambiente --exe    so imprime um interpretador que serve
"""
import importlib.util
import os
import shutil
import subprocess
import sys

# (modulo, pacote para instalar, para que serve, essencial?)
DEPENDENCIAS = [
    ("numpy",          "numpy",          "toda a matematica", True),
    ("shapely",        "shapely",        "eixos, cortes e buffers", True),
    ("geopandas",      "geopandas",      "ler a base da ANA e gravar GeoJSON", True),
    ("pyproj",         "pyproj",         "converter coordenadas", True),
    ("rasterio",       "rasterio",       "ler o terreno", True),
    ("pandas",         "pandas",         "tabelas das secoes", True),
    ("scipy",          "scipy",          "vies local no QA/QC do SIG-SC", False),
    ("ras_commander",  "ras-commander",  "montar secoes, corrigir e rodar o "
                                         "HEC-RAS", False),
    ("h5py",           "h5py",           "ler o log do solver no .p01.hdf", False),
]

# candidatos a interpretador, na ordem em que sao tentados
CANDIDATOS = [
    r"C:\Users\haas\miniforge3\envs\hecras-qc\python.exe",
    r"C:\Users\haas\miniforge3\envs\vale\python.exe",
    r"C:\Users\haas\miniforge3\python.exe",
    sys.executable,
]


def falta(exe=None):
    """O que falta neste interpretador (ou no interpretador dado)."""
    if exe and os.path.abspath(exe) != os.path.abspath(sys.executable):
        codigo = ("import importlib.util as u;"
                  "print(','.join(m for m in %r if not u.find_spec(m)))"
                  % [d[0] for d in DEPENDENCIAS])
        try:
            r = subprocess.run([exe, "-c", codigo], capture_output=True,
                               text=True, timeout=90)
            ausentes = set(x for x in r.stdout.strip().split(",") if x)
        except Exception:                                    # noqa: BLE001
            return None
    else:
        ausentes = {m for m, _, _, _ in DEPENDENCIAS
                    if importlib.util.find_spec(m) is None}
    return [d for d in DEPENDENCIAS if d[0] in ausentes]


def versao(exe):
    try:
        r = subprocess.run([exe, "-c",
                            "import sys;print('.'.join(map(str,sys.version_info[:3])))"],
                           capture_output=True, text=True, timeout=60)
        return r.stdout.strip()
    except Exception:                                        # noqa: BLE001
        return "?"


def melhor_interpretador():
    """O primeiro candidato que existe e tem menos dependencias faltando."""
    achados = []
    for exe in CANDIDATOS:
        if not exe or not os.path.exists(exe):
            continue
        f = falta(exe)
        if f is None:
            continue
        essenciais = sum(1 for d in f if d[3])
        achados.append((essenciais, len(f), exe, f))
    if not achados:
        return None, None
    achados.sort()
    _, _, exe, f = achados[0]
    return exe, f


def comando_instalacao(ausentes, exe):
    """Como instalar, preferindo conda quando o interpretador e do miniforge."""
    pacotes = [d[1] for d in ausentes]
    if not pacotes:
        return None
    # MAMBA SEMPRE: resolve o ambiente em segundos onde o conda leva
    # minutos, e e o solver padrao do Miniforge.
    gerente = shutil.which("mamba") or shutil.which("conda")
    base = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
    e_conda = "miniforge" in exe.lower() or "conda" in exe.lower() or \
              os.path.exists(os.path.join(base, "conda-meta"))
    # ras-commander so existe no PyPI
    pip_only = {"ras-commander"}
    via_conda = [p for p in pacotes if p not in pip_only]
    via_pip = [p for p in pacotes if p in pip_only]
    linhas = []
    if e_conda and gerente and via_conda:
        nome = os.path.basename(gerente).split(".")[0]
        linhas.append(f"{nome} install -c conda-forge -y " + " ".join(via_conda))
        if via_pip:
            linhas.append(f'"{exe}" -m pip install ' + " ".join(via_pip))
    else:
        linhas.append(f'"{exe}" -m pip install ' + " ".join(pacotes))
    return linhas


def relatorio(log=print):
    exe, ausentes = melhor_interpretador()
    log("=" * 68)
    log("AMBIENTE PARA O MODELO DO VALE")
    log("=" * 68)
    if exe is None:
        log("  Nenhum interpretador Python utilizavel foi encontrado.")
        log("  Procurei em:")
        for c in CANDIDATOS:
            log(f"     {c}")
        return 2, None

    log(f"  interpretador : {exe}")
    log(f"  versao        : Python {versao(exe)}")
    presentes = [d for d in DEPENDENCIAS if d not in ausentes]
    log(f"  bibliotecas   : {len(presentes)} de {len(DEPENDENCIAS)} presentes")
    log("")
    if ausentes:
        log("  FALTANDO:")
        for m, pacote, para, essencial in ausentes:
            marca = "ESSENCIAL" if essencial else "opcional "
            log(f"     [{marca}] {m:<16} {para}")
        log("")
        cmds = comando_instalacao(ausentes, exe)
        log("  Para instalar, rode no Prompt (Miniforge):")
        for c in cmds:
            log(f"     {c}")
        log("")
        essenciais = [d for d in ausentes if d[3]]
        if essenciais:
            log("  Sem as ESSENCIAIS o programa nao roda de jeito nenhum.")
            return 1, exe
        log("  As que faltam sao opcionais:")
        log("     sem scipy         -- o QA/QC usa vies global em vez de local")
        log("     sem ras-commander -- nao monta secao pelo HEC-RAS, nao roda")
        log("                          o solver; a formatacao propria assume")
        log("     sem h5py          -- nao le o log do solver de dentro do HDF")
        return 0, exe
    log("  Tudo presente. Pode rodar:")
    log(f'     "{exe}" -m vale tudo --auto fonte=copernicus')
    return 0, exe


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--exe" in argv:
        exe, _ = melhor_interpretador()
        print(exe or "")
        return 0 if exe else 2
    codigo, _ = relatorio()
    return codigo


if __name__ == "__main__":
    raise SystemExit(main())
