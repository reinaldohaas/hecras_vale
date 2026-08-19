# -*- coding: utf-8 -*-
"""
Checagem do modelo ANTES de simular -- o que a GUI do HEC-RAS faz no menu.

    python -m vale.checagem                        confere o projeto padrao
    python -m vale.checagem projeto=vale_v2        outro projeto
    python -m vale.checagem projeto=acu_v3 --rios  detalha rio a rio

POR QUE ISTO EXISTE. A checagem ja estava no passo 8 e rodava CEGA. Ela
reportava `FlowType.GEOMETRY_ONLY` e "No cross section data found in geometry
HDF", e mesmo assim seguia -- porque o `.g01.hdf` nao existia ainda. Ele so
nasce dentro da pasta isolada, no passo 9, DEPOIS da checagem. Conferir uma
geometria compilada que ninguem compilou nao acusa nada, e "nao acusou nada"
foi lido como "esta bom" a noite inteira.

O QUE COMPILA O .g01.hdf. Nao e a biblioteca: `GeomMesh.compile_geometry`
recusa a tarefa com todas as letras -- "ras-commander cannot generate .g##.hdf
from .g## text with RasMapperLib". E `RasProcess.compute_geometry` COMPLETA um
HDF existente, nao cria. Quem cria e o pre-processador do proprio HEC-RAS, e
ele roda como parte de um plano. Entao aqui o plano e rodado com

    Run HTab=-1     pre-processa a geometria (e o que queremos)
    Run UNet=0      NAO simula
    Run PostProcess=0
    Run RASMapper=0

que e o equivalente sem tela de abrir a geometria na GUI e mandar recalcular as
tabelas hidraulicas. Sai barato perto de uma simulacao e produz exatamente o
arquivo que a checagem precisa ler.

E POR QUE ANTES DE INTEGRAR. Cada rio deve passar sozinho antes de os doze
rodarem juntos. Corrigir na integracao e o que vinha sendo feito, e o preco
apareceu: um conserto no condicionamento do perfil, feito para o Itajai do
Norte, derrubou o Itajai-Acu -- que estava rodando -- em 50 segundos com 193%
de erro de volume. Componente compartilhado se mede nos doze, nao no que
quebrou.
"""
import argparse
import os
import re
import shutil
import sys

from .config import Opcoes


def _plano_so_htab(pasta, projeto, log=print):
    """Deixa o .p01 em modo pre-processar-e-parar. Devolve o texto original."""
    caminho = os.path.join(pasta, f"{projeto}.p01")
    original = open(caminho, encoding="latin-1", errors="replace").read()
    t = original
    for chave, valor in (("Run HTab", "-1"), ("Run UNet", "0"),
                         ("Run PostProcess", "0"), ("Run RASMapper", "0")):
        if re.search(rf"^{chave}=", t, re.M):
            t = re.sub(rf"^{chave}=.*$", f"{chave}={valor}", t, flags=re.M)
        else:
            t += f"{chave}={valor}\n"
    open(caminho, "w", encoding="latin-1", errors="replace").write(t)
    log("      plano ajustado: Run HTab=-1, Run UNet=0 (pre-processa, nao simula)")
    return original


def compilar(op, log=print):
    """Roda o pre-processador do HEC-RAS e devolve o caminho do .g01.hdf.

    Em pasta ISOLADA, pelos mesmos motivos do vale/executar.py: uma pasta com
    varios projetos faz o '01' resolver para outro, e HDF de execucao anterior
    e lido como se fosse desta.
    """
    from ras_commander import RasCmdr, init_ras_project

    from .executar import isolar

    prj = op.caminho(f"{op.projeto}.prj")
    if not os.path.exists(prj):
        raise SystemExit(f"nao ha {prj}; rode o passo 7 antes")
    novo = isolar(prj, os.path.join(
        os.environ.get("TEMP", "."), "vale_checagem", op.projeto))
    log(f"      isolado em {novo.parent}")
    _plano_so_htab(str(novo.parent), op.projeto, log)

    ras = init_ras_project(str(novo), op.ras_exe)
    r = RasCmdr.compute_plan("01", ras_object=ras, force_rerun=True,
                             clear_geompre=True)
    hdf = novo.parent / f"{op.projeto}.g01.hdf"
    log(f"      {r}")
    if not hdf.exists():
        raise SystemExit(
            "o pre-processador nao produziu o .g01.hdf. Sem ele a checagem "
            "roda cega, e checagem cega nao serve para liberar nada.\n"
            f"procure {hdf}")
    log(f"      geometria compilada: {hdf.name} "
        f"({hdf.stat().st_size / 1e6:.0f} MB)")
    return str(novo), str(hdf)


def conferir(prj, log=print):
    """RasCheck.run_all sobre a geometria JA compilada. Devolve o resultado.

    O PRIMEIRO ARGUMENTO E O PLANO, e nao o projeto. A assinatura diz "plan
    number (e.g. '01') or path to plan HDF file", e passar o .prj faz o HdfXsec
    tentar abrir o .prj COMO HDF -- "Unable to synchronously open file (file
    signature not found)". Dali ele conclui "No cross section data found in
    geometry HDF" e cai para geometry_only. O sintoma parece falta de geometria
    compilada e e argumento errado; foi assim que a checagem passou a noite
    inteira cega mesmo com o .g01.hdf de 96 MB ao lado.
    """
    from ras_commander import RasCheck, init_ras_project

    ras = init_ras_project(prj, Opcoes().ras_exe)
    r = RasCheck.run_all("01", ras_object=ras)
    # geometry_only e o ESPERADO aqui: nao simulamos, entao nao ha resultado
    # para as checagens de balanco, pico e estabilidade. O que precisa ser
    # conferido e outra coisa -- se ela LEU as secoes. Quando nao le, o
    # proprio RasCheck emite um erro de SYSTEM dizendo isso, e e esse o sinal.
    cega = False
    try:
        df = r.to_dataframe()
        if df is not None and len(df) and "message" in df.columns:
            cega = df["message"].astype(str).str.contains(
                "No cross section data", case=False).any()
    except Exception:                                        # noqa: BLE001
        pass
    if cega:
        log("      ATENCAO: a checagem NAO leu as secoes -- o resultado "
            "abaixo nao libera nada")
    return r


def relatorio(r, log=print, por_rio=False, limite=20):
    """Imprime o resultado agrupado, e devolve o numero de ERROS."""
    n_err = int(getattr(r, "get_error_count", lambda: 0)())
    n_avi = int(getattr(r, "get_warning_count", lambda: 0)())
    log("")
    log(f"   tipo de escoamento : {getattr(r, 'flow_type', '?')}")
    log(f"   erros              : {n_err}")
    log(f"   avisos             : {n_avi}")
    try:
        df = r.to_dataframe()
    except Exception:                                        # noqa: BLE001
        return n_err
    if df is None or not len(df):
        return n_err
    cols = [c for c in ("severity", "check_type", "river", "reach", "station",
                        "message") if c in df.columns]
    sev = "severity" if "severity" in df.columns else None
    if sev is not None:
        log("")
        for s, n in df[sev].value_counts().items():
            log(f"   {str(s):<10} {n}")
    if por_rio and "river" in df.columns:
        log("")
        log(f"   {'rio':<16}{'erros':>7}{'avisos':>8}")
        for rio, g in df.groupby("river"):
            e = int((g[sev].astype(str).str.upper() == "ERROR").sum()) if sev else 0
            a = len(g) - e
            log(f"   {str(rio):<16}{e:>7}{a:>8}")
    graves = df
    if sev is not None:
        graves = df[df[sev].astype(str).str.upper() == "ERROR"]
        if not len(graves):
            graves = df
    log("")
    for _, x in graves.head(limite).iterrows():
        log("   " + " | ".join(f"{c}={x[c]}" for c in cols))
    if len(graves) > limite:
        log(f"   ... e mais {len(graves) - limite}")
    return n_err


def main(argv=None):
    p = argparse.ArgumentParser(
        description="confere o modelo antes de simular, como a GUI do HEC-RAS")
    p.add_argument("pares", nargs="*", help="chave=valor (ex: projeto=vale_v2)")
    p.add_argument("--rios", action="store_true", help="detalha rio a rio")
    p.add_argument("--manter", action="store_true",
                   help="nao apaga a pasta de checagem no fim")
    a = p.parse_args(argv)

    op = Opcoes()
    op.aplicar([x for x in a.pares if "=" in x])
    print("=" * 74)
    print(f"CHECAGEM DO MODELO -- projeto {op.projeto}")
    print("=" * 74)
    prj, hdf = compilar(op, print)
    r = conferir(prj, print)
    n = relatorio(r, print, por_rio=a.rios)
    if not a.manter:
        shutil.rmtree(os.path.dirname(prj), ignore_errors=True)
    print("")
    if n:
        print(f"REPROVADO: {n} erro(s). Nao integre antes de zerar isto.")
        return 1
    print("APROVADO: nenhum erro. O modelo pode ser simulado.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
