# -*- coding: utf-8 -*-
"""
Roda o plano e le o log do solver.

Existe por dois motivos, os dois pagos com sessoes inteiras de depuracao:

PROJETO ERRADO. RasCmdr.compute_plan('01') resolve o plano dentro da PASTA, e
uma pasta com varios projetos faz ele computar outro e devolver SUCCESS. Aqui o
projeto e copiado sozinho para uma pasta isolada antes de rodar, entao '01' so
pode ser o dele.

RESULTADO VELHO. Copiando a pasta inteira vao junto o .p01.hdf, o .u01.hdf e o
.g01.hdf de execucoes anteriores. O solver le deles, e o resumo sai identico ao
da rodada passada -- a conclusao seria que a correcao nao fez efeito, quando o
que houve foi ler o log de antes. Os gerados ficam de fora da copia.

O log nao esta em computeMsgs.txt: o Compute via COM nao escreve esse arquivo,
so a GUI escreve. Ele esta dentro do .p01.hdf, e e de la que se le.
"""
import os
import re
import shutil
import pathlib

GERADOS = (".p01.hdf", ".u01.hdf", ".g01.hdf", ".O01", ".O02", ".r01",
           ".x01", ".bco01", ".ic.o01", ".dss", ".b01")


def isolar(prj, destino=None):
    # ABSOLUTO, e nao como veio. O .bat passa "modelo\projeto.prj", relativo, e
    # o symlink do terreno herdava esse caminho: o link em
    # %TEMP%\vale_runs\<proj>\Terrain apontava para "modelo\Terrain", que
    # resolve DENTRO da pasta temporaria e nao existe. O link existia, tinha 14
    # bytes e nao levava a lugar nenhum -- e o RAS Mapper abria sem terreno,
    # sem mancha de inundacao e sem dizer por que.
    prj = pathlib.Path(prj).resolve()
    raiz, nome = prj.parent, prj.stem
    destino = pathlib.Path(destino or (pathlib.Path(
        os.environ.get("TEMP", ".")) / "vale_runs" / nome))
    if destino.exists():
        shutil.rmtree(destino, ignore_errors=True)
    destino.mkdir(parents=True, exist_ok=True)
    for f in raiz.glob(f"{nome}.*"):
        if any(f.name.lower().endswith(e.lower()) for e in GERADOS):
            continue
        shutil.copy2(f, destino / f.name)
    terreno = raiz / "Terrain"
    if terreno.is_dir():
        # o terreno pode ter varios GB: liga-se por junction em vez de copiar
        alvo = destino / "Terrain"
        # CONFERIR que o link RESOLVE, e nao so que a chamada nao levantou. No
        # Windows sem modo desenvolvedor o symlink pode ser criado e nao levar
        # a lugar nenhum -- e ai o RAS Mapper abre sem terreno, sem mancha de
        # inundacao e sem dizer por que. Nao levantar excecao nao e prova de
        # que funcionou.
        ok = False
        try:
            os.symlink(terreno, alvo, target_is_directory=True)
            ok = any(alvo.glob("*.hdf"))
        except (OSError, NotImplementedError):
            ok = False
        if not ok:
            if alvo.exists() or alvo.is_symlink():
                try:
                    alvo.unlink()
                except OSError:
                    shutil.rmtree(alvo, ignore_errors=True)
            # so o terreno DESTE projeto: a pasta guarda o de todas as rodadas
            alvo.mkdir(parents=True, exist_ok=True)
            for f in terreno.glob(f"{prj.stem}*"):
                shutil.copy2(f, alvo / f.name)
    return destino / prj.name


def rodar(prj, ras_exe, log=print):
    from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan
    novo = isolar(prj)
    log(f"      isolado em {novo.parent}")
    p = init_ras_project(str(novo), ras_exe)
    r = RasCmdr.compute_plan("01", ras_object=p, force_rerun=True,
                             clear_geompre=True)
    hdf = novo.with_suffix(".p01.hdf")
    try:
        msgs = str(HdfResultsPlan.get_compute_messages(hdf))
    except Exception as e:                                   # noqa: BLE001
        msgs = f"(sem mensagens: {e})"
    return r, msgs, str(novo)


def erros_de_dado(pasta):
    """Le o *.data_errors.txt -- onde o HEC-RAS escreve por que NAO rodou.

    Este arquivo tinha a resposta e o programa o ignorava. No Benedito isolado
    dizia, em duas linhas: "Boundary at River: Benedito Reach: R1 RS: 75.00 /
    Stage(s) in time series data are below the cross section minimum" -- a mare
    de 0,3 m contra uma foz a 50 m de altitude. Em vez disso a rodada terminou
    com "NENHUM PROBLEMA DETECTADO".
    """
    import glob
    saida = []
    for c in glob.glob(os.path.join(pasta or ".", "*.data_errors.txt")):
        try:
            t = open(c, encoding="latin-1", errors="replace").read().strip()
        except OSError:
            continue
        if t:
            saida.append(t)
    return "\n".join(saida)


def rodou_de_fato(log_txt):
    """O solver produziu alguma coisa?

    Sem isto, ausencia de log era lida como sucesso: resumir() dizia
    "instavel em: nao (completou)" quando nao havia mensagem NENHUMA, e o
    passo 9 terminava em 2 segundos declarando exito. Um solver que nao roda e
    a pior das falhas para reportar como sucesso, porque todo o resto do
    relatorio -- auditoria, figuras, resumo -- continua saindo bonito.
    """
    if not (log_txt or "").strip():
        return False
    return bool(re.search(r"Unsteady Flow Simulation|Finished|"
                          r"Volume Accounting|went unstable", log_txt))


def diagnostico(log_txt, pasta=None):
    """O essencial em forma de dicionario, para a checagem automatica ver.

    Sem isto o passo 9 podia abortar a simulacao e a rodada terminar com
    "NENHUM PROBLEMA DETECTADO": o unico caso que falhava era o solver nao
    rodar. Simulacao que roda e ABORTA no meio e falha igual -- so que mais
    silenciosa, porque produz log, produz HDF e produz figura.
    """
    if not rodou_de_fato(log_txt):
        return {"rodou": False, "completou": False, "volume": None,
                "abortou_em": None, "dados": erros_de_dado(pasta)}
    fim = re.search(r"went unstable at:\s*(\S+\s+\S+)", log_txt)
    onde = re.search(r"Minimum error exceeds allowable tolerance at\s+(\S+)\s+"
                     r"(\S+)\s*\n\s*\n(\S+)\s+(\S+)\s+([\d.]+)", log_txt)
    vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", log_txt)
    return {"rodou": True,
            "completou": fim is None and onde is None,
            "volume": float(vol.group(1)) if vol else None,
            "abortou_em": (f"{onde.group(3)} {onde.group(4)} RS {onde.group(5)}"
                           if onde else None),
            "instavel_em": fim.group(1) if fim else None,
            "dados": ""}


def resumir(log_txt, pasta=None):
    """O essencial: ate onde chegou, onde doeu, quanto de erro de volume."""
    if not rodou_de_fato(log_txt):
        L = ["O SOLVER NAO RODOU -- nao ha mensagem de computacao nenhuma.",
             "Isto NAO e sucesso: nao houve simulacao para dar certo."]
        dados = erros_de_dado(pasta)
        if dados:
            L += ["", "o HEC-RAS recusou os dados:"]
            L += ["   " + l for l in dados.splitlines()]
        else:
            L += ["", "sem *.data_errors.txt; procure o .p01.hdf e o log do RAS"]
        return "\n".join(L)
    fim = re.search(r"went unstable at:\s*(\S+\s+\S+)", log_txt)
    onde = re.search(r"Minimum error exceeds allowable tolerance at\s+(\S+)\s+"
                     r"(\S+)\s*\n\s*\n(\S+)\s+(\S+)\s+([\d.]+)", log_txt)
    vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", log_txt)
    piores, ultimo = [], None
    for m in re.finditer(r"^(\d{2}\w{3}\d{4}\s+\d{2}:\d{2}:\d{2})\s+(\S+)\s+"
                         r"(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$",
                         log_txt, re.M):
        ultimo = m.group(1)
        # o grupo 5 e a COTA DA LAMINA, e ficava de fora. Sem ela a linha
        # diz so "erro 11.18 no RS 29574" e nao da para saber que a agua
        # estava 2,5 cm acima do leito -- que era a causa. Custou duas rodadas.
        piores.append((float(m.group(6)), m.group(2), m.group(3),
                       m.group(4), m.group(1), int(m.group(7)),
                       float(m.group(5))))
    L = [f"instavel em: {fim.group(1) if fim else 'nao (completou)'}",
         f"ultimo passo com log: {ultimo}",
         f"erro de volume: {vol.group(1) + '%' if vol else 'n/d'}"]
    if onde:
        L.append(f"ABORTOU EM: {onde.group(3)} {onde.group(4)} "
                 f"RS {onde.group(5)}  ({onde.group(1)} {onde.group(2)})")
    L += ["", "maiores erros de nivel por passo:"]
    for e, rio, rea, rs, t, it, ws in sorted(piores, reverse=True)[:12]:
        L.append(f"   {e:8.2f} m  {rio:<16}{rea:<4}RS {rs:>10}  "
                 f"WSEL {ws:9.2f}  {t}  it={it}")
    ex = re.search(r"Extrapolated above Cross Section Table at:\s*\*+\s*\n(.*?)\n\n",
                   log_txt, re.S)
    if ex:
        L += ["", "extrapolou acima da tabela hidraulica:"]
        L += ["   " + l.strip() for l in ex.group(1).strip().splitlines()[:12]]
    return "\n".join(L)


def pontos_criticos(log_txt, n=4):
    """Onde o solver sofreu, com a lamina que ele calculou.

    Devolve (rio, rs, wsel, motivo) para virar figura. A lamina e o que
    distingue "geometria ruim" de "modelo rodando seco" -- e sem ela as duas
    coisas produzem exatamente o mesmo log.
    """
    saida, vistos = [], set()
    for m in re.finditer(r"^(\d{2}\w{3}\d{4}\s+\d{2}:\d{2}:\d{2})\s+(\S+)\s+"
                         r"(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$",
                         log_txt, re.M):
        saida.append((float(m.group(6)), m.group(2), float(m.group(4)),
                      float(m.group(5)),
                      f"solver: erro de nivel {float(m.group(6)):.2f} m em "
                      f"{int(m.group(7))} iteracoes ({m.group(1)})"))
    saida.sort(reverse=True)
    fim = []
    for _e, rio, rs, ws, motivo in saida:
        if (rio, rs) in vistos:
            continue
        vistos.add((rio, rs))
        fim.append((rio, rs, ws, motivo))
        if len(fim) >= n:
            break
    onde = re.search(r"Minimum error exceeds allowable tolerance at\s+(\S+)\s+"
                     r"(\S+)\s*\n\s*\n(\S+)\s+(\S+)\s+([\d.]+)", log_txt)
    if onde:
        par = (onde.group(3), float(onde.group(5)))
        if par not in vistos:
            fim.append((par[0], par[1], None,
                        "solver: erro excede a tolerancia -- ABORTOU AQUI"))
    return fim


def resultados(hdf):
    """Serie de nivel e vazao por secao, para a visualizacao."""
    from ras_commander import HdfResultsPlan
    return {"wse": HdfResultsPlan.get_wse(hdf)}
