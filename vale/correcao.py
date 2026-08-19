# -*- coding: utf-8 -*-
"""
Correcao de erros: as ferramentas do proprio HEC-RAS, mais a auditoria local.

Preferencia deliberada pelas ferramentas do RAS. Elas conhecem o formato e as
regras internas melhor que qualquer coisa que eu escreva, e o que elas
consertam nao volta a quebrar na proxima versao do programa. O que nao existe
la -- auditar a geometria contra o TERRENO -- fica por conta do hecras_qc.

  RasFixit.fix_bank_stations         margem que nao casa com um sta da tabela
  RasFixit.fix_ineffective_flow      escoamento inefetivo nas secoes largas
  RasFixit.fix_htab_starting_elev    tabela hidraulica comecando errado
  GeomHtabUtils.calculate_optimal_xs_htab
                                     as HTab por secao. Sem isso todas ficam no
                                     padrao, e o log enche de "Extrapolated
                                     above Cross Section Table" -- que e a agua
                                     passando do TOPO DA TABELA, nao do topo da
                                     secao
  RasCheck.run_all                   a checagem do RAS antes de computar

Cada uma e opcional (ver Opcoes.usar_*), e o que falhar e reportado sem
derrubar o resto: uma versao do ras-commander sem alguma delas nao pode
impedir o modelo de ser gerado.
"""


import glob
import os


def _tentar(nome, fn, log):
    try:
        r = fn()
        log(f"      {nome}: ok" + (f" -- {r}" if isinstance(r, str) else ""))
        return r
    except (ImportError, AttributeError) as e:
        log(f"      {nome}: indisponivel nesta versao do ras-commander ({e})")
    except Exception as e:                                   # noqa: BLE001
        log(f"      {nome}: FALHOU -- {e}")
    return None


def aplicar(op, g01, log=print):
    """Roda as ferramentas de correcao do HEC-RAS sobre a geometria escrita.

    backup=False em todas, e UMA copia nossa antes de comecar. O padrao do
    RasFixit e backup=True, e ele grava o .g01 INTEIRO a cada secao editada:
    com 2.077 secoes e tres passes sairam 1.735 arquivos de 12 MB, 21 GB, e o
    disco chegou a 100% -- o passo seguinte morreu com "Espaco insuficiente no
    disco". Uma copia basta: a geometria se regenera pelo passo 7 em ~90 s.
    """
    import shutil
    resultados = {}
    copia = str(g01) + ".antes_da_correcao"
    shutil.copy2(g01, copia)
    log(f"      copia unica antes da correcao: {os.path.basename(copia)}")

    if op.usar_fixit and getattr(op, "corrigir_margens", False):
        from ras_commander import RasFixit
        resultados["bank"] = _tentar(
            "RasFixit.fix_bank_stations",
            lambda: RasFixit.fix_bank_stations(g01, backup=False), log)
    elif op.usar_fixit:
        log("      fix_bank_stations pulado: o build_cross_section ja insere "
            "as estacas das margens")
        log("         (17,5 min e 26 GB na rodada anterior para mudar zero "
            "secoes; corrigir_margens=true religa)")

    if op.usar_fixit:
        from ras_commander import RasFixit
        resultados["htab_el"] = _tentar(
            "RasFixit.fix_htab_starting_elevations",
            lambda: RasFixit.fix_htab_starting_elevations(g01, backup=False),
            log)

    if op.usar_ineffective:
        from ras_commander import RasFixit
        # Secao larga demais atravessa meandro do proprio rio e conta o
        # escoamento duas vezes. A area inefetiva e o remedio previsto: a agua
        # ocupa o volume mas nao conduz. Foi o que faltou quando alargar a
        # secao derrubou a simulacao de 30 passos para 2.
        resultados["ineffective"] = _tentar(
            "RasFixit.fix_ineffective_flow",
            lambda: RasFixit.fix_ineffective_flow(g01, backup=False), log)

    # As ferramentas gravam .bak POR SECAO mesmo com backup=False -- o
    # parametro nao alcanca todos os caminhos internos. Com 2.077 secoes
    # sairam 1.811 copias de 12 MB, 21 GB, e o disco foi a 100% duas vezes.
    # Limpar depois funciona qualquer que seja a ferramenta culpada.
    lixo = [f for f in glob.glob(str(g01) + ".bak*")]
    if lixo:
        n = sum(os.path.getsize(f) for f in lixo)
        for f in lixo:
            try:
                os.remove(f)
            except OSError:
                pass
        log(f"      removidos {len(lixo)} backups automaticos "
            f"({n/1e9:.1f} GB); a copia unica acima e suficiente")
    return resultados


def checar(op, prj, log=print):
    """RasCheck.run_all -- a checagem do proprio HEC-RAS, antes de computar.

    DUAS COISAS QUE ESTAVAM ERRADAS AQUI, e juntas produziram um "ok" numa
    geometria que o solver recusou minutos depois.

    NINGUEM OLHAVA O RESULTADO. `run_all` devolve um CheckResults com
    flow_type, get_error_count() e as mensagens. O codigo antigo so verificava
    se a chamada nao levantou excecao e logava "ok". Uma checagem que encontra
    erros e retorna normalmente era lida como aprovacao.

    ELA RODAVA CEGA. Sem a geometria compilada o HdfXsec nao abre o .g01.hdf
    ("file signature not found"), o RasCheck cai para "geometry_only" e as
    checagens de secao nao acontecem. O log dizia, uma linha antes do "ok":
    "No results found in plan HDF, running geometry-only checks". Compilar
    antes custa segundos e e o que da a ela o que conferir.

    Agora: compila, roda, LE o resultado, e checagem cega nunca vira "ok".
    """
    if not op.usar_check:
        return None
    from ras_commander import RasCheck

    g01hdf = op.caminho(f"{op.projeto}.g01.hdf")
    if not os.path.exists(g01hdf):
        from ras_commander import init_ras_project
        from ras_commander.geom import GeomMesh
        _tentar("compilar a geometria antes de checar",
                lambda: GeomMesh.compile_geometry(
                    "01", ras_object=init_ras_project(prj, op.ras_exe)), log)

    r = _tentar("RasCheck.run_all", lambda: RasCheck.run_all(prj), log)
    if r is None:
        log("      AVISO: RasCheck nao rodou -- a geometria segue SEM conferir")
        return None

    n_err = int(getattr(r, "get_error_count", lambda: 0)())
    n_avi = int(getattr(r, "get_warning_count", lambda: 0)())
    tipo = str(getattr(r, "flow_type", "?"))
    cega = not os.path.exists(g01hdf) or "geometry_only" in tipo.lower()
    log(f"      RasCheck: {tipo}, {n_err} erro(s), {n_avi} aviso(s)"
        + ("  -- CEGA, sem a geometria compilada" if cega else ""))

    if n_err or cega:
        try:
            df = r.to_dataframe()
            for _, x in df.head(15).iterrows():
                log("         " + " | ".join(f"{k}={x[k]}" for k in df.columns
                                             if k in ("severity", "check_type",
                                                      "river", "station",
                                                      "message"))[:160])
        except Exception:                                    # noqa: BLE001
            pass
    return r


def conferir_g01(g01, log=print):
    """Le de volta o .g01 GRAVADO e confere o que o solver recusa.

    Nao substitui o RasCheck: e a rede que pega o caso em que ele nao roda,
    roda cego ou nao cobre. Confere exatamente as tres coisas que ja pararam
    este modelo, no texto que o HEC-RAS vai ler:

        estaca repetida       "Station and elevation data contains duplicate
                              points" -- recusa o modelo INTEIRO por uma secao
        estaca fora de ordem  perfil que volta sobre si
        Bank Sta solta        margem que nao coincide com nenhuma estaca

    Devolve lista de (rio, rs, motivo). Ler de volta o arquivo e o unico jeito
    de conferir o que foi gravado: tudo que se checa antes checa a intencao.
    """
    import re
    achados = []
    txt = open(g01, encoding="latin-1", errors="replace").read()
    rio = reach = "?"
    for bloco in re.split(r"(?=River Reach=|Type RM Length)", txt):
        if bloco.startswith("River Reach="):
            p = bloco.split("\n", 1)[0][len("River Reach="):].split(",")
            rio, reach = p[0].strip(), (p[1].strip() if len(p) > 1 else "")
            continue
        if not bloco.startswith("Type RM Length"):
            continue
        m = re.match(r"Type RM Length L Ch R = 1\s*,\s*([\d.]+)", bloco)
        if not m:
            continue
        rs = m.group(1)
        q = re.search(r"#Sta/Elev=\s*(\d+)\s*\n(.*?)(?=\n[#A-Z])", bloco, re.S)
        if not q:
            continue
        n = int(q.group(1))
        cru = q.group(2).replace("\n", "")
        campos = [cru[k:k + 8] for k in range(0, 2 * n * 8, 8)]
        if len(campos) < 2 * n:
            achados.append((rio, rs, f"bloco #Sta/Elev truncado "
                                     f"({len(campos)//2} de {n} pontos)"))
            continue
        sta = campos[::2]
        rep = [i for i in range(1, len(sta)) if sta[i] == sta[i - 1]]
        if rep:
            achados.append((rio, rs, f"estaca repetida no(s) ponto(s) "
                                     f"{','.join(str(i+1) for i in rep[:6])}"))
        try:
            v = [float(s) for s in sta]
            fora = [i for i in range(1, len(v)) if v[i] < v[i - 1]]
            if fora:
                achados.append((rio, rs, f"estaca fora de ordem no ponto "
                                         f"{fora[0]+1}"))
        except ValueError:
            achados.append((rio, rs, "estaca ilegivel no #Sta/Elev"))
        b = re.search(r"Bank Sta=([\d.\-]+),([\d.\-]+)", bloco)
        if b:
            # os dois lados SEM espaco. Comparar "  277.15" formatado em 8
            # colunas contra um conjunto ja sem espaco acusou 22.502 margens
            # soltas num arquivo em que nenhuma estava -- checagem errada e
            # pior que checagem nenhuma, porque parece que conferiu.
            conj = {s.strip() for s in sta}
            for lado, val in (("esquerda", b.group(1)), ("direita", b.group(2))):
                if f"{float(val):.2f}" not in conj:
                    achados.append((rio, rs, f"Bank Sta {lado} ({val}) nao "
                                             f"coincide com nenhuma estaca"))
    achados += _conferir_contra_a_biblioteca(g01, log)
    if achados:
        log(f"      conferencia do .g01: {len(achados)} problema(s)")
        for r_, rs_, motivo in achados[:12]:
            log(f"         {r_} RS {rs_}: {motivo}")
    else:
        log("      conferencia do .g01: nenhuma estaca repetida, fora de "
            "ordem ou margem solta")
    return achados


def _conferir_contra_a_biblioteca(g01, log=print, n=25, semente=20260818):
    """A BIBLIOTECA COMO ORACULO: confere a nossa leitura numa amostra.

    O contrato do ras-commander (`ras_commander/geom/AGENTS.md`) diz, com
    todas as letras: "use existing parser and formatter helpers rather than
    hand-rolling string slicing". A varredura acima e exatamente isso -- corta
    o arquivo em campos de 8 caracteres a mao --, e a primeira versao dela
    acusou 22.502 margens soltas num arquivo onde nao havia nenhuma, por
    comparar "  277.15" com espacos contra valores sem espaco.

    Mas o helper e POR SECAO e reanalisa os 18,5 MB a cada chamada: medido em
    0,103 s, o que da 19 minutos para as 11.251 secoes -- a mesma patologia
    que fez o passo 8 queimar 44 minutos em Python de thread unica.
    `get_cross_sections`, que le o arquivo inteiro numa passada, leva 0,5 s.

    Entao: varredura nossa para percorrer tudo, e o helper conferindo uma
    AMOSTRA. Custa ~2,6 s e transforma "confio no meu fatiador" em "meu
    fatiador concorda com o parser da biblioteca nestas 25 secoes". Divergiu,
    e problema meu, e sai dito assim -- porque a referencia e ela.
    """
    try:
        import random

        import numpy as np

        from ras_commander.geom import GeomCrossSection as G
    except ImportError:
        return []
    try:
        meta = G.get_cross_sections(g01)
    except Exception as e:                                   # noqa: BLE001
        log(f"      (oraculo indisponivel: {e})")
        return []
    if not len(meta):
        return []
    linhas = meta.sample(min(n, len(meta)),
                         random_state=semente).itertuples()
    txt = open(g01, encoding="latin-1", errors="replace").read()
    fora = []
    for r in linhas:
        rio, reach, rs = str(r.River).strip(), str(r.Reach).strip(), str(r.RS).strip()
        try:
            df = G.get_station_elevation(g01, rio, reach, rs)
        except Exception:                                    # noqa: BLE001
            continue
        if df is None or not len(df):
            continue
        nosso = _estacas_da_secao(txt, rio, rs)
        if nosso is None:
            fora.append((rio, rs, "a varredura nao achou a secao que a "
                                  "biblioteca acha"))
            continue
        deles = [f"{v:.2f}" for v in df.iloc[:, 0].astype(float)]
        if len(nosso) != len(deles):
            fora.append((rio, rs, f"a varredura leu {len(nosso)} pontos e a "
                                  f"biblioteca {len(deles)}"))
        elif any(a != b for a, b in zip(nosso, deles)):
            k = next(i for i, (a, b) in enumerate(zip(nosso, deles)) if a != b)
            fora.append((rio, rs, f"estaca {k+1} difere da biblioteca: "
                                  f"{nosso[k]} contra {deles[k]}"))
    if fora:
        log(f"      ORACULO: a varredura DIVERGE da biblioteca em "
            f"{len(fora)} de {min(n, len(meta))} secoes -- o defeito e nosso")
    else:
        log(f"      oraculo: varredura confere com a biblioteca em "
            f"{min(n, len(meta))} secoes")
    return fora


def _estacas_da_secao(txt, rio, rs):
    """As estacas de UMA secao, pela mesma varredura usada em conferir_g01.

    O RIO TAMBEM, e nao so a RS. Sem isso a busca devolvia a primeira secao do
    ARQUIVO com aquela river station, e RS repete entre rios -- varios tem
    secao na RS 75. O oraculo pegou isto na primeira vez que rodou, acusando
    Itajai_Norte RS 75.00 com 280 pontos contra 59 da biblioteca: eram duas
    secoes de rios diferentes sendo comparadas uma com a outra.
    """
    import re
    atual = None
    for bloco in re.split(r"(?=River Reach=|Type RM Length)", txt):
        if bloco.startswith("River Reach="):
            atual = bloco.split("\n", 1)[0][len("River Reach="):].split(",")[0].strip()
            continue
        if not bloco.startswith("Type RM Length") or atual != str(rio).strip():
            continue
        m = re.match(r"Type RM Length L Ch R = 1\s*,\s*([\d.]+)", bloco)
        if not m or abs(float(m.group(1)) - float(rs)) > 0.005:
            continue
        q = re.search(r"#Sta/Elev=\s*(\d+)\s*\n(.*?)(?=\n[#A-Z])", bloco, re.S)
        if not q:
            return None
        n = int(q.group(1))
        cru = q.group(2).replace("\n", "")
        return [cru[k:k + 8].strip() for k in range(0, 2 * n * 8, 16)]
    return None


def auditar_terreno(op, g01, geotiffs, log=print):
    """Confere o que foi GRAVADO contra o terreno, na mesma linha.

    O resto da auditoria olha a cutline; esta olha as COTAS. Sem ela um erro de
    escavacao rebaixou 1.232 secoes a uma cota unica e todos os testes
    continuaram verdes -- uma bacia chata com parede vertical passa em "salto
    de area" e em "secao rasa".
    """
    try:
        from hecras_qc.modelo import rodar
    except ImportError:
        log("      hecras_qc.modelo indisponivel; auditoria de terreno pulada")
        return None
    return rodar(g01, geotiffs[0], f"EPSG:{op and 31982}")


def auditar_secoes(op, g01, geotiff, eixo_geojson, log=print):
    """QC geometrico das cutlines contra o terreno (hecras_qc)."""
    try:
        from hecras_qc import cross_sections, qc as _qc
        from hecras_qc.dem import DEM
        from hecras_qc.ras_geometry import exportar
        from hecras_qc.river_axis import EixoRio
    except ImportError:
        log("      hecras_qc indisponivel; QC das secoes pulado")
        return None
    _, secoes_geojson = exportar(g01, "EPSG:31982")
    dem = DEM(geotiff)
    eixo = EixoRio.ler(eixo_geojson, dem.crs_metrico)
    sec = cross_sections.carregar(secoes_geojson, dem.crs_metrico, eixo)
    lim = _qc.Limiares(espacamento=max(op.res_corredor * 2.0, 2.0))
    for s in sec:
        s.extrair(dem, lim.espacamento, eixo, lim.proeminencia_min)
    _qc.avaliar_todas(sec, lim)
    c = _qc.contagem(sec)
    log(f"      QC: OK {c[_qc.OK]} | atencao {c[_qc.ATENCAO]} | "
        f"incerto {c[_qc.INCERTO]} | critica {c[_qc.CRITICA]}")
    return sec
