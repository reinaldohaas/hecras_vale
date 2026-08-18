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
    """RasCheck.run_all -- a checagem do proprio HEC-RAS, antes de computar."""
    if not op.usar_check:
        return None
    from ras_commander import RasCheck
    return _tentar("RasCheck.run_all", lambda: RasCheck.run_all(prj), log)


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
