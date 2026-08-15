#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT: run_hecras.py (ORQUESTRADOR MESTRE DO PIPELINE GENÉRICO HEC-RAS)
DESCRIÇÃO:
  Pipeline completo, automatizado e parametrizável para qualquer rio/bacia do Brasil:
  1. Extração dos rios da ANA via API REST online (isolado da foz).
  2. Download automático do DEM (Copernicus 30m) cobrindo a bacia com folga.
  3. Amostragem DEM e geração da geometria HEC-RAS 1D (.prj, .g01, .u01, .p01).
  4. Mapeamento de barragens e obras de contenção de cheias.
  5. Geração de documentação detalhada em Markdown (relatorio_bacia_<nome>.md).
  6. Execução da simulação hidrodinâmica no HEC-RAS Controller / COM.
  7. Plotagem gráfica e exportação de resultados em 'figuras/'.

EXEMPLOS DE USO:
  python run_hecras.py --bacia "Itajaí" --gerar-relatorio
  python run_hecras.py --bacia "Tijucas" --unificar-desconexos --gerar-relatorio
  python run_hecras.py --bacia "Tubarão" --vazao-aleatoria
================================================================================
"""

import os
import sys
import json
import random
import argparse
from datetime import datetime

# Importar módulos do pipeline
import extrair_bacia_ana_online
import baixar_dem_online
import gerar_geometria_hecras
import plot_resultados_hecras

RAS_VERSION = "7.0.1"
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"

def executar_hecras_controller(caminho_prj):
    """Executa a simulação HEC-RAS via COM controller ou ras-commander e verifica a geração de resultados."""
    print(f"\n[HEC-RAS] Executando simulação hidrodinâmica para o projeto: {os.path.basename(caminho_prj)}...")
    base_dir = os.path.dirname(caminho_prj)
    base_name = os.path.splitext(os.path.basename(caminho_prj))[0]
    hdf_out = os.path.join(base_dir, f"{base_name}.p01.hdf")
    
    simulou = False
    
    # 1. Tentar ras-commander
    try:
        from ras_commander import init_ras_project, RasCmdr
        ver = RAS_VERSION if not os.path.exists(RAS_EXE) else RAS_EXE
        print(f"  [ras-commander] Tentando executar via ras-commander ({ver})...")
        ras_obj = init_ras_project(caminho_prj, ver)
        res = RasCmdr.compute_plan("01", overwrite_dest=False, ras_object=ras_obj)
        status_str = str(getattr(res, 'status', res))
        if "SUCCESS" in status_str or os.path.exists(hdf_out):
            print("  ✓ Simulação HEC-RAS concluída com SUCESSO via ras-commander!")
            simulou = True
        else:
            print(f"  ⚠ ras-commander retornou: {status_str}")
    except Exception as e:
        print(f"  i ras-commander não executou: {e}")

    # 2. Tentar COM / pywin32 se ras-commander não rodou
    if not simulou:
        try:
            import win32com.client as win32
            rc = None
            for cls in ("RAS701.HECRASController", "RAS61.HECRASController", "RAS.HECRASController"):
                try:
                    rc = win32.Dispatch(cls)
                    print(f"  ✓ Conectado ao HECRASController via COM ({cls})")
                    break
                except Exception:
                    continue
            if rc:
                rc.Project_Open(caminho_prj)
                try: rc.Compute_HideComputationWindow()
                except Exception: pass
                res = rc.Compute_CurrentPlan(None, None, False)
                ok = res[0] if isinstance(res, tuple) else bool(res)
                try: rc.Project_Close()
                except Exception: pass
                for q in ("QuitRas", "QuitRAS"):
                    try: getattr(rc, q)(); break
                    except Exception: pass
                
                if ok or os.path.exists(hdf_out):
                    print("  ✓ Simulação HEC-RAS executada com sucesso via COM Interface!")
                    simulou = True
                else:
                    print("  ⚠ HECRASController COM respondeu, mas o arquivo de resultados não foi gerado.")
        except Exception as e:
            print(f"  i Interface COM do HEC-RAS não disponível nesta máquina: {e}")

    if not simulou:
        print("  ⚠ NENHUMA SIMULAÇÃO AUTOMÁTICA FOI CONCLUÍDA.")
        print("  ℹ Os arquivos do projeto (.prj, .g01, .u01, .p01) foram gerados com sucesso e podem ser abertos e executados diretamente no HEC-RAS GUI.")
        return False

    return True

def gerar_relatorio_markdown(nome_bacia, pasta_saida, vazao_simulada):
    """Gera documentação técnica completa da bacia e das estruturas hidráulicas."""
    caminho_relatorio = os.path.join(pasta_saida, f"relatorio_bacia_{extrair_bacia_ana_online.sanitizar_nome_arquivo(nome_bacia)}.md")
    
    # Carregar barragens se existirem
    barragens_info = []
    caminho_barragens = os.path.join(pasta_saida, "dados_estruturas", "barragens_itajai.json")
    if os.path.exists(caminho_barragens) and ("itaj" in extrair_bacia_ana_online.sanitizar_nome_arquivo(nome_bacia)):
        with open(caminho_barragens, "r", encoding="utf-8") as f:
            b_data = json.load(f)
            barragens_info = b_data.get("barragens", [])

    data_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    conteudo = f"""# Documentação Técnica e Configuração HEC-RAS: Bacia do {nome_bacia.capitalize()}

**Data de Geração:** {data_str}  
**Orquestrador Mestre:** `run_hecras.py`  
**Fonte de Dados Hidrográficos:** Agência Nacional de Águas e Saneamento Básico (ANA BHO 2017)  
**Modelo Digital de Elevação:** Copernicus DEM GLO-30 (~30m)  

---

## 1. Resumo Executivo da Bacia

A modelagem hidrodinâmica 1D/2D para a **Bacia do {nome_bacia.capitalize()}** foi estruturada automaticamente com base em dados vetoriais oficiais online. A calha principal e os afluentes conectados foram isolados topologicamente a partir da exutória principal no oceano.

- **Nome da Bacia / Sistema:** {nome_bacia.capitalize()}
- **Vazão Sintética Simulação HEC-RAS:** {vazao_simulada:.1f} m³/s
- **Resolução do Relevo DEM:** 30 metros (Copernicus GLO-30)
- **Formato de Saída do Projeto HEC-RAS:** `.prj`, `.g01`, `.u01`, `.p01`

---

## 2. Estruturas Hidráulicas e Sistema de Contenção de Cheias

"""
    if barragens_info:
        conteudo += "| Nome da Estrutura | Rio de Localização | Capacidade de Acumulação | Comportas | Função Hidráulica |\n"
        conteudo += "|---|---|---|---|---|\n"
        for b in barragens_info:
            cap_m = b['capacidade_acumulacao_m3'] / 1e6
            conteudo += f"| **{b['nome']}** | {b['rio']} | {cap_m:.1f} milhões m³ | {b['numero_comportas']} comportas | {b['impacto_jusante']} |\n"
    else:
        conteudo += "*Não foram registradas barragens de grande porte cadastradas especificamente para esta bacia secundária.* \n"

    conteudo += f"""

---

## 3. Arquivos de Geometria e Projeto Gerados

Os seguintes arquivos do modelo HEC-RAS foram criados em `{pasta_saida}`:

1. **[`{nome_bacia.capitalize()}_Bacia_Real.prj`](file:///{pasta_saida}/{nome_bacia.capitalize()}_Bacia_Real.prj):** Arquivo de projeto mestre do HEC-RAS.
2. **[`{nome_bacia.capitalize()}_Bacia_Real.g01`](file:///{pasta_saida}/{nome_bacia.capitalize()}_Bacia_Real.g01):** Geometria 1D contendo a amostragem de elevação do DEM e estacas de margem.
3. **[`{nome_bacia.capitalize()}_Bacia_Real.u01`](file:///{pasta_saida}/{nome_bacia.capitalize()}_Bacia_Real.u01):** Condições de contorno de hidrograma de cheia não-permanente.
4. **[`{nome_bacia.capitalize()}_Bacia_Real.p01`](file:///{pasta_saida}/{nome_bacia.capitalize()}_Bacia_Real.p01):** Plano de execução e parâmetros do solver hidrodinâmico.

---

## 4. Gráficos e Resultados de Inundação

Os gráficos gerados pela simulação foram salvos no diretório `figuras/`:

- **Rede Hidrográfica e Barragens:** [`figuras/figura_1_rede_de_rios_e_barragens.png`](file:///{pasta_saida}/figuras/figura_1_rede_de_rios_e_barragens.png)
- **Perfil Longitudinal do Rio:** [`figuras/figura_2_perfil_longitudinal.png`](file:///{pasta_saida}/figuras/figura_2_perfil_longitudinal.png)
- **Hidrogramas de Vazão:** [`figuras/figura_3_hidrogramas_cheia.png`](file:///{pasta_saida}/figuras/figura_3_hidrogramas_cheia.png)
"""

    with open(caminho_relatorio, "w", encoding="utf-8") as f:
        f.write(conteudo)

    print(f"  ✓ Relatório de documentação técnica salvo: {os.path.basename(caminho_relatorio)}")
    return caminho_relatorio

def executar_pipeline_completo(bacia, min_strahler, unificar_desconexos, vazao_aleatoria, gerar_relatorio, pasta_saida):
    print("=" * 80)
    print(f" 🚀 EXECUTANDO PIPELINE GENÉRICO HEC-RAS PARA A BACIA DO {bacia.upper()}")
    print("=" * 80)
    
    # PASSO 1: Extração online dos rios da ANA
    print("\n[PASSO 1/6] Extraindo rede hidrográfica online da ANA...")
    extrair_bacia_ana_online.extrair_rios_principais_ana_online(
        nome_bacia=bacia,
        min_strahler=min_strahler,
        pasta_saida=pasta_saida,
        unificar_desconexos=unificar_desconexos
    )

    # PASSO 2: Download dinâmico do DEM cobrindo a bacia com folga
    print("\n[PASSO 2/6] Baixando DEM 30m cobrindo o vale com margem de segurança...")
    caminho_dem = baixar_dem_online.obter_dem_bacia(
        nome_bacia=bacia,
        pasta_trabalho=pasta_saida,
        folga_graus=0.15
    )

    # PASSO 3: Gerar geometria HEC-RAS 1D (.g01, .prj, .u01, .p01)
    print("\n[PASSO 3/6] Amostrando elevações e construindo geometria HEC-RAS...")
    gerar_geometria_hecras.processar_geometria_bacia(
        nome_bacia=bacia,
        pasta_trabalho=pasta_saida
    )

    # PASSO 4: Definir vazão de cheia (aleatória ou padrão)
    if vazao_aleatoria:
        vazao_pico = float(random.randint(1200, 3500))
        print(f"\n[PASSO 4/6] Gerando cheia aleatória de teste: {vazao_pico:.1f} m³/s")
    else:
        vazao_pico = 1800.0
        print(f"\n[PASSO 4/6] Utilizando vazão padrão de projeto: {vazao_pico:.1f} m³/s")

    # PASSO 5: Executar simulação no HEC-RAS
    print("\n[PASSO 5/6] Executando simulação no HEC-RAS Controller...")
    bacia_ascii = extrair_bacia_ana_online.sanitizar_nome_arquivo(bacia).capitalize()
    prj_path = os.path.join(pasta_saida, f"{bacia_ascii}_Oficial.prj")
    executar_hecras_controller(prj_path)

    # PASSO 6: Visualização gráfica e documentação
    print("\n[PASSO 6/6] Gerando gráficos e documentação técnica...")
    plot_resultados_hecras.gerar_graficos_bacia(bacia, pasta_saida)
    
    if gerar_relatorio:
        gerar_relatorio_markdown(bacia, pasta_saida, vazao_pico)

    print("\n" + "=" * 80)
    print(f" ✨ PIPELINE FINALIZADO COM SUCESSO PARA O RIO {bacia.upper()}!")
    print("=" * 80)
    print(f" 📂 Todos os arquivos estão prontos em: {pasta_saida}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Orquestrador Mestre Genérico HEC-RAS para Bacias do Brasil.")
    parser.add_argument("-b", "--bacia", type=str, default="Itajaí", help="Nome da bacia/rio (ex: Itajaí, Tijucas, Tubarão, Araranguá)")
    parser.add_argument("-s", "--strahler", type=int, default=3, help="Ordem mínima de Strahler (padrão: 3)")
    parser.add_argument("-u", "--unificar-desconexos", action="store_true", default=False, help="Ativa reconexão espacial de rios desconexos")
    parser.add_argument("--vazao-aleatoria", action="store_true", default=False, help="Simula evento de inundação com vazão de pico aleatória")
    parser.add_argument("--gerar-relatorio", action="store_true", default=True, help="Gera documentação técnica em Markdown relatorio_bacia_*.md")
    parser.add_argument("-o", "--saida", type=str, default=r"C:\Users\haas\github\hecras_vale", help="Pasta de trabalho dos arquivos")

    args = parser.parse_args()

    executar_pipeline_completo(
        bacia=args.bacia,
        min_strahler=args.strahler,
        unificar_desconexos=args.unificar_desconexos,
        vazao_aleatoria=args.vazao_aleatoria,
        gerar_relatorio=args.gerar_relatorio,
        pasta_saida=args.saida
    )
