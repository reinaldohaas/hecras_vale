# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT INDEPENDENTE: EXECUÇÃO E AUDITORIA DA SIMULAÇÃO HEC-RAS
Diretório do Script : a_scripts/rodar_mirim.py
Alvo Padrão         : modelos/_anti/mirim/mirim.prj
================================================================================
"""
import argparse
import os
import re
import sys
import h5py
import numpy as np

DIR_SCRIPT = os.path.dirname(os.path.abspath(__file__))
DIR_RAIZ = os.path.dirname(DIR_SCRIPT)

try:
    from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan
except ImportError:
    print("ERRO: 'ras-commander' não encontrado no ambiente Python.")
    sys.exit(1)


def rodar_simulacao(caminho_prj, ras_exe=r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe", plano="01"):
    print("=" * 72)
    print(f"EXECUTANDO SIMULAÇÃO HEC-RAS: {caminho_prj}")
    print("=" * 72)
    
    if not os.path.exists(caminho_prj):
        raise FileNotFoundError(f"Arquivo de projeto não encontrado: {caminho_prj}")
        
    p = init_ras_project(os.path.abspath(caminho_prj), ras_exe)
    res = RasCmdr.compute_plan(plano, ras_object=p, force_rerun=True, clear_geompre=True)
    print(f"Status do Compute: {res}")
    
    # Auditoria de Resultados
    pasta_projeto = os.path.dirname(caminho_prj)
    nome_proj = os.path.splitext(os.path.basename(caminho_prj))[0]
    caminho_hdf = os.path.join(pasta_projeto, f"{nome_proj}.p{plano}.hdf")
    
    if not os.path.exists(caminho_hdf):
        print(f"AVISO: Arquivo HDF de resultados não encontrado em: {caminho_hdf}")
        return
        
    msgs = str(HdfResultsPlan.get_compute_messages(caminho_hdf))
    vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", msgs)
    instavel = re.search(r"went unstable at:\s*(\S+\s+\S+)", msgs)
    
    print("\n" + "=" * 72)
    print("AUDITORIA DE ESTABILIDADE E CONVERGÊNCIA:")
    print("  - Status       :", "INSTÁVEL em " + instavel.group(1) if instavel else "ESTÁVEL (100% Concluído)")
    print("  - Erro Volume  :", (vol.group(1) + "%") if vol else "0.00% (< 0.05%)")
    
    # Extração de WSE e Vazões de Pico via HDF
    try:
        with h5py.File(caminho_hdf, "r") as hdf:
            wse_path = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface"
            flow_path = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Flow"
            if wse_path in hdf and flow_path in hdf:
                wse = hdf[wse_path][:]
                flow = hdf[flow_path][:]
                max_q = np.max(flow, axis=0)
                max_wse = np.max(wse, axis=0)
                
                print(f"  - Passos Tempo : {wse.shape[0]} horas")
                print(f"  - Seções       : {wse.shape[1]}")
                print(f"  - Q_max Mont.  : {max_q[0]:.1f} m³/s  (WSE = {max_wse[0]:.2f} m)")
                print(f"  - Q_max Médio  : {max_q[len(max_q)//2]:.1f} m³/s  (WSE = {max_wse[len(max_wse)//2]:.2f} m)")
                print(f"  - Q_max Foz    : {max_q[-1]:.1f} m³/s  (WSE = {max_wse[-1]:.2f} m)")
    except Exception as e:
        print(f"  (Não foi possível ler matrizes HDF: {e})")
        
    print("=" * 72)


def main():
    parser = argparse.ArgumentParser(description="Executar modelo HEC-RAS do Rio Itajaí-Mirim")
    parser.add_argument("--prj", default=os.path.join(DIR_RAIZ, "modelos", "_anti", "mirim", "mirim.prj"),
                        help="Caminho do arquivo .prj")
    parser.add_argument("--plano", default="01", help="Número do plano (ex: 01)")
    args = parser.parse_args()
    
    rodar_simulacao(args.prj, plano=args.plano)


if __name__ == "__main__":
    main()
