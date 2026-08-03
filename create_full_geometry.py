import os
import math
import numpy as np

def format_16_num(val):
    return f"{val:>16.4f}"

def format_8(val):
    return f"{val:>8.2f}"

def create_full_basin_hecras():
    project_name = "Itajai_Bacia_Completa"
    prj_file = f"{project_name}.prj"
    geom_file = f"{project_name}.g01"
    flow_file = f"{project_name}.u01"
    plan_file = f"{project_name}.p01"

    print("Gerando projeto HEC-RAS da Bacia Completa com formato .u01 limpo e compatível...")

    # 1. PRJ File
    with open(prj_file, "w") as f:
        f.write(f"Proj Title={project_name}\n")
        f.write("Current Plan=p01\n")
        f.write("Default Exp/Contr=0.3,0.1\n")
        f.write("SI Units\n")
        f.write("Geom File=g01\n")
        f.write("Unsteady File=u01\n")
        f.write("Plan File=p01\n")
        f.write("Y Axis Title=Elevation\n")
        f.write("X Axis Title(PR)=Distance\n")
        f.write("X Axis Title(CS)=Station\n")

    # 2. Geometry File (.g01)
    with open(geom_file, "w") as f:
        f.write("Geom Title=Geometria Completa Vale do Itajai com 3 Barragens\n")
        f.write("Program Version=7.01\n")

        # --- TRECHO 1: Rio Itajaí do Sul ---
        f.write("River Reach= Itajai_do_Sul,Trecho_Sul\n")
        st_sul = np.arange(100000, -1000, -1000)
        f.write(f"Reach XY= 2 \n")
        f.write(format_16_num(650000) + format_16_num(6970000) + format_16_num(660000) + format_16_num(7000000) + "\n")
        
        for st in st_sul:
            z_bottom = 350.0 + (st / 100000.0) * 100.0
            st_str = f"{st:.2f}"
            rl = 1000.0 if st > 0 else 0.0
            
            f.write(f"Type RM Length L Ch R = 1 , {st_str:>8} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
            f.write("Node Last Edited Time= Aug/02/2026 00:00:00\n")
            
            if st == 50000:
                f.write("Inline Structure= 1 , 50000.00 , 0 \n")
                f.write("Type= Dam\n")
                f.write("Inline Structure Title= Barragem Sul (Ituporanga)\n")
                f.write("Inline Structure Node Last Edited Time= Aug/02/2026 00:00:00\n")
                f.write("Inline Structure Weir= 1 \n")
                f.write("Weir Crest Elev= 390.0 \n")
                f.write("Weir Coeff= 1.6 \n")
                f.write("Inline Structure Gate= 5 \n")
                f.write("Gate Group= Comportas Sul\n")
                f.write("Gate Shape= Rectangular\n")
                f.write("Gate Width= 2.5 \n")
                f.write("Gate Height= 2.5 \n")
                f.write("Gate Invert= 360.0 \n")
            
            f.write("#Sta/Elev= 5 \n")
            pts = [(-80, z_bottom + 8), (-50, z_bottom), (0, z_bottom), (50, z_bottom), (80, z_bottom + 8)]
            line = "".join([format_8(px) + format_8(py) for px, py in pts])
            f.write(line + "\n")
            f.write("Bank Sta=-50, 50\n")
            f.write("#Mann= 3 , -1 , 0 \n")
            f.write(format_8(-80) + format_8(0.05) + format_8(0) + format_8(-50) + format_8(0.035) + format_8(0) + format_8(50) + format_8(0.05) + format_8(0) + "\n")

        # --- TRECHO 2: Rio Itajaí do Oeste ---
        f.write("River Reach= Itajai_do_Oeste,Trecho_Oeste\n")
        st_oeste = np.arange(100000, -1000, -1000)
        f.write(f"Reach XY= 2 \n")
        f.write(format_16_num(600000) + format_16_num(7010000) + format_16_num(660000) + format_16_num(7000000) + "\n")

        for st in st_oeste:
            z_bottom = 340.0 + (st / 100000.0) * 110.0
            st_str = f"{st:.2f}"
            rl = 1000.0 if st > 0 else 0.0
            
            f.write(f"Type RM Length L Ch R = 1 , {st_str:>8} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
            f.write("Node Last Edited Time= Aug/02/2026 00:00:00\n")
            
            if st == 50000:
                f.write("Inline Structure= 1 , 50000.00 , 0 \n")
                f.write("Type= Dam\n")
                f.write("Inline Structure Title= Barragem Oeste (Taió)\n")
                f.write("Inline Structure Node Last Edited Time= Aug/02/2026 00:00:00\n")
                f.write("Inline Structure Weir= 1 \n")
                f.write("Weir Crest Elev= 360.0 \n")
                f.write("Weir Coeff= 1.6 \n")
                f.write("Inline Structure Gate= 7 \n")
                f.write("Gate Group= Comportas Oeste\n")
                f.write("Gate Shape= Rectangular\n")
                f.write("Gate Width= 2.0 \n")
                f.write("Gate Height= 2.0 \n")
                f.write("Gate Invert= 345.0 \n")

            f.write("#Sta/Elev= 5 \n")
            pts = [(-80, z_bottom + 8), (-50, z_bottom), (0, z_bottom), (50, z_bottom), (80, z_bottom + 8)]
            line = "".join([format_8(px) + format_8(py) for px, py in pts])
            f.write(line + "\n")
            f.write("Bank Sta=-50, 50\n")
            f.write("#Mann= 3 , -1 , 0 \n")
            f.write(format_8(-80) + format_8(0.05) + format_8(0) + format_8(-50) + format_8(0.035) + format_8(0) + format_8(50) + format_8(0.05) + format_8(0) + "\n")

        # --- JUNÇÃO 1: Rio do Sul ---
        f.write("Junction= Junc_Rio_do_Sul, , 660000, 7000000\n")
        f.write("Upstream Reach= Itajai_do_Sul,Trecho_Sul\n")
        f.write("Upstream Reach= Itajai_do_Oeste,Trecho_Oeste\n")
        f.write("Downstream Reach= Itajai_Acu,Trecho_Principal\n")

        # --- TRECHO 3: Rio Itajaí-Açu ---
        f.write("River Reach= Itajai_Acu,Trecho_Principal\n")
        st_acu = np.arange(150000, -1000, -1000)
        f.write(f"Reach XY= 2 \n")
        f.write(format_16_num(660000) + format_16_num(7000000) + format_16_num(730000) + format_16_num(7020000) + "\n")

        for st in st_acu:
            z_bottom = -15.0 + (st / 150000.0) * 355.0
            st_str = f"{st:.2f}"
            rl = 1000.0 if st > 0 else 0.0

            f.write(f"Type RM Length L Ch R = 1 , {st_str:>8} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
            f.write("Node Last Edited Time= Aug/02/2026 00:00:00\n")

            f.write("#Sta/Elev= 5 \n")
            pts = [(-100, z_bottom + 12), (-75, z_bottom), (0, z_bottom), (75, z_bottom), (100, z_bottom + 12)]
            line = "".join([format_8(px) + format_8(py) for px, py in pts])
            f.write(line + "\n")
            f.write("Bank Sta=-75, 75\n")
            f.write("#Mann= 3 , -1 , 0 \n")
            f.write(format_8(-100) + format_8(0.05) + format_8(0) + format_8(-75) + format_8(0.035) + format_8(0) + format_8(75) + format_8(0.05) + format_8(0) + "\n")

        # --- TRECHO 4: Rio Itajaí do Norte ---
        f.write("River Reach= Itajai_do_Norte,Trecho_Norte\n")
        st_norte = np.arange(80000, -1000, -1000)
        f.write(f"Reach XY= 2 \n")
        f.write(format_16_num(640000) + format_16_num(7050000) + format_16_num(680000) + format_16_num(7010000) + "\n")

        for st in st_norte:
            z_bottom = 200.0 + (st / 80000.0) * 150.0
            st_str = f"{st:.2f}"
            rl = 1000.0 if st > 0 else 0.0

            f.write(f"Type RM Length L Ch R = 1 , {st_str:>8} , {rl:.2f},{rl:.2f},{rl:.2f}\n")
            f.write("Node Last Edited Time= Aug/02/2026 00:00:00\n")

            if st == 40000:
                f.write("Inline Structure= 1 , 40000.00 , 0 \n")
                f.write("Type= Dam\n")
                f.write("Inline Structure Title= Barragem Norte (Jose Boiteux)\n")
                f.write("Inline Structure Node Last Edited Time= Aug/02/2026 00:00:00\n")
                f.write("Inline Structure Weir= 1 \n")
                f.write("Weir Crest Elev= 300.0 \n")
                f.write("Weir Coeff= 1.6 \n")
                f.write("Inline Structure Gate= 2 \n")
                f.write("Gate Group= Comportas Norte\n")
                f.write("Gate Shape= Rectangular\n")
                f.write("Gate Width= 4.0 \n")
                f.write("Gate Height= 4.0 \n")
                f.write("Gate Invert= 270.0 \n")

            f.write("#Sta/Elev= 5 \n")
            pts = [(-80, z_bottom + 8), (-50, z_bottom), (0, z_bottom), (50, z_bottom), (80, z_bottom + 8)]
            line = "".join([format_8(px) + format_8(py) for px, py in pts])
            f.write(line + "\n")
            f.write("Bank Sta=-50, 50\n")
            f.write("#Mann= 3 , -1 , 0 \n")
            f.write(format_8(-80) + format_8(0.05) + format_8(0) + format_8(-50) + format_8(0.035) + format_8(0) + format_8(50) + format_8(0.05) + format_8(0) + "\n")

    # 3. Unsteady Flow File (.u01) - Sintaxe estritamente limpa idêntica ao Itajai_Blumenau.u01
    with open(flow_file, "w") as f:
        f.write("Flow Title=Cenario_Previsao_Bacia_Completa\n")
        f.write("Program Version=7.01\n")

        # Entrada Montante: Rio Itajaí do Sul
        f.write("Boundary Location= Itajai_do_Sul,Trecho_Sul, 100000.00 \n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 24 \n")
        f.write(" ".join(["1200"] * 24) + "\n")

        # Entrada Montante: Rio Itajaí do Oeste
        f.write("Boundary Location= Itajai_do_Oeste,Trecho_Oeste, 100000.00 \n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 24 \n")
        f.write(" ".join(["1500"] * 24) + "\n")

        # Entrada Montante: Rio Itajaí do Norte
        f.write("Boundary Location= Itajai_do_Norte,Trecho_Norte, 80000.00 \n")
        f.write("Interval= 1HOUR\n")
        f.write("Flow Hydrograph= 24 \n")
        f.write(" ".join(["2000"] * 24) + "\n")

        # Jusante: Foz do Rio Itajaí-Açu
        f.write("Boundary Location= Itajai_Acu,Trecho_Principal, 0.00 \n")
        f.write("Interval= 1HOUR\n")
        f.write("Stage Hydrograph= 24 \n")
        f.write(" ".join(["0"] * 12 + ["1"] * 12) + "\n")

        f.write("Initial Stage= 2 \n")
        f.write("Initial Flow= 1000 \n")

    # 4. Plan File (.p01)
    with open(plan_file, "w") as f:
        f.write("Plan Title=Simulacao_Bacia_Completa\n")
        f.write("Program Version=7.01\n")
        f.write("Short Identifier=001\n")
        f.write("Simulation Date=01AUG2026,00,02AUG2026,00\n")
        f.write("Geom File=g01\n")
        f.write("Flow File=u01\n")
        f.write("Subcritical Flow\n")
        f.write("Computation Interval=1MIN\n")
        f.write("Output Interval=1HOUR\n")
        f.write("Instantaneous Interval=1HOUR\n")
        f.write("Mapping Interval=1HOUR\n")
        f.write("Run HTab=-1\n")
        f.write("Run UNet=-1\n")
        f.write("Run PostProcess=-1\n")
        f.write("Run RASMapper=-1\n")

    print("Projeto HEC-RAS com sintaxe limpa de .u01 gerado com sucesso!")

if __name__ == "__main__":
    create_full_basin_hecras()
