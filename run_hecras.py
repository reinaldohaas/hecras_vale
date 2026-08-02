import os
import win32com.client
import time

def get_hecras_controller():
    """Conecta ao HEC-RAS via COM Interface."""
    com_classes = [
        "RAS701.HECRASController",
        "RAS610.HECRASController",
        "RAS507.HECRASController",
        "RAS.HECRASController"
    ]
    
    rc = None
    for com_cls in com_classes:
        try:
            print(f"Tentando conectar ao {com_cls}...")
            rc = win32com.client.Dispatch(com_cls)
            print(f"Sucesso ao conectar: {com_cls}")
            break
        except Exception:
            continue
            
    if rc is None:
        raise Exception("Não foi possível conectar ao HEC-RAS Controller.")
    return rc

def run_hecras_simulation():
    current_dir = os.path.abspath(os.path.dirname(__file__))
    
    # Tenta primeiro o projeto da bacia completa, se existir, senão o projeto simples
    full_project = os.path.join(current_dir, "Itajai_Bacia_Completa.prj")
    simple_project = os.path.join(current_dir, "Itajai_Blumenau.prj")
    
    project_file = full_project if os.path.exists(full_project) else simple_project
    
    print(f"Abrindo projeto HEC-RAS: {project_file}")
    rc = get_hecras_controller()
    
    rc.Project_Open(project_file)
    print("Iniciando simulação hidráulica da bacia (Compute_CurrentPlan)...")
    
    try:
        res = rc.Compute_CurrentPlan(None, None, True)
        if isinstance(res, tuple):
            success = res[0]
            msg = res[1] if len(res) > 1 else ""
        else:
            success = res
            msg = ""
    except Exception as e:
        success = False
        msg = str(e)
        
    if success:
        print("\n==========================================")
        print("Simulação da Bacia Completa finalizada com SUCESSO!")
        print("Resultados disponíveis para extração de nível e vazão.")
        print("==========================================")
    else:
        print(f"\nAviso de Execução (Msg/Código): {msg}")
        print("Para compilar a nova geometria multi-trecho com as 3 barragens pela 1ª vez:")
        print("1. Abra o HEC-RAS GUI e carregue 'Itajai_Bacia_Completa.prj'.")
        print("2. Abra a Janela Geometric Data e clique em File > Save Geometric Data.")
        print("3. Vá em Run > Unsteady Flow Analysis e clique em COMPUTE.")

    rc.Project_Close()
    rc.QuitRAS()

if __name__ == "__main__":
    run_hecras_simulation()
