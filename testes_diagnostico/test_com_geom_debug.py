import win32com.client as win32
import os

prj_path = r"C:\Users\haas\github\hecras_vale\Itajai_Bacia_Real.prj"
print(f"Testando abertura de {prj_path} no HECRASController via COM...")

rc = win32.Dispatch("RAS701.HECRASController")
rc.Project_Open(prj_path)
print("Projeto aberto no Controller!")

# Tentar ler se há erros geométricos
n_err = rc.Geom_GetNumCrossSections()
print(f"Total de seções lidas pelo HEC-RAS: {n_err}")

rc.QuitRAS()
