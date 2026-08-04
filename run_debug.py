import sys
import traceback

with open("debug_output.txt", "w", encoding="utf-8") as log:
    try:
        log.write("Iniciando 03_gerar_geometria.py...\n")
        log.flush()
        import subprocess
        # Run 03_gerar_geometria.py and capture output
        res = subprocess.run([r"C:\Users\haas\miniforge3\python.exe", "03_gerar_geometria.py"],
                             capture_output=True, text=True, cwd=r"C:\Users\haas\github\hecras_vale")
        log.write(f"STDOUT:\n{res.stdout}\n")
        log.write(f"STDERR:\n{res.stderr}\n")
        log.write(f"RETURNCODE: {res.returncode}\n")
    except Exception as e:
        log.write(f"EXCEPTION:\n{traceback.format_exc()}\n")
print("Done debug run.")
