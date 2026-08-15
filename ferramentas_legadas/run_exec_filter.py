import subprocess

with open("filter_log.txt", "w", encoding="utf-8") as f:
    res = subprocess.run([r"C:\Users\haas\miniforge3\python.exe", "filter_basin.py"], capture_output=True, text=True, cwd=r"C:\Users\haas\github\hecras_vale")
    f.write(f"STDOUT:\n{res.stdout}\n")
    f.write(f"STDERR:\n{res.stderr}\n")
    f.write(f"EXIT CODE: {res.returncode}\n")
print("Done.")
