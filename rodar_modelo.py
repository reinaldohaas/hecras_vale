#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Roda o modelo hidrodinamico da Bacia do Itajai de ponta a ponta.

    gerar_rede_hecras   ->  .g01/.u01/.p01/.prj
    HEC-RAS 7.0.1       ->  .p01.hdf
    gerar_mancha_hecras ->  app/manchas_inundacao_hecras.geojson

Uso:
    python rodar_modelo.py                # cheia sintetica (projeto Itajai_Rede)
    python rodar_modelo.py 1983           # evento historico com chuva real
    python rodar_modelo.py 1983 --sem-barragens
    python rodar_modelo.py --todos        # 1983, 2008, 2011 e 2023 em sequencia
    python rodar_modelo.py 1983 --so-mancha    # pula a simulacao

Eventos disponiveis: os arquivos itajai_flood_model/data/rainfall_events/
chuva_real_<EVENTO>.csv. Hoje: 1983, 2008, 2011, 2023, 2008_epagri_ana.

Codigo de saida: 0 se tudo correu bem, 1 caso contrario -- serve para
agendador de tarefas / CI.
"""
import os
import sys
import time
import subprocess

AQUI = os.path.abspath(os.path.dirname(__file__))
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"
EVENTOS_PADRAO = ["1983", "2008", "2011", "2023"]

# Referencia observada por evento: (local, RS em km da foz, Q de pico m3/s)
OBSERVADO = {
    "1983": [("Rio do Sul", 181.5, 3900), ("Blumenau", 67.7, 5850), ("Foz", 2.8, 6900)],
    "2008": [("Rio do Sul", 181.5, 650),  ("Blumenau", 67.7, 4200), ("Foz", 2.8, 5700)],
    "2011": [("Rio do Sul", 181.5, 3200), ("Blumenau", 67.7, 4650), ("Foz", 2.8, 5400)],
    "2023": [("Rio do Sul", 181.5, 2850), ("Blumenau", 67.7, 3950), ("Foz", 2.8, 4800)],
}


def log(msg, nivel=""):
    marca = {"ok": "  [OK] ", "erro": "  [ERRO] ", "": "  "}.get(nivel, "  ")
    print(f"{marca}{msg}", flush=True)


def titulo(t):
    print("\n" + "=" * 70, flush=True)
    print(t, flush=True)
    print("=" * 70, flush=True)


# ------------------------------------------------------------------ etapas
def gerar(evento, barragens):
    """Gera .g01/.u01/.p01/.prj. Devolve o nome do projeto."""
    import gerar_rede_hecras as G
    G.EVENTO = evento
    G.BARRAGENS = barragens
    G.PROJECT = "Itajai_Rede"           # main() acrescenta _<EVENTO> se houver
    G.main()
    return f"Itajai_Rede_{evento}" if evento else "Itajai_Rede"


def status_hdf(projeto):
    """Status VERDADEIRO da simulacao.

    O COM do HEC-RAS devolve 'Computations Completed' mesmo quando a
    simulacao foi interrompida por instabilidade; a verdade esta no HDF.
    """
    hdf = os.path.join(AQUI, f"{projeto}.p01.hdf")
    if not os.path.exists(hdf):
        return None, "sem arquivo de resultados", 0
    try:
        import h5py
        with h5py.File(hdf, "r") as f:
            s = f["Results/Unsteady/Summary"]
            sol = s.attrs.get("Solution", b"?")
            sol = sol.decode() if isinstance(sol, bytes) else str(sol)
            g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
                  "Unsteady Time Series/Cross Sections"]
            n = g["Water Surface"].shape[0]
        return ("Success" in sol or "Successfully" in sol), sol, n
    except Exception as e:
        return None, f"erro lendo o HDF: {e}", 0


def simular(projeto):
    """Roda o HEC-RAS pela interface COM."""
    try:
        import win32com.client as win32
    except ImportError:
        log("pywin32 nao instalado: pip install pywin32", "erro")
        return False
    prj = os.path.join(AQUI, f"{projeto}.prj")    # caminho Windows (barras \)
    rc = None
    for cls in ("RAS701.HECRASController", "RAS61.HECRASController",
                "RAS.HECRASController"):
        try:
            rc = win32.Dispatch(cls)
            break
        except Exception:
            continue
    if rc is None:
        log("HECRASController nao registrado. O HEC-RAS 7.0.1 esta instalado?", "erro")
        return False
    t0 = time.time()
    rc.Project_Open(prj)
    try:
        rc.Compute_CurrentPlan(None, None, True)
    except Exception as e:
        log(f"excecao no compute: {e}", "erro")
    for fechar in ("Project_Close",):
        try:
            getattr(rc, fechar)()
        except Exception:
            pass
    for sair in ("QuitRas", "QuitRAS"):
        try:
            getattr(rc, sair)()
            break
        except Exception:
            pass
    ok, sol, n = status_hdf(projeto)
    log(f"{sol}  ({n} passos, {time.time()-t0:.0f}s)", "ok" if ok else "erro")
    return bool(ok)


def validar(projeto, evento):
    """Compara o pico simulado com o observado, quando houver referencia."""
    ref = OBSERVADO.get(str(evento))
    if not ref:
        return
    try:
        import h5py
        import numpy as np
    except ImportError:
        return
    hdf = os.path.join(AQUI, f"{projeto}.p01.hdf")
    with h5py.File(hdf, "r") as f:
        g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
              "Unsteady Time Series/Cross Sections"]
        q = g["Flow"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = [x["River"].decode().strip() for x in at]
        rs = np.array([float(x["RS"].decode()) for x in at])
    acu = [k for k, x in enumerate(riv) if x == "Itajai_Acu"]
    if not acu:
        return
    print(f"\n  {'local':<12}{'simulado':>10}{'observado':>11}{'erro':>9}")
    for nome, km, obs in ref:
        i = min(acu, key=lambda k: abs(rs[k] - km * 1000))
        sim = float(q[:, i].max())
        print(f"  {nome:<12}{sim:>10.0f}{obs:>11}{100*(sim-obs)/obs:>8.1f}%")


def kml(projeto):
    r = subprocess.run([sys.executable, os.path.join(AQUI, "exportar_kml.py"),
                        projeto], cwd=AQUI)
    return r.returncode == 0


def secoes_app(projeto):
    r = subprocess.run([sys.executable, os.path.join(AQUI, "gerar_secoes_app.py"),
                        projeto], cwd=AQUI)
    return r.returncode == 0


def mancha(projeto):
    r = subprocess.run([sys.executable, os.path.join(AQUI, "gerar_mancha_hecras.py"),
                        projeto], cwd=AQUI)
    return r.returncode == 0


# ------------------------------------------------------------------ pipeline
def rodar(evento, barragens=True, so_mancha=False):
    titulo(f"BACIA DO ITAJAI  |  {'evento ' + evento if evento else 'cheia sintetica'}"
           f"{'' if barragens else '  |  SEM barragens'}")
    os.chdir(AQUI)
    if so_mancha:
        projeto = f"Itajai_Rede_{evento}" if evento else "Itajai_Rede"
    else:
        log("[1/5] gerando geometria e condicoes de contorno...")
        projeto = gerar(evento, barragens)
        log("[2/5] simulando no HEC-RAS...")
        if not simular(projeto):
            log("simulacao nao concluiu; a mancha nao sera gerada", "erro")
            return False
        validar(projeto, evento)
    log("[3/5] gerando a mancha de inundacao...")
    if not mancha(projeto):
        log("falha ao gerar a mancha", "erro")
        return False
    log("[4/5] exportando as secoes para o app...")
    secoes_app(projeto)
    log("[5/5] exportando KMZ para o Google Earth...")
    kml(projeto)
    log(f"concluido: {projeto}", "ok")
    return True


def main():
    args = [a for a in sys.argv[1:]]
    barragens = "--sem-barragens" not in args
    so_mancha = "--so-mancha" in args
    todos = "--todos" in args
    eventos = [a for a in args if not a.startswith("--")]

    if todos:
        alvos = EVENTOS_PADRAO
    elif eventos:
        alvos = eventos
    else:
        alvos = [None]

    falhas = []
    for ev in alvos:
        if not rodar(ev, barragens, so_mancha):
            falhas.append(ev or "sintetico")

    titulo("RESUMO")
    for ev in alvos:
        nome = ev or "sintetico"
        log(f"{nome:<18}{'FALHOU' if nome in falhas else 'ok'}",
            "erro" if nome in falhas else "ok")
    if falhas:
        log(f"{len(falhas)} de {len(alvos)} falharam", "erro")
        return 1
    print("\n  Para ver no navegador:")
    print("    python -m http.server 8050 --directory app")
    print("    http://localhost:8050/index.html  ->  aba Mapa  ->  camada 'Mancha HEC-RAS'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
