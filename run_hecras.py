"""
Roda a simulacao HEC-RAS e mostra um resumo (nivel/vazao maximos por trecho).

Uso:
    python run_hecras.py                 # auto: modelo real, senao o idealizado
    python run_hecras.py Itajai_Bacia_Real
    python run_hecras.py Itajai_Bacia_Completa

Backends (tenta nesta ordem):
    1) ras-commander  (pip install ras-commander ; Python >= 3.11, roda por linha
       de comando, NAO precisa de pywin32)
    2) COM / pywin32  (pip install pywin32 ; usa o RAS701.HECRASController)

Resumo dos resultados requer h5py (pip install h5py).
"""
import os
import sys

HERE = os.path.abspath(os.path.dirname(__file__))
RAS_VERSION = "7.0.1"
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"


def pick_project(arg):
    if arg:
        base = arg[:-4] if arg.lower().endswith(".prj") else arg
    else:
        base = ("Itajai_Bacia_Real"
                if os.path.exists(os.path.join(HERE, "Itajai_Bacia_Real.prj"))
                else "Itajai_Bacia_Completa")
    prj = os.path.join(HERE, base + ".prj")
    if not os.path.exists(prj):
        sys.exit(f"Projeto nao encontrado: {prj}")
    return base, prj


def run_rascommander(prj):
    from ras_commander import init_ras_project, RasCmdr
    ver = RAS_VERSION if not os.path.exists(RAS_EXE) else RAS_EXE
    print(f"[ras-commander] init_ras_project({os.path.basename(prj)}, {ver})")
    init_ras_project(prj, ver)
    res = RasCmdr.compute_plan("01", overwrite_dest=False)
    ok = res[0] if isinstance(res, tuple) else bool(res)
    return ok


def run_com(prj):
    import win32com.client as win32
    rc = None
    for cls in ("RAS701.HECRASController", "RAS61.HECRASController",
                "RAS.HECRASController"):
        try:
            rc = win32.Dispatch(cls); print(f"[COM] {cls}"); break
        except Exception:
            continue
    if rc is None:
        raise RuntimeError("Sem HECRASController COM (pip install pywin32).")
    rc.Project_Open(prj)
    res = rc.Compute_CurrentPlan(None, None, True)
    ok = res[0] if isinstance(res, tuple) else bool(res)
    try: rc.Project_Close()
    except Exception: pass
    for q in ("QuitRas", "QuitRAS"):
        try: getattr(rc, q)(); break
        except Exception: pass
    return ok


def summarize(base):
    try:
        import h5py
    except ImportError:
        print("(pip install h5py para ver o resumo dos resultados)"); return
    hdf = os.path.join(HERE, base + ".p01.hdf")
    if not os.path.exists(hdf):
        print("(sem arquivo de resultados .p01.hdf)"); return
    with h5py.File(hdf, "r") as f:
        g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
              "Unsteady Time Series/Cross Sections"]
        ws, q = g["Water Surface"][:], g["Flow"][:]
        attr = f["Geometry/Cross Sections/Attributes"][:]
        reach = [r["Reach"].decode().strip() for r in attr]
        river = [r["River"].decode().strip() for r in attr]
        key = [f"{rv} / {rc}" for rv, rc in zip(river, reach)]
        print(f"\n{ws.shape[0]} passos x {ws.shape[1]} secoes\n")
        print(f"{'Rio / Trecho':<34}{'Q max (m3/s)':>13}{'Nivel max (m)':>15}")
        for k in dict.fromkeys(key):
            idx = [i for i, kk in enumerate(key) if kk == k]
            print(f"{k:<34}{q[:, idx].max():>13.1f}{ws[:, idx].max():>15.2f}")
        print(f"\nVazao maxima global: {q.max():.1f} m3/s")


def main():
    base, prj = pick_project(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"Projeto: {prj}")
    ok = None
    for backend in (run_rascommander, run_com):
        try:
            ok = backend(prj)
            break
        except ImportError as e:
            print(f"  ({backend.__name__} indisponivel: {e})")
        except Exception as e:
            print(f"  ({backend.__name__} falhou: {e})")
    if ok is None:
        sys.exit("Nenhum backend disponivel. Instale: pip install ras-commander  (ou pywin32)")
    if not ok:
        print("\n*** Simulacao NAO concluiu com sucesso. "
              "Veja o arquivo .p01.computeMsgs.txt ***"); return
    print("\n=== Simulacao concluida com SUCESSO ===")
    summarize(base)


if __name__ == "__main__":
    main()
