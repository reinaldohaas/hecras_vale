# -*- coding: utf-8 -*-
"""Constroi DO ZERO o modelo calibrado de julho/1983 do taha_ai raiz.

    python scripts/construir_1983.py            # constroi g01/u01/p01
    python scripts/construir_1983.py --rodar    # e roda + confere reguas
    python scripts/construir_1983.py --ate 6    # para apos a etapa 6

E a receita INTEIRA da campanha de 26-27/08/2026, executavel: parte de
`taha_ai.g01.antes_do_reparo_1983` + `taha_ai.u01.antes_do_observado` +
os dados publicos ja no repo (doc/ana_1983, doc/larguras_sigsc) e chega
ao estado que fechou os 31 dias com 0,0008% de erro de volume
(`taha_ai.g01.estado_mes_completo`).

ETAPAS (cada uma e um programa proprio, com conferencia propria):

  g01 -----------------------------------------------------------------
   1 rebaixar_foz --minimo 0.1        fozes penduradas descem ao receptor
   2 afogar_soleiras                  soleiras adversas afogadas
   3 nivelar_bancadas                 canal na encosta desce ao talvegue
   4 encolher_canal --fator 2 (SEM FBDS, poupar Mirim<10km)
                                      largura real: lamina SIG-SC x2
   5 interpolar @2m: Acu R1 139-167 / Norte R2 5-30 / Mirim 80-105
                                      secoes curtas nas serras
   6 amputar Mirim>40k, Benedito>23k, Cedros>12k
                                      contorno entra NA regua da ANA
   7 construir_barragem Sul RS 33500  vertedouro 399, fenda 1,3 (JICA)
   8 interpolar pe da barragem @0.25 (Sul 29-33,5) e @0.15 (Sul 31-33,5)
   9 interpolar canion @1m (Acu R1 139-167, respeita a estrutura)
  10 esticar_htab 60 m               tabela nao pode ser extrapolada

  u01 -----------------------------------------------------------------
   a contorno_observado               observado ANA por sistema, fracao
                                      de cabeceira por AREA, mare M2
   b converter_inicial                Initial RS= fantasma -> Flow Loc
   (as amputacoes da etapa 6 movem os contornos junto)
   c dividir_lateral Sul em 33500 (46% no pool = bacia real 1150 km2)

  p01 -----------------------------------------------------------------
   d Simulation Date 01-31JUL1983 (contorno faz), Interval=1MIN,
     ZTol/ZSATol 0.02, MxIter 40, LPI 0.8/4 -- conferidos, CRLF sempre

ACEITE (--rodar): "Finished Unsteady Flow Simulation", erro de volume
< 0,01%, e os picos dentro das faixas da campanha (Blumenau -9+-4%,
Taio/Benedito/Arrozeira/Brusque/Barra +-5%).
"""
import os
import shutil
import subprocess
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
SCRIPTS = os.path.join(RAIZ, "scripts")

# sem Library\bin no PATH o numpy do miniforge crasha CALADO (exit 127);
# arma aqui para este processo e para todos os filhos
_LIB = os.path.join(os.path.dirname(PY), "Library", "bin")
if os.path.isdir(_LIB) and _LIB not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _LIB + os.pathsep + os.environ.get("PATH", "")


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def rodar(cmd, resumo):
    print(f"\n=== {resumo}")
    print("    " + " ".join(os.path.basename(c) for c in cmd[:3])
          + " " + " ".join(cmd[3:]))
    r = subprocess.run(cmd, cwd=RAIZ, capture_output=True, text=True)
    cauda = (r.stdout + r.stderr).strip().split("\n")
    for l in cauda[-4:]:
        print("    | " + l)
    if r.returncode != 0:
        raise SystemExit(f"ETAPA FALHOU ({resumo}) -- veja acima")


def prog(nome, *args):
    return [PY, os.path.join(SCRIPTS, nome)] + list(args)


def g(ext):
    return os.path.join(RAIZ, f"taha_ai.{ext}")


def promover(ext):
    # etapa sem nada a fazer nao grava saida (ex.: "nenhuma soleira
    # adversa") -- o g01 corrente segue valendo
    if not os.path.exists(g(ext)):
        print("    | (sem saida: etapa nao tinha nada a fazer)")
        return
    shutil.copy2(g(ext), g("g01"))
    os.remove(g(ext))


def main(argv):
    ate = _arg(argv, "--ate", 99, int)
    quer_rodar = "--rodar" in argv
    # --sem-amputar: mantem os rios INTEIROS (no relevo SIG-SC a serra e
    # real; a amputacao era remedio para a serra sintetica)
    sem_amputar = "--sem-amputar" in argv
    # --slot N: grava a variante como gNN/uNN/pNN registrados no .prj,
    # SEM tocar g01/u01/p01 ([[nunca-sobrescrever-o-projeto]])
    slot = _arg(argv, "--slot", 0, int)
    # --base troca a geometria crua de partida (ex.: taha_ai.r00, o
    # relevo do MDT 1 m gerado por relevo_nas_secoes.py)
    base_g = _arg(argv, "--base", "taha_ai.g01.antes_do_reparo_1983")
    guardados = {}
    if slot:
        for ext in ("g01", "u01", "p01"):
            arq = g(ext)
            if os.path.exists(arq):
                guardados[ext] = arq + ".antes_do_slot"
                shutil.copy2(arq, guardados[ext])

    for arq in (base_g, "taha_ai.u01.antes_do_observado"):
        if not os.path.exists(os.path.join(RAIZ, arq)):
            raise SystemExit(f"falta o insumo {arq}")
    print(f"base de geometria: {base_g}")

    # ------------------------------------------------------------- u01/p01
    shutil.copy2(os.path.join(RAIZ, "taha_ai.u01.antes_do_observado"),
                 g("u01"))
    rodar(prog("contorno_observado.py", "taha_ai",
               "--series", "doc/ana_1983"),
          "a. contorno observado de 1983 (fracao de cabeceira por area)")
    rodar(prog("converter_inicial.py", "taha_ai.u01"),
          "b. Initial RS= fantasma -> Initial Flow Loc=")

    # ------------------------------------------------------------- g01
    shutil.copy2(os.path.join(RAIZ, base_g), g("g01"))
    passos = [
        (1, prog("rebaixar_foz.py", "taha_ai.g01", "--saida", "z01",
                 "--minimo", "0.1"), "z01",
         "1. fozes penduradas rebaixadas"),
        (2, prog("afogar_soleiras.py", "taha_ai.g01", "--saida", "z02"),
         "z02", "2. soleiras adversas afogadas"),
        (3, prog("nivelar_bancadas.py", "taha_ai.g01", "--saida", "z03"),
         "z03", "3. bancadas de encosta niveladas"),
        # medidas CONGELADAS (doc/larguras_g95, do commit 32eb1e9): a
        # receita e sensivel ao dado -- com os CSVs pos-teto o modelo
        # reconstruido morre em 13,8 h; com estes fecha o mes
        (4, prog("encolher_canal.py", "taha_ai.g01", "--saida", "z04",
                 "--fator", "1.3", "--fbds", "NAO",
                 "--medidas", "doc/larguras_g95",
                 "--poupar", "Itajai_Mirim:10"), "z04",
         "4. canais encolhidos a lamina SIG-SC x1,3 (insumo congelado)"),
        (5, prog("interpolar_secoes.py", "taha_ai.g01", "--saida", "z05",
                 "--queda-max", "2.0",
                 "--trecho", "Itajai_Acu,R1,139000,167000",
                 "--trecho", "Itajai_Norte,R2,5000,30000",
                 "--trecho", "Itajai_Mirim,R1,80000,105000"), "z05",
         "5. serras interpoladas a 2 m de queda"),
    ]
    for num, cmd, ext, resumo in passos:
        if num > ate:
            return
        rodar(cmd, resumo)
        promover(ext)

    if ate >= 6 and sem_amputar:
        print("\n=== 6. amputacoes PULADAS (--sem-amputar: rios inteiros)")
    if ate >= 6 and not sem_amputar:
        for reach, corte in [("Itajai_Mirim,R1", "40000"),
                             ("Rio_Benedito,R1", "23000"),
                             ("Rio_dos_Cedros,R1", "12000")]:
            rodar(prog("amputar_cabeceira.py", "taha_ai.g01",
                       "--reach", reach, "--rs-corte", corte,
                       "--saida", "z06"),
                  f"6. amputado {reach} > {corte} (contorno na regua)")
            promover("z06")
    if ate >= 7:
        rodar(prog("construir_barragem.py", "taha_ai.g01",
                   "--saida", "z07", "--rio", "Itajai_Sul",
                   "--reach", "R1", "--rs", "33500",
                   "--crista", "399.0", "--topo", "402.0",
                   "--larg-vertedouro", "100", "--fenda", "1.3",
                   "--nome", "Barragem Sul (Ituporanga)"),
              "7. Barragem Sul construida (JICA)")
        promover("z07")
    if ate >= 8:
        rodar(prog("interpolar_secoes.py", "taha_ai.g01", "--saida", "z08",
                   "--queda-max", "0.25",
                   "--trecho", "Itajai_Sul,R1,29000,33499"),
              "8a. pe da barragem a 0,25 m de queda")
        promover("z08")
    if ate >= 9:
        rodar(prog("interpolar_secoes.py", "taha_ai.g01", "--saida",
                   "z09", "--queda-max", "1.0",
                   "--trecho", "Itajai_Acu,R1,139000,167000"),
              "9. canion do Acu a 1 m de queda")
        promover("z09")
        rodar(prog("interpolar_secoes.py", "taha_ai.g01", "--saida",
                   "z10", "--queda-max", "0.15",
                   "--trecho", "Itajai_Sul,R1,31000,33499"),
              "8b. pe da barragem a 0,15 m de queda")
        promover("z10")
    if ate >= 10:
        rodar(prog("esticar_htab.py", "taha_ai.g01", "--saida", "z11",
                   "--alcance", "60"),
              "10. tabelas hidraulicas esticadas a 60 m")
        promover("z11")
        rodar(prog("dividir_lateral.py", "taha_ai.u01",
                   "--rio", "Itajai_Sul", "--reach", "R1",
                   "--em", "33500", "--fracao-acima", "0.46"),
              "c. lateral do Sul repartida na barragem")

    # p01: garantias
    p01 = g("p01")
    t = open(p01, encoding="latin-1").read().replace("\r\n", "\n")
    for chave, valor in [("Computation Interval=", "1MIN"),
                         ("UNET ZTol= ", "0.02 ")]:
        if chave + valor not in t:
            print(f"    AVISO p01: esperava '{chave}{valor}'")
    open(p01, "w", encoding="latin-1", newline="\r\n").write(t)
    print("\nd. p01 conferido (CRLF garantido)")

    plano = "01"
    if slot:
        # empacota a variante em gNN/uNN/pNN e RESTAURA g01/u01/p01
        gs, us, ps = f"g{slot:02d}", f"u{slot:02d}", f"p{slot:02d}"
        plano = f"{slot:02d}"
        shutil.copy2(g("g01"), g(gs))
        shutil.copy2(g("u01"), g(us))
        tg = open(g(gs), encoding="latin-1").read().replace("\r\n", "\n")
        tg = tg.replace("Geom Title=taha_ai - eixo do relevo Copernicus",
                        f"Geom Title=taha_ai - variante slot {slot}", 1)
        open(g(gs), "w", encoding="latin-1", newline="\r\n").write(tg)
        import re as _re
        tp = open(g("p01"), encoding="latin-1").read().replace("\r\n",
                                                               "\n")
        tp = tp.replace("Geom File=g01", f"Geom File={gs}", 1)
        tp = tp.replace("Flow File=u01", f"Flow File={us}", 1)
        tp = _re.sub(r"Plan Title=.*",
                     f"Plan Title=1983 variante slot {slot}", tp, count=1)
        tp = _re.sub(r"Short Identifier=.*",
                     f"Short Identifier=1983_slot{slot}", tp, count=1)
        open(g(ps), "w", encoding="latin-1", newline="\r\n").write(tp)
        for ext, backup in guardados.items():
            shutil.copy2(backup, g(ext))
        tj = open(g("prj"), encoding="latin-1").read().replace("\r\n",
                                                               "\n")
        if f"Geom File={gs}" not in tj:
            tj = tj.replace("Geom File=g01",
                            f"Geom File=g01\nGeom File={gs}", 1)
        if f"Plan File={ps}" not in tj:
            tj = tj.replace("Plan File=p01",
                            f"Plan File=p01\nPlan File={ps}", 1)
        open(g("prj"), "w", encoding="latin-1", newline="\r\n").write(tj)
        print(f"\nCONSTRUIDO no slot: {gs}/{us}/{ps} registrados; "
              f"g01/u01/p01 restaurados intactos.")
    else:
        print("\nCONSTRUIDO. g01/u01/p01 prontos para o plano 01.")

    if not quer_rodar:
        print("(use --rodar para computar e conferir as reguas)")
        return
    lixos = ["taha_ai.p01.data_errors.txt",
             f"taha_ai.p{plano}.data_errors.txt",
             f"taha_ai.g{slot:02d}.hdf" if slot else "taha_ai.g01.hdf"]
    for lixo in lixos:
        if os.path.exists(os.path.join(RAIZ, lixo)):
            os.remove(os.path.join(RAIZ, lixo))
    sys.path.insert(0, RAIZ)
    from ras_commander import init_ras_project, RasCmdr
    from vale.terreno import HECRAS_DIR
    p = init_ras_project(os.path.join(RAIZ, "taha_ai.prj"),
                         os.path.join(HECRAS_DIR, "Ras.exe"))
    print(f"\ncomputando o plano {plano} (31 dias, ~30-40 min)...")
    RasCmdr.compute_plan(plano, ras_object=p, force_rerun=True)
    import h5py
    with h5py.File(os.path.join(RAIZ, f"taha_ai.p{plano}.hdf"),
                   "r") as f:
        txt = bytes(f["Results/Summary/Compute Messages (text)"][()]) \
            .decode("utf-8", "replace")
    ok = "Finished Unsteady Flow Simulation" in txt
    vol = [l for l in txt.split("\n") if "percentage" in l]
    print("\nACEITE:")
    print(f"   terminou o mes: {ok}   (tem de ser True)")
    print(f"   {vol[-1].strip() if vol else 'sem linha de volume'}")
    rodar(prog("comparar_com_observado.py", f"taha_ai.p{plano}.hdf"),
          "reguas simulado x observado")


if __name__ == "__main__":
    main(sys.argv[1:])
