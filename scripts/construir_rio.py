# -*- coding: utf-8 -*-
"""UM comando: do relevo ate um modelo HEC-RAS validado, para qualquer rio.

    python scripts/construir_rio.py Itajai_Acu
    python scripts/construir_rio.py --todos

Faz as quatro etapas em ordem, e a quarta decide se as tres anteriores valem:

    1. geometria do MDT SIG-SC 1 m            rio_do_relevo.py
    2. projeto com projecao e contorno        projeto_rio_avulso.py
    3. VALIDACAO SEM RODAR O SOLVER           ler_erros_geometria.py
    4. pedido de batimetria                   batimetria.py

E SE AINDA HOUVER ERRO, ELE APERTA SOZINHO. O que controla a geometria e a
TAXA de variacao da largura da secao, nao a largura: a edge line liga as
pontas das secoes, e se a meia-largura anda D metros para o lado enquanto o
rio anda dx para a frente, ela faz angulo atan(D/dx) e a partir de certo ponto
dobra sobre a vizinha. Medido nos seis rios da bacia:

      Norte  12 m por secao (0,08 m/m)  ->    2 erros
      Acu    32 m por secao (0,21)      ->  190 erros
      Oeste  88 m por secao (0,59)      ->  303 erros

Entao o laco comeca em `--taxa` e vai apertando enquanto o validador do
proprio HEC-RAS acusar mais que `--limite` mensagens. Apertar SO ENCOLHE a
secao -- nenhuma fica mais larga do que o terreno mediu --, e o preco esta no
relatorio: secao menor contem menos cheia.

POR QUE NAO RODA O SOLVER

  Validar geometria pelo preprocessador geometrico custa cerca de um minuto;
  descobrir o mesmo rodando a simulacao custa dez ou mais. Rodar so faz
  sentido depois que a geometria esta boa, e essa e uma etapa separada.

CADA RIO EM UMA PASTA SO, com o nome do rio. Nada de sufixos de tentativa.
"""
import argparse
import os
import re
import subprocess
import sys

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
PY = os.path.join(RAIZ, ".venv", "Scripts", "python.exe")
if not os.path.exists(PY):
    PY = sys.executable

RIOS = ["Itajai_Acu", "Itajai_Mirim", "Itajai_Norte", "Itajai_Sul",
        "Itajai_Oeste", "Rio_Benedito"]
TAXAS = [0.15, 0.10, 0.07, 0.05]     # do mais folgado ao mais apertado


def roda(args, mostrar=()):
    p = subprocess.run([PY] + args, cwd=RAIZ, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    saida = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0:
        print(saida[-1500:])
        raise SystemExit(f"falhou: {' '.join(args)}")
    for l in saida.split("\n"):
        if any(k in l for k in mostrar):
            print("   " + l.strip())
    return saida


def erros_de(saida):
    m = re.search(r"colecoes de erro:\s*(\d+)\s+mensagens:\s*(\d+)", saida)
    if not m:
        return None, None, saida
    fat = re.search(r"Fatal (\d+)", saida)
    return int(m.group(2)), int(fat.group(1)) if fat else 0, saida


def construir(rio, pasta, limite, taxas, dx, cada):
    print(f"\n{'='*68}\n{rio}\n{'='*68}")
    g = os.path.join(pasta, os.path.basename(pasta) + ".g01")
    melhor = None
    for taxa in taxas:
        print(f"\n-- taxa de largura {taxa:g} m/m")
        roda(["scripts/rio_do_relevo.py", "--rio", rio, "--saida", pasta,
              "--dx", str(dx), "--taxa", str(taxa), "--monotono"],
             mostrar=("secoes :", "descartadas", "calha ", "secao ",
                      "talvegue:", "edge line:"))
        roda(["scripts/projeto_rio_avulso.py", g, "--rio-fonte", rio],
             mostrar=("hidrograma:", "jusante   :"))
        hdf = g + ".hdf"
        if os.path.exists(hdf):
            os.remove(hdf)
        n, fat, saida = erros_de(roda(["scripts/ler_erros_geometria.py", g]))
        if n is None:
            print("   nao consegui ler o validador"); break
        print(f"   VALIDADOR: {n} mensagens, {fat} Fatal")
        if melhor is None or n < melhor[0]:
            melhor = (n, fat, taxa)
        if n <= limite:
            break
    if melhor is None:
        return None
    n, fat, taxa = melhor
    if taxa != taxas[-1] and n > limite:
        pass
    print(f"\n   melhor: {n} mensagens ({fat} Fatal) com taxa {taxa:g} m/m")
    ped = os.path.join("doc", f"batimetria_{os.path.basename(pasta)}.csv")
    roda(["scripts/batimetria.py", "pedir", g, "--cada", str(cada),
          "--saida", ped], mostrar=("pontos    :",))
    return {"rio": rio, "pasta": pasta, "erros": n, "fatal": fat,
            "taxa": taxa, "pedido": ped}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rio", nargs="?")
    ap.add_argument("--todos", action="store_true")
    ap.add_argument("--saida", default=None)
    ap.add_argument("--limite", type=int, default=6,
                    help="mensagens do validador aceitas antes de apertar")
    ap.add_argument("--taxa", type=float, default=None,
                    help="fixa a taxa em vez de deixar o laco procurar")
    ap.add_argument("--dx", type=float, default=150.0)
    ap.add_argument("--cada", type=float, default=2000.0)
    a = ap.parse_args()
    if not a.rio and not a.todos:
        raise SystemExit("informe um rio ou --todos.  Rios: " + ", ".join(RIOS))
    taxas = [a.taxa] if a.taxa else TAXAS
    alvos = RIOS if a.todos else [a.rio]

    res = []
    for rio in alvos:
        pasta = a.saida or os.path.join("modelo", rio.lower())
        r = construir(rio, pasta, a.limite, taxas, a.dx, a.cada)
        if r:
            res.append(r)

    print(f"\n{'='*68}\nRESUMO\n{'='*68}")
    print(f"{'rio':<16}{'pasta':<26}{'erros':>7}{'Fatal':>7}{'taxa':>8}")
    for r in res:
        print(f"{r['rio']:<16}{r['pasta']:<26}{r['erros']:>7}{r['fatal']:>7}"
              f"{r['taxa']:>8.2f}")
    ruins = [r for r in res if r["erros"] > a.limite]
    if ruins:
        print(f"\nacima do limite de {a.limite}: "
              + ", ".join(r["rio"] for r in ruins))
        print("   nesses a varzea e plana e a secao vai ao teto; apertar mais")
        print("   a taxa reduz o erro e tambem a capacidade de conter a cheia")
    return res


if __name__ == "__main__":
    main()
