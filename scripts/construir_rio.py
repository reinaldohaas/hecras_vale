# -*- coding: utf-8 -*-
"""UM comando: do relevo ate um modelo HEC-RAS validado, para qualquer rio.

    python scripts/construir_rio.py Itajai_Acu
    python scripts/construir_rio.py --todos

Faz TUDO, em ordem, e a validacao decide se as etapas anteriores valem:

    por rio
      1. geometria do MDT SIG-SC 1 m          rio_do_relevo.py
      2. projeto com projecao e contorno      projeto_rio_avulso.py
      3. VALIDACAO SEM RODAR O SOLVER         ler_erros_geometria.py
      4. pedido de batimetria                 batimetria.py
    uma vez, no fim
      5. terreno sobre a UNIAO dos rios       terreno_30m.py
      6. limpeza do vazio negativo            limpar_vazio_negativo.py
      7. religa o terreno em cada projeto     projeto_rio_avulso.py
      8. confere a edge line NO HDF DO RAS    conferir_edge_lines.py

O TERRENO E DA BACIA, E NAO DE UM RIO

  Havia um terreno de 30 m e ele cobria 100% do Mirim, 60% do Acu, 5% do
  Benedito e ZERO do Norte, do Sul e do Oeste -- fora feito sobre o dominio de
  um rio so. Quem abrisse qualquer um dos outros no RAS Mapper via o modelo
  sem relevo nenhum. Aqui ele sai da uniao dos seis, de uma vez, e entra no
  `.rasmap` de todos. Vem DEPOIS das geometrias porque o dominio vem delas.

  `--sem-terreno` pula as etapas 5 a 7, que sao as caras: 765 folhas do
  SIG-SC a 1 m sobre 150 x 137 km. Sem elas o modelo abre sem relevo.

POR QUE A CONFERENCIA DA EDGE LINE E SEPARADA DA VALIDACAO

  O "Validate Geometry" NAO conta a auto-interseccao das edge lines na sua
  lista: o Oeste marcava zero mensagens e o RAS Mapper avisava assim mesmo.
  A etapa 8 le `/Geometry/River Edge Lines` do HDF -- o traco que o RAS usa
  para montar a superficie de interpolacao -- e mede nele.

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


def terreno(pastas, nome="vale30"):
    """Terreno unico sobre a uniao das geometrias, limpo do vazio negativo.

    Devolve o caminho do `.hdf`, ou None se nao deu.
    """
    print(f"\n{'='*68}\nTERRENO DA BACIA\n{'='*68}")
    geoms = [os.path.join(p, os.path.basename(p) + ".g01") for p in pastas]
    geoms = [g for g in geoms if os.path.exists(g)]
    if not geoms:
        print("   nenhuma geometria -- nada a fazer")
        return None
    roda(["scripts/terreno_30m.py"] + geoms + ["--nome", nome],
         mostrar=("dominio   :", "fonte     :", "raster:", "cobertura",
                  "celulas entre", "OK    ", "FALTA"))
    raiz = os.path.dirname(os.path.dirname(os.path.abspath(geoms[0])))
    pt = os.path.join(raiz, "Terrain")
    tif = os.path.join(pt, "MDT_SIGSC_30m.tif")
    vrt = os.path.join(pt, nome + "_sigsc_1m.vrt")
    # A LIMPEZA NAO E OPCIONAL. Reduzir 1 m -> 30 m e uma media de 900 pixels,
    # e as folhas que gravam o vazio como numero negativo grande contaminam a
    # media: medidas 879 celulas assim na bacia, e 471 delas com media
    # POSITIVA -- plausiveis, erradas, e invisiveis a olho.
    if os.path.exists(tif) and os.path.exists(vrt):
        roda(["scripts/limpar_vazio_negativo.py", tif, "--vrt", vrt,
              "--nome", nome],
             mostrar=("com minimo negativo", "media POSITIVA", "raster limpo:",
                      "cota ", "OK    ", "FALTA"))
    h = os.path.join(pt, nome + "_Terreno.hdf")
    return h if os.path.exists(h) else None


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
    ap.add_argument("--sem-terreno", action="store_true",
                    help="pula as etapas 5 a 7; o modelo abre sem relevo")
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

    # ---- 5 a 7: terreno da bacia, depois das geometrias
    if res and not a.sem_terreno:
        thdf = terreno([r["pasta"] for r in res])
        if thdf:
            print("\n   religando o terreno em cada projeto")
            for r in res:
                g = os.path.join(r["pasta"],
                                 os.path.basename(r["pasta"]) + ".g01")
                roda(["scripts/projeto_rio_avulso.py", g, "--rio-fonte",
                      r["rio"], "--terreno", thdf], mostrar=("terreno   :",))

    # ---- 8: a edge line, medida NA QUE O RAS CONSTRUIU
    for r in res:
        h = os.path.join(r["pasta"],
                         os.path.basename(r["pasta"]) + ".g01.hdf")
        r["dobras"] = None
        if os.path.exists(h):
            m = re.search(r"TOTAL: (\d+)",
                          roda(["scripts/conferir_edge_lines.py", h]))
            if m:
                r["dobras"] = int(m.group(1))

    print(f"\n{'='*68}\nRESUMO\n{'='*68}")
    print(f"{'rio':<16}{'pasta':<26}{'erros':>7}{'Fatal':>7}{'taxa':>8}"
          f"{'dobras':>8}")
    for r in res:
        d = r.get("dobras")
        print(f"{r['rio']:<16}{r['pasta']:<26}{r['erros']:>7}{r['fatal']:>7}"
              f"{r['taxa']:>8.2f}{('?' if d is None else str(d)):>8}")
    ruins = [r for r in res if r["erros"] > a.limite]
    if ruins:
        print(f"\nacima do limite de {a.limite}: "
              + ", ".join(r["rio"] for r in ruins))
        print("   nesses a varzea e plana e a secao vai ao teto; apertar mais")
        print("   a taxa reduz o erro e tambem a capacidade de conter a cheia")
    return res


if __name__ == "__main__":
    main()
