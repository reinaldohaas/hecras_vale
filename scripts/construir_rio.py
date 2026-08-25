# -*- coding: utf-8 -*-
"""UM comando: do relevo ate um modelo HEC-RAS validado, para qualquer rio.

    python scripts/construir_rio.py Itajai_Acu
    python scripts/construir_rio.py --todos

Faz TUDO, em ordem, e a validacao decide se as etapas anteriores valem:

    por rio
      1. geometria do MDT SIG-SC 1 m          rio_do_relevo.py
      2. projeto com projecao e contorno      projeto_rio_avulso.py
      3. VALIDACAO SEM RODAR O SOLVER         ler_erros_geometria.py
      4. pedido de batimetria                 batimetria.py pedir
      5. BATIMETRIA DO LEGADO -> g02          batimetria_do_legado.py +
         com PORTEIRO DE EIXO ALTO            batimetria.py aplicar
    uma vez, no fim
      6. terreno sobre a UNIAO dos rios       terreno_30m.py
      7. limpeza do vazio negativo            limpar_vazio_negativo.py
      8. religa o terreno em cada projeto     projeto_rio_avulso.py
      9. confere a edge line NO HDF DO RAS    conferir_edge_lines.py

O PORTEIRO DE EIXO ALTO (passo 5)

  Antes de ancorar a batimetria, o pipeline compara o talvegue lido do MDT com
  o fundo levantado de 1983. Onde o eixo esquematico corre pela ENCOSTA de um
  vale encaixado, a lamina do MDT fica dezenas a centenas de metros acima do
  fundo, e ancorar ali cavaria um canion -- foi o que deixou o Benedito
  inviavel (rebaixamento mediano 104 m) e mandou o solver do Acu a
  instabilidade (55 m no R1). Rio com trecho assim NAO ganha g02: o pipeline
  imprime a faixa em km, gera a figura (diagnostico_eixo_alto.py) e marca na
  tabela final "REFAZER EIXO". O g02 velho, se existir, e removido -- rodar o
  pipeline duas vezes tem de dar o mesmo estado.

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


LIMIAR_EIXO = 25.0     # m; rebaixamento acima disto = eixo pela encosta
TRECHO_EIXO = 1.0      # km; menos que isto e blip de uma secao, nao eixo


def eixo_alto(g01, rio, limiar=LIMIAR_EIXO, trecho_km=TRECHO_EIXO):
    """(km_ini, km_fim, reb_max) do trecho onde ancorar cavaria um canion.

    Compara o talvegue do MDT (lamina) com o fundo levantado do legado, ao
    longo do rio. Rebaixamento maior que `limiar` nao e batimetria: a calha
    real destes rios tem ~11 m, e um pedido de 25+ m diz que o eixo corre
    pela encosta e a secao pegou o vale errado.

    O criterio e o TRECHO CONTIGUO, nao o ponto. Medido nos seis rios: os que
    rodam bem tem no maximo blips isolados de 0,0-0,2 km um pouco acima do
    limiar (Mirim 25,1 m numa secao; Norte 26,8 m em 200 m), enquanto os
    quebrados tem 20,6 km (Acu R1, ate 117 m) e 24,6 km (Benedito, ate 284 m).
    Trecho menor que `trecho_km` nao barra. Devolve None se o rio e sadio.
    """
    import numpy as np
    from qc_secoes import ler_secoes
    from batimetria_do_legado import secoes_levantadas, LEGADO
    S = ler_secoes(g01)
    S.sort(key=lambda d: -d["rs"])
    rs = np.array([d["rs"] for d in S])
    z = np.array([float(np.asarray(d["z"], float).min()) for d in S])
    ch = np.array([float(d["len_ch"]) for d in S])
    x = np.r_[0.0, np.cumsum(ch[:-1])]
    L = secoes_levantadas(LEGADO, rio)
    if L is None:
        return None
    o = np.argsort(-L[:, 0])
    fundo = np.interp(x, np.interp(-L[o, 0], -rs, x), L[o, 3])
    reb = z - fundo
    m = reb > limiar
    if not m.any():
        return None
    # trechos contiguos de m; so conta o que tiver extensao >= trecho_km
    corta = np.flatnonzero(np.diff(m.astype(int)))
    ini = np.r_[0, corta + 1]
    fim = np.r_[corta, len(m) - 1]
    ruins = [(x[i], x[f]) for i, f in zip(ini, fim)
             if m[i] and (x[f] - x[i]) / 1000 >= trecho_km]
    if not ruins:
        return None
    x0 = min(r[0] for r in ruins)
    x1 = max(r[1] for r in ruins)
    dentro = m & (x >= x0) & (x <= x1)
    return float(x0 / 1000), float(x1 / 1000), float(reb[dentro].max())


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

    # ---- 5. batimetria do legado -> g02, atras do porteiro de eixo alto
    roda(["scripts/batimetria_do_legado.py", ped, "--rio", rio],
         mostrar=("casados", "REBAIXAMENTO", "ATENCAO"))
    g2 = os.path.join(pasta, os.path.basename(pasta) + ".g02")
    alto = eixo_alto(g, rio)
    if alto is None:
        roda(["scripts/batimetria.py", "aplicar", g, "--pontos", ped,
              "--saida", "g02"],
             mostrar=("contradeclives", "declividade", "leito bate"))
        roda(["scripts/projeto_rio_avulso.py", g2, "--rio-fonte", rio],
             mostrar=("jusante   :",))
        print(f"   batimetria aplicada -> {g2}")
    else:
        km0, km1, reb = alto
        fig = os.path.join("doc", "figuras",
                           f"eixo_alto_{os.path.basename(pasta)}.png")
        roda(["scripts/diagnostico_eixo_alto.py", "--rios", rio,
              "--saida", fig])
        # determinismo: sem g02 valido, nao pode sobrar um g02 velho no lugar
        if os.path.exists(g2):
            os.remove(g2)
        print(f"   EIXO ALTO de {km0:.0f} a {km1:.0f} km (rebaixamento ate "
              f"{reb:.0f} m): g02 NAO aplicado -- refazer o eixo pelo "
              f"talvegue do MDT. Figura: {fig}")
    return {"rio": rio, "pasta": pasta, "erros": n, "fatal": fat,
            "taxa": taxa, "pedido": ped, "eixo_alto": alto,
            "g02": alto is None}


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

    # ---- 6 a 8: terreno da bacia, depois das geometrias
    if res and not a.sem_terreno:
        thdf = terreno([r["pasta"] for r in res])
        if thdf:
            print("\n   religando o terreno em cada projeto")
            for r in res:
                # A GEOMETRIA EM USO manda: religar com o g01 quando o rio tem
                # g02 reescreveria `Geom File=g01` no projeto, e o plano
                # voltaria a rodar SEM batimetria -- o mesmo defeito, calado,
                # que ja custou uma rodada inteira (item 3 do RETOMAR).
                base = os.path.join(r["pasta"], os.path.basename(r["pasta"]))
                g = base + (".g02" if os.path.exists(base + ".g02")
                            else ".g01")
                roda(["scripts/projeto_rio_avulso.py", g, "--rio-fonte",
                      r["rio"], "--terreno", thdf], mostrar=("terreno   :",))

    # ---- 9: a edge line, medida NA QUE O RAS CONSTRUIU
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
    print(f"{'rio':<16}{'erros':>7}{'Fatal':>7}{'taxa':>7}{'dobras':>8}"
          f"   batimetria")
    for r in res:
        d = r.get("dobras")
        if r.get("g02"):
            bat = "g02 aplicado"
        elif r.get("eixo_alto"):
            k0, k1, reb = r["eixo_alto"]
            bat = (f"REFAZER EIXO ({k0:.0f}-{k1:.0f} km, "
                   f"reb ate {reb:.0f} m)")
        else:
            bat = "?"
        print(f"{r['rio']:<16}{r['erros']:>7}{r['fatal']:>7}"
              f"{r['taxa']:>7.2f}{('?' if d is None else str(d)):>8}   {bat}")
    ruins = [r for r in res if r["erros"] > a.limite]
    if ruins:
        print(f"\nacima do limite de {a.limite}: "
              + ", ".join(r["rio"] for r in ruins))
        print("   nesses a varzea e plana e a secao vai ao teto; apertar mais")
        print("   a taxa reduz o erro e tambem a capacidade de conter a cheia")
    return res


if __name__ == "__main__":
    main()
