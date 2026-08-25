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

O LEGADO SINTETICO (passo 5)

  Nas cabeceiras do Acu (RS ~143-164 km) e do Benedito (RS 18-44 km) o "fundo
  levantado" do legado e uma RETA DESENHADA: declive exatamente 8,00 m/km com
  residuo rms de 1-2 MILIMETROS por dezenas de secoes -- nos dois rios, a
  mesma constante ("rede real ANA + relevo DEM", diz o proprio titulo).
  Ancorar nisso pediria cavar ate 284 m. O diagnostico anterior ("eixo pela
  encosta") estava ERRADO: procurado o caminho de menor cota num corredor de
  +-1500 m, nao ha vale mais baixo -- o eixo esta certo, o legado e que e
  ficcao ali. O detector no batimetria_do_legado.py descarta essas ancoras
  (reta local com rms < 5 cm), o aplicar nao interpola por cima do vao, e o
  MDT lidar fica valendo no trecho. `eixo_alto` segue medindo, agora so para
  RELATAR a faixa na tabela e na figura.

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


def roda(args, mostrar=(), aceitar_falha=False):
    p = subprocess.run([PY] + args, cwd=RAIZ, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    saida = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0 and not aceitar_falha:
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
    """(km_ini, km_fim, reb_max) do trecho onde o legado diverge do MDT.

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



def alvos_de_reparo(g, csv_erros):
    """RS a descartar, lidos DOS DEFEITOS MEDIDOS -- e nao de palpite.

    Duas fontes, as duas do proprio HEC-RAS:
      - as XS dos Fatal "XS intersects > N banklines" (RS no proprio erro);
      - para cada DOBRA da edge line lida do `.g01.hdf` (a linha densificada
        que o RAS constroi -- a licao do Mirim: medir nela, nunca no proxy),
        a secao MAIS ESTREITA entre as 3 mais proximas do ponto da dobra.
    """
    import csv as _csv
    import h5py
    import numpy as np
    from qc_secoes import ler_secoes
    from conferir_edge_lines import cruzamentos
    alvos = set()
    if os.path.exists(csv_erros):
        for r in _csv.DictReader(open(csv_erros, encoding="utf-8"),
                                 delimiter=";"):
            if "intersects >" in r.get("mensagem", ""):
                m = re.search(r"\(([\d.]+)\)", r.get("onde", ""))
                if m:
                    alvos.add(round(float(m.group(1)), 2))
    h = g + ".hdf"
    if os.path.exists(h):
        S = ler_secoes(g)
        C = np.array([np.asarray(d["cut"], float).mean(0) for d in S])
        W = np.array([float(d["sta"][-1] - d["sta"][0]) for d in S])
        RS = np.array([d["rs"] for d in S])
        with h5py.File(h, "r") as f:
            cam = "/Geometry/River Edge Lines"
            if cam + "/Polyline Info" in f:
                info = f[cam + "/Polyline Info"][:]
                pts = f[cam + "/Polyline Points"][:]
                for l in info:
                    P = pts[int(l[0]):int(l[0]) + int(l[1])]
                    for i, _j in cruzamentos(P):
                        pm = 0.5 * (P[i] + P[i + 1])
                        d2 = np.hypot(C[:, 0] - pm[0], C[:, 1] - pm[1])
                        viz = np.argsort(d2)[:3]
                        k = int(viz[int(np.argmin(W[viz]))])
                        alvos.add(round(float(RS[k]), 2))
    return sorted(alvos)


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
    print(f"\n   melhor: {n} mensagens ({fat} Fatal) com taxa {taxa:g} m/m")

    # ================= O CONTRATO DE ACEITE, do usuario, na letra =========
    #   1. validador com Fatal -> NAO gera nem aponta g02;
    #   2. as linhas sao medidas no HDF da geometria EM USO;
    #   3. depois de escrever g02, o g02.hdf velho e removido -- HDF stale
    #      nao vale como prova;
    #   4. GRAVES do qc_perfis REPROVAM, nao so imprimem;
    #   5. aprovado = 0 Fatal + 0 GRAVES + TOTAL 0 de linhas, no g01 E no g02.
    # Reprovado fica reprovado NA TABELA, com o motivo -- e o projeto volta a
    # apontar para o g01, para nao sobrar plano mirando geometria invalida.
    base = os.path.basename(pasta)
    g2 = os.path.join(pasta, base + ".g02")

    def _reprova(motivo):
        for fx in (g2, g2 + ".hdf"):
            if os.path.exists(fx):
                os.remove(fx)
        roda(["scripts/projeto_rio_avulso.py", g, "--rio-fonte", rio],
             mostrar=())
        print(f"   REPROVADO: {motivo}")
        return {"rio": rio, "pasta": pasta, "erros": n, "fatal": fat,
                "taxa": taxa, "status": "REPROVADO: " + motivo}

    def _portoes(geom, rotulo):
        """0 Fatal, 0 GRAVES e TOTAL 0 nas linhas -- ou o motivo da reprova."""
        n_, fat_, _ = erros_de(roda(["scripts/ler_erros_geometria.py", geom]))
        print(f"   [{rotulo}] validador: {n_} mensagens, {fat_} Fatal")
        if fat_:
            return f"{fat_} Fatal no validador ({rotulo})"
        qs = roda(["scripts/qc_perfis.py", geom], mostrar=("GRAVE",))
        mq = re.search(r"GRAVES (\d+)", qs)
        graves = int(mq.group(1)) if mq else -1
        if graves != 0:
            return f"{graves} GRAVES no qc_perfis ({rotulo})"
        cs = roda(["scripts/conferir_edge_lines.py", geom + ".hdf"],
                  mostrar=("bank line", "edge line"),
                  aceitar_falha=True)
        mt = re.search(r"TOTAL: (\d+)", cs)
        tot = int(mt.group(1)) if mt else -1
        if tot != 0:
            return f"{tot} defeito(s) de edge/bank line ({rotulo})"
        return None

    # ---- LACO DE REPARO: o que eu (o assistente) vinha fazendo na mao --
    # rodar, medir a dobra no HDF, escolher a secao, rodar de novo -- vira
    # software: reprova por linha => le os defeitos medidos => descarta as
    # participantes (CSV deterministico, como o rede_descartes) => reconstroi
    # e re-mede, ate limpar ou esgotar 4 voltas.
    rep_csv = os.path.join("doc", f"reparo_{base}.csv")
    if os.path.exists(rep_csv):
        os.remove(rep_csv)
    descartadas = []
    motivo = _portoes(g, "g01")
    volta = 0
    while motivo is not None and volta < 4:
        volta += 1
        novos = [r_ for r_ in alvos_de_reparo(
                     g, os.path.join(pasta, base + "_erros.csv"))
                 if r_ not in descartadas]
        if not novos:
            return _reprova(motivo + " -- reparo sem alvo novo")
        descartadas += novos
        import csv as _csv
        with open(rep_csv, "w", newline="", encoding="utf-8") as f_:
            w_ = _csv.writer(f_, delimiter=";")
            w_.writerow(["rs"])
            for x_ in descartadas:
                w_.writerow([x_])
        print(f"   REPARO {volta}: +{len(novos)} secao(oes) participante(s) "
              f"de defeito de linha descartada(s) (total "
              f"{len(descartadas)}) -> {rep_csv}")
        roda(["scripts/rio_do_relevo.py", "--rio", rio, "--saida", pasta,
              "--dx", str(dx), "--taxa", str(taxa), "--monotono",
              "--excluir", rep_csv], mostrar=("reparo:", "secoes :"))
        roda(["scripts/projeto_rio_avulso.py", g, "--rio-fonte", rio],
             mostrar=())
        motivo = _portoes(g, "g01")
    if motivo is not None:
        return _reprova(motivo + f" -- apos {volta} volta(s) de reparo")

    # ---- 5. batimetria do legado -> g02 (so chega aqui com g01 limpo)
    ped = os.path.join("doc", f"batimetria_{base}.csv")
    roda(["scripts/batimetria.py", "pedir", g, "--cada", str(cada),
          "--saida", ped], mostrar=("pontos    :",))
    roda(["scripts/batimetria_do_legado.py", ped, "--rio", rio],
         mostrar=("casados", "REBAIXAMENTO", "ATENCAO", "SINTETICO"))
    alto = eixo_alto(g, rio)
    if alto is not None:
        km0, km1, reb = alto
        fig = os.path.join("doc", "figuras", f"eixo_alto_{base}.png")
        roda(["scripts/diagnostico_eixo_alto.py", "--rios", rio,
              "--saida", fig])
        print(f"   legado sintetico de {km0:.0f} a {km1:.0f} km "
              f"(lamina-fundo ate {reb:.0f} m): ancoras descartadas, "
              f"MDT mantido. Figura: {fig}")
    roda(["scripts/batimetria.py", "aplicar", g, "--pontos", ped,
          "--saida", "g02"],
         mostrar=("contradeclives", "declividade", "leito bate", "VAO"))
    roda(["scripts/projeto_rio_avulso.py", g2, "--rio-fonte", rio],
         mostrar=("jusante   :",))
    if os.path.exists(g2 + ".hdf"):
        os.remove(g2 + ".hdf")        # item 3 do contrato: stale nao prova
    motivo = _portoes(g2, "g02")
    if motivo:
        return _reprova(motivo)
    print("   APROVADO: g02 com 0 Fatal, 0 GRAVES e 0 defeitos de linha")
    return {"rio": rio, "pasta": pasta, "erros": n, "fatal": fat,
            "taxa": taxa, "pedido": ped, "eixo_alto": alto,
            "status": "APROVADO g02"}


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

    # ---- IDEMPOTENCIA: refazer 765 folhas so quando o dominio CRESCEU.
    # A decisao "precisa reconstruir?" e do software, nao de quem chama com
    # flag: se o raster existente ainda contem todas as geometrias com pelo
    # menos 500 m de folga, nada a refazer -- o religar (passo 8) continua.
    # A folga de construcao e 2000 m; exigir 500 m aqui tolera a geometria
    # crescer ate 1500 m para dentro da folga antiga sem disparar o warp.
    raiz0 = os.path.dirname(os.path.dirname(os.path.abspath(geoms[0])))
    pt0 = os.path.join(raiz0, "Terrain")
    tif0 = os.path.join(pt0, "MDT_SIGSC_30m.tif")
    h0 = os.path.join(pt0, nome + "_Terreno.hdf")
    if os.path.exists(tif0) and os.path.exists(h0):
        try:
            import rasterio
            from terreno_30m import dominio
            (x0, y0, x1, y1), _ = dominio(
                geoms if len(geoms) > 1 else geoms[0], 500.0)
            with rasterio.open(tif0) as s_:
                b = s_.bounds
            if (b.left <= x0 and b.bottom <= y0
                    and b.right >= x1 and b.top >= y1):
                print(f"   terreno existente cobre o dominio "
                      f"({(x1-x0)/1000:.0f} x {(y1-y0)/1000:.0f} km dentro "
                      f"de {(b.right-b.left)/1000:.0f} x "
                      f"{(b.top-b.bottom)/1000:.0f} km) -- nada a refazer")
                return h0
            print("   o dominio cresceu alem do terreno existente -- "
                  "refazendo")
        except Exception as e:                   # noqa: BLE001
            print(f"   nao deu para conferir o terreno existente ({e}) -- "
                  "refazendo")
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

    # ---- 9: as linhas, medidas no HDF da geometria EM USO (item 2 do
    # contrato de aceite): g02.hdf quando o rio foi aprovado com g02, g01.hdf
    # quando reprovado (o projeto voltou a apontar o g01)
    for r in res:
        base_ = os.path.join(r["pasta"], os.path.basename(r["pasta"]))
        h = base_ + (".g02.hdf" if os.path.exists(base_ + ".g02.hdf")
                     else ".g01.hdf")
        r["dobras"] = None
        if os.path.exists(h):
            m = re.search(r"TOTAL: (\d+)",
                          roda(["scripts/conferir_edge_lines.py", h],
                               aceitar_falha=True))
            if m:
                r["dobras"] = int(m.group(1))

    print(f"\n{'='*68}\nRESUMO\n{'='*68}")
    print(f"{'rio':<16}{'erros':>7}{'Fatal':>7}{'taxa':>7}{'dobras':>8}"
          f"   veredito")
    for r in res:
        d = r.get("dobras")
        print(f"{r['rio']:<16}{r['erros']:>7}{r['fatal']:>7}"
              f"{r['taxa']:>7.2f}{('?' if d is None else str(d)):>8}   "
              f"{r.get('status', '?')}")
    ruins = [r for r in res if r["erros"] > a.limite]
    if ruins:
        print(f"\nacima do limite de {a.limite}: "
              + ", ".join(r["rio"] for r in ruins))
        print("   nesses a varzea e plana e a secao vai ao teto; apertar mais")
        print("   a taxa reduz o erro e tambem a capacidade de conter a cheia")
    return res


if __name__ == "__main__":
    main()
