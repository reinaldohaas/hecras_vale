# -*- coding: utf-8 -*-
"""
Orquestrador: roda o modelo do Vale um passo por vez.

    python -m vale                     lista os passos e o que ja foi feito
    python -m vale opcoes              mostra todas as opcoes e os valores
    python -m vale 1                   roda so o passo 1
    python -m vale 1-4                 roda do 1 ao 4, confirmando cada um
    python -m vale tudo --sim          roda tudo sem perguntar
    python -m vale 4 selecao=todos corredor_m=500

Cada passo grava o resultado em disco e o proximo le de la. Isso permite parar,
conferir, mudar uma opcao e retomar do meio -- que e como este modelo precisa
ser construido, porque quase toda decisao so se avalia olhando o resultado da
anterior.

Nada e feito por conta propria: com --sim o programa avisa o que vai fazer e
faz; sem --sim, pergunta. As opcoes ficam todas em vale/config.py, visiveis.
"""
import os
import pickle
import sys
import time

from .config import Opcoes

ESTADO = "estado.pkl"


# --------------------------------------------------------------- utilitario
def _log(msg=""):
    print(msg, flush=True)


def _logs_para_stdout():
    """Manda as mensagens de biblioteca para stdout, e nao para stderr.

    O ras-commander instala um StreamHandler em stderr, e o .bat canaliza a
    saida por "2>&1 | Tee-Object". O PowerShell 5.1 embrulha CADA linha de
    stderr de um executavel nativo num NativeCommandError, com cabecalho,
    CategoryInfo e FullyQualifiedErrorId. O resultado e que um INFO banal como
    "Checking HTAB starting elevations for 2077 cross sections" ocupa seis
    linhas e parece falha -- e este projeto ja perdeu horas com erro de verdade
    escondido no meio de ruido que parecia erro. Mandando para stdout, o
    Tee-Object recebe texto limpo e o que sobrar em stderr e erro mesmo.
    """
    import logging
    for h in list(logging.root.handlers):
        if isinstance(h, logging.StreamHandler) and h.stream is sys.stderr:
            h.setStream(sys.stdout)
    for nome in list(logging.root.manager.loggerDict):
        for h in list(getattr(logging.getLogger(nome), "handlers", [])):
            if isinstance(h, logging.StreamHandler) and h.stream is sys.stderr:
                h.setStream(sys.stdout)


def carregar(op, chave, obrigatorio=True):
    caminho = op.caminho(ESTADO)
    if not os.path.exists(caminho):
        if obrigatorio:
            raise SystemExit(
                f"passo anterior nao rodou: nao ha {caminho}.\n"
                f"Rode 'python -m vale 1' primeiro.")
        return None
    with open(caminho, "rb") as f:
        d = pickle.load(f)
    if chave not in d and obrigatorio:
        raise SystemExit(
            f"falta '{chave}' no estado. Rode o passo que o produz "
            f"(veja 'python -m vale').")
    return d.get(chave)


def salvar(op, **kw):
    caminho = op.caminho(ESTADO)
    os.makedirs(op.saida, exist_ok=True)
    d = {}
    if os.path.exists(caminho):
        with open(caminho, "rb") as f:
            d = pickle.load(f)
    d.update(kw)
    d["_quando"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(caminho, "wb") as f:
        pickle.dump(d, f)
    return caminho


# -------------------------------------------------------------------- passos
def passo_rios(op):
    from .rios import catalogo, selecionar, imprimir
    cat = catalogo(op.area_minima, op.bho)
    esc = selecionar(cat, op.selecao)
    imprimir(cat, esc)
    salvar(op, catalogo=cat, selecao=esc)
    return f"{len(esc)} rios selecionados de {len(cat)} no catalogo"


def passo_eixos(op):
    from .eixos import montar, gravar
    eixos, principal, perdidos = montar(op.selecao, op.area_minima, op.bho)
    _log(f"   principal: {principal['ras']}")
    for d in sorted(eixos, key=lambda x: -x["area"]):
        _log(f"   {d['ras']:<16}{d['area']:>8.0f} km2  {d['km_eixo']:>7.1f} km"
             f"   -> {d.get('receptor') or 'foz'}")
    if perdidos:
        _log(f"   SEM EIXO: {', '.join(d['ras'] for d in perdidos)}")
    gravar(eixos, op.caminho(f"{op.projeto}_eixos.geojson"))
    salvar(op, eixos=eixos)
    return f"{len(eixos)} eixos"


def passo_terreno(op):
    from . import terreno
    eixos = carregar(op, "eixos")
    _log(f"   fonte: {op.fonte}")
    terreno.estimativa(op, eixos, _log)
    tifs = terreno.preparar(op, eixos, _log)
    hdf = None
    try:
        if not op.terreno_hdf:
            _log("   terreno .hdf pulado (terreno_hdf=false); o modelo roda "
                 "igual, so o RAS Mapper fica sem relevo")
            raise RuntimeError("pulado por opcao")
        hdf = terreno.hdf(op, tifs, _log)
    except Exception as e:                                   # noqa: BLE001
        if "pulado por opcao" not in str(e):
            _log(f"   terreno .hdf NAO gerado ({e}); o modelo roda mesmo "
                 f"assim, mas o RAS Mapper nao mostra relevo")
    salvar(op, geotiffs=tifs, terreno_hdf=hdf)
    return f"{len(tifs)} GeoTIFF" + (" + .hdf" if hdf else "")


def passo_secoes(op):
    from . import secoes, terreno
    from .progresso import Progresso
    eixos = carregar(op, "eixos")
    tifs = carregar(op, "geotiffs")
    # total estimado pelo comprimento dos eixos: o espacamento e adaptativo,
    # entao o numero exato so se sabe cortando -- mas para dar tempo restante
    # uma estimativa serve, e o contador se corrige sozinho no fim
    est = sum(max(int(d["linha"].length / op.espacamento_min), 1) for d in eixos)
    prog = Progresso(est, "secoes", log=None)
    am = terreno.Amostrador(tifs)
    try:
        xs = {d["ras"]: secoes.cortar_rio(d, am, op, _log, prog) for d in eixos}
    finally:
        am.fechar()
        prog.fim()
    salvar(op, xs=xs)
    return f"{sum(len(v) for v in xs.values())} secoes"


def passo_perfil(op):
    """Condiciona o perfil a partir das secoes CRUAS do passo 4.

    Copia antes de mexer. Sem isso o passo nao e reexecutavel: ele apara a
    cabeceira e grava a lista aparada, entao rodar de novo apara de novo --
    o Itajai do Oeste caiu de 156 para 127 secoes so por ter sido rodado tres
    vezes. E o passo 6 sobrescreve d["z"] com a calha, de modo que o passo 5
    passaria a ler o talvegue da CALHA em vez do terreno, afundando o rio a
    cada volta.
    """
    import copy
    from . import perfil
    eixos = carregar(op, "eixos")
    xs = copy.deepcopy(carregar(op, "xs"))
    por_ras = {d["ras"]: d for d in eixos}
    # do maior para o menor: o afluente ancora no leito de quem o recebe, e o
    # receptor precisa estar condicionado antes
    for d in sorted(eixos, key=lambda x: -x["area"]):
        ras = d["ras"]
        if ras not in xs:
            continue
        # primeira confluencia deste rio, para o aparo nao engoli-la
        filhos = [m for m in eixos if m.get("receptor") == ras
                  and m.get("recebe_em") is not None]
        rs_lim = None
        if filhos:
            rs_lim = d["linha"].length - min(m["recebe_em"] for m in filhos)
        xs[ras] = perfil.condicionar(xs[ras], op, rs_lim, ras, _log)
        rec = d.get("receptor")
        if rec and rec in xs and d.get("recebe_em") is not None:
            alvo_rs = por_ras[rec]["linha"].length - d["recebe_em"]
            leito = min(xs[rec], key=lambda x: abs(x["rs"] - alvo_rs))
            perfil.ancorar(xs[ras], perfil.cota(leito), _log, ras, op)
        perfil.manning(xs[ras], op)
    salvar(op, xs_cond=xs)
    return f"perfis condicionados ({sum(len(v) for v in xs.values())} secoes)"


def passo_calha(op):
    """Escava a calha sobre o perfil condicionado, sempre em copia nova."""
    import copy
    from . import calha
    xs = copy.deepcopy(carregar(op, "xs_cond"))
    for ras, v in xs.items():
        calha.escavar_rio(v, op, _log, ras)
    salvar(op, xs_pronto=xs)
    return ("calha escavada" if op.escavar else "calha imposta sem batimetria")


def passo_escrever(op):
    from . import hidrologia, projeto, rede
    eixos = carregar(op, "eixos")
    xs = carregar(op, "xs_pronto")
    trechos, juncoes = rede.montar(eixos, xs, _log)
    if not trechos:
        raise SystemExit("nenhum trecho montado -- confira os passos 4 a 6")

    cab, lat = hidrologia.series(op, eixos, xs, None, _log)

    por_rio = {}
    for t in trechos:
        por_rio.setdefault(t["rio"], []).append(t)

    # o contorno de cabeceira vai no trecho mais de MONTANTE de cada rio; a
    # lateral, no trecho que contem a faixa
    precisa = rede.contornos_necessarios(trechos, juncoes)
    _log(f"      {len(precisa)} trechos precisam de contorno de montante")
    cabs = []
    for c in cab:
        alvos = [t for t in precisa if t["rio"] == c["rio"]]
        if not alvos:
            continue
        t = min(alvos, key=lambda t: t["a"])
        cabs.append({**c, "reach": t["reach"], "xs": t["xs"]})
    lats = []
    for l in lat:
        v = por_rio.get(l["rio"]) or []
        if not v:
            continue
        t = max(v, key=lambda t: len(t["xs"]))
        rss = sorted((x["rs"] for x in t["xs"]), reverse=True)
        if len(rss) < 4:
            continue
        lats.append({**l, "reach": t["reach"], "rs_hi": rss[1], "rs_lo": rss[-2]})

    # A vazao inicial so pode ser calculada AGORA: ela depende de qual trecho
    # ficou com o contorno de cabeceira e qual ficou com a lateral, e isso e
    # decidido logo acima. Antes disso o valor era por RIO, e o mesmo numero ia
    # para todos os trechos -- ver hidrologia.inicial_por_trecho.
    q0 = hidrologia.inicial_por_trecho(trechos, juncoes, cabs, lats, _log)
    for t in trechos:
        t["q_base"] = max(q0.get((t["rio"], t["reach"]), 0.0), 0.05)

    principal = max(eixos, key=lambda d: d["area"])["ras"]
    saida = por_rio[principal][-1]

    # Mare so se a foz estiver no mar. Rodando um rio isolado a foz dele vira a
    # saida do modelo, e impor mare de 0,3 m a uma secao com fundo em 50 m faz
    # o HEC-RAS recusar os dados antes de computar.
    import numpy as _np
    ult = saida["xs"][-1]
    cota_saida = float(_np.min(ult["z"]))
    if cota_saida > op.cota_mare:
        mare = None
        v = saida["xs"][-4:]
        dz = float(_np.min(v[0]["z"]) - _np.min(v[-1]["z"]))
        dx = float(v[0]["rs"] - v[-1]["rs"])
        decl = dz / dx if dx > 0 else op.decl_minima
        _log(f"      foz em {cota_saida:.1f} m (acima de {op.cota_mare:.0f} m): "
             f"contorno de jusante por profundidade normal, "
             f"declividade {100*max(decl, op.decl_minima):.3f}%")
    else:
        mare = hidrologia.mare(op.horas)
        decl = None
    g01 = projeto.geometria(op, trechos, juncoes)
    u01 = projeto.fluxo(op, trechos, cabs, saida, mare, lats,
                        hidrologia.INICIO, decl)
    p01, prj = projeto.plano(op, hidrologia.INICIO)
    rmap = projeto.rasmap(op, carregar(op, "terreno_hdf", False))
    for c in (g01, u01, p01, prj, rmap):
        _log(f"   {os.path.basename(c)}")
    salvar(op, trechos=trechos, juncoes=juncoes, g01=g01, prj=prj)
    return f"{len(trechos)} trechos, {len(juncoes)} juncoes"


def passo_corrigir(op):
    from . import correcao
    g01 = carregar(op, "g01")
    prj = carregar(op, "prj")
    tifs = carregar(op, "geotiffs")
    correcao.aplicar(op, g01, _log)
    correcao.checar(op, prj, _log)
    eixo_geojson = op.caminho(f"{op.projeto}_eixos.geojson")
    if os.path.exists(eixo_geojson) and tifs:
        correcao.auditar_secoes(op, g01, tifs[0], eixo_geojson, _log)
        correcao.auditar_terreno(op, g01, tifs, _log)
    return "correcao e auditoria aplicadas"


def _impressao_digital(op):
    """Diz QUAL geometria vai rodar, no proprio log da rodada.

    Em 18/08/2026 uma simulacao foi lida como "a correcao nao teve efeito"
    quando na verdade rodara sobre um .g01 de duas horas antes: o log saiu
    identico ao anterior, byte por byte, e log identico nunca e resultado de
    correcao sem efeito -- e sinal de que nada foi recomputado. Sem carimbo,
    descobrir isso exigiu abrir o .g01 e medir o fundo de uma secao. Com
    carimbo, basta comparar duas linhas.
    """
    import hashlib
    for nome in (f"{op.projeto}.g01", f"{op.projeto}.u01"):
        c = op.caminho(nome)
        if not os.path.exists(c):
            continue
        b = open(c, "rb").read()
        _log(f"   {nome}: {hashlib.sha1(b).hexdigest()[:10]}  "
             f"{len(b)/1e6:.1f} MB  "
             f"{time.strftime('%d/%m %H:%M', time.localtime(os.path.getmtime(c)))}")


def passo_rodar(op):
    from . import executar
    prj = carregar(op, "prj")
    _impressao_digital(op)
    r, msgs, pasta = executar.rodar(prj, op.ras_exe, _log)
    _log(f"   {r}")
    _log("")
    _log(executar.resumir(msgs, pasta))
    with open(op.caminho("compute.log"), "w", encoding="utf-8") as f:
        f.write(msgs)
    # Falhar em rodar tem de FALHAR o passo. Sem isto o passo 9 terminava em
    # 2 s, o passo 10 gerava a pagina de uma cheia que nunca foi computada, e a
    # rodada fechava com "NENHUM PROBLEMA DETECTADO".
    if not executar.rodou_de_fato(msgs):
        dados = executar.erros_de_dado(pasta)
        raise SystemExit(
            "o solver nao rodou -- nenhuma mensagem de computacao.\n"
            + (dados or "sem *.data_errors.txt para explicar"))
    # Uma figura por ponto que o solver acusou, COM A LAMINA que ele calculou.
    # E o que separa "geometria ruim" de "modelo rodando seco": as duas coisas
    # produzem o mesmo log, e so a figura mostra 2 cm de agua num entalhe.
    try:
        from . import figura
        figura.do_solver(_estado(op), op.caminho("figuras"),
                         executar.pontos_criticos(msgs), _log)
    except Exception as e:                                   # noqa: BLE001
        _log(f"   figuras do solver nao geradas: {e}")
    salvar(op, run_pasta=pasta,
           solver=executar.diagnostico(msgs, pasta))
    return "rodado (log em compute.log)"


def passo_visual(op):
    from . import visual
    trechos = carregar(op, "trechos")
    destino = visual.gerar(op, trechos)
    _log(f"   {destino}")
    return destino


PASSOS = [
    ("rios", "catalogo da ANA e selecao", passo_rios),
    ("eixos", "eixos dos rios, das linhas da ANA", passo_eixos),
    ("terreno", "mosaico SIG-SC (corredor 1 m + fundo) e o .hdf", passo_terreno),
    ("secoes", "cortar as secoes transversais", passo_secoes),
    ("perfil", "condicionar o perfil longitudinal", passo_perfil),
    ("calha", "escavar a calha e o pilot channel", passo_calha),
    ("escrever", "gravar geometria, fluxo e plano do HEC-RAS", passo_escrever),
    ("corrigir", "ferramentas de correcao do HEC-RAS e auditoria", passo_corrigir),
    ("rodar", "computar e ler o log do solver", passo_rodar),
    ("visual", "pagina interativa da cheia", passo_visual),
]


# ---------------------------------------------------------------------- CLI
def _faixa(txt):
    txt = str(txt).strip().lower()
    if txt in ("tudo", "todos", ""):
        return list(range(1, len(PASSOS) + 1))
    saida = []
    for parte in txt.split(","):
        parte = parte.strip()
        if "-" in parte:
            a, b = (int(x) for x in parte.split("-"))
            saida += list(range(min(a, b), max(a, b) + 1))
        elif parte.isdigit():
            saida.append(int(parte))
        else:
            nomes = [i for i, p in enumerate(PASSOS, 1) if p[0] == parte]
            if not nomes:
                raise SystemExit(f"passo desconhecido: {parte!r}")
            saida += nomes
    return [n for n in saida if 1 <= n <= len(PASSOS)]


def listar(op):
    caminho = op.caminho(ESTADO)
    feito = {}
    if os.path.exists(caminho):
        with open(caminho, "rb") as f:
            feito = pickle.load(f)
    chaves = {1: "selecao", 2: "eixos", 3: "geotiffs", 4: "xs",
              5: "xs_cond", 6: "xs_pronto", 7: "g01", 9: "run_pasta"}
    print(f"projeto '{op.projeto}'   pasta '{op.saida}'   "
          f"selecao '{op.selecao}'\n")
    for i, (nome, desc, _) in enumerate(PASSOS, 1):
        marca = "ok" if chaves.get(i) in feito else "--"
        print(f"  [{marca}] {i:>2}  {nome:<10} {desc}")
    print(f"\nultimo passo gravado: {feito.get('_quando', '(nada ainda)')}")
    print("\n  python -m vale 1        roda um passo")
    print("  python -m vale 1-4      roda uma faixa")
    print("  python -m vale tudo --sim      roda tudo sem perguntar")
    print("  python -m vale tudo --auto     roda tudo, checando e corrigindo")
    print("  python -m vale.qaqc            qualifica o SIG-SC vs Copernicus")
    print("  python -m vale opcoes")


def _estado(op):
    cam = op.caminho(ESTADO)
    if not os.path.exists(cam):
        return {}
    with open(cam, "rb") as f:
        return pickle.load(f)


def _auto(op, passos, tentativas=3):
    """Roda os passos EM SERIE, checando e corrigindo a cada um.

    E o que a construcao a mao nao consegue. Cada erro desta reconstrucao foi
    achado olhando log, e o seguinte so aparecia depois que o anterior saia da
    frente -- varios se mascaravam entre si. O passo nao reexecutavel escondia
    a contaminacao da entrada, que por sua vez escondia o valor absurdo no
    terreno: foram tres medicoes so para separar os tres, e cada uma custou
    uma rodada inteira.

    Aqui a checagem roda sempre, na ordem, e a correcao entra na hora em que o
    problema aparece. Quando uma correcao muda uma opcao, o passo afetado e
    REEXECUTADO -- ate 'tentativas' vezes. Correcao que nao muda nada nao
    reexecuta, e e assim que o laco termina em vez de girar.
    """
    from . import validacao
    from .progresso import Etapas

    peso_terreno = {"sigsc": 60.0, "misto": 25.0, "copernicus": 1.0}
    pesos = {"rios": 0.2, "eixos": 0.5,
             "terreno": peso_terreno.get(op.fonte, 10.0),
             "secoes": 10.0, "perfil": 1.0, "calha": 1.0, "escrever": 3.0,
             "corrigir": 2.0, "rodar": 5.0, "visual": 0.5}
    et = Etapas({PASSOS[n - 1][0]: pesos.get(PASSOS[n - 1][0], 1.0)
                 for n in passos}, _log)

    relatorio = []
    i = 0
    while i < len(passos):
        n = passos[i]
        nome, desc, fn = PASSOS[n - 1]
        print()
        print("=" * 74)
        print(f"[{n}] {nome} -- {desc}")
        print("=" * 74)
        et.inicia(nome)
        for tentativa in range(1, tentativas + 1):
            t0 = time.time()
            # depois de cada passo: as bibliotecas sao importadas por dentro
            # dos passos, entao o handler em stderr so existe depois do
            # primeiro import -- reaplicar e barato e pega todos
            _logs_para_stdout()
            try:
                saida = fn(op)
                _logs_para_stdout()
            except SystemExit:
                raise
            except Exception as e:                           # noqa: BLE001
                import traceback
                traceback.print_exc()
                relatorio.append((n, nome, f"FALHOU: {e}"))
                print()
                print(f"FALHOU no passo {n} ({nome}): {e}")
                print()
                print(_resumo_auto(relatorio))
                return 1
            print()
            print(f"-> {saida}   [{time.time() - t0:.0f} s]")

            probs, refazer = validacao.checar_passo(n, _estado(op), op, _log)
            for x in probs:
                relatorio.append((n, nome, str(x)))
            # REGRA: todo problema detectado vira imagem. Texto e tabela
            # escondem geometria quebrada -- o pilot channel que achatou 1.232
            # secoes passou por toda a auditoria numerica sem acusar nada.
            if probs:
                try:
                    from . import figura
                    figura.dos_problemas(_estado(op), probs,
                                         op.caminho("figuras"), _log)
                except Exception as e:                       # noqa: BLE001
                    _log(f"   figuras nao geradas: {e}")
            if not refazer:
                break
            if tentativa >= tentativas:
                print(f"   [check] {tentativas} tentativas no passo {n}; "
                      f"seguindo com o problema anotado")
                break
            alvo_refaz = refazer if refazer in passos else n
            print(f"   [check] reexecutando o passo {alvo_refaz} "
                  f"(tentativa {tentativa + 1} de {tentativas})")
            i = passos.index(alvo_refaz)
            n = passos[i]
            nome, desc, fn = PASSOS[n - 1]
        et.termina()
        i += 1

    print()
    print(_resumo_auto(relatorio))
    print()
    print(et.resumo())
    return 0


def _resumo_auto(relatorio):
    if not relatorio:
        return "NENHUM PROBLEMA DETECTADO pelas checagens automaticas."
    L = [f"PROBLEMAS DETECTADOS: {len(relatorio)}", ""]
    atual = None
    for n, nome, txt in relatorio:
        if (n, nome) != atual:
            atual = (n, nome)
            L.append(f"  passo {n} ({nome}):")
        L.append(f"     {txt}")
    L += ["", "  Os marcados 'sem correcao automatica' precisam de decisao."]
    return chr(10).join(L)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    sim = "--sim" in argv
    auto = "--auto" in argv
    argv = [a for a in argv if a not in ("--sim", "--auto")]
    conf = None
    for a in list(argv):
        if a.startswith("--config="):
            conf = a.split("=", 1)[1]
            argv.remove(a)
    op = Opcoes.ler(conf) if conf else Opcoes()

    alvo = argv[0] if argv and not argv[0].count("=") else None
    pares = [a for a in argv if "=" in a]
    op.aplicar(pares)
    if sim or auto:
        op.confirmar = False

    if alvo == "custo":
        from . import terreno
        op.coerir(_log)
        eixos = carregar(op, "eixos", False)
        if not eixos:
            print("rode 'python -m vale 2' primeiro (os eixos definem o "
                  "corredor).")
            return 1
        terreno.estimativa(op, eixos, print)
        return 0
    if alvo == "opcoes":
        for k, v in op.dict().items():
            print(f"  {k:<20} {v!r}")
        return 0
    if alvo is None:
        listar(op)
        return 0

    os.makedirs(op.saida, exist_ok=True)
    op.gravar()
    op.coerir(_log)
    passos = _faixa(alvo)
    if auto:
        return _auto(op, passos)

    for n in passos:
        nome, desc, fn = PASSOS[n - 1]
        print(f"\n{'='*74}\n[{n}] {nome} -- {desc}\n{'='*74}")
        if op.confirmar:
            r = input("rodar este passo? [S/n/q] ").strip().lower()
            if r == "q":
                print("parado.")
                return 0
            if r == "n":
                print("pulado.")
                continue
        t0 = time.time()
        try:
            saida = fn(op)
        except SystemExit:
            raise
        except Exception as e:                               # noqa: BLE001
            import traceback
            traceback.print_exc()
            print(f"\nFALHOU no passo {n} ({nome}): {e}")
            print("As opcoes estao em", op.caminho("opcoes.json"))
            return 1
        print(f"\n-> {saida}   [{time.time()-t0:.0f} s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
