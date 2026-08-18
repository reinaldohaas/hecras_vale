# -*- coding: utf-8 -*-
"""
Erros conhecidos: como DETECTAR e como CORRIGIR cada um, automaticamente.

Este arquivo existe porque a alternativa foi tentada e nao funciona. Cada erro
abaixo foi encontrado a mao, um por vez, olhando log: rodar, ver quebrar,
diagnosticar, corrigir, rodar de novo. Sao horas por erro -- e o erro seguinte
so aparece depois que o anterior sai da frente, porque varios deles mascaravam
uns aos outros.

Aqui cada um vira um par (detectar, corrigir) que o batch aplica sozinho, na
ordem, assim que o problema aparece. O que se ganha nao e velocidade: e que a
correcao deixa de depender de alguem estar olhando o log na hora certa.

DUAS REGRAS, e as duas foram aprendidas errando:

  detectar SEMPRE MEDE, nunca supoe. Metade dos diagnosticos desta reconstrucao
  foram falsos porque a medicao estava contaminada -- area comparada em cota
  absoluta num rio em declive, largura de calha medida com o espacamento do
  lado errado, mediana usada onde o problema estava na cauda.

  corrigir devolve True SO SE mudou alguma coisa. E isso que faz o batch
  reexecutar o passo em vez de entrar em laco.
"""
import glob
import os

import numpy as np

# ---------------------------------------------------------------- limiares
ESCAVACAO_ABSURDA = 40.0     # m; acima disso nao e condicionamento, e invencao
ACHATAMENTO_MAX = 0.60       # fracao da secao numa cota so
COTA_MINIMA = -20.0          # m; abaixo disso nao ha leito no Vale do Itajai
BAK_LIMITE_GB = 1.0          # backups automaticos acumulados


class Problema:
    def __init__(self, tipo, alvo, valor, detalhe=""):
        self.tipo = tipo
        self.alvo = alvo
        self.valor = valor
        self.detalhe = detalhe

    def __str__(self):
        v = f"{self.valor:.2f}" if isinstance(self.valor, float) else self.valor
        extra = f"  ({self.detalhe})" if self.detalhe else ""
        return f"{self.tipo}: {self.alvo} = {v}{extra}"


# ============================================================== DETECCOES
def secoes_sem_terreno(estado, op):
    """Secao cujo talvegue nao pode ser terreno.

    Um unico ponto errado destroi o rio inteiro: a imposicao de monotonia
    obriga cada secao a ficar abaixo da anterior, entao um 0,02 m no meio de
    vizinhas de 387 m arrasta os 94 km seguintes com ele. No Itajai do Oeste
    isso deu escavacao mediana de 344 m -- por causa de DUAS secoes em 156.

    Criterio robusto (mediana movel e desvio absoluto mediano), nao limiar
    fixo: serve para rio de planicie e de serra sem ajuste.
    """
    achados = []
    for rio, v in (estado.get("xs") or {}).items():
        if len(v) < 7:
            continue
        z = np.array([float(x["z_terreno"]) for x in v], float)
        med = np.array([np.median(z[max(0, i - 3):i + 4]) for i in range(len(z))])
        mad = float(np.median(np.abs(z - med))) or 1.0
        fora = np.abs(z - med) > max(8.0 * mad, 10.0)
        for i in np.flatnonzero(fora):
            achados.append(Problema(
                "talvegue impossivel", f"{rio} RS {v[i]['rs']:.0f}",
                float(z[i]), f"vizinhas em {med[i]:.2f} m"))
    return achados


def escavacao_excessiva(estado, op):
    """Perfil condicionado fundo demais em relacao ao terreno."""
    achados = []
    for rio, v in (estado.get("xs_cond") or {}).items():
        if not v:
            continue
        e = np.array([float(x.get("z_terreno", 0.0)) - float(x.get("z_alvo", 0.0))
                      for x in v], float)
        m = float(np.median(e))
        if m > ESCAVACAO_ABSURDA:
            achados.append(Problema("escavacao excessiva", rio, m,
                                    f"mediana de {len(v)} secoes"))
    return achados


def cota_impossivel(estado, op):
    """Leito abaixo do nivel do mar longe da foz."""
    achados = []
    for rio, v in (estado.get("xs_cond") or {}).items():
        for x in v:
            z = float(x.get("z_alvo", 0.0))
            if z < COTA_MINIMA:
                achados.append(Problema(
                    "leito abaixo do mar", f"{rio} RS {x['rs']:.0f}", z))
                break
    return achados


def contrapendente(estado, op):
    """Leito subindo rio abaixo -- uma barragem dentro do modelo."""
    achados = []
    for rio, v in (estado.get("xs_cond") or {}).items():
        w = sorted(v, key=lambda x: -x["rs"])
        for a, b in zip(w, w[1:]):
            za = float(a.get("z_alvo", 0.0))
            zb = float(b.get("z_alvo", 0.0))
            if zb > za + 0.01:
                achados.append(Problema(
                    "contrapendente", f"{rio} RS {b['rs']:.0f}", zb - za,
                    "sobe rio abaixo"))
                break
    return achados


def _plana(x):
    z = np.asarray(x["z"], float)
    _, c = np.unique(np.round(z, 2), return_counts=True)
    return float(c.max() / len(z))


def secoes_achatadas(estado, op):
    """Secao que NOS achatamos -- comparada com o terreno dela, nao sozinha.

    Nao ha teste de area nem de altura que pegue isto: uma bacia chata com
    parede vertical passa nos dois. Foi o que um pilot channel mal escrito fez
    com as 1.232 secoes de um modelo, sem que nenhuma auditoria acusasse --
    quem viu foi o usuario, olhando uma figura.

    CONTRA O TERRENO, e nao contra um limiar absoluto. Medindo so o modelo, a
    checagem acusava dois lugares onde o vale E plano: a foz do Acu (RS 75, o
    estuario) e o Itajai do Oeste RS 83.278, onde o terreno tem 0,0 m de
    desnivel em 575 m. Nos dois o modelo estava MENOS plano que o terreno --
    0,64 contra 0,71 e 0,75 contra 1,00 --, ou seja, a calha nao achatou nada;
    o vale e que e plano. Acusar isso como erro gasta tentativa de correcao
    (gastou tres, encolhendo o pilot de 25 para 5 m sem efeito) e, pior,
    ensina a ignorar a checagem. Vale plano sem desnivel e outro problema, com
    outra correcao: ver desnivel_minimo em vale/secoes.py.
    """
    achados = []
    cru = estado.get("xs") or {}
    for rio, v in (estado.get("xs_pronto") or {}).items():
        terreno = {round(x["rs"], 2): x for x in (cru.get(rio) or [])}
        pior, alvo, ft = 0.0, None, 0.0
        for x in v:
            f = _plana(x)
            t = terreno.get(round(x["rs"], 2))
            base = _plana(t) if t is not None else 0.0
            if f - base > pior:                 # o quanto NOS achatamos
                pior, alvo, ft = f - base, x, f
        if ft > ACHATAMENTO_MAX and pior > 0.15 and alvo is not None:
            achados.append(Problema(
                "secao achatada", f"{rio} RS {alvo['rs']:.0f}", ft,
                f"{100*pior:.0f} pontos percentuais mais plana que o terreno"))
    return achados


def juncao_invalida(estado, op):
    """Juncao com um trecho entrando e um saindo -- o HEC-RAS recusa."""
    achados = []
    for j in (estado.get("juncoes") or []):
        if len(j.get("up") or []) < 2:
            achados.append(Problema(
                "juncao invalida", j.get("nome", "?"), len(j.get("up") or []),
                "precisa de 2 ou mais a montante"))
    return achados


def rio_desconectado(estado, op):
    """Rio que nao entra em juncao nenhuma e nao e o principal.

    A agua dele nao chega a lugar nenhum. Aconteceu com Itajai_Sul e
    Itajai_Oeste, que desaguam na quilometragem ZERO do Acu -- o Acu nasce da
    juncao dos dois, e a regra geral (dividir o receptor) nao se aplicava.
    """
    achados = []
    eixos = estado.get("eixos") or []
    juncoes = estado.get("juncoes") or []
    if not juncoes or not eixos:
        return achados
    tem_up = {r for j in juncoes for r, _ in (j.get("up") or [])}
    principal = max(eixos, key=lambda d: d["area"])["ras"]
    for d in eixos:
        if d["ras"] == principal or d["ras"] in tem_up:
            continue
        achados.append(Problema(
            "rio desconectado", d["ras"], float(d.get("area", 0.0)),
            f"desagua em {d.get('receptor')} mas sem juncao"))
    return achados


def simulacao_incompleta(estado, op):
    """A simulacao chegou ao fim? Abortar no meio e falha, nao aviso.

    A rodada de 18/08/2026 abortou em 01AUG 01:13 com 18,76% de erro de volume
    e terminou anunciando "NENHUM PROBLEMA DETECTADO pelas checagens
    automaticas". O unico caso coberto era o solver nao rodar; simulacao que
    roda e para no meio produz log, produz HDF e produz figura -- e por isso
    passa despercebida com mais facilidade que a que nao roda.
    """
    s = estado.get("solver") or {}
    if not s:
        return []
    achados = []
    if not s.get("rodou"):
        achados.append(Problema("solver nao rodou", "plano", 0.0,
                                s.get("dados") or "sem *.data_errors.txt"))
    elif not s.get("completou"):
        achados.append(Problema(
            "simulacao abortou", s.get("abortou_em") or "?",
            float(s.get("volume") or 0.0),
            f"% de erro de volume; instavel em {s.get('instavel_em')}"))
    elif s.get("volume") is not None and abs(s["volume"]) > 1.0:
        achados.append(Problema("erro de volume alto", "plano",
                                float(s["volume"]), "% -- acima de 1%"))
    return achados


def backups_ocupando_disco(estado, op):
    """Backups automaticos do RasFixit.

    Ele grava o .g01 INTEIRO a cada secao editada mesmo com backup=False -- o
    parametro nao alcanca todos os caminhos internos. Com 2.077 secoes sairam
    1.811 copias de 12 MB: 21 GB, e o disco foi a 100% duas vezes, matando o
    passo seguinte com "Espaco insuficiente no disco".
    """
    g01 = estado.get("g01")
    if not g01:
        return []
    lixo = glob.glob(str(g01) + ".bak*")
    b = sum(os.path.getsize(f) for f in lixo if os.path.exists(f))
    if b / 1e9 > BAK_LIMITE_GB:
        return [Problema("backups automaticos", f"{len(lixo)} arquivos",
                         b / 1e9, "GB ocupados")]
    return []


def htab_ausente(estado, op):
    """Tabela hidraulica no padrao do HEC-RAS.

    Sem HTab por secao o log enche de "Extrapolated above Cross Section
    Table" -- que e a agua passando do topo da TABELA, nao do topo da SECAO.
    Confundir as duas custou uma tarde perseguindo altura de secao.
    """
    g01 = estado.get("g01")
    if not g01 or not os.path.exists(g01):
        return []
    with open(g01, encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    n_xs = txt.count("#Sta/Elev=")
    n_ht = txt.count("XS HTab Starting El and Incr=")
    if n_xs and n_ht < 0.9 * n_xs:
        return [Problema("htab ausente", f"{n_ht} de {n_xs} secoes",
                         float(n_ht), "tabela hidraulica no padrao")]
    return []


def sem_contorno(estado, op):
    """Trecho que precisa de contorno de montante e nao tem."""
    g01, u01 = estado.get("g01"), None
    if not g01:
        return []
    u01 = str(g01)[:-4] + ".u01"
    if not os.path.exists(u01):
        return []
    with open(u01, encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    if txt.count("Flow Hydrograph=") == 0:
        return [Problema("sem contorno", "u01", 0.0,
                         "nenhum hidrograma de cabeceira gravado")]
    return []


# =============================================================== CORRECOES
def corrigir_talvegue(estado, op, probs, log=print):
    """A substituicao pela mediana local ja esta em perfil.rejeitar_absurdos.

    O que falta e reexecutar o passo 5 para que ela aconteca -- e o passo 5 so
    da o resultado certo se partir do terreno, nao do z_alvo anterior.
    """
    log(f"      -> {len(probs)} secoes serao substituidas pela mediana local "
        f"ao reexecutar o condicionamento")
    return True


def corrigir_escavacao(estado, op, probs, log=print):
    """Aperta o piso de escavacao e recondiciona."""
    antes = op.escavacao_maxima
    novo = max(6.0, antes * 0.6)
    if abs(novo - antes) < 0.1:
        log("      -> piso ja no minimo; nao ha o que apertar")
        return False
    op.escavacao_maxima = novo
    log(f"      -> escavacao_maxima {antes:.1f} -> {novo:.1f} m")
    return True


def _faixa_plana(x):
    """Largura, em metros, da maior faixa CONTIGUA numa cota so."""
    z = np.round(np.asarray(x["z"], float), 2)
    sta = np.asarray(x["sta"], float)
    v, c = np.unique(z, return_counts=True)
    m = z == v[c.argmax()]
    melhor, i = 0.0, 0
    while i < len(m):
        if m[i]:
            j = i
            while j < len(m) and m[j]:
                j += 1
            melhor = max(melhor, float(sta[j - 1] - sta[i]))
            i = j
        else:
            i += 1
    return melhor


def corrigir_achatamento(estado, op, probs, log=print):
    """Reduz o pilot channel -- SO SE a faixa plana couber na largura dele.

    MEDIR ANTES DE MEXER. A versao anterior encolhia o pilot sem olhar onde
    estava o plano, e na execucao de 18/08/2026 gastou as tres tentativas
    baixando pilot_largura de 25 para 12, 6 e 5 m enquanto a metrica nao saia
    do lugar (0,75 -> 0,75 -> 0,75 no Oeste; 0,75 -> 0,77 -> 0,77 no Acu). Nao
    podia sair: a faixa plana do Acu RS 75 tem 1.122 m contiguos na cota 0,00
    -- e a lamina do estuario gravada no Copernicus, que e um MDS -- e a do
    Oeste RS 83.278 tem 360 m na cota 349,50, que e o espelho da Barragem
    Oeste. Nenhum entalhe de 25 m alcanca isso. Corretor que nao mede o que vai
    corrigir queima tentativa e ainda sugere que a causa era outra.
    """
    plano = 0.0
    onde = ""
    for p in probs:
        rio = p.alvo.split(" RS ")[0]
        try:
            rs = float(p.alvo.split(" RS ")[1])
        except (IndexError, ValueError):
            continue
        v = (estado.get("xs_pronto") or {}).get(rio) or []
        if not v:
            continue
        x = min(v, key=lambda a: abs(a["rs"] - rs))
        larg = _faixa_plana(x)
        if larg > plano:
            plano, onde = larg, p.alvo
    if plano > op.pilot_largura * 1.5:
        log(f"      -> faixa plana de {plano:.0f} m em {onde}, contra um pilot "
            f"de {op.pilot_largura:.0f} m: o plano esta no TERRENO, nao no "
            f"entalhe")
        log("         (Copernicus e MDS: lamina d'agua de estuario e de "
            "represa entram como plano na cota do espelho)")
        log("         sem correcao automatica -- so troca de fonte de terreno "
            "resolve")
        return False
    if op.pilot_largura <= 5.0:
        log("      -> pilot ja no minimo; problema nao e o entalhe")
        return False
    antes = op.pilot_largura
    op.pilot_largura = max(5.0, antes * 0.5)
    log(f"      -> pilot_largura {antes:.0f} -> {op.pilot_largura:.0f} m")
    return True


def corrigir_backups(estado, op, probs, log=print):
    """Apaga e segue: nao ha o que reexecutar, e so limpeza."""
    g01 = estado.get("g01")
    lixo = glob.glob(str(g01) + ".bak*")
    b = 0
    for f in lixo:
        try:
            b += os.path.getsize(f)
            os.remove(f)
        except OSError:
            pass
    log(f"      -> removidos {len(lixo)} backups ({b/1e9:.1f} GB)")
    return False


def corrigir_htab(estado, op, probs, log=print):
    if op.usar_htab:
        log("      -> usar_htab ja esta ligado; a chamada falhou. "
            "Confira o log do passo que escreve a geometria.")
        return False
    op.usar_htab = True
    log("      -> usar_htab ligado")
    return True


def corrigir_juncao(estado, op, probs, log=print):
    """Afrouxa o minimo de secoes por trecho, que e o que impede a divisao."""
    from . import rede
    if rede.MIN_SECOES_TRECHO <= 2:
        log("      -> MIN_SECOES_TRECHO ja no minimo")
        return False
    rede.MIN_SECOES_TRECHO -= 1
    log(f"      -> MIN_SECOES_TRECHO -> {rede.MIN_SECOES_TRECHO}")
    return True


# ================================================================ REGISTRO
# (passo apos o qual checar, nome, detectar, corrigir, passo a reexecutar)
CHECAGENS = [
    (4, "talvegue impossivel", secoes_sem_terreno, corrigir_talvegue, 5),
    (5, "escavacao excessiva", escavacao_excessiva, corrigir_escavacao, 5),
    (5, "leito abaixo do mar", cota_impossivel, corrigir_escavacao, 5),
    (5, "contrapendente", contrapendente, None, None),
    (6, "secao achatada", secoes_achatadas, corrigir_achatamento, 6),
    (7, "juncao invalida", juncao_invalida, corrigir_juncao, 7),
    (7, "rio desconectado", rio_desconectado, None, None),
    (7, "htab ausente", htab_ausente, corrigir_htab, 7),
    (7, "sem contorno", sem_contorno, None, None),
    (8, "backups automaticos", backups_ocupando_disco, corrigir_backups, None),
    (9, "simulacao incompleta", simulacao_incompleta, None, None),
]


def checar_passo(n, estado, op, log=print):
    """Checagens do passo n. Devolve (problemas, passo_a_reexecutar_ou_None)."""
    achados, refazer = [], None
    for passo, nome, detectar, corrigir, refaz in CHECAGENS:
        if passo != n:
            continue
        try:
            probs = detectar(estado, op) or []
        except Exception as e:                                # noqa: BLE001
            log(f"   [check] {nome}: NAO PODE SER MEDIDO ({e})")
            continue
        if not probs:
            continue
        log(f"   [check] {nome}: {len(probs)} ocorrencia(s)")
        for x in probs[:3]:
            log(f"      {x}")
        if len(probs) > 3:
            log(f"      ... e mais {len(probs) - 3}")
        achados += probs
        if corrigir is None:
            log("      -> sem correcao automatica; anotado para revisao")
            continue
        if corrigir(estado, op, probs, log) and refaz:
            refazer = min(refazer, refaz) if refazer else refaz
    return achados, refazer
