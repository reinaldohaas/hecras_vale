# -*- coding: utf-8 -*-
"""
Escavacao da calha e o pilot channel. Aplicada UMA vez, no fim.

Por que escavar agora faz sentido, e antes nao fazia. Com MDS de 30 m
(Copernicus) o dado JA CONTINHA a lamina d'agua, gravada como um plano na cota
do espelho: escavar em cima disso contava a profundidade duas vezes, o modelo
partia seco e o assentamento tinha de encher 382 hm3 de canal inventado. Com
MDT (SIG-SC, solo exposto) o espelho d'agua nao esta la -- o leito submerso
esta genuinamente ausente do dado, e a batimetria sintetica passa a ser a
melhor aproximacao disponivel em vez de uma duplicacao.

Fica registrado o que ela e: relacao regional de Leopold & Maddock aplicada a
area de drenagem, NAO batimetria medida. Se houver levantamento do Itajai-Acu,
e aqui que ele entra.

CALHA IMPOSTA, nao subtraida. Subtrair uma profundidade constante preserva a
forma do terreno -- vale em V da calha em V, vale em U da calha em U -- e a
area molhada salta entre vizinhas: eram 121 pares com mais de 3x de diferenca
a 1 m de lamina, contra 31 a 8 m, ou seja o defeito estava na CALHA e nao na
planicie. Imposta como trapezio de fundo plano, a conducao de estiagem passa a
ser continua por construcao.

O LIMITE TEM DE VOLTAR AO TERRENO. Escrito como
    z = min(z, base + sobe * PROFUNDIDADE)
o limite satura no fundo longe do canal, e min(z, fundo) rebaixa a secao
INTEIRA. Foi assim que as 1.232 secoes de um modelo viraram bacias chatas com
um entalhe no meio -- e a auditoria continuou dizendo "0 saltos de area, 0
secoes rasas", porque uma bacia chata com parede vertical passa nos dois.
"""
import numpy as np


def canal(area_km2, op):
    """(profundidade, largura) da calha para a area de drenagem."""
    a = max(float(area_km2), 1.0)
    prof = op.canal_kh * a ** op.canal_eh if op.escavar else 0.0
    return float(prof), float(op.canal_kw * a ** op.canal_ew)


def vazao_projeto(area_km2, k=2.6, exp=0.80):
    """Vazao de pico de referencia, so para DIMENSIONAR a secao.

    Ancorada nos ~5.700 m3/s de 1983 na foz, com 15.000 km2. A vazao que roda
    vem da hidrologia; esta serve para saber que altura a secao precisa ter.
    """
    return k * max(float(area_km2), 1.0) ** exp


def altura_para_vazao(sta, z, n, S, Q, h_max=30.0):
    """Menor lamina sobre o talvegue que conduz Q, por Manning na propria secao.

    O criterio anterior era um numero fixo (12 m de parede para todos), igual
    para o Itajai-Acu e para um afluente de 30 km2. O solver listava
    nominalmente as secoes que passavam do topo da tabela, e TODAS estavam
    travadas nesse valor -- o piso, nao o terreno.
    """
    z = np.asarray(z, float)
    z0 = float(np.min(z))
    S = max(float(S or 0.0), 1e-4)
    n = max(float(n or 0.035), 0.02)
    for h in np.arange(1.0, h_max + 0.01, 0.5):
        prof = np.clip(z0 + h - z, 0.0, None)
        A = float(np.trapezoid(prof, sta)) if hasattr(np, "trapezoid") \
            else float(np.trapz(prof, sta))
        if A <= 0.0:
            continue
        molh = (prof[:-1] + prof[1:]) > 0.0
        P = float(np.sum(np.hypot(np.diff(sta), np.diff(z))[molh]))
        if P > 0.0 and (A / n) * (A / P) ** (2.0 / 3.0) * S ** 0.5 >= Q:
            return float(h)
    return float(h_max)


def margens(sta, z, i0, prof_canal, folga=3.0):
    """Estacas das margens: do talvegue ate o topo da calha mais a folga.

    Medir a folga a partir do TALVEGUE poe a margem dentro do canal escavado e
    o modelo passa a achar que tudo extravasa. A margem e o topo da calha.
    """
    lim = z[i0] + prof_canal + folga
    e = i0
    while e > 0 and z[e - 1] < lim:
        e -= 1
    d = i0
    while d < len(z) - 1 and z[d + 1] < lim:
        d += 1
    e = min(max(e, 1), len(sta) - 3)
    d = max(min(d, len(sta) - 2), e + 1)
    return round(float(sta[e]), 2), round(float(sta[d]), 2)


def largura_pilot(d, op, larg_calha):
    """Largura do entalhe: a que da a profundidade-alvo com a VAZAO DE BASE.

    O entalhe existe para dar um caminho bem condicionado a lamina baixa. Com
    largura constante ele nao faz isso nos dois extremos do rio ao mesmo tempo:
    25 m sao um entalhe de verdade numa secao de 3 km na foz do Acu, e sao a
    calha INTEIRA na cabeceira do Benedito, onde a largura de margens plenas
    por Leopold (75 km2) da 28 m. Resultado medido la: 0,58 m3/s espalhados por
    25 m a 5% de declividade dao 6 CENTIMETROS de lamina, e o solver nao
    resolve isso -- foi o que abortou a rodada com a vazao correta.

    Manning para canal largo (R ~ h):  b = Q*n / (h^(5/3) * sqrt(S)).

    Nao adianta compensar com vazao: para a mesma lamina seria preciso
    multiplicar a vazao de base por quinze, o que da 115 L/s/km2 -- cheia, e
    nao escoamento de base. A largura e a variavel certa.
    """
    from .hidrologia import AREA_REF_FOZ, Q_REF_FOZ
    q = Q_REF_FOZ * (max(float(d.get("area_km2", 1.0)), 1.0) / AREA_REF_FOZ)
    q *= op.base_frac
    S = max(float(d.get("S_terreno") or 0.0), op.decl_minima)
    n = float(d.get("n") or 0.05)
    b = q * n / (op.pilot_prof_alvo ** (5.0 / 3.0) * np.sqrt(S))
    teto = min(op.pilot_largura, max(larg_calha, op.pilot_largura_min))
    return float(np.clip(b, op.pilot_largura_min, teto))


def escavar(d, op, altura_minima=12.0, folga_altura=1.4, altura_max=30.0):
    """Aplica a calha, o pilot channel e -- so onde precisa -- a parede."""
    sta = d["sta"]
    z = np.array(d["z"], float)
    i0 = d["i_thal"]
    prof, larg = canal(d["area_km2"], op)
    alvo = float(d.get("z_alvo", z[i0] - prof))
    dist = np.abs(sta - sta[i0])

    # 1) trapezio de fundo plano em 'alvo', subindo ate reencontrar o terreno
    meia = larg / 2.0
    talude = max(larg * 0.25, 30.0)
    sobe = np.clip((dist - meia) / talude, 0.0, 1.0)
    z = np.minimum(z, alvo + sobe * np.maximum(z - alvo, 0.0))

    # 2) pilot channel. Com o fundo chato de ~100 m, em lamina baixa a area
    #    molhada e o raio hidraulico ficam mal definidos -- alguns centimetros
    #    espalhados por 100 m -- e a conducao oscila desde a primeira iteracao.
    #    O entalhe da um caminho continuo e bem condicionado sem alterar a
    #    capacidade de cheia (25 m x 1,5 m sao 37 m2 num rio de milhares).
    if op.pilot_prof > 0:
        base_p = alvo - op.pilot_prof
        larg_p = largura_pilot(d, op, larg)
        d["pilot_largura"] = larg_p
        meia_p = larg_p / 2.0
        talude_p = max(larg_p * 0.6, 10.0)
        sobe_p = np.clip((dist - meia_p) / talude_p, 0.0, 1.0)
        z = np.minimum(z, base_p + sobe_p * np.maximum(z - base_p, 0.0))
    else:
        base_p = alvo

    # 2b) NADA ABAIXO DO FUNDO. A escavacao so sabe BAIXAR, entao um buraco do
    #     terreno mais fundo que o alvo sobrevive -- e o HEC-RAS toma o MINIMO
    #     da secao como cota de fundo, nao o z_alvo que o perfil escolheu. O
    #     perfil sai monotonico e a geometria nao: com o leito perto do terreno
    #     apareceram 59 contrapendentes na secao (49 so no Iraputa, a pior de
    #     8,31 m) com z_alvo perfeitamente monotonico. Enquanto a escavacao era
    #     de 13 m isso ficava escondido, porque nenhum buraco chegava tao fundo.
    #     Os buracos sao artefato do MDS (clareira lida como depressao); levanta
    #     -los e mais honesto que deixar o solver achar um poco no meio do rio.
    z = np.maximum(z, base_p)

    # 2c) O TALVEGUE E EXATAMENTE base_p, e nao "o que calhou". O entalhe e
    #     cortado por `minimum` sobre os PONTOS da secao, e o piso dele so e
    #     alcancado por um ponto que caia DENTRO da meia-largura. Com o entalhe
    #     no minimo de 3 m (meia-largura 1,5 m) e os pontos a cada 5 m
    #     (espacamento_pontos), nenhum ponto cai la: a secao nunca chega a
    #     base_p e o fundo dela fica onde a amostragem permitiu.
    #
    #     Isso quebra a monotonia que o passo 5 garantiu. O HEC-RAS toma o
    #     MINIMO da secao como fundo, e o minimo passa a variar de secao para
    #     secao: no Trombudo RS 39950 o alvo desce 0,25 m e o leito SOBE 0,29 m,
    #     porque uma secao foi entalhada 1,50 m e a vizinha 0,96 m. Eram 11
    #     contrapendentes nos 12 rios, todos nascidos aqui.
    #
    #     Com o ponto do talvegue cravado em base_p o fundo de toda secao passa
    #     a ser alvo - pilot_prof, que e monotonico por construcao, porque alvo
    #     e monotonico e a profundidade e a mesma.
    z[i0] = base_p

    # 3) parede vertical, so onde o vale e plano demais para conter a cheia.
    #    E o recurso padrao do proprio HEC-RAS ("glass wall"). A altura vem da
    #    VAZAO que a secao conduz, nao de um numero igual para todos.
    precisa = altura_para_vazao(sta, z, d.get("n"), d.get("S_terreno"),
                                vazao_projeto(d["area_km2"]), altura_max)
    minima = float(np.clip(folga_altura * precisa, altura_minima, altura_max))
    util = float(z.max() - z[i0])
    if util < minima:
        topo = z[i0] + minima
        z[0] = max(z[0], topo)
        z[-1] = max(z[-1], topo)
        d["parede"] = round(minima - util, 2)

    d["z"] = z
    d["h_precisa"] = round(precisa, 2)
    # a profundidade que a margem enxerga e a ESCAVACAO -- terreno menos fundo
    # --, nao a diferenca para o ponto vizinho. Com fundo plano o vizinho
    # tambem esta em 'alvo', isso da zero, e as margens sao procuradas DENTRO
    # da propria calha: houve rio saindo com area de calha 0,00 em toda secao.
    d["lb"], d["rb"] = margens(sta, z, i0,
                               max(d.get("z_terreno", z[i0]) - alvo, 0.5))
    return d


def escavar_rio(xs, op, log=print, rotulo=""):
    for d in xs:
        escavar(d, op)
    if xs:
        p = np.array([d.get("parede", 0.0) for d in xs])
        e = np.array([d.get("z_terreno", 0.0) - d.get("z_alvo", 0.0) for d in xs])
        log(f"      {rotulo}: escavacao {np.median(e):.2f} m (mediana), "
            f"parede em {int((p > 0).sum())} de {len(xs)} secoes")
    return xs
