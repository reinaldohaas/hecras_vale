# -*- coding: utf-8 -*-
"""
Perfil longitudinal: a cota do talvegue ao longo do rio.

Trabalha sobre um ESCALAR por secao (z_alvo), nunca sobre a geometria. A secao
so muda em vale/calha.py, no fim, quando o perfil ja esta resolvido. Enquanto
o condicionamento roda, a secao continua sendo o terreno cru.

Tres armadilhas que este arquivo evita, todas medidas:

ESCADA DE POCOS E QUEDAS. Impor declividade minima descendo o rio e maxima
subindo sao dois clamps RIGIDOS, e entre os dois nao sobra valor intermediario:
o perfil sai alternando exatamente os extremos. No Itajai-Mirim, com secoes de
150 m, dava -2,40 / -0,02 / -0,02 / -2,02 / -2,41 -- dezenas de degraus
seguidos, cada um um salto transcritico. As iteracoes do solver batiam o teto
em TODA linha desde o aquecimento.

JOELHO. Limitar a declividade nao impede a mudanca BRUSCA dela: o perfil pode
respeitar os dois limites e ainda passar de 0,07% para 0,70% num vao de 150 m.
Foi o que sobrou no Acu na entrada da garganta do Salto Pilao. Por isso o
alisamento acontece tambem no dominio da declividade.

MONOTONIA POR ULTIMO. O piso de escavacao (o talvegue nao se afasta mais que X
metros do terreno) LEVANTA o leito onde o terreno tem um alto local. Aplicado
depois da monotonia, ele desfaz o que ela garantiu: no Itajai_Sul RS 41.340 o
leito ficou em 375,09 m entre vizinhas em 366,18 e 365,36 -- um corcovo de 9 m
em 300 m, uma barragem dentro do modelo. Entre respeitar o piso e nao ter
contrapendente, nao ter contrapendente ganha: escavar 15 m num ponto e uma
aproximacao, um degrau ao contrario e um erro de fisica.
"""
import numpy as np


def cota(d):
    """Cota alvo do talvegue -- o escalar que carrega a decisao."""
    if "z_alvo" in d:
        return float(d["z_alvo"])
    return float(d["z"][d["i_thal"]])


def teto_declividade(xs, op):
    """Teto de declividade DO RIO, tirado do proprio terreno.

    Um limite unico nao serve: com 0,8% para todos, um afluente de serra que
    desce 9% no terreno recebe um clamp que cava um canion -- houve 279 m de
    escavacao sob o Benedito. O teto sai do percentil 90 da declividade real,
    limitado pela validade de Jarrett.
    """
    S = []
    for a, b in zip(xs, xs[1:]):
        dx = a["rs"] - b["rs"]
        if dx > 0:
            S.append(abs(cota(a) - cota(b)) / dx)
    if not S:
        return op.decl_maxima
    return float(np.clip(np.percentile(S, 90) * 1.2,
                         op.decl_maxima, op.decl_teto))


def aparar_cabeceira(xs, dmax, op, rs_limite=None, rotulo="", log=print):
    """Corta o trecho de torrente no alto do rio.

    Um rio que desce 10% nao e escoamento gradualmente variado, e um 1D nao o
    representa de jeito nenhum. Cortar nao perde agua: a area de drenagem
    aparada continua entrando como vazao lateral nas secoes que sobram.

    O corte procura a secao mais A JUSANTE que ainda esta acima do teto, numa
    janela de algumas secoes -- parar no primeiro par abaixo do teto deixava o
    resto da torrente dentro do modelo, porque um patamar isolado interrompia
    a busca.

    rs_limite: RS da primeira confluencia. O corte nao pode passar dela, senao
    o trecho acima da juncao some inteiro, a juncao fica com um trecho
    entrando e um saindo, e o HEC-RAS recusa a geometria antes de computar.
    """
    if len(xs) < 8:
        return xs
    rs = np.array([d["rs"] for d in xs], float)
    zt = np.array([cota(d) for d in xs], float)
    lim = max(int(op.corte_max_fracao * len(xs)), 0)
    if rs_limite is not None:
        acima = int(np.searchsorted(-rs, -float(rs_limite)))
        lim = min(lim, max(acima - op.secoes_acima_juncao, 0))
    jan = 3
    corte = 0
    for i in range(min(lim, len(xs) - jan - 1)):
        dx = rs[i] - rs[i + jan]
        if dx > 0 and (zt[i] - zt[i + jan]) / dx > dmax:
            corte = i + 1
    if corte:
        log(f"      {rotulo}: aparadas {corte} secoes de cabeceira "
            f"({(xs[0]['rs'] - xs[corte]['rs'])/1000:.1f} km acima de "
            f"{100*dmax:.1f}%)")
    return xs[corte:]


def isotonica(y, w=None):
    """O perfil NAO-CRESCENTE mais proximo do medido (Pool Adjacent Violators).

    Devolve a sequencia nao-crescente que minimiza o desvio quadratico em
    relacao a y. Onde a medicao viola a monotonia, substitui o trecho inteiro
    pela media dele -- nao cava nem aterra sistematicamente, que e a diferenca
    que interessa aqui.
    """
    y = np.asarray(y, float)
    w = np.ones_like(y) if w is None else np.asarray(w, float)
    val, pes, tam = [], [], []
    for i in range(len(y)):
        val.append(y[i])
        pes.append(w[i])
        tam.append(1)
        while len(val) > 1 and val[-2] < val[-1]:       # violou: tem de cair
            v2, p2, t2 = val.pop(), pes.pop(), tam.pop()
            v1, p1, t1 = val.pop(), pes.pop(), tam.pop()
            val.append((v1 * p1 + v2 * p2) / (p1 + p2))
            pes.append(p1 + p2)
            tam.append(t1 + t2)
    out = np.empty_like(y)
    k = 0
    for v, n in zip(val, tam):
        out[k:k + n] = v
        k += n
    return out


def alisar(xs, dmax, op, rotulo="", log=print, n_iter=400):
    """Leva o talvegue medido ao perfil monotonico mais proximo dele.

    ISOTONICA, E NAO ALISAR-E-ACUMULAR. Ate 18/08/2026 este passo era um laco
    de 400 iteracoes que alternava media movel de 3 pontos com
    minimum.accumulate. As duas operacoes so sabem BAIXAR: a media derruba os
    picos e o acumulado toma o minimo corrente, e 400 ciclos disso afundam o
    leito. Medido nos doze rios, o leito acabava 4,76 m abaixo do talvegue na
    mediana e 24,43 m no pior caso -- no Iraputa, 13,5 m de escavacao e 102 m
    de largura para um rio cujo pico e 10,65 m3/s. Um canion, nao um canal.

    A causa e o talvegue de MDS ser serrilhado: no Iraputa ele SOBE rio abaixo
    em 67 dos 151 vaos, porque ao amostrar o fundo do vale a linha pula entre
    copa de mata e clareira. A monotonia tinha de eliminar as 67 subidas, e o
    unico jeito que o laco conhecia era baixar tudo que vinha depois.

    A regressao isotonica resolve o mesmo problema pelo criterio certo: o
    perfil nao-crescente MAIS PROXIMO do medido. Onde o terreno viola a
    monotonia ela promedia o trecho em vez de rebaixar o rio inteiro. Mesmos
    doze rios: mediana 4,76 -> 0,00 m, p90 11,40 -> 0,78 m, max 24,43 -> 9,11 m.

    Depois dela vem, nesta ordem: declividade minima (a isotonica deixa
    patamares exatamente planos onde agrupou, e patamar plano nao escoa), teto
    de declividade, alisamento da DECLIVIDADE (que e o que tira o joelho da
    entrada de garganta) e a monotonia por ultimo.
    """
    if len(xs) < 5:
        return
    rs = np.array([d["rs"] for d in xs], float)
    t = rs[0] - rs                              # distancia rio abaixo
    z = np.array([cota(d) for d in xs], float)
    # z_terreno, e nao z[i_thal]: em secao cujo talvegue foi rejeitado por
    # absurdo, z_terreno ja e o valor corrigido
    terreno = np.array([float(d.get("z_terreno", d["z"][d["i_thal"]]))
                        for d in xs], float)
    piso = terreno - op.escavacao_maxima
    dmin = op.decl_minima

    # Os dois limites de declividade sao a MESMA isotonica em coordenada
    # deslocada, e e isso que evita rebaixar o rio:
    #     declividade >= dmin   <=>   (z + dmin*t) nao-crescente
    #     declividade <= dmax   <=>   (z + dmax*t) nao-decrescente
    # O minimum.accumulate que estava aqui devolve o maior perfil VIAVEL abaixo
    # do medido -- sempre baixando, e a queda se acumula: com dmin de 0,01% sao
    # 4,7 m ao longo do Iraputa e 18,7 m ao longo do Acu, so de impor
    # declividade minima. A isotonica devolve o perfil viavel MAIS PROXIMO, que
    # inclina o patamar em torno da media dele em vez de ancorar no topo.
    # (a primeira ja garante monotonia: se z + dmin*t nao cresce e t cresce,
    #  entao z decresce)
    z = isotonica(z + dmin * t) - dmin * t                     # decl. minima
    z = -isotonica(-(z + dmax * t)) - dmax * t                 # teto
    z = np.maximum(z, piso)

    dt = np.diff(t)
    if (dt > 0).sum() >= 5:
        decl = np.zeros_like(dt)
        ok = dt > 0
        decl[ok] = (z[:-1] - z[1:])[ok] / dt[ok]
        jan = 7
        k = np.ones(jan) / jan
        peso = np.convolve(np.pad(dt, jan // 2, mode="edge"), k, "valid")
        suave = np.convolve(np.pad(decl * dt, jan // 2, mode="edge"), k, "valid")
        decl = np.clip(np.divide(suave, peso, out=decl.copy(), where=peso > 0),
                       dmin, dmax)
        # A QUEDA TOTAL TEM DE SER A MESMA -- e a docstring sempre disse que
        # era, mas nao era: o clip em [dmin, dmax] muda a media das
        # declividades, e como a reintegracao anda da foz para a cabeceira o
        # vies se acumula secao a secao. Em rio curto nao aparece; no Itajai do
        # Norte, com 340 secoes em 134 km, empurrava o leito 6 m para baixo do
        # terreno, e no Acu 3 m em 188 km. Reescalar devolve exatamente a queda
        # que a isotonica escolheu e deixa o alisamento fazer so o que promete:
        # REDISTRIBUIR a queda, tirando o joelho, sem mudar o total.
        queda = float(z[0] - z[-1])
        soma = float(np.sum(decl * dt))
        if soma > 1e-9 and queda > 1e-9:
            decl = np.clip(decl * (queda / soma), dmin, dmax)
        # reintegra a partir da FOZ que a isotonica devolveu, e nao mais de um
        # z_foz capturado antes de tudo -- era esse pino que punha a ultima
        # secao acima da vizinha e criava contrapendente em sete rios
        zn = np.empty_like(z)
        zn[-1] = z[-1]
        for i in range(len(z) - 2, -1, -1):
            zn[i] = zn[i + 1] + decl[i] * dt[i]
        z = np.maximum(np.minimum(zn, terreno), piso)

    # Declividade minima por ULTIMO, e sem reaplicar o piso depois dela -- o
    # piso levanta o leito onde o terreno tem um alto local, e aplicado depois
    # desfaz o que a monotonia garantiu (foi assim que nasceu o corcovo de 9 m
    # no Itajai_Sul RS 41.340).
    #
    # A FOZ ENTRA NA MONOTONIA. Ate 18/08/2026 havia aqui um "z[-1] = z_foz"
    # que devolvia a foz ao valor CRU do terreno, capturado antes de qualquer
    # alisamento, depois que a monotonia ja tinha baixado o resto do rio. O
    # ultimo vao herdava toda a diferenca e o leito SUBIA rio abaixo: 4,40 m no
    # Acu, 3,39 m no Taio, 1,22 m no Iraputa -- em sete dos doze rios, sempre
    # na ultima secao (RS 75), e sempre exatamente na secao onde o solver
    # depois acusava o maior erro de nivel (18,17 m no Iraputa RS 75).
    # Os dois limites, ALTERNADOS ate valerem juntos. Aplicar um depois do
    # outro uma vez so nao basta: a projecao de declividade minima agrupa
    # patamares e, ao faze-lo, pode reabrir um degrau acima do teto. Medido
    # antes de alternar: Benedito com 11,43% contra teto de 5,00%, Taio com
    # 9,90% contra 3,74%. Sao conjuntos convexos, entao alternar converge.
    for _ in range(30):
        z = isotonica(z + dmin * t) - dmin * t
        z = -isotonica(-(z + dmax * t)) - dmax * t
        dt_ = t[1:] - t[:-1]
        s_ = np.where(dt_ > 0, (z[:-1] - z[1:]) / np.maximum(dt_, 1e-9), 0.0)
        if s_.size == 0 or (s_.max() <= dmax * 1.001
                            and s_.min() >= dmin * 0.999):
            break
    z = isotonica(z + dmin * t) - dmin * t      # monotonia garantida no fim
    fundo = int(np.sum(z < piso - 0.01))
    if fundo:
        log(f"      {rotulo}: {fundo} secoes abaixo do piso de "
            f"{op.escavacao_maxima:.0f} m, para nao criar contrapendente")
    for d, zi in zip(xs, z):
        d["z_alvo"] = float(zi)


def rejeitar_absurdos(xs, rotulo="", log=print, janela=7, fator=8.0):
    """Tira do talvegue os pontos que nao podem ser terreno.

    UM ponto ruim destroi o rio inteiro, e nao e exagero: no Itajai do Oeste
    uma unica secao leu 0,02 m entre vizinhas de 390 m: um vazio do MDT que a
    reamostragem deixou perto de zero em vez de NoData. A imposicao de
    monotonia obriga cada secao a ficar abaixo da anterior, entao esse 0,02
    arrastou os 94 km seguintes para baixo dele -- o rio saiu com leito de
    -9,11 m na foz e escavacao mediana de 344 m. Duas secoes de 156 bastaram.

    O criterio e ROBUSTO, nao um limiar fixo: compara-se com a mediana movel e
    com a dispersao tipica do proprio rio (desvio absoluto mediano). Assim vale
    tanto para um rio de planicie quanto para um de serra, sem ajuste.
    """
    if len(xs) < janela:
        return xs
    z = np.array([cota(d) for d in xs], float)
    k = janela // 2
    pad = np.pad(z, k, mode="edge")
    med = np.array([np.median(pad[i:i + janela]) for i in range(len(z))])
    desvio = np.abs(z - med)
    mad = float(np.median(desvio)) or 1.0
    fora = desvio > max(fator * mad, 10.0)
    if fora.any():
        for i in np.flatnonzero(fora):
            d = xs[i]
            d["z_rejeitado"] = float(z[i])
            # os TRES lugares que dependem do talvegue, senao a correcao se
            # desfaz: alisar() rele o terreno cru e o usa como TETO
            # (z = clip(m, piso, terreno)), reintroduzindo o valor absurdo
            # tres linhas depois de ele ter sido removido; e a propria secao
            # fica com um furo de um ponto.
            d["z_alvo"] = float(med[i])
            d["z_terreno"] = float(med[i])
            d["z"] = np.maximum(np.asarray(d["z"], float), float(med[i]))
        pior = int(np.argmax(desvio))
        log(f"      {rotulo}: {int(fora.sum())} secoes com talvegue impossivel "
            f"substituidas pela mediana local (a pior: RS {xs[pior]['rs']:.0f}, "
            f"{z[pior]:.2f} m entre vizinhas de {med[pior]:.2f} m)")
    return xs


def condicionar(xs, op, rs_limite=None, rotulo="", log=print):
    """Talvegue monotonico, com declividade e curvatura sob controle."""
    if len(xs) < 3:
        return xs
    # SEMPRE do terreno, nunca do z_alvo anterior. Com setdefault, rodar o
    # passo de novo partia do resultado da vez passada: um perfil ja arruinado
    # entrava liso no filtro de absurdos, que nao achava nada para rejeitar --
    # e a correcao parecia nao funcionar. Num programa feito para parar,
    # conferir e retomar, cada passo tem de dar o mesmo resultado a partir da
    # mesma entrada.
    for d in xs:
        d["z_alvo"] = float(d["z"][d["i_thal"]])
        d.pop("z_rejeitado", None)
    # ANTES de qualquer clamp: um valor absurdo propaga por todo o rio
    xs = rejeitar_absurdos(xs, rotulo, log)
    dmax = teto_declividade(xs, op)
    if dmax > op.decl_maxima * 1.05:
        log(f"      {rotulo}: rio de serra, teto {100*dmax:.2f}% "
            f"(padrao {100*op.decl_maxima:.1f}%)")
    xs = aparar_cabeceira(xs, dmax, op, rs_limite, rotulo, log)
    alisar(xs, dmax, op, rotulo, log)
    return xs


def ancorar(xs, cota_foz, log=print, rotulo="", op=None):
    """Casa a foz do afluente com o leito do receptor.

    So a CALHA se move, e o deslocamento e afunilado -- cheio na foz, zero na
    cabeceira. Deslocar a secao inteira preserva o relevo relativo e afunda o
    ABSOLUTO: no Mirim isso punha o ponto mais alto de secoes de 1.474 m de
    largura a 2 m ABAIXO do nivel do mar, e o HEC-RAS nao conseguia sequer
    estabelecer a lamina inicial.
    """
    if len(xs) < 2 or cota_foz is None:
        return xs
    delta = float(cota_foz) - cota(xs[-1])
    if abs(delta) < 0.01:
        return xs
    rs = np.array([d["rs"] for d in xs], float)
    peso = (rs[0] - rs) / max(rs[0] - rs[-1], 1e-6)      # 0 no topo, 1 na foz
    for d, w in zip(xs, peso):
        d["z_alvo"] = cota(d) + delta * float(w)

    # O AFUNILAMENTO PODE INVERTER A DECLIVIDADE. O deslocamento cresce rio
    # abaixo; somado a um perfil que ja descia pouco, ele levanta a parte de
    # jusante acima da de montante. Enquanto o leito era escavado 13 m isso nao
    # aparecia -- a foz ficava tao abaixo do receptor que o delta era grande e
    # negativo. Com o leito perto do terreno, o Iraputa passou a ter 16
    # contrapendentes so por causa desta funcao. Reprojeta impondo declividade
    # minima com a FOZ PRESA (peso alto no ultimo ponto): e a mesma isotonica,
    # e devolve o perfil mais proximo que respeita as duas coisas.
    if op is not None and len(xs) >= 3:
        t = rs[0] - rs
        z = np.array([cota(d) for d in xs], float)
        w = np.ones(len(z))
        w[-1] = 1e6                       # a foz nao se move: e o no da rede
        z = isotonica(z + op.decl_minima * t, w) - op.decl_minima * t
        for d, zi in zip(xs, z):
            d["z_alvo"] = float(zi)
    log(f"      {rotulo}: foz ancorada em {cota_foz:.2f} m "
        f"(deslocamento {delta:+.2f} m, afunilado ate a cabeceira)")
    return xs


def manning(xs, op, razao_planicie=1.8):
    """Rugosidade por Jarrett (1984): n = 0,39 * S^0,38 * R^-0,16.

    Usa a declividade do TERRENO, guardada no corte, e nao a do perfil
    condicionado -- esta ultima e a do clamp, e subestima a rugosidade
    justamente nos afluentes de serra.
    Fora da faixa de validade (0,002 a 0,052) cai para valores tabelados.
    """
    for d in xs:
        S = float(d.get("S_terreno") or 0.0)
        R = max(0.3 * (max(d["area_km2"], 1.0) ** 0.15), 0.5)
        if 0.002 <= S <= 0.052:
            n = 0.39 * S ** 0.38 * R ** -0.16
        elif S > 0.052:
            n = 0.075
        else:
            n = 0.035
        d["n"] = round(float(np.clip(n, 0.028, 0.090)), 4)
        d["n_planicie"] = round(d["n"] * razao_planicie, 3)
    return xs
