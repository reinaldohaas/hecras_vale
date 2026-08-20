# -*- coding: utf-8 -*-
"""
Escava a calha NO TERRENO, pelo mecanismo do proprio HEC-RAS.

O QUE MUDA. Ate aqui a calha era imposta secao a secao (vale/calha.py): cada
secao recebia um trapezio e um entalhe piloto, desenhados nos pontos dela. Isso
poe a calha na GEOMETRIA e deixa o TERRENO intacto, e as duas coisas passam a
discordar -- o RAS Mapper mostra um vale sem rio, e a mancha de inundacao sai
de uma lamina que nao tem leito por baixo.

Pior que isso, e medido: o entalhe so alcanca o fundo se algum PONTO da secao
cair dentro da meia-largura dele. Com entalhe de 3 m e pontos a cada 5 m,
nenhum cai. Em 10 das 8.956 secoes o fundo ficou onde a amostragem permitiu, e
como o HEC-RAS toma o MINIMO da secao como cota de fundo, essas dez viraram
contrapendente num perfil que o passo 5 entregou perfeitamente monotonico.

`RasTerrainModWriter.add_channel_modification` resolve pela raiz: escreve um
canal trapezoidal no HDF do TERRENO, em modo TakeLower (so rebaixa onde o canal
fica abaixo do chao), e registra no .rasmap. Quem corta a secao depois ja
encontra o canal la -- nao ha o que entalhar, nao ha ponto para cravar, e a
consistencia entre terreno e geometria e por construcao.

A ORDEM PASSA A IMPORTAR, e e o preco:

    1. amostrar o talvegue ao longo do EIXO         (terreno original)
    2. condicionar o perfil                          (vale/perfil.py)
    3. escavar o canal NO TERRENO                    (este modulo)
    4. cortar as secoes do terreno JA MODIFICADO     (vale/secoes.py)

Inverter 3 e 4 nao adianta: as secoes leriam o terreno sem canal.

PROFUNDIDADE CONSTANTE POR MODIFICACAO. `add_channel_modification` recebe UM
`depth` e UMA `width` para a polilinha inteira. O nosso perfil condicionado
varia os dois ao longo do rio -- a profundidade e `z_terreno - z_alvo`, que
muda seção a seção, e a largura vem de Leopold. Por isso o rio e partido em
SEGMENTOS onde os dois sao aproximadamente constantes, e cada segmento vira uma
modificacao. Segmento demais custa tempo de escrita; segmento de menos achata o
perfil. `tolerancia_prof` controla o corte.
"""
import os

import numpy as np


def segmentar(xs, op, log=print):
    """Parte o rio em trechos de profundidade e largura quase constantes.

    Devolve lista de dicts com `s0`, `s1` (chainage), `prof`, `larg` e os
    indices das secoes. O criterio e a PROFUNDIDADE, porque e ela que decide a
    cota do leito; a largura acompanha porque as duas crescem com a area.
    """
    tol = float(getattr(op, "tolerancia_prof", 0.5))
    w = sorted(xs, key=lambda d: -d["rs"])           # cabeceira -> foz
    prof = np.array([max(float(d.get("z_terreno", 0.0))
                         - float(d.get("z_alvo", 0.0)), 0.0) for d in w])
    # Leopold, os MESMOS coeficientes que a calha por secao usava
    larg = np.array([op.canal_kw * max(float(d.get("area_km2", 1.0)), 1.0)
                     ** op.canal_ew for d in w], float)
    # QUANTIZAR, e nao cortar quando a faixa estoura. A primeira versao fechava
    # o segmento assim que a amplitude passava da tolerancia e recomecava ali:
    # com a profundidade oscilando em torno de um valor, isso picota o rio em
    # segmentos de duas secoes -- 104 so no Cedros, o que daria mais de mil
    # modificacoes nos doze rios. Quantizando, secoes vizinhas caem na mesma
    # faixa e o segmento cresce naturalmente.
    faixa = np.round(prof / max(tol, 1e-6)).astype(int)
    faixa = np.round(np.convolve(faixa, np.ones(5) / 5.0, "same")).astype(int)
    corte = np.flatnonzero(np.diff(faixa) != 0) + 1
    segs = []
    for ini, fim in zip(np.r_[0, corte], np.r_[corte, len(w)]):
        if fim - ini < 2:
            continue
        segs.append({"i0": int(ini), "i1": int(fim - 1),
                     "prof": float(np.median(prof[ini:fim])),
                     "larg": float(np.median(larg[ini:fim])),
                     "xs": w[ini:fim]})
    log(f"      {len(segs)} segmento(s) de calha "
        f"(profundidade {prof.min():.2f} a {prof.max():.2f} m, "
        f"tolerancia {tol:.2f} m)")
    return segs


def _polilinha(eixo, s0, s1, passo=50.0):
    """Pontos do eixo entre duas chainages, para a polilinha da modificacao."""
    n = max(int((s1 - s0) / passo) + 1, 2)
    return np.array([[p.x, p.y] for p in
                     (eixo.interpolate(float(s))
                      for s in np.linspace(s0, s1, n))], float)


def _corredor(eixo, xs, op):
    """Poligono do canal e a cota do LEITO em cada vertice do contorno.

    O anel desce pela margem esquerda, da cabeceira para a foz, e volta pela
    direita. A cota de cada vertice e o `z_alvo` da secao correspondente -- a
    mesma dos dois lados, porque o fundo do canal e plano na transversal e a
    inclinacao transversal fica por conta do que o terreno ja tem fora dele.

    A LARGURA VARIA ao longo do rio (Leopold), e por isso o poligono e
    construido ponto a ponto em vez de ser um buffer de largura unica.
    """
    L = eixo.length
    w = sorted(xs, key=lambda d: -d["rs"])           # cabeceira -> foz
    esq, dir_, z_esq, z_dir = [], [], [], []
    for d in w:
        s = L - float(d["rs"])
        s = float(np.clip(s, 0.0, L))
        meia = 0.5 * max(op.canal_kw * max(float(d.get("area_km2", 1.0)), 1.0)
                         ** op.canal_ew, op.pilot_largura_min)
        a = eixo.interpolate(max(s - 25.0, 0.0))
        b = eixo.interpolate(min(s + 25.0, L))
        tx, ty = b.x - a.x, b.y - a.y
        n = float(np.hypot(tx, ty)) or 1.0
        rx, ry = -ty / n, tx / n                      # normal a esquerda
        p = eixo.interpolate(s)
        z = float(d.get("z_alvo", d.get("z_terreno", 0.0)))
        esq.append((p.x + meia * rx, p.y + meia * ry))
        dir_.append((p.x - meia * rx, p.y - meia * ry))
        z_esq.append(z)
        z_dir.append(z)
    # ANEL FECHADO: o primeiro ponto repetido no fim, e a cota dele junto. A
    # funcao valida `boundary_elevations` contra a contagem do poligono JA
    # FECHADO, entao entregar o anel aberto falha com "must be a 1D array
    # matching the closed polygon boundary point count".
    anel = esq + dir_[::-1] + [esq[0]]
    cotas = z_esq + z_dir[::-1] + [z_esq[0]]
    return np.array(anel, float), np.array(cotas, float)


def escavar_por_poligono(op, eixos, xs_por_rio, terreno_hdf, rasmap,
                         log=print):
    """Escava a calha com as COTAS do perfil condicionado, nao com profundidade.

    `add_channel_modification` escava profundidade constante ABAIXO DO TERRENO:
    o fundo acompanha as ondulacoes do DEM e o perfil monotonico que o passo 5
    calculou se perde -- e e justamente ele que impede leito subindo rio
    abaixo. Segmentar nao resolve, porque dentro do segmento o fundo continua
    seguindo o chao.

    `add_modification_polygon` com `elevation_method='provided'` recebe a cota
    de cada vertice do contorno. Entao o canal sai EXATAMENTE em `z_alvo`, e o
    condicionamento sobrevive.

    `mode='take_lower'` de proposito: onde o terreno ja e mais fundo que o
    perfil, ele fica -- escavar so baixa, nunca aterra. Aterro dentro do rio
    foi um dos defeitos que esta reconstrucao passou o dia removendo.

    Um poligono por RIO. Sao centenas de vertices num anel longo e estreito, e
    e uso fora do caso tipico da funcao (acude, area umida), mas a estrutura e
    a mesma: contorno com cotas.
    """
    from ras_commander.terrain import RasTerrainModWriter as W

    if not os.path.exists(terreno_hdf):
        raise SystemExit(
            f"nao ha terreno em {terreno_hdf}. A calha e escrita NO TERRENO "
            f"agora; sem ele nao ha o que escavar.")
    feitas = []
    for d in eixos:
        ras = d["ras"]
        xs = xs_por_rio.get(ras) or []
        if len(xs) < 3:
            continue
        anel, cotas = _corredor(d["linha"], xs, op)
        nome = f"Calha_{ras}"[:32]
        W.add_modification_polygon(
            terrain_hdf_path=terreno_hdf, name=nome, polygon_coords=anel,
            elevation_method="provided", boundary_elevations=cotas,
            mode="take_lower", rasmap_path=rasmap, group_name="Calha")
        feitas.append({"rio": ras, "nome": nome, "vertices": len(anel),
                       "cota_min": float(cotas.min()),
                       "cota_max": float(cotas.max())})
        log(f"      {ras:<16} poligono de {len(anel)} vertices, "
            f"leito de {cotas.max():.1f} a {cotas.min():.1f} m")
    log(f"      {len(feitas)} calha(s) escritas em "
        f"{os.path.basename(terreno_hdf)}")
    return feitas


def escavar(op, eixos, xs_por_rio, terreno_hdf, rasmap, log=print):
    """Escreve as modificacoes de canal no HDF do terreno.

    Devolve a lista do que foi escrito, para a auditoria. NAO engole excecao:
    se a escrita falhar, o terreno fica sem canal e todo o resto do modelo
    seria construido sobre um vale sem rio -- e isso ja aconteceu nesta
    reconstrucao com o entalhe, silenciosamente.
    """
    from ras_commander.terrain import RasTerrainModWriter as W

    if not os.path.exists(terreno_hdf):
        raise SystemExit(
            f"nao ha terreno em {terreno_hdf}. A calha e escrita NO TERRENO "
            f"agora; sem ele nao ha o que escavar. Rode o passo 3 com "
            f"terreno_hdf=true.")

    feitas = []
    for d in eixos:
        ras = d["ras"]
        xs = xs_por_rio.get(ras) or []
        if len(xs) < 3:
            continue
        eixo = d["linha"]
        L = eixo.length
        for k, seg in enumerate(segmentar(xs, op, log), 1):
            w = seg["xs"]
            s0 = L - float(w[0]["rs"])
            s1 = L - float(w[-1]["rs"])
            if s1 - s0 < op.espacamento_piso:
                continue
            pts = _polilinha(eixo, s0, s1, passo=op.passo_polilinha)
            nome = f"{ras}_{k:02d}"[:32]
            W.add_channel_modification(
                terrain_hdf_path=terreno_hdf, rasmap_path=rasmap, name=nome,
                polyline_points=pts,
                width=max(seg["larg"], op.pilot_largura_min),
                depth=max(seg["prof"], op.prof_minima_canal),
                left_slope=op.talude_canal, right_slope=op.talude_canal,
                max_extent=max(4.0 * seg["larg"], 200.0),
                group_name="Calha")
            feitas.append({"rio": ras, "nome": nome, "km0": s0 / 1000.0,
                           "km1": s1 / 1000.0, "prof": seg["prof"],
                           "larg": seg["larg"], "pontos": len(pts)})
        log(f"      {ras:<16} {sum(1 for f in feitas if f['rio'] == ras)} "
            f"modificacao(oes) de canal")
    log(f"      total: {len(feitas)} modificacoes escritas em "
        f"{os.path.basename(terreno_hdf)}")
    return feitas


def conferir(op, eixos, terreno_hdf, rasmap, log=print, n=3):
    """Le o terreno de volta e mede o canal que foi escrito.

    A biblioteca tem `compare_before_after_profiles` para exatamente isto --
    quem confere e ela, e nao a nossa aritmetica. Sem esta leitura, "escrevi a
    modificacao" e afirmacao, nao medida: este projeto ja teve um symlink de
    terreno que existia, tinha 14 bytes e nao levava a lugar nenhum.
    """
    from ras_commander.terrain import RasTerrainModWriter as W

    try:
        mods = W.list_modifications(terreno_hdf)
    except Exception as e:                                   # noqa: BLE001
        log(f"      nao consegui listar as modificacoes: {e}")
        return []
    log(f"      {len(mods)} modificacao(oes) no terreno")
    saida = []
    for d in eixos[:n]:
        eixo = d["linha"]
        pts = _polilinha(eixo, 0.0, eixo.length,
                         passo=max(eixo.length / 200.0, 50.0))
        try:
            cmp = W.compare_before_after_profiles(
                terrain_hdf_path=terreno_hdf, rasmap_path=rasmap,
                x_coords=list(pts[:, 0]), y_coords=list(pts[:, 1]))
            saida.append((d["ras"], cmp))
            log(f"      {d['ras']:<16} perfil conferido "
                f"({len(cmp)} pontos)")
        except Exception as e:                               # noqa: BLE001
            log(f"      {d['ras']:<16} conferencia falhou: {e}")
    return saida
