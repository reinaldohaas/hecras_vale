# -*- coding: utf-8 -*-
"""
Correcao GEOMETRICA da secao. Nunca por cota.

Ponto fundamental do programa: quando uma secao esta ruim, a tentacao e mexer
nas elevacoes -- rebaixar um ponto, cavar um canal, levantar as pontas. Isso
apaga o problema do relatorio e o mantem no modelo. Aqui a cota vem sempre do
terreno; o que se altera e ONDE se corta.

Duas estrategias, nesta ordem:

  ESTENDER   alonga a secao para os dois lados e reamostra. Resolve o caso mais
             comum -- secao curta demais, com o canal caindo perto da ponta --
             preservando a posicao e a orientacao originais, que muitas vezes
             foram escolhidas por alguem que conhecia o rio.

  PERPENDICULAR   se estender nao resolve, o problema nao e comprimento: e
             posicao ou angulo. Ai se acha o cruzamento com o eixo, toma-se a
             direcao local do rio e gera-se uma secao nova, centrada no eixo e
             perpendicular a ele.

Nenhuma das duas escreve nada. Elas devolvem uma proposta com o proprio QC, e
quem aceita e o usuario.
"""
import numpy as np

from . import qc as _qc

PASSOS_EXTENSAO = (10.0, 20.0, 30.0, 50.0, 100.0)


def _prolongar(linha, dist_esq, dist_dir):
    """Estende a LineString nas duas pontas, mantendo a direcao das pontas."""
    c = list(linha.coords)
    from shapely.geometry import LineString
    if len(c) < 2:
        return linha
    (x0, y0), (x1, y1) = c[0], c[1]
    n = float(np.hypot(x1 - x0, y1 - y0)) or 1.0
    a = (x0 - dist_esq * (x1 - x0) / n, y0 - dist_esq * (y1 - y0) / n)
    (xa, ya), (xb, yb) = c[-2], c[-1]
    n = float(np.hypot(xb - xa, yb - ya)) or 1.0
    b = (xb + dist_dir * (xb - xa) / n, yb + dist_dir * (yb - ya) / n)
    return LineString([a] + c + [b])


def estender(secao, dem, eixo, lim, passos=PASSOS_EXTENSAO, simetrico=True):
    """Tenta resolver alongando. Devolve a primeira proposta aceitavel.

    Se nenhuma resolver, devolve a melhor delas mesmo assim -- com o QC dela,
    para o usuario ver que nao resolveu, em vez de nao receber nada.
    """
    tentativas = []
    for d in passos:
        de, dd = (d, d) if simetrico else _assimetrico(secao, d)
        nova = secao.copia_com_geometria(_prolongar(secao.geom, de, dd),
                                         f"estendida +{de:.0f}/+{dd:.0f} m")
        nova.extrair(dem, lim.espacamento, eixo, lim.proeminencia_min)
        _qc.avaliar(nova, (), lim)
        tentativas.append(nova)
        if nova.qc.status == _qc.OK:
            return nova, tentativas
    melhor = max(tentativas, key=lambda s: s.qc.nota) if tentativas else None
    return melhor, tentativas


def _assimetrico(secao, d):
    """Estende so do lado em que o talvegue esta encostado."""
    p = secao.posicao_relativa
    if not np.isfinite(p):
        return d, d
    return (2.0 * d, 0.0) if p < 0.5 else (0.0, 2.0 * d)


def perpendicular(secao, dem, eixo, lim, largura=None):
    """Gera uma secao nova, centrada no eixo e perpendicular a ele."""
    if eixo is None:
        return None
    p, dist, linha = eixo.ponto_no_eixo(secao.geom)
    meia = 0.5 * float(largura if largura else secao.largura)
    geom = eixo.perpendicular(p, meia, meia, linha)
    nova = secao.copia_com_geometria(
        geom, f"perpendicular ao eixo (largura {2*meia:.0f} m)")
    nova.extrair(dem, lim.espacamento, eixo, lim.proeminencia_min)
    _qc.avaliar(nova, (), lim)
    return nova


def propor(secao, dem, eixo, lim, passos=PASSOS_EXTENSAO):
    """A melhor proposta disponivel, com todas as candidatas para comparar.

    Ordem deliberada: estender antes de recriar. Recriar perpendicular joga
    fora a posicao original da secao, que pode ter sido escolhida por um motivo
    que o programa nao ve -- uma ponte, um estreitamento, uma estacao
    fluviometrica.
    """
    candidatas = []
    melhor_ext, tentativas = estender(secao, dem, eixo, lim, passos)
    candidatas.extend(tentativas)
    if melhor_ext is not None and melhor_ext.qc.status == _qc.OK:
        return melhor_ext, candidatas
    perp = perpendicular(secao, dem, eixo, lim)
    if perp is not None:
        candidatas.append(perp)
        # perpendicular mais larga, se a de mesma largura ainda nao resolveu
        if perp.qc.status != _qc.OK:
            p2 = perpendicular(secao, dem, eixo, lim,
                               largura=1.5 * secao.largura)
            if p2 is not None:
                candidatas.append(p2)
    if not candidatas:
        return None, []
    melhor = max(candidatas, key=lambda s: s.qc.nota)
    return (melhor if melhor.qc.nota > (secao.qc.nota if secao.qc else 0.0)
            else None), candidatas


def comparar(original, proposta):
    """Tabela lado a lado das duas, para a tela de comparacao."""
    def linha(s):
        return {"origem": s.origem,
                "largura": round(s.largura, 1),
                "talvegue_pct": (round(100 * s.posicao_relativa, 1)
                                 if np.isfinite(s.posicao_relativa) else None),
                "z_talvegue": (round(s.z_talvegue, 2)
                               if np.isfinite(s.z_talvegue) else None),
                "profundidade": (round(s.profundidade_relativa, 2)
                                 if np.isfinite(s.profundidade_relativa) else None),
                "orientacao": (round(s.azimute, 1) if s.azimute is not None
                               else None),
                "status": s.qc.status if s.qc else "?",
                "qc": s.qc.nota if s.qc else 0.0}
    return {"original": linha(original), "proposta": linha(proposta)}
