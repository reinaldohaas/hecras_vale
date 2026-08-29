# -*- coding: utf-8 -*-
"""
Controle de qualidade geometrico das secoes.

Cada teste responde a uma pergunta unica e diz por que reprovou, em numeros.
Um "CRITICA" sem motivo mensuravel nao ajuda ninguem a decidir.

Os testes:

  A  TALVEGUE NA EXTREMIDADE   posicao relativa do talvegue na secao. E o
                               criterio principal: canal na borda significa
                               que a secao nao cobre a planicie de um dos
                               lados, e a lamina encontra uma parede.
  B  CANAL NAO IDENTIFICADO    profundidade relativa (media das extremidades
                               menos o talvegue). Pequena demais: ou a secao
                               esta fora do vale, ou o DEM nao resolve o canal.
  C  SALTO DE TALVEGUE         cota do talvegue contra as vizinhas. Um degrau
                               isolado e quase sempre secao mal posicionada,
                               nao relevo.
  D  LARGURA ANORMAL           largura contra a mediana das vizinhas.
  E  ORIENTACAO                angulo com o eixo do rio. Secao obliqua mede
                               largura maior que a real.

O programa NAO existe para deixar tudo verde. Quando o terreno nao permite
concluir, o resultado e INCERTO -- revisao manual --, que e diferente de OK e
diferente de CRITICA.
"""
from dataclasses import dataclass, field, asdict

import numpy as np

OK = "OK"
ATENCAO = "ATENCAO"
CRITICA = "CRITICA"
INCERTO = "INCERTO"

ORDEM = {OK: 0, ATENCAO: 1, INCERTO: 2, CRITICA: 3}


@dataclass
class Limiares:
    """Tudo que o usuario pode ajustar na interface."""
    pos_ok_min: float = 0.20          # talvegue entre 20% e 80% -> OK
    pos_ok_max: float = 0.80
    pos_atencao_min: float = 0.10     # entre 10% e 20% (ou 80-90%) -> ATENCAO
    pos_atencao_max: float = 0.90
    profundidade_min: float = 1.0     # m, abaixo disso: canal nao identificado
    proeminencia_min: float = 0.5     # m, limiar de deteccao do talvegue
    salto_talvegue: float = 3.0       # m contra a mediana das vizinhas
    razao_largura: float = 2.0        # x a mediana das vizinhas
    desvio_ortogonal: float = 25.0    # graus fora de 90
    espacamento: float = 2.0          # m, amostragem do DEM
    n_vizinhas: int = 4               # quantas de cada lado nos testes C e D

    def dict(self):
        return asdict(self)


@dataclass
class Resultado:
    status: str = OK
    nota: float = 100.0
    motivos: list = field(default_factory=list)
    testes: dict = field(default_factory=dict)

    @property
    def resumo(self):
        return "; ".join(self.motivos) if self.motivos else "sem problemas"


def _pior(a, b):
    return a if ORDEM[a] >= ORDEM[b] else b


def teste_a_posicao(s, lim):
    """Talvegue na extremidade."""
    p = s.posicao_relativa
    if not np.isfinite(p):
        return INCERTO, "posicao do talvegue indefinida", 0.0
    if p < lim.pos_atencao_min or p > lim.pos_atencao_max:
        return CRITICA, f"talvegue a {100*p:.1f}% da largura", 0.0
    if p < lim.pos_ok_min or p > lim.pos_ok_max:
        return ATENCAO, f"talvegue a {100*p:.1f}% da largura", 0.5
    # nota cheia no centro, caindo ate a faixa de atencao
    d = abs(p - 0.5) / max(0.5 - lim.pos_ok_min, 1e-6)
    return OK, "", float(np.clip(1.0 - 0.4 * d, 0.0, 1.0))


def teste_b_canal(s, lim):
    """Canal identificavel?"""
    pr = s.profundidade_relativa
    inc = (s.talvegue or {}).get("incerto", False)
    if not np.isfinite(pr):
        return INCERTO, "profundidade relativa indefinida", 0.0
    if pr < lim.profundidade_min:
        return INCERTO, (f"CANAL NAO IDENTIFICADO (profundidade relativa "
                         f"{pr:.2f} m)"), 0.0
    if inc:
        return INCERTO, "canal raso para o limiar de proeminencia", 0.3
    return OK, "", float(np.clip(pr / (4.0 * lim.profundidade_min), 0.3, 1.0))


def teste_c_salto(s, vizinhas, lim):
    """Salto de cota do talvegue contra a TENDENCIA local, nao a mediana.

    Comparar com a mediana das vizinhas so funciona em rio plano. Num trecho
    de serra o leito cai varios metros entre secoes vizinhas por fisica, e o
    teste acusava quase toda secao: no Itajai-Mirim, com 150 m de espacamento
    e 1,6% de declividade, a queda legitima entre vizinhas ja e 2,4 m.

    O que se procura e o ponto FORA DA LINHA -- aquele que destoa da propria
    tendencia do trecho. Entao ajusta-se uma reta pelas vizinhas (cota contra
    River Station) e mede-se o residuo desta secao em relacao a ela.
    """
    pares = [(v.rs, v.z_talvegue) for v in vizinhas
             if v is not s and np.isfinite(v.z_talvegue)
             and isinstance(v.rs, float) and np.isfinite(v.rs)]
    if len(pares) < 3 or not np.isfinite(s.z_talvegue) or \
            not isinstance(s.rs, float) or not np.isfinite(s.rs):
        return OK, "", 1.0
    x = np.array([p[0] for p in pares], float)
    y = np.array([p[1] for p in pares], float)
    # Reta por minimos quadrados na forma fechada, de proposito: np.polyfit
    # chama o LAPACK, e nesta instalacao (numpy sobre a mesma stack que o
    # GDAL do rasterio) ele derruba o processo em codigo nativo, sem traceback
    # -- o programa morria com exit 127 logo depois de extrair os perfis. Com
    # dois parametros a forma fechada e exata e nao depende de biblioteca
    # externa nenhuma.
    xm, ym = float(x.mean()), float(y.mean())
    sxx = float(((x - xm) ** 2).sum())
    if sxx < 1e-9:
        esperado = float(np.median(y))
    else:
        a = float(((x - xm) * (y - ym)).sum() / sxx)
        esperado = float(a * (s.rs - xm) + ym)
    # Alem do limiar absoluto, exige-se que o residuo destoe da DISPERSAO das
    # vizinhas. Sobre DEM de 30 m o talvegue oscila alguns metros por ruido de
    # amostragem, e so o limiar absoluto marcava 21% das secoes do Mirim --
    # ruido apresentado como anomalia. O teste C pergunta "qual delas nao se
    # encaixa", e isso e desvio relativo ao espalhamento, nao valor fixo.
    # MAD em vez de desvio padrao para o proprio ponto suspeito nao inflar a
    # referencia contra a qual esta sendo julgado.
    res = y - (a * (x - xm) + ym) if sxx >= 1e-9 else y - esperado
    mad = float(np.median(np.abs(res - np.median(res))))
    limite = max(lim.salto_talvegue, 4.0 * 1.4826 * mad)
    d = abs(s.z_talvegue - esperado)
    if d > limite:
        return CRITICA, (f"talvegue {d:.1f} m fora da tendencia do trecho "
                         f"(esperado {esperado:.2f} m, medido "
                         f"{s.z_talvegue:.2f} m)"), 0.0
    return OK, "", float(np.clip(1.0 - d / limite, 0.0, 1.0))


def teste_d_largura(s, vizinhas, lim):
    """Largura anormal contra as vizinhas."""
    ls = [v.largura for v in vizinhas if v is not s and v.largura > 0]
    if len(ls) < 2 or s.largura <= 0:
        return OK, "", 1.0
    med = float(np.median(ls))
    r = max(s.largura / med, med / s.largura)
    if r > lim.razao_largura:
        return ATENCAO, (f"largura {s.largura:.0f} m contra {med:.0f} m das "
                         f"vizinhas ({r:.1f}x)"), 0.2
    return OK, "", float(np.clip(1.0 - (r - 1.0) / (lim.razao_largura - 1.0),
                                 0.0, 1.0))


def teste_e_orientacao(s, lim):
    """Perpendicularidade em relacao ao eixo."""
    if s.azimute is None or not np.isfinite(s.azimute):
        return OK, "", 1.0
    desvio = abs(90.0 - float(s.azimute))
    if desvio > lim.desvio_ortogonal:
        return ATENCAO, f"secao a {s.azimute:.0f} graus do eixo (ideal 90)", 0.2
    return OK, "", float(np.clip(1.0 - desvio / lim.desvio_ortogonal, 0.0, 1.0))


PESOS = {"A": 0.40, "B": 0.25, "C": 0.15, "D": 0.10, "E": 0.10}


def avaliar(s, vizinhas=(), lim=None):
    """Roda os cinco testes e consolida status e nota de 0 a 100."""
    lim = lim or Limiares()
    if not s.valida:
        r = Resultado(INCERTO, 0.0, ["sem terreno valido na secao (NoData)"],
                      {})
        s.qc = r
        return r

    brutos = {
        "A": teste_a_posicao(s, lim),
        "B": teste_b_canal(s, lim),
        "C": teste_c_salto(s, vizinhas, lim),
        "D": teste_d_largura(s, vizinhas, lim),
        "E": teste_e_orientacao(s, lim),
    }
    status, motivos, nota = OK, [], 0.0
    testes = {}
    for k, (st, msg, sub) in brutos.items():
        status = _pior(status, st)
        if msg:
            motivos.append(msg)
        nota += PESOS[k] * float(sub)
        testes[k] = {"status": st, "motivo": msg, "nota": round(100 * sub, 1)}
    r = Resultado(status, round(100.0 * nota, 1), motivos, testes)
    s.qc = r
    return r


def avaliar_todas(secoes, lim=None):
    """Analisar todas: cada secao contra suas vizinhas NO MESMO trecho.

    O agrupamento por (rio, trecho) nao e detalhe. Com varios rios na mesma
    lista, as vizinhas de uma secao do Benedito acabariam sendo secoes do
    Mirim, e os testes C (salto de talvegue) e D (largura anormal) passariam a
    comparar rios diferentes -- um afluente de serra contra um rio de planicie.
    Toda secao viraria anomalia, e o relatorio, ruido.
    """
    lim = lim or Limiares()
    n = lim.n_vizinhas
    grupos = {}
    for s in secoes:
        grupos.setdefault((s.rio, s.reach), []).append(s)
    for g in grupos.values():
        for i, s in enumerate(g):
            viz = g[max(0, i - n):i] + g[i + 1:i + 1 + n]
            avaliar(s, viz, lim)
    return secoes


def contagem(secoes):
    c = {OK: 0, ATENCAO: 0, INCERTO: 0, CRITICA: 0}
    for s in secoes:
        if s.qc:
            c[s.qc.status] = c.get(s.qc.status, 0) + 1
    return c
