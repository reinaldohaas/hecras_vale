# -*- coding: utf-8 -*-
"""
Onde esta o canal dentro do perfil.

O menor ponto do perfil NAO e necessariamente o talvegue. Um buraco de uma
celula na margem -- ruido de DEM, sombra de radar, um valor NoData mal
preenchido -- ganha do canal de verdade se o criterio for so "cota minima":

        _________
                  \\_
                    \\______        <- canal real
                    canal
                           \\__ *    <- buraco de 1 celula, mais fundo

Aqui a decisao usa PROEMINENCIA, que e quanto uma depressao afunda em relacao
ao entorno antes de encontrar um caminho para outra depressao mais funda. E a
grandeza que separa "vale" de "cova": o buraco de uma celula tem proeminencia
de centimetros; o canal, de metros. Ao lado dela entram a proximidade ao eixo
do rio -- que e informacao independente do DEM -- e uma penalidade por estar
encostado na borda.

Sao devolvidos os TRES pontos, sempre, para o usuario poder discordar:
  - minimo absoluto do perfil;
  - depressao principal (maior proeminencia);
  - talvegue provavel (o escolhido).

Quando nem a maior proeminencia passa do limiar, o resultado e INCERTO. Nao ha
chute: o programa existe para achar secao ruim, nao para maquiar.
"""
import numpy as np
from scipy.signal import find_peaks, savgol_filter

PROEMINENCIA_MIN = 0.5     # m; abaixo disto nao e canal, e rugosidade
JANELA_SUAVIZA = 7         # pontos; so para DETECTAR, nunca para exportar
PESO_PROEMINENCIA = 1.0
PESO_EIXO = 1.0
PESO_BORDA = 0.6


def suavizar(z, janela=JANELA_SUAVIZA):
    """Savitzky-Golay curto: tira o degrau de celula sem mover o vale.

    Media movel arredondaria o fundo do canal e deslocaria o minimo; o
    polinomio local preserva a posicao do extremo, que e o que interessa.
    """
    z = np.asarray(z, float)
    ok = np.isfinite(z)
    if ok.sum() < 5:
        return z.copy()
    zz = z.copy()
    if not ok.all():                       # interpola so para filtrar
        zz[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), z[ok])
    j = int(min(max(janela | 1, 5), (len(zz) // 2) * 2 - 1))
    if j < 5:
        return z.copy()
    s = savgol_filter(zz, j, 2)
    s[~ok] = np.nan                        # o buraco continua sendo buraco
    return s


def candidatos(sta, z, proeminencia_min=PROEMINENCIA_MIN):
    """Minimos locais com proeminencia, do mais proeminente para o menos."""
    z = np.asarray(z, float)
    ok = np.isfinite(z)
    if ok.sum() < 3:
        return []
    zz = z.copy()
    if not ok.all():
        # tapa os buracos por CIMA (maximo local) para nao criar vale falso
        zz[~ok] = np.nanmax(z[ok])
    idx, prop = find_peaks(-zz, prominence=float(proeminencia_min))
    saida = []
    for i, p in zip(idx, prop["prominences"]):
        saida.append({"i": int(i), "sta": float(sta[i]), "z": float(zz[i]),
                      "proeminencia": float(p)})
    # extremidades: find_peaks nunca marca o primeiro nem o ultimo ponto, e um
    # canal cortado ao meio pela borda aparece exatamente ali. Sem isto uma
    # secao mal posicionada -- justamente o alvo do programa -- fica sem
    # candidato nenhum e cai em INCERTO por engano.
    for i in (0, len(zz) - 1):
        if ok[i]:
            viz = zz[1:min(len(zz), 25)] if i == 0 else zz[max(0, len(zz) - 25):-1]
            if len(viz) and zz[i] < np.nanmin(viz) + 1e-9:
                saida.append({"i": int(i), "sta": float(sta[i]),
                              "z": float(zz[i]),
                              "proeminencia": float(np.nanmax(viz) - zz[i]),
                              "borda": True})
    saida.sort(key=lambda c: -c["proeminencia"])
    return saida


def detectar(sta, z, sta_eixo=None, proeminencia_min=PROEMINENCIA_MIN,
             usar_suavizado=True):
    """Talvegue provavel, depressao principal e minimo absoluto.

    sta_eixo: estacao, no perfil, em que a secao cruza o eixo do rio. E a
    informacao mais valiosa que existe aqui, porque nao vem do DEM: mesmo com
    terreno ruim, o canal esta perto do eixo. Se for None, so o DEM decide.
    """
    sta = np.asarray(sta, float)
    z = np.asarray(z, float)
    ok = np.isfinite(z)
    fora = {"i_talvegue": None, "i_min_abs": None, "i_principal": None,
            "incerto": True, "motivo": "perfil sem cota valida",
            "candidatos": [], "proeminencia": 0.0, "profundidade_relativa": 0.0}
    if ok.sum() < 3:
        return fora

    z_det = suavizar(z) if usar_suavizado else z
    i_min_abs = int(np.nanargmin(z))
    cand = candidatos(sta, z_det, proeminencia_min)
    L = float(sta[-1] - sta[0]) or 1.0

    if not cand:
        return dict(fora, i_min_abs=i_min_abs, i_talvegue=i_min_abs,
                    i_principal=i_min_abs, incerto=True,
                    motivo=f"nenhuma depressao com mais de "
                           f"{proeminencia_min:.2f} m de proeminencia",
                    profundidade_relativa=profundidade_relativa(z, i_min_abs))

    i_principal = cand[0]["i"]
    prom_max = cand[0]["proeminencia"]

    # pontuacao: proeminencia relativa + proximidade ao eixo - penalidade de
    # borda. Sem o eixo, a proeminencia decide praticamente sozinha.
    for c in cand:
        nota = PESO_PROEMINENCIA * (c["proeminencia"] / prom_max)
        if sta_eixo is not None:
            d = abs(c["sta"] - float(sta_eixo))
            nota += PESO_EIXO * float(np.exp(-(d / (0.15 * L)) ** 2))
        pos = (c["sta"] - sta[0]) / L
        nota -= PESO_BORDA * max(0.0, 1.0 - min(pos, 1.0 - pos) / 0.10) * 0.5
        c["nota"] = float(nota)
    melhor = max(cand, key=lambda c: c["nota"])

    return {"i_talvegue": int(melhor["i"]),
            "i_min_abs": i_min_abs,
            "i_principal": int(i_principal),
            "proeminencia": float(melhor["proeminencia"]),
            "profundidade_relativa": profundidade_relativa(z, int(melhor["i"])),
            "candidatos": cand[:8],
            "incerto": bool(melhor["proeminencia"] < proeminencia_min),
            "motivo": ("" if melhor["proeminencia"] >= proeminencia_min else
                       "canal raso demais para ser identificado com seguranca")}


def profundidade_relativa(z, i, n_borda=5):
    """Cota media das extremidades menos a cota do talvegue.

    Mede se existe canal, nao onde ele esta. Uma secao inteiramente em rampa
    tem minimo (a ponta mais baixa) mas nao tem canal, e sai com valor
    pequeno ou negativo.
    """
    z = np.asarray(z, float)
    if i is None or not np.isfinite(z[i]):
        return 0.0
    e = z[:n_borda][np.isfinite(z[:n_borda])]
    d = z[-n_borda:][np.isfinite(z[-n_borda:])]
    if not len(e) or not len(d):
        return 0.0
    return float(0.5 * (e.mean() + d.mean()) - z[i])
