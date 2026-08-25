# -*- coding: utf-8 -*-
"""Geometria de um rio a partir do RELEVO, sem nada esculpido.

    python scripts/rio_do_relevo.py --rio Rio_Benedito --saida modelo/benedito

NAO IMPORTA NADA DA CADEIA DE CORRECOES. Usa so `mdt_sigsc` (acesso ao MDT) e
`ras_io` (gravar em CRLF, que o HEC-RAS exige). Tudo o mais e calculado aqui, a
partir do terreno.

O QUE ESTE GERADOR NAO FAZ, e por que ele existe

  `gerar_mirim_do_zero.py` produz um modelo que converge em 2 iteracoes e tem 2
  erros de geometria -- mas o talvegue dele sao OITO NUMEROS escritos no
  codigo, interpolados por PCHIP, e a calha e uma parabola:

      z[calha] = z_alvo + (cota_margem - z_alvo) * (dist_norm ** 2)
      z_lob    = max(z_lob, z_alvo + 2.5)
      z[esq]   = np.maximum(z[esq], z_lob)
      z[0]     = max(z[0], 4.50)

  Medido contra o MDT em 64 secoes: dentro da calha a mediana e -0,50 m mas o
  p90 chega a +39 m; na planicie a mediana e zero e o p90 e +25 m, que e o
  `np.maximum` levantando o terreno onde ele desce. Converge porque e liso por
  construcao, e nao porque descreve o rio.

  Aqui NADA e esculpido: o perfil e o que o MDT da, ponto a ponto.

O QUE E MEDIDO, E COMO

  talvegue    o ponto mais baixo do MDT perto do eixo. RESSALVA QUE NAO SE
              resolve: o MDT ve a LAMINA D'AGUA, nao o fundo. Onde ha agua, o
              "talvegue" e a superficie livre. Sem batimetria nao ha como
              saber o fundo, e inventa-lo e o que este gerador recusa fazer.

  margens     andando para fora do talvegue, o primeiro ponto de cada lado que
              sobe `FOLGA_CALHA` acima dele. E o topo do encaixe, medido.

  meia-largura  continua para fora ate subir `ALVO_SECAO` acima do talvegue,
              com teto em `MEIA_MAX`. Onde a varzea e plana o teto manda, e o
              relatorio diz em quantas secoes isso aconteceu -- ali a secao
              1D nao contem a cheia, e isso pede armazenamento ou 2D.

  Manning     `N_CALHA` e `N_PLANICIE`, constantes e declaradas. Nao ha dado de
              rugosidade nesta bacia; fingir que ha seria o mesmo erro.

O PERFIL LONGITUDINAL SAI CRU

  Sem suavizacao, sem declividade minima forcada, sem monotonicidade imposta.
  O relatorio mede quantos contradeclives e quantos degraus o terreno traz, e
  a decisao de tratar isso e de quem le -- com `--monotono` disponivel, que
  aplica regressao isotonica: ela ajusta uma curva nao-crescente aos valores
  MEDIDOS, sem sair da faixa deles.
"""
import argparse
import json
import os
import sys

import numpy as np
from shapely.geometry import LineString, Point

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402
from ras_io import escrever                            # noqa: E402

EIXOS = "eixos_do_relevo.geojson"
DX = 150.0            # m entre secoes
PASSO = 4.0           # m entre pontos amostrados na cutline
JANELA = 60.0         # m para cada lado, ao medir a tangente do eixo
MEIA_MAX = 400.0      # m; teto da meia-largura
MEIA_MIN = 60.0       # m; piso, para a secao nunca degenerar
FOLGA_CALHA = 1.5     # m acima do talvegue = topo da margem
DESCE_MAX = 0.30
DESCE_FORA = 0.30     # m abaixo do fundo da calha; abaixo disso e outro vale
JANELA_MARGEM = 3      # secoes para cada lado, na mediana movel da margem
ALVO_SECAO = 8.0      # m acima do talvegue = onde a secao pode parar
TAXA_LARGURA = 0.15   # quanto a meia-largura pode variar por metro de rio
K_CURV = 0.50         # fracao do raio de curvatura admitida do lado de dentro
EXTRA_EIXO = 150.0   # m de eixo alem da primeira e da ultima secao; nunca
                     # menos que UM espacamento -- com 40 m o Oeste levava
                     # 'XS must intersect exactly one Reach' na primeira
                     # secao, com a interseccao medindo 40,00 m do inicio
                     # em todo teste geometrico que montei. Nao descobri a
                     # tolerancia que o RAS usa; medi que 150 m limpa.
BUSCA = 500.0         # m; ate onde se procura terreno alto
N_CALHA, N_PLANICIE = 0.032, 0.055
WKT = ('PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",'
       'DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],'
       'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
       'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
       'PARAMETER["False_Northing",10000000.0],'
       'PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],'
       'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')


def eixo_do_rio(nome, caminho=EIXOS):
    d = json.load(open(caminho, encoding="utf-8"))
    for f in d["features"]:
        if f["properties"].get("nome") == nome:
            return LineString(np.asarray(f["geometry"]["coordinates"], float))
    raise SystemExit(f"'{nome}' nao esta em {caminho}. Ha: "
                     f"{[f['properties'].get('nome') for f in d['features']]}")


def curvatura(eixo, s, janela):
    """Raio de curvatura e sinal do giro no ponto `s` do eixo.

    Tres pontos a `janela` de distancia definem um circulo; o raio dele e o
    raio de curvatura local, e o sinal do produto vetorial diz para que lado o
    eixo vira. Janela curta demais mede ruido da polilinha, longa demais mede a
    corda -- por isso ela vem proporcional ao espacamento das secoes.
    """
    a = max(s - janela, 0.0)
    b = min(s + janela, eixo.length)
    P0 = np.array(eixo.interpolate(a).coords[0])
    P1 = np.array(eixo.interpolate(0.5 * (a + b)).coords[0])
    P2 = np.array(eixo.interpolate(b).coords[0])
    v1, v2 = P1 - P0, P2 - P1
    cr = v1[0] * v2[1] - v1[1] * v2[0]
    d01 = float(np.hypot(*(P1 - P0)))
    d12 = float(np.hypot(*(P2 - P1)))
    d02 = float(np.hypot(*(P2 - P0)))
    if abs(cr) < 1e-9 or d01 * d12 * d02 == 0:
        return np.inf, 0.0
    return d01 * d12 * d02 / (2.0 * abs(cr)), float(np.sign(cr))


def isotonica(z):
    """Maior curva NAO-CRESCENTE que melhor ajusta z (pool adjacent violators).

    Nao inventa valor fora da faixa medida: cada patamar da saida e a MEDIA de
    um bloco de valores de entrada.
    """
    v = [float(z[0])]
    w = [1.0]
    for x in z[1:]:
        v.append(float(x))
        w.append(1.0)
        while len(v) > 1 and v[-2] < v[-1]:      # violou o nao-crescente
            x2 = v.pop()
            w2 = w.pop()
            x1 = v.pop()
            w1 = w.pop()
            v.append((x1 * w1 + x2 * w2) / (w1 + w2))
            w.append(w1 + w2)
    saida = []
    for x, k in zip(v, w):
        saida += [x] * int(k)
    return np.array(saida)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rio", required=True)
    ap.add_argument("--saida", required=True)
    ap.add_argument("--reach", default="R1")
    ap.add_argument("--dx", type=float, default=DX)
    ap.add_argument("--taxa", type=float, default=TAXA_LARGURA,
                    help="quanto a meia-largura pode variar por metro de rio")
    ap.add_argument("--monotono", action="store_true",
                    help="ajusta o talvegue por regressao isotonica")
    a = ap.parse_args()

    eixo = eixo_do_rio(a.rio)
    L = eixo.length
    est = np.arange(0.0, L, a.dx)
    if L - est[-1] > 20.0:
        est = np.append(est, L)
    print(f"rio    : {a.rio}   eixo {L/1000:.2f} km   {len(est)} secoes "
          f"a cada {a.dx:g} m")

    # ---- geometria das cutlines e amostragem do MDT
    off = np.arange(-BUSCA, BUSCA + PASSO / 2, PASSO)
    base, normais, pts = [], [], []
    for s in est:
        p = np.array(eixo.interpolate(s).coords[0])
        q0 = np.array(eixo.interpolate(max(s - JANELA, 0.0)).coords[0])
        q1 = np.array(eixo.interpolate(min(s + JANELA, L)).coords[0])
        t = q1 - q0
        nt = float(np.hypot(*t))
        if nt < 1e-6:
            continue
        t /= nt
        n = np.array([-t[1], t[0]])
        base.append((float(s), p))
        normais.append(n)
        for o in off:
            pts.append(p + o * n)
    pts = np.array(pts)
    bb = (pts[:, 0].min() - 60, pts[:, 1].min() - 60,
          pts[:, 0].max() + 60, pts[:, 1].max() + 60)
    tiles = tiles_do_dominio(bb)
    print(f"MDT    : {len(pts)} pontos a {PASSO:g} m sobre {len(tiles)} folhas")
    Z = MosaicoSigsc(tiles=tiles).cota(pts[:, 0], pts[:, 1]) \
        .reshape(len(base), len(off))

    # ---- onde o MDT e VAZIO, a secao LEVANTADA do legado entra INTEIRA.
    #
    # No SIG-SC a lamina d'agua e 0.0 = nodata, entao o estuario inteiro sai
    # "sem MDT utilizavel": no Acu eram 109 secoes descartadas e um VAO DE
    # 1500 m colado no contorno de mare. O legado tem exatamente esse trecho
    # levantado (Acu R4, RS 5534 a 75; canal do porto com 2.116 m ENTRE AS
    # MARGENS). A primeira versao reamostrava esse perfil no meu grid de
    # +-BUSCA m -- e cortava o canal pela metade: a secao adotada nao chegava
    # ao outro lado do rio, como se viu no RAS Mapper. Agora a secao entra
    # INTEIRA, com a cutline, as margens e o perfil DELA, uma unica vez (na
    # estacao de RS mais proximo -- adotar a mesma em varias estacoes
    # duplicaria a cutline). Ela nao passa pelos filtros de largura, que
    # existem para ruido do MDT, nao para levantamento; e entra na lista so
    # depois deles. O casamento e por RS (imune ao pescoco de meandro) com
    # folga de meio espacamento do legado, e distancia 2D larga de sanidade.
    from batimetria_do_legado import secoes_levantadas, LEGADO
    LEG_RS = 800.0          # m de RS; ~metade do maior espacamento do legado
    LEG_XY = 1200.0         # m; sanidade contra casamento grosseiramente errado
    _leg = secoes_levantadas(LEGADO, a.rio, completas=True) or []
    _leg_rs = np.array([d["rs"] for d in _leg]) if _leg else np.array([])
    _s_base = np.array([b[0] for b in base])
    _dono = (np.array([int(np.argmin(np.abs((L - _s_base) - r)))
                       for r in _leg_rs]) if len(_leg_rs) else np.array([]))

    def adotar(k, s, p, fundo_max=None):
        """A secao levantada cuja estacao-dona e `k`, inteira, ou None.

        Com `fundo_max`, so adota se o fundo levantado esta ABAIXO dele --
        e o criterio da zona de mare: o lidar mede a lamina (~0 m) e o
        levantamento sabe que o canal desce a -10,85; se o legado nao for
        mais fundo que o que o MDT ja ve, nao ha razao para trocar.
        """
        for j in np.flatnonzero(_dono == k) if len(_dono) else []:
            d = _leg[int(j)]
            if abs(d["rs"] - (L - s)) > LEG_RS:
                continue
            if float(np.hypot(d["x"] - p[0], d["y"] - p[1])) > LEG_XY:
                continue
            m = (d["sta"] >= d["lb"]) & (d["sta"] <= d["rb"])
            zt_leg = float(d["z"][m].min() if m.any() else d["z"].min())
            if fundo_max is not None and zt_leg >= fundo_max:
                continue
            return {"s": s, "rs": round(float(L - s), 2), "pronta": True,
                    "sta": np.asarray(d["sta"], float),
                    "z": np.asarray(d["z"], float),
                    "lb": float(d["lb"]), "rb": float(d["rb"]),
                    "cut": (np.asarray(d["cut"][0], float),
                            np.asarray(d["cut"][1], float)),
                    "zt": zt_leg}
        return None

    # ---- cada secao, medida
    secoes, prontas, no_teto, sem_dado = [], [], 0, 0
    for k, ((s, p), n) in enumerate(zip(base, normais)):
        z = Z[k]
        centro = np.abs(off) <= 30.0
        # O criterio e COBERTURA, nao presenca. No estuario a lamina e vazio
        # mas sobram respingos de borda d'agua: meia duzia de pixels no
        # centro passavam no teste "tem dado", viravam secao-lixo (76 m de
        # largura, 5 cm de fundo, na boca do porto) e morriam adiante no
        # filtro de faixa vazia -- 19 assim, e a foz ficava com vao de
        # 1500 m. Centro com menos de metade das amostras finitas nao e
        # medida: vai para a adocao do legado, ou fora.
        if float(np.isfinite(z[centro]).mean()) < 0.5:
            d = adotar(k, s, p)
            if d is not None:
                prontas.append(d)
            else:
                sem_dado += 1
            continue
        # NA ZONA DE MARE O LIDAR VE A LAMINA, NAO O CANAL. Na boca do porto
        # a estacao tinha cobertura boa (a restinga) e virava secao de 224 m
        # com fundo 0,17 m -- enquanto o levantamento sabe que ali o canal
        # tem 2.116 m entre margens e fundo -10,85. Onde o talvegue do MDT
        # esta abaixo de 2 m e o legado conhece fundo mais de 2 m abaixo
        # dele, a secao levantada vale mais que o respingo.
        zt_c = float(np.nanmin(np.where(centro, z, np.nan)))
        if zt_c < 2.0:
            d = adotar(k, s, p, fundo_max=zt_c - 2.0)
            if d is not None:
                prontas.append(d)
                continue
        i0 = int(np.nanargmin(np.where(centro, z, np.nan)))
        zt = float(z[i0])

        def anda(sinal, alvo):
            """Anda para fora ate subir `alvo`, ou ate CAIR abaixo do talvegue.

            A segunda parada e o que faltava: a busca ia ate 500 m e podia
            entrar em OUTRO curso d'agua mais baixo. Medido no Benedito, 63 das
            294 secoes ficavam com o ponto mais baixo FORA das proprias
            margens, ate 8,64 m abaixo do invert da calha, e o HEC-RAS acusou
            62 delas com "hTab starting values below the XS invert". Terreno
            abaixo do talvegue deste rio nao e planicie deste rio.
            """
            i = i0
            ult = i0
            while 0 < i < len(off) - 1:
                i += sinal
                if not np.isfinite(z[i]):
                    continue
                if z[i] < zt - DESCE_MAX:
                    return ult
                ult = i
                if z[i] >= zt + alvo:
                    return i
            return None

        # QUAL LADO E A ESQUERDA. `n` e a tangente girada 90 graus no sentido
        # anti-horario, e num mapa (x para leste, y para norte) isso e a
        # margem ESQUERDA olhando para jusante. Como o ponto de offset `o` foi
        # amostrado em `p + o*n`, offset POSITIVO e a esquerda.
        #
        # E o HEC-RAS exige a cutline percorrida da ESQUERDA para a DIREITA:
        # a estaca 0 fica no offset MAIS POSITIVO e cresce descendo o offset.
        # Ter escrito ao contrario custou "XS is reversed" em 293 das 294
        # secoes, e bank line atravessando o rio em todo o trecho.
        i_esq_calha = anda(+1, FOLGA_CALHA)
        i_dir_calha = anda(-1, FOLGA_CALHA)
        i_esq_sec = anda(+1, ALVO_SECAO)
        i_dir_sec = anda(-1, ALVO_SECAO)
        if i_esq_sec is None or i_dir_sec is None:
            no_teto += 1
        me = float(np.clip(abs(off[i_esq_sec]) if i_esq_sec is not None
                           else MEIA_MAX, MEIA_MIN, MEIA_MAX))
        md = float(np.clip(abs(off[i_dir_sec]) if i_dir_sec is not None
                           else MEIA_MAX, MEIA_MIN, MEIA_MAX))
        # LIMITE DE CURVATURA, do lado de DENTRO da curva.
        # Numa curva de raio R as secoes vizinhas convergem pelo lado interno e
        # se encontram no centro de curvatura: passar de R ali faz a secao
        # cruzar a vizinha e dobrar sobre o proprio rio, que e o que produz a
        # auto-interseccao das edge lines. Medido nos seis rios, a relacao e
        # direta -- Norte com secao de 128 m tem 2 erros, Sul com 208 m tem 15,
        # Acu com 264 m tem 190 e Oeste com 460 m tem 303, e em ambos os piores
        # 76% a 86% dos erros sao auto-interseccao de edge line.
        R, giro = curvatura(eixo, s, max(2.0 * a.dx, 200.0))
        if np.isfinite(R):
            lim = K_CURV * R
            if giro > 0:            # eixo vira a esquerda: lado interno e o +
                me = min(me, max(lim, MEIA_MIN))
            elif giro < 0:
                md = min(md, max(lim, MEIA_MIN))
        secoes.append({"s": s, "rs": round(float(L - s), 2),
                       "p": p, "n": n, "z_perfil": z, "zt_bruto": zt,
                       "me": me, "md": md, "off_t": float(off[i0]),
                       "d_esq": (float(off[i_esq_calha] - off[i0])
                                 if i_esq_calha is not None else np.nan),
                       "d_dir": (float(off[i0] - off[i_dir_calha])
                                 if i_dir_calha is not None else np.nan),
                       # fundo DA CALHA, e nao da secao: e a referencia para
                       # decidir o que esta abaixo do rio e portanto fora dele
                       "z_calha": float(zt)})

    print(f"secoes : {len(secoes)}   sem MDT utilizavel: {sem_dado}   "
          f"adotadas do legado (inteiras): {len(prontas)}   "
          f"pararam no teto de {MEIA_MAX:g} m: {no_teto}")

    # ---- a MARGEM E MEDIDA, mas a medida e ruidosa
    # Detectar o topo do encaixe secao a secao da calha de 41 a 371 m, com
    # SALTO DE ATE 752 m entre vizinhas -- onde a varzea e plana o ponto que
    # sobe 1,5 m acima do talvegue pode estar muito longe. A bank line liga
    # esses pontos e vira um zigue-zague que cruza o rio 82 vezes.
    # A mediana movel nao inventa nada: cada valor de saida E UMA DAS MEDIDAS
    # da vizinhanca. Filtra o salto e preserva a variacao real do rio.
    # A MEIA-LARGURA DA SECAO TAMBEM E RUIDOSA, e nao so a da calha. A edge
    # line liga as PONTAS das secoes: onde a largura salta de uma para a
    # outra, ela faz gancho e cruza. Medido nos rios ja gerados, o salto de
    # largura entre vizinhas separa os limpos dos sujos com clareza --
    # Norte 12 m de mediana e 2 erros, Acu 32 m e 190, Oeste 88 m e 303.
    for lado in ("me", "md", "d_esq", "d_dir"):
        v = np.array([s[lado] for s in secoes], float)
        bom = np.isfinite(v)
        if bom.sum() < 3:
            v[:] = np.nanmedian(v) if bom.any() else 20.0
        else:
            v[~bom] = np.interp(np.flatnonzero(~bom), np.flatnonzero(bom),
                                v[bom])
        sv = v.copy()
        # NAO chamar estes de `a` e `b`: `a` e o namespace dos argumentos, e
        # sobrescreve-lo faz `a.monotono` estourar mais adiante.
        for i in range(len(v)):
            j0 = max(0, i - JANELA_MARGEM)
            j1 = min(len(v), i + JANELA_MARGEM + 1)
            sv[i] = float(np.median(v[j0:j1]))
        bruto = np.abs(np.diff(v))
        filt = np.abs(np.diff(sv))
        print(f"   {lado}: salto entre vizinhas   bruto mediana "
              f"{np.median(bruto):.0f} m / max {bruto.max():.0f}   ->   "
              f"filtrado {np.median(filt):.0f} / {filt.max():.0f}")
        # LIMITE DE TAXA, que e o que generaliza para qualquer rio.
        # A mediana movel tira o pico isolado mas deixa a largura mudar de uma
        # secao para a outra o quanto o terreno quiser. A edge line liga as
        # PONTAS: se a meia-largura anda `D` metros para o lado enquanto o rio
        # anda `dx` para a frente, a linha faz angulo `atan(D/dx)` com a secao,
        # e a partir de certo ponto ela dobra e cruza a vizinha.
        #
        # Medido nos rios ja gerados, essa taxa separa os limpos dos sujos
        # melhor do que qualquer largura absoluta:
        #     Norte  12 m por secao (0,08)  ->    2 erros
        #     Acu    32 m por secao (0,21)  ->  190 erros
        #     Oeste  88 m por secao (0,59)  ->  303 erros
        #
        # O limite e aplicado em duas passadas, para frente e para tras, e SO
        # ENCOLHE: nenhuma secao fica mais larga do que o terreno mediu. Por
        # isso ele nunca inventa planicie -- no maximo deixa de usar parte da
        # que existe, e o relatorio diz quanto.
        teto = a.taxa * a.dx
        lim = sv.copy()
        for i in range(1, len(lim)):
            lim[i] = min(lim[i], lim[i - 1] + teto)
        for i in range(len(lim) - 2, -1, -1):
            lim[i] = min(lim[i], lim[i + 1] + teto)
        cortado = float(np.sum(sv - lim))
        if cortado > 1.0:
            print(f"      limite de taxa ({a.taxa:g} m/m): retirou "
                  f"{cortado/max(len(lim),1):.1f} m por secao em media")
        for s, x in zip(secoes, lim):
            s[lado] = float(x)

    # LIMITE DE CURVATURA, aplicado DEPOIS do filtro e nao antes.
    # Aplicado antes, a mediana movel o desfazia: bastava uma vizinha larga
    # para o valor limitado voltar a subir. Medido no Acu, nos 14 segmentos em
    # que a edge line se auto-intersecta o raio mediano e 545 m contra 1.420 m
    # do rio inteiro, e a meia-largura e 234 m contra 92 m -- e secao larga em
    # curva fechada, e nao largura sozinha nem curva sozinha.
    #
    # A ponta da secao descreve uma curva paralela ao eixo, deslocada de `w`.
    # Pelo lado de dentro essa paralela encolhe, e em `w = R` ela colapsa no
    # centro de curvatura: dali em diante ela anda para tras e cruza a si
    # mesma. `K_CURV` fica abaixo de 1 com folga, porque as secoes sao
    # discretas e a colisao chega antes do limite continuo.
    n_curv = 0
    for s in secoes:
        R, giro = curvatura(eixo, s["s"], max(2.0 * a.dx, 200.0))
        if not np.isfinite(R):
            continue
        lim = max(K_CURV * R, MEIA_MIN)
        if giro > 0 and s["me"] > lim:
            s["me"] = lim
            n_curv += 1
        elif giro < 0 and s["md"] > lim:
            s["md"] = lim
            n_curv += 1
    if n_curv:
        print(f"   limite de curvatura (K={K_CURV:g}): apertou {n_curv} secoes")

    def montar(s, me, md):
        """Perfil, cutline e margens de UMA secao, para as larguras dadas.

        Isolado em funcao porque a passada anti-dobra, mais abaixo, chama de
        novo -- com largura menor e so nas secoes culpadas.
        """
        z = s["z_perfil"]
        # A SECAO PARA NO DIVISOR, E NAO NA LARGURA PEDIDA.
        #
        # Andando do talvegue para fora, o terreno sobe ate um divisor e depois
        # DESCE -- para o vale vizinho, para uma cava, para o proprio rio noutro
        # meandro. O que esta depois do divisor nao troca agua com esta secao
        # nesta cota, e o HEC-RAS nao sabe disso: ele ve um ponto baixo dentro
        # da secao, enche primeiro, e o solver passa a procurar nivel numa
        # secao com dois fundos. Medido no Mirim: 309 das 765 secoes tinham o
        # ponto mais baixo FORA da calha, uma delas 27,43 m abaixo do fundo --
        # e a rodada batia nas 40 iteracoes em todo instante e terminava com
        # "Solution Solver Failed" e 92,38% de erro de volume.
        #
        # Cortar no divisor NAO altera cota nenhuma: so encurta a secao. O que
        # ficou de fora nao vira cota inventada -- vira area de armazenamento,
        # se for o caso, e isso e decisao de modelagem e nao de amostragem.
        inv = s.get("z_calha")
        if inv is not None:
            i0 = int(np.argmin(np.abs(off - s["off_t"])))
            for sinal in (+1, -1):
                lim = me if sinal > 0 else md
                j = i0
                while True:
                    j += sinal
                    if j < 0 or j >= len(off) or abs(off[j] - s["off_t"]) > lim:
                        break
                    if np.isfinite(z[j]) and z[j] < inv - DESCE_FORA:
                        novo = abs(off[j - sinal] - s["off_t"])
                        if sinal > 0:
                            me = min(me, s["off_t"] + novo)
                        else:
                            md = min(md, novo - s["off_t"])
                        break
        dentro = (off <= me) & (off >= -md) & np.isfinite(z)
        # A SECAO TEM DE ALCANCAR O PROPRIO RIO, dos dois lados. Onde o MDT
        # falta perto do eixo -- no Acu sao 109 secoes sem dado utilizavel, e o
        # estuario tem buracos -- os pontos validos podem ficar todos de um
        # lado: a cutline passa ao largo e o HEC-RAS acusa "XS doesn't
        # intersect the associated Reach" junto com "XS intersects < 2
        # banklines". Secao que nao chega no rio nao e secao: sai.
        tem_esq = bool((dentro & (off > 0)).any())
        tem_dir = bool((dentro & (off < 0)).any())
        if dentro.sum() < 8 or not (tem_esq and tem_dir):
            s["sta"] = None
            return False
        ordem = np.argsort(-off[dentro])          # estaca cresce esq -> dir
        # A CUTLINE NASCE DO PERFIL, e nao da meia-largura pedida. Os pontos
        # sao amostrados de 4 em 4 m e o primeiro e o ultimo raramente caem
        # exatamente em `me` e `-md`: montar a cutline com os valores pedidos
        # faz a polilinha ficar mais longa que a amplitude das estacas, e o
        # HEC-RAS acusa "XS Profile length is different than Polyline length".
        o_esq = float(off[dentro].max())
        o_dir = float(off[dentro].min())
        s["o_esq"], s["o_dir"] = o_esq, o_dir
        s["sta"] = np.round((o_esq - off[dentro])[ordem], 2)
        s["z"] = np.round(z[dentro][ordem], 2)
        s["zt"] = float(s["z"].min())
        s["cut"] = (s["p"] + o_esq * s["n"], s["p"] + o_dir * s["n"])
        sta = s["sta"]
        # A CALHA TEM DE CONTER O EIXO. As margens sao medidas a partir do
        # talvegue, que nem sempre cai sobre o eixo; onde ele esta todo de um
        # lado, a margem "esquerda" acaba a DIREITA do eixo e a bank line
        # atravessa o rio. Obrigar o intervalo a cruzar o offset zero e o que
        # impede isso por construcao, e nao por conserto depois.
        o_l = max(s["off_t"] + s["d_esq"], PASSO)
        o_r = min(s["off_t"] - s["d_dir"], -PASSO)
        # A ESTACA DE UM OFFSET E `o_esq - o`, e nao `me - o`: a estaca zero
        # esta na PONTA DO PERFIL, nao na meia-largura pedida. Enquanto usei
        # `me` aqui, toda bank station saiu deslocada de `me - o_esq` metros
        # para fora -- exatamente o descolamento entre bank station e bank line
        # que o validador vinha acusando.
        # O SNAP NAO PODE DEVOLVER A MARGEM AO LADO ERRADO. `argmin(|sta-alvo|)`
        # escolhia a estaca mais proxima DOS DOIS LADOS: com a margem grampeada
        # a PASSO do eixo, "mais proxima" caia do outro lado do zero e o grampo
        # evaporava -- medido no Acu, 12 secoes de 1086 com a calha fora do
        # eixo (ex.: RS 32300, eixo na estaca 400,0 e lb=408) e 22 cruzamentos
        # de bank line com o rio, em pares, todos no cinturao de meandros. A
        # esquerda arredonda PARA FORA (estaca <= alvo) e a direita idem
        # (estaca >= alvo): o eixo fica dentro por construcao.
        cand = sta[sta <= o_esq - o_l]
        lb = float(cand.max()) if len(cand) else float(sta[0])
        cand = sta[sta >= o_esq - o_r]
        rb = float(cand.min()) if len(cand) else float(sta[-1])
        if rb <= lb:
            # canal degenerado: abre para as estacas vizinhas DO EIXO (estaca
            # `o_esq`, o offset zero), e nao do talvegue -- em volta do
            # talvegue as duas podiam cair do MESMO lado do eixo, e era dai
            # que a bank line atravessava o rio.
            j = int(np.argmin(np.abs(sta - o_esq)))
            lb = float(sta[max(j - 1, 0)])
            rb = float(sta[min(j + 1, len(sta) - 1)])
        s["lb"], s["rb"] = lb, rb
        return True

    for s in secoes:
        montar(s, s["me"], s["md"])
    antes = len(secoes)
    secoes = [s for s in secoes if s.get("sta") is not None]
    if len(secoes) < antes:
        print(f"   descartadas por faixa vazia depois do filtro: "
              f"{antes - len(secoes)}")

    # ---- A EDGE LINE NAO PODE DOBRAR SOBRE SI MESMA
    #
    # O HEC-RAS liga as pontas esquerdas de todas as secoes numa polilinha (a
    # "edge line"), faz o mesmo do lado direito, e usa as duas para montar a
    # superficie de interpolacao entre secoes. Se uma delas se cruza, o RAS
    # Mapper avisa "The generated edge lines have self intersections, the
    # interpolation surface may not generate correctly" e a superficie sai
    # errada -- e nada disso aparece na contagem do Validate Geometry.
    #
    # Tentei DEDUZIR a causa e nao fecha. Suspeitei da curvatura do eixo
    # (secao larga na volta fechada dobra sobre a vizinha) e limitei a
    # meia-largura a K*R; sobrou dobra em secao com raio de 10.909 m, ou seja,
    # em reta. Suspeitei do salto de largura e ja havia mediana movel mais
    # limite de taxa; sobrou. A causa e composta e nao vale mais palpite.
    #
    # Entao a condicao passa a ser IMPOSTA E VERIFICADA, e nao inferida:
    # mede-se o cruzamento com a propria geometria e encolhe-se quem participa
    # dele, ate nao haver nenhum. Encolher SO TIRA largura -- nenhuma secao
    # fica maior do que o terreno mediu -- e o preco esta no relatorio.
    from shapely.strtree import STRtree

    def dobras(pts):
        """Indices dos vertices envolvidos em cruzamento da polilinha."""
        seg = [LineString([pts[i], pts[i + 1]]) for i in range(len(pts) - 1)]
        if len(seg) < 3:
            return set()
        arv = STRtree(seg)
        maus = set()
        for i, g in enumerate(seg):
            for j in arv.query(g):
                j = int(j)
                if abs(i - j) <= 1:
                    continue
                if g.intersects(seg[j]):
                    maus.update((i, i + 1, j, j + 1))
        return maus

    def recuos(pts, k):
        """Pontas que ANDAM PARA TRAS ao longo do eixo.

        A poligonal das pontas nao e a edge line que o HEC-RAS constroi: no
        Oeste ele monta 862 vertices para 380 secoes, densificando o traco, e
        o cruzamento aparece SO na versao densificada -- meu teste de corda
        dava zero enquanto o validador acusava dois pontos de auto-
        interseccao. Lida a edge line do proprio HDF, a dobra estava onde as
        pontas direitas de RS 30548 e 30398 avancam cem metros para fora e a
        linha volta por cima do caminho de ida.

        O que descreve isso e a PROJECAO da ponta sobre o eixo: se ela recua,
        a edge line anda para tras, e densificada ela se cruza. Condicao mais
        apertada que a corda, e que converge -- encolher puxa a ponta para o
        eixo, e sobre o eixo a ordem e monotona por construcao.
        """
        pr = np.array([eixo.project(Point(q)) for q in pts])
        maus = set()
        for i in np.where(np.diff(pr) <= 0)[0]:
            i = int(i)
            # encolhe quem esta mais para fora; empate, os dois.
            a = abs(secoes[i]["o_esq" if k == 0 else "o_dir"])
            b = abs(secoes[i + 1]["o_esq" if k == 0 else "o_dir"])
            if a >= b:
                maus.add(i)
            if b >= a:
                maus.add(i + 1)
        return maus

    PISO = 3.0 * PASSO           # nenhuma secao encolhe abaixo disto
    n_enc, voltas, restou = 0, 0, 0
    while voltas < 40:
        voltas += 1
        maus = set()
        for lado, k in (("me", 0), ("md", 1)):
            pts = [tuple(s["cut"][k]) for s in secoes]
            for i in dobras(pts) | recuos(pts, k):
                maus.add((i, lado))
        if not maus:
            break
        mexeu = False
        for i, lado in maus:
            s = secoes[i]
            # encolhe a partir da largura EFETIVA (`o_esq`/`o_dir`), e nao da
            # pedida: onde o MDT falta, a pedida ja e maior que a real e
            # multiplica-la nao mexeria na cutline -- laco eterno.
            real = s["o_esq"] if lado == "me" else -s["o_dir"]
            novo_v = max(min(s[lado], real) * 0.85, PISO)
            if novo_v < s[lado] - 1e-6:
                s[lado] = novo_v
                mexeu = True
                n_enc += 1
            montar(s, s["me"], s["md"])
        secoes = [s for s in secoes if s.get("sta") is not None]
        if not mexeu:
            restou = len(maus)
            break
    else:
        restou = len(maus)
    pts_e = [tuple(s["cut"][0]) for s in secoes]
    pts_d = [tuple(s["cut"][1]) for s in secoes]
    sobra = (len(dobras(pts_e)) + len(dobras(pts_d))
             + len(recuos(pts_e, 0)) + len(recuos(pts_d, 1)))
    print(f"edge line: {n_enc} encolhimentos em {voltas} passadas   "
          f"cruzamentos+recuos restantes {sobra}"
          + ("   (no piso de %.0f m)" % PISO if sobra else ""))

    # ---- talvegue: cru, ou isotonico
    zt = np.array([s["zt"] for s in secoes])
    cru = zt.copy()
    if a.monotono:
        novo = isotonica(zt)
        for s, z0, z1 in zip(secoes, cru, novo):
            d = z1 - z0
            st = s["sta"]
            m = (st >= s["lb"]) & (st <= s["rb"])
            if m.any():
                prof = s["z"][m].max() - s["z"][m]
                pmax = prof.max()
                peso = prof / pmax if pmax > 1e-9 else np.zeros_like(prof)
                s["z"][m] = np.round(s["z"][m] + d * peso, 2)
            s["zt"] = float(s["z"].min())
        zt = np.array([s["zt"] for s in secoes])
        print(f"talvegue: regressao isotonica aplicada   "
              f"ajuste mediano {np.median(np.abs(novo-cru)):.2f} m   "
              f"max {np.abs(novo-cru).max():.2f} m")

    # ---- as ADOTADAS entram agora, inteiras, na ordem da estacao. Depois dos
    # filtros de largura e da isotonica de proposito: largura levantada nao e
    # ruido a filtrar, e batimetria levantada nao se "corrige" por regressao.
    if prontas:
        secoes = sorted(secoes + prontas, key=lambda q: q["s"])

    # ---- NA ZONA ADOTADA, O PENTE DE LAMINA CAI. Entre duas secoes
    # levantadas vizinhas o MDT ainda gerava secoes de lamina (fundo ~0,0)
    # alternando com as dragadas (fundo -10,8): uma falsa soleira de 10 m a
    # cada 150 m -- um pente que nenhum solver engole. Secao do MDT cujo
    # fundo esta mais de 2 m ACIMA de ambas as levantadas vizinhas (a menos
    # de 2,5 km) e lamina sobre agua funda, nao leito: cai. O mesmo vale para
    # o rabo alem da ultima levantada (era a restinga de 224 m recebendo a
    # mare no lugar do canal do porto).
    if prontas:
        pr = sorted((q["s"], q["zt"]) for q in secoes if q.get("pronta"))
        pr_s = np.array([q[0] for q in pr])
        pr_z = np.array([q[1] for q in pr])
        mantem, n_pente = [], 0
        for q in secoes:
            if q.get("pronta"):
                mantem.append(q)
                continue
            # so na ZONA DE MARE: rio acima a secao de lamina carrega a
            # varzea real do lidar, e o fundo dela e ancorado depois pelo
            # `batimetria.py aplicar` -- cortar la jogaria fora a planicie.
            if q["zt"] >= 2.0:
                mantem.append(q)
                continue
            i = int(np.searchsorted(pr_s, q["s"]))
            cima = i - 1 if i > 0 else None
            baixo = i if i < len(pr_s) else None
            cai = False
            if cima is not None and baixo is not None:
                if (pr_s[baixo] - pr_s[cima] <= 2500.0
                        and q["zt"] > max(pr_z[cima], pr_z[baixo]) + 2.0):
                    cai = True
            elif cima is not None:
                if (q["s"] - pr_s[cima] <= 1200.0
                        and q["zt"] > pr_z[cima] + 2.0):
                    cai = True
            elif baixo is not None:
                if (pr_s[baixo] - q["s"] <= 1200.0
                        and q["zt"] > pr_z[baixo] + 2.0):
                    cai = True
            if cai:
                n_pente += 1
            else:
                mantem.append(q)
        if n_pente:
            print(f"   pente de lamina na zona adotada: {n_pente} secao(oes) "
                  "do MDT descartada(s) entre/junto a secoes levantadas")
        secoes = mantem

    # ---- PONTA DEGENERADA NAO RECEBE CONTORNO. Na boca do porto sobrava uma
    # secao de 76 m de largura e 9 cm de fundo -- um respingo de borda d'agua
    # com cobertura suficiente para passar no criterio, mas sem nenhuma
    # representatividade -- e era NELA que a mare entraria. Secao de ponta com
    # menos de 30% da largura mediana das 10 vizinhas cai; as internas ficam,
    # que estreitamento no meio do rio pode ser real.
    def _ext(q):
        return float(q["sta"][-1] - q["sta"][0])
    for lado_, idx_ in (("jusante", -1), ("montante", 0)):
        while len(secoes) > 12:
            viz = (secoes[-11:-1] if idx_ == -1 else secoes[1:11])
            med = float(np.median([_ext(q) for q in viz]))
            if _ext(secoes[idx_]) >= 0.3 * med:
                break
            s_ = secoes.pop(idx_)
            print(f"   ponta de {lado_} degenerada descartada: RS "
                  f"{s_['rs']:.1f} ({_ext(s_):.0f} m de largura, vizinhas "
                  f"{med:.0f} m)")

    # ---- escreve
    os.makedirs(a.saida, exist_ok=True)
    nome = os.path.basename(a.saida.rstrip("/\\"))
    g = os.path.join(a.saida, f"{nome}.g01")
    l = [f"Geom Title={nome}", "Program Version=7.01"]
    P = np.vstack([np.vstack(s["cut"]) for s in secoes])
    l.append("Viewing Rectangle= %.2f , %.2f , %.2f , %.2f "
             % (P[:, 0].min(), P[:, 0].max(), P[:, 1].max(), P[:, 1].min()))
    l.append("Spatial Reference System=" + WKT)
    l.append("")
    l.append(f"River Reach={a.rio:<16.16},{a.reach:<16.16}")
    # O EIXO SEGUE ALEM DA PRIMEIRA E DA ULTIMA SECAO. Elas caem exatamente
    # sobre as pontas da polilinha e ali a cutline TOCA o eixo em vez de
    # cruza-lo: o HEC-RAS recusa com "XS doesn't intersect the associated
    # Reach", e junto vem "XS intersects < 2 banklines" e "XS must intersect
    # exactly one Reach" -- seis Fatal vindos de duas secoes. Prolongar o eixo
    # nao move secao nenhuma e nao muda comprimento de trecho, que sai de
    # `Length Ch`.
    cc = np.asarray(eixo.coords, float)
    d0 = cc[0] - cc[1]
    d0 /= max(float(np.hypot(*d0)), 1e-9)
    d1 = cc[-1] - cc[-2]
    d1 /= max(float(np.hypot(*d1)), 1e-9)
    folga = max(EXTRA_EIXO, a.dx)
    cc = np.vstack([cc[0] + folga * d0, cc, cc[-1] + folga * d1])
    c = [tuple(q) for q in cc]
    l.append(f"Reach XY= {len(c)} ")
    ss = [f"{x:16.4f}{y:16.4f}" for x, y in c]
    for k in range(0, len(ss), 2):
        l.append("".join(ss[k:k + 2]))
    l.append("Rch Text X Y=0,0,0,0")
    l.append("")
    for i, s in enumerate(secoes):
        d = (round(float(secoes[i + 1]["s"] - s["s"]), 2)
             if i + 1 < len(secoes) else 0.0)
        l.append(f"Type RM Length L Ch R = 1 ,{s['rs']:.2f},"
                 f"{d:8.2f},{d:8.2f},{d:8.2f}")
        l.append(f"Bank Sta={s['lb']:.2f},{s['rb']:.2f}")
        l.append("XS GIS Cut Line= 2")
        l.append("".join(f"{q[0]:16.2f}{q[1]:16.2f}" for q in s["cut"]))
        l.append(f"#Sta/Elev= {len(s['sta'])} ")
        pf = [f"{x:8.2f}{y:8.2f}" for x, y in zip(s["sta"], s["z"])]
        for k in range(0, len(pf), 5):
            l.append("".join(pf[k:k + 5]))
        l.append("#Mann= 3 , 0 , 0 ")
        l.append(f"{s['sta'][0]:8.2f}{N_PLANICIE:8.3f}{0:8d}"
                 f"{s['lb']:8.2f}{N_CALHA:8.3f}{0:8d}"
                 f"{s['rb']:8.2f}{N_PLANICIE:8.3f}{0:8d}")
        # HTAB ANCORADO NO INVERT DA CALHA, e nao no minimo da secao. O
        # HEC-RAS compara com o ponto mais baixo ENTRE AS MARGENS; ancorar no
        # minimo da secao poe o HTab abaixo do invert sempre que houver ponto
        # mais fundo na planicie, e ele reseta a tabela sozinho.
        mcal = (s["sta"] >= s["lb"]) & (s["sta"] <= s["rb"])
        z_inv = float(s["z"][mcal].min()) if mcal.any() else s["zt"]
        l.append(f"XS HTab Starting El and Incr={z_inv+0.02:.2f},0.100, 500 ")
        l.append("XS HTab Horizontal Distribution=-1,-1,-1")
        l.append("XS Rating Curve= 0 ,0")
        l.append("Exp/Cntr=0.3,0.1")
        l.append("")
    escrever(g, "\n".join(l))

    # ---- o que o terreno entregou, sem maquiagem
    lc = np.array([s["rb"] - s["lb"] for s in secoes])
    ls = np.array([s["sta"][-1] for s in secoes])
    npt = np.array([len(s["sta"]) for s in secoes])
    # recalcula do `secoes` FINAL: as adotadas do legado entraram depois da
    # isotonica, e medir dz num array de antes do merge quebrava o relatorio
    zt = np.array([s["zt"] for s in secoes])
    dz = np.diff(zt)
    ch = np.array([secoes[i + 1]["s"] - secoes[i]["s"]
                   for i in range(len(secoes) - 1)])
    decl = np.abs(dz) / np.maximum(ch, 1e-9)
    print(f"\ngeometria: {g}")
    print(f"   talvegue     : {zt.min():.2f} a {zt.max():.2f} m")
    print(f"   sobem p/ jusante: {int((dz > 1e-9).sum())} de {len(dz)}"
          f"   (o terreno traz isto; --monotono trata)")
    print(f"   pares de cota igual: {int((np.abs(dz) < 0.005).sum())}")
    print(f"   declividade  : mediana {np.median(decl):.5f}   "
          f">2% em {int((decl > 0.02).sum())} trechos   max {decl.max():.3f}")
    print(f"   calha        : mediana {np.median(lc):.0f} m   "
          f"p10 {np.percentile(lc,10):.0f}   p90 {np.percentile(lc,90):.0f}")
    print(f"   secao        : mediana {np.median(ls):.0f} m   "
          f"max {ls.max():.0f} m")
    print(f"   pontos/secao : mediana {np.median(npt):.0f}   max {npt.max()}"
          "   (limite do HEC-RAS: 500)")
    return g


if __name__ == "__main__":
    main()
