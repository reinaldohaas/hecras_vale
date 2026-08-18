# -*- coding: utf-8 -*-
"""
Todo erro detectado vira uma imagem.

POR QUE ISTO EXISTE. Os erros que mais custaram tempo neste projeto eram
invisiveis em numero e obvios em figura. O Pilot Channel que achatou as 1.232
secoes de um modelo passou por toda a auditoria sem acusar nada -- area, altura
e largura continuavam plausiveis; quem viu foi o usuario, olhando um perfil. A
lamina de 2 cm que travou o solver no Iraputa estava em duas linhas de log que
diziam so "WSEL 501.50": para saber que aquilo eram 2 cm acima do leito era
preciso ir buscar o leito em outro arquivo.

POR QUE SVG ESCRITO A MAO, e nao matplotlib. Porque matplotlib quebrou. No
miniforge desta maquina qualquer desenho derruba o processo com 0xc06d007f
(delay-load de DLL falhou), sem traceback, exit 127 -- o mesmo modo de falhar
do rasterio.sample() que ja custou um dia neste projeto. Uma regra do tipo
"todo erro vira figura" nao pode depender de uma biblioteca que some. SVG e
texto, o projeto ja escreve .g01 e .rasmap a mao, e assim a figura sai em
qualquer ambiente que rode o resto do programa. Abre em qualquer navegador.

A figura de secao tem tres paineis, e os tres respondem perguntas diferentes:

  PERFIL LONGITUDINAL: o problema e local ou e do rio inteiro? Um degrau
  isolado e uma secao ruim; uma escada e o condicionamento.

  SECAO INTEIRA: o que foi feito com esta secao? E onde aparecem calha
  inventada, parede vertical, achatamento e escavacao absurda.

  RECORTE NO ENTALHE: a lamina que importa costuma ter centimetros numa secao
  de centenas de metros. No grafico inteiro ela e uma linha em cima de outra.
"""
import os

import numpy as np

# ---------------------------------------------------------------- desenho
L_TERRENO = "#9a8f80"
L_MODELO = "#1f4e79"
L_ALVO = "#c0392b"
L_AGUA = "#2e86c1"
L_EIXO = "#5b6770"
L_GRADE = "#e3e6e8"


def _esc(v, v0, v1, p0, p1):
    """Mapeia dado -> pixel. Sem faixa, cai no meio (evita divisao por zero)."""
    if v1 - v0 == 0:
        return (p0 + p1) / 2.0
    return p0 + (np.asarray(v, float) - v0) * (p1 - p0) / (v1 - v0)


def _ticks(v0, v1, n=5):
    """Valores redondos dentro da faixa."""
    if v1 <= v0:
        return [v0]
    bruto = (v1 - v0) / max(n, 1)
    mag = 10.0 ** np.floor(np.log10(bruto))
    passo = min([m * mag for m in (1, 2, 2.5, 5, 10)],
                key=lambda p: abs(p - bruto))
    ini = np.ceil(v0 / passo) * passo
    return [ini + k * passo for k in range(int((v1 - ini) / passo) + 1)]


def _txt(x, y, s, tam=11, cor="#222", anc="start", peso="normal"):
    s = str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return ('<text x="%.1f" y="%.1f" font-size="%s" fill="%s" '
            'text-anchor="%s" font-weight="%s" '
            'font-family="DejaVu Sans, Segoe UI, sans-serif">%s</text>'
            % (x, y, tam, cor, anc, peso, s))


def _linha(xs, ys, cor, largura=1.5, tracejado=None):
    if len(xs) == 0:
        return ""
    d = " ".join(("M" if i == 0 else "L") + "%.1f,%.1f" % (x, y)
                 for i, (x, y) in enumerate(zip(xs, ys)))
    ds = ' stroke-dasharray="%s"' % tracejado if tracejado else ""
    return ('<path d="%s" fill="none" stroke="%s" stroke-width="%s"%s '
            'stroke-linejoin="round"/>' % (d, cor, largura, ds))


def _area(xs, y_sup, y_inf, cor, opac=0.16):
    """Faixa preenchida entre duas curvas.

    Existe porque a versao anterior desenhava as duas linhas e mais nada: o
    usuario olhou a figura do Iraputa e disse "nao vejo nada de errado com esse
    perfil" -- com o leito 13,5 m abaixo do talvegue bem ali, em duas linhas
    que, sem o vao pintado e sem o numero, so pareciam um perfil liso embaixo
    de um terreno serrilhado. Duas linhas separadas nao leem como um erro; uma
    faixa vermelha de 13 m le.
    """
    if len(xs) == 0:
        return ""
    ida = " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L", x, y)
                   for i, (x, y) in enumerate(zip(xs, y_sup)))
    volta = " ".join("L%.1f,%.1f" % (x, y)
                     for x, y in zip(xs[::-1], y_inf[::-1]))
    return ('<path d="%s %s Z" fill="%s" fill-opacity="%s" stroke="none"/>'
            % (ida, volta, cor, opac))


class _Painel:
    """Uma caixa com eixos, em pixels, e os dados que cabem nela.

    O conteudo e RECORTADO na caixa. Sem isso a margem alta de uma secao sai do
    painel e risca os vizinhos: na primeira versao eram 202 pontos fora da area
    util numa figura so, e o recorte do entalhe aparecia atravessado por linhas
    do perfil longitudinal.
    """

    _n = 0

    def __init__(self, x0, y0, x1, y1, titulo, rot_x, rot_y):
        _Painel._n += 1
        self.id = "cp%d" % _Painel._n
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        self.titulo, self.rot_x, self.rot_y = titulo, rot_x, rot_y
        self.dx0 = self.dx1 = self.dy0 = self.dy1 = None
        self.corpo = []

    def faixa(self, xs, ys, folga=0.06):
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        xs, ys = xs[np.isfinite(xs)], ys[np.isfinite(ys)]
        if len(xs) == 0 or len(ys) == 0:
            return
        mx = (xs.max() - xs.min()) * folga or 1.0
        my = (ys.max() - ys.min()) * folga or 1.0
        a, b, c, d = xs.min() - mx, xs.max() + mx, ys.min() - my, ys.max() + my
        self.dx0 = a if self.dx0 is None else min(self.dx0, a)
        self.dx1 = b if self.dx1 is None else max(self.dx1, b)
        self.dy0 = c if self.dy0 is None else min(self.dy0, c)
        self.dy1 = d if self.dy1 is None else max(self.dy1, d)

    def px(self, x):
        return _esc(x, self.dx0, self.dx1, self.x0, self.x1)

    def py(self, y):
        return _esc(y, self.dy0, self.dy1, self.y1, self.y0)   # y p/ cima

    def linha(self, x, y, cor, lw=1.5, tracejado=None):
        self.corpo.append(_linha(self.px(np.asarray(x, float)),
                                 self.py(np.asarray(y, float)), cor, lw,
                                 tracejado))

    def hline(self, y, cor, lw=1.2, tracejado=None):
        p = self.py(y)
        self.corpo.append(_linha([self.x0, self.x1], [p, p], cor, lw,
                                 tracejado))

    def vline(self, x, cor, lw=1.2, tracejado=None):
        p = self.px(x)
        self.corpo.append(_linha([p, p], [self.y0, self.y1], cor, lw,
                                 tracejado))

    def ponto(self, x, y, cor, r=4.5):
        self.corpo.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s"/>'
                          % (self.px(x), self.py(y), r, cor))

    def nota(self, x, y, s, cor="#222", dx=6, dy=-6, tam=10.5):
        self.corpo.append(_txt(self.px(x) + dx, self.py(y) + dy, s, tam, cor))

    def area(self, x, y_sup, y_inf, cor, opac=0.16):
        self.corpo.append(_area(self.px(np.asarray(x, float)),
                                self.py(np.asarray(y_sup, float)),
                                self.py(np.asarray(y_inf, float)), cor, opac))

    def svg(self):
        out = ['<rect x="%s" y="%s" width="%s" height="%s" fill="#fff" '
               'stroke="%s" stroke-width="0.8"/>'
               % (self.x0, self.y0, self.x1 - self.x0, self.y1 - self.y0,
                  L_EIXO)]
        if self.dx0 is None:
            return "".join(out)
        for v in _ticks(self.dx0, self.dx1):
            p = self.px(v)
            if not (self.x0 - 1 <= p <= self.x1 + 1):
                continue
            out.append(_linha([p, p], [self.y0, self.y1], L_GRADE, 0.7))
            out.append(_txt(p, self.y1 + 15, "%g" % v, 10, L_EIXO, "middle"))
        for v in _ticks(self.dy0, self.dy1):
            p = self.py(v)
            if not (self.y0 - 1 <= p <= self.y1 + 1):
                continue
            out.append(_linha([self.x0, self.x1], [p, p], L_GRADE, 0.7))
            out.append(_txt(self.x0 - 6, p + 3.5, "%g" % v, 10, L_EIXO, "end"))
        out.append('<clipPath id="%s"><rect x="%s" y="%s" width="%s" '
                   'height="%s"/></clipPath>'
                   % (self.id, self.x0, self.y0, self.x1 - self.x0,
                      self.y1 - self.y0))
        out.append('<g clip-path="url(#%s)">' % self.id)
        out += self.corpo
        out.append("</g>")
        out.append(_txt(self.x0, self.y0 - 8, self.titulo, 12, "#222",
                        peso="600"))
        out.append(_txt((self.x0 + self.x1) / 2, self.y1 + 32, self.rot_x,
                        10.5, L_EIXO, "middle"))
        out.append('<g transform="translate(%s,%s) rotate(-90)">%s</g>'
                   % (self.x0 - 40, (self.y0 + self.y1) / 2,
                      _txt(0, 0, self.rot_y, 10.5, L_EIXO, "middle")))
        return "".join(out)


def _legenda(x, y, itens):
    out = []
    for i, (cor, rot) in enumerate(itens):
        yy = y + i * 16
        out.append(_linha([x, x + 22], [yy, yy], cor, 2.2))
        out.append(_txt(x + 28, yy + 3.5, rot, 10.5, "#333"))
    return "".join(out)


def _gravar(partes, destino):
    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    with open(destino, "w", encoding="utf-8") as f:
        f.write("\n".join(partes))
    return destino


# ------------------------------------------------------------ figura de secao
def _serie(v, chave):
    """Perfil longitudinal de um rio, do topo para a foz."""
    w = sorted(v, key=lambda x: -x["rs"])
    rs = np.array([x["rs"] for x in w], float)
    if chave == "leito":
        z = np.array([float(np.min(x["z"])) for x in w], float)
    else:
        z = np.array([float(x.get("z_terreno", np.min(x["z"]))) for x in w],
                     float)
    return rs, z


def desenhar(estado, rio, rs_alvo, destino, motivo="", wsel=None, log=print):
    """Tres paineis para um problema em rio/RS.

    wsel: cota da lamina, quando o solver informou. E o que transforma
    "o modelo abortou aqui" em "o modelo abortou com 2 cm de agua aqui".
    """
    cru = (estado.get("xs") or {}).get(rio) or []
    pronto = (estado.get("xs_pronto") or {}).get(rio) or []
    if not pronto:
        return None
    x = min(pronto, key=lambda a: abs(a["rs"] - rs_alvo))
    xc = min(cru, key=lambda a: abs(a["rs"] - x["rs"])) if cru else None

    W, H = 1380, 620
    partes = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
              'viewBox="0 0 %d %d">' % (W, H, W, H),
              '<rect width="%d" height="%d" fill="#fbfbfa"/>' % (W, H),
              _txt(24, 30, "%s   RS %.2f" % (rio, x["rs"]), 17, "#111",
                   peso="700")]
    if motivo:
        partes.append(_txt(24, 50, motivo, 12, L_ALVO))

    # ------------------------------------------------ 1. perfil longitudinal
    p1 = _Painel(90, 80, 700, 540, "%s -- perfil longitudinal" % rio,
                 "RS (km acima da foz)", "cota (m)")
    rs_t, z_t = _serie(cru or pronto, "terreno")
    rs_m, z_m = _serie(pronto, "leito")
    p1.faixa(rs_m / 1000.0, z_m)
    p1.faixa(rs_t / 1000.0, z_t)
    # o VAO entre o terreno e o leito, pintado. E ele que diz se o modelo esta
    # representando o rio ou inventando um canion; sem preenchimento e sem
    # numero, duas linhas separadas passam por "perfil liso".
    n = min(len(rs_t), len(rs_m))
    if n >= 2 and np.allclose(rs_t[-n:], rs_m[-n:]):
        gap = z_t[-n:] - z_m[-n:]
        p1.area(rs_t[-n:] / 1000.0, z_t[-n:], z_m[-n:], L_ALVO, 0.15)
        p1.corpo.append(_txt(p1.x0 + 12, p1.y0 + 52,
                             "leito %.1f m abaixo do talvegue (mediana do "
                             "rio; max %.1f m)"
                             % (np.median(gap), gap.max()), 11.5, L_ALVO,
                             peso="600"))
    p1.linha(rs_t / 1000.0, z_t, L_TERRENO, 1.2)
    p1.linha(rs_m / 1000.0, z_m, L_MODELO, 1.8)
    p1.vline(x["rs"] / 1000.0, L_ALVO, 1.1, "4,3")
    zi = float(np.min(x["z"]))
    p1.ponto(x["rs"] / 1000.0, zi, L_ALVO)
    i = int(np.argmin(np.abs(rs_m - x["rs"])))
    if 0 < i < len(rs_m) - 1:
        dxm = rs_m[i - 1] - rs_m[i + 1]
        s = 100.0 * (z_m[i - 1] - z_m[i + 1]) / dxm if dxm > 0 else 0.0
        p1.nota(x["rs"] / 1000.0, zi, "declividade local %.3f%%" % s, L_ALVO,
                dx=8, dy=-10)
    partes.append(p1.svg())
    partes.append(_legenda(110, 100, [(L_TERRENO, "terreno (talvegue)"),
                                      (L_MODELO, "leito do modelo")]))

    # ------------------------------------------------------ 2. secao inteira
    sta = np.asarray(x["sta"], float)
    z = np.asarray(x["z"], float)
    p2 = _Painel(800, 80, 1330, 300, "secao RS %.2f -- inteira" % x["rs"],
                 "estaca (m)", "cota (m)")
    p2.faixa(sta, z)
    if xc is not None:
        p2.faixa(np.asarray(xc["sta"], float), np.asarray(xc["z"], float))
        p2.linha(xc["sta"], xc["z"], L_TERRENO, 1.2)
    p2.linha(sta, z, L_MODELO, 1.6)
    for b in (x.get("lb"), x.get("rb")):
        if b is not None:
            p2.vline(float(b), "#95a5a6", 0.8, "2,3")
    if wsel is not None:
        p2.hline(float(wsel), L_AGUA, 1.3)
    # o tamanho do que foi cavado, em numero: profundidade abaixo do terreno,
    # largura no topo da escavacao e largura do fundo
    if xc is not None:
        zc_ = np.asarray(xc["z"], float)
        sc_ = np.asarray(xc["sta"], float)
        if len(zc_) == len(z):
            p2.area(sta, zc_, np.minimum(z, zc_), L_ALVO, 0.15)
            cav = z < zc_ - 0.01
            if cav.any():
                fun = z <= float(np.min(z)) + 0.05
                p2.corpo.append(_txt(
                    p2.x0 + 10, p2.y0 + 18,
                    "escavado %.1f m abaixo do terreno   |   topo %.0f m   |  "
                    " fundo %.0f m"
                    % ((zc_ - z).max(), sc_[cav].max() - sc_[cav].min(),
                       sta[fun].max() - sta[fun].min()), 11, L_ALVO,
                    peso="600"))
    partes.append(p2.svg())

    # -------------------------------------------------- 3. recorte no entalhe
    leito = float(np.min(z))
    fundo = z <= leito + 0.05
    p3 = _Painel(800, 370, 1330, 540, "recorte no entalhe", "estaca (m)",
                 "cota (m)")
    if fundo.any():
        c0, c1 = float(sta[fundo].min()), float(sta[fundo].max())
        larg = max(c1 - c0, 1.0)
        alt = max(float(wsel) - leito, 0.0) if wsel is not None else 0.0
        jan = (sta >= c0 - 3 * larg) & (sta <= c1 + 3 * larg)
        p3.faixa(sta[jan], np.clip(z[jan], leito - 0.5,
                                   leito + max(2.5, 8 * alt)))
        if xc is not None:
            sc = np.asarray(xc["sta"], float)
            zc = np.asarray(xc["z"], float)
            jc = (sc >= c0 - 3 * larg) & (sc <= c1 + 3 * larg)
            p3.linha(sc[jc], zc[jc], L_TERRENO, 1.2)
        p3.linha(sta[jan], z[jan], L_MODELO, 1.8)
        if wsel is not None:
            p3.hline(float(wsel), L_AGUA, 1.6)
            p3.corpo.append(_txt(p3.x0 + 8, p3.py(float(wsel)) - 6,
                                 "lamina %.1f cm" % (100 * alt), 12, L_AGUA,
                                 peso="600"))
        p3.corpo.append(_txt(p3.x1 - 8, p3.y0 + 16,
                             "fundo de %.0f m   |   secao de %.0f m"
                             % (larg, sta[-1]), 10.5, L_EIXO, "end"))
    partes.append(p3.svg())
    partes.append("</svg>")
    return _gravar(partes, destino)


# -------------------------------------------------------- figura de contorno
def hidrograma(rio, serie, q_inicial, destino, intervalo_h=1.0, motivo=""):
    """Hidrograma de contorno com a condicao inicial marcada.

    Existe por causa de um erro que nenhuma checagem de geometria pegaria: a
    condicao inicial e o inicio do hidrograma sao escritos por funcoes
    diferentes (hidrologia.inicial e hidrologia.series) e podem discordar. Na
    rodada de 18/08/2026 discordavam por 19 a 38 vezes em TODOS os rios -- o
    modelo era inicializado cheio e esvaziado no primeiro passo por um contorno
    vinte vezes menor. Em numero sao duas listas em arquivos distintos; em
    figura e um degrau na origem que nao da para nao ver.
    """
    q = np.asarray(serie, float)
    t = np.arange(len(q)) * float(intervalo_h)
    W, H = 900, 470
    partes = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
              'viewBox="0 0 %d %d">' % (W, H, W, H),
              '<rect width="%d" height="%d" fill="#fbfbfa"/>' % (W, H),
              _txt(24, 30, "%s -- hidrograma de contorno" % rio, 17, "#111",
                   peso="700")]
    if motivo:
        partes.append(_txt(24, 50, motivo, 12, L_ALVO))
    p = _Painel(95, 85, 860, 390, "vazao no contorno de montante",
                "tempo (h)", "vazao (m3/s)")
    p.faixa(t, q)
    p.faixa([0.0, 0.0], [0.0, float(q_inicial) * 1.15])
    p.linha(t, q, L_MODELO, 1.8)
    p.hline(float(q_inicial), L_ALVO, 1.4, "5,4")
    p.ponto(0.0, float(q[0]), L_AGUA, 5)
    razao = float(q_inicial) / max(float(q[0]), 1e-9)
    p.corpo.append(_txt(p.px(t[-1] * 0.32), p.py(float(q_inicial)) - 8,
                        "condicao inicial do .u01:  %.1f m3/s" % q_inicial,
                        12, L_ALVO, peso="600"))
    p.corpo.append(_txt(p.px(t[-1] * 0.32), p.py(float(q[0])) + 18,
                        "hidrograma em t=0:  %.2f m3/s   (%.0f vezes menor)"
                        % (q[0], razao), 12, L_AGUA, peso="600"))
    partes.append(p.svg())
    partes.append(_legenda(120, 105, [(L_MODELO, "hidrograma do contorno"),
                                      (L_ALVO, "condicao inicial")]))
    partes.append("</svg>")
    return _gravar(partes, destino)


# --------------------------------------------------------------- em lote
def _rio_rs(alvo):
    """Separa 'Itajai_Acu RS 75' em ('Itajai_Acu', 75.0)."""
    if " RS " not in str(alvo):
        return str(alvo), None
    rio, rs = str(alvo).split(" RS ", 1)
    try:
        return rio, float(rs)
    except ValueError:
        return rio, None


def dos_problemas(estado, problemas, pasta, log=print, limite=12):
    """Uma figura por problema que aponte para um RS.

    Problema de rio inteiro (escavacao mediana, rio desconectado) nao vira
    figura de secao: nao ha RS para marcar, e figura sem alvo e ruido.
    """
    feitas = []
    for p in problemas[:limite]:
        rio, rs = _rio_rs(getattr(p, "alvo", ""))
        if rs is None or rio not in (estado.get("xs_pronto") or {}):
            continue
        nome = "erro_%s_%s_%.0f.svg" % (
            str(getattr(p, "tipo", "problema")).replace(" ", "_"), rio, rs)
        c = desenhar(estado, rio, rs, os.path.join(pasta, nome),
                     motivo=str(p), log=log)
        if c:
            feitas.append(c)
    if feitas:
        log("      %d figura(s) do erro em %s" % (len(feitas), pasta))
    return feitas


def do_solver(estado, pasta, avisos, log=print):
    """Figuras dos pontos que o SOLVER acusou, com a lamina que ele calculou.

    avisos: lista de (rio, rs, wsel, motivo). Sai do log de computacao -- e o
    unico lugar onde existe a cota da lamina, que e o que revela um modelo
    rodando seco.
    """
    feitas = []
    for rio, rs, wsel, motivo in avisos:
        nome = "solver_%s_%.0f.svg" % (rio, float(rs))
        c = desenhar(estado, rio, float(rs), os.path.join(pasta, nome),
                     motivo=motivo, wsel=wsel, log=log)
        if c:
            feitas.append(c)
    if feitas:
        log("      %d figura(s) do solver em %s" % (len(feitas), pasta))
    return feitas


# ------------------------------------------------------------ planta do 2D
def planta_2d(est, destino, titulo="", log=print):
    """Vista em planta da area 2D, do eixo e dos contornos.

    E a unica figura em que a ESCALA DOS DOIS EIXOS TEM DE SER IGUAL. Um
    painel que estica x e y de forma independente mostra o buffer de 1.500 m
    como uma faixa de larguras diferentes ao longo do rio, e o erro que se
    procura numa planta -- tampa torta, contorno em cima do outro, area
    partida em dois -- e justamente de forma.
    """
    pol, eixo, bcs = est["area"], est["eixo"], est["bcs"]
    px, py = np.array(pol.exterior.coords).T
    W, H = 1180, 900
    p = ["<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' "
         "viewBox='0 0 %d %d'><rect width='%d' height='%d' fill='#fff'/>"
         % (W, H, W, H, W, H)]
    p.append(_txt(70, 34, titulo or ("area 2D %s" % est.get("nome_area", "")),
                  17, "#111", peso="700"))
    p.append(_txt(70, 56, "%.0f km2   |   %d vertices no perimetro   |   "
                  "%d contornos   |   eixo %.0f km"
                  % (pol.area / 1e6, len(px), len(bcs), eixo.length / 1000.0),
                  12, "#555"))

    a = _Painel(80, 100, W - 60, H - 90, "", "leste (m, UTM 22S)",
                "norte (m, UTM 22S)")
    a.faixa(px, py)
    # ESCALA IGUAL: alarga a faixa mais curta ate bater a razao do painel
    razao = (a.x1 - a.x0) / (a.y1 - a.y0)
    dx, dy = a.dx1 - a.dx0, a.dy1 - a.dy0
    if dx / dy < razao:
        f = (razao * dy - dx) / 2.0
        a.dx0, a.dx1 = a.dx0 - f, a.dx1 + f
    else:
        f = (dx / razao - dy) / 2.0
        a.dy0, a.dy1 = a.dy0 - f, a.dy1 + f

    a.area(px, py, np.full_like(py, a.dy0), "#3b82c4", 0.10)
    a.linha(px, py, "#2c6fbb", 1.4)
    ex, ey = np.array(eixo.coords).T
    a.linha(ex, ey, "#1b5e20", 1.8)
    cor = {"montante": "#c0392b", "foz": "#8e44ad", "lateral": "#e67e22"}
    for b in bcs:
        bx, by = np.array(b["coords"]).T
        a.linha(bx, by, cor[b["tipo"]], 3.4)
    a.ponto(ex[0], ey[0], "#1b5e20", 5.0)
    a.nota(ex[0], ey[0], "nascente", "#1b5e20")
    a.ponto(ex[-1], ey[-1], "#8e44ad", 5.0)
    a.nota(ex[-1], ey[-1], "foz", "#8e44ad")
    p.append(a.svg())
    p.append(_legenda(W - 300, 120, [
        ("#2c6fbb", "perimetro da area 2D"), ("#1b5e20", "eixo do rio"),
        ("#c0392b", "contorno de montante"), ("#8e44ad", "foz (mare)"),
        ("#e67e22", "contornos laterais")]))
    q = sum(float(b["serie"].max()) for b in bcs if "serie" in b)
    p.append(_txt(70, H - 22, "vazao de pico somada nos contornos: %.0f m3/s"
                  "   |   escala igual nos dois eixos" % q, 10, "#888"))
    p.append("</svg>")
    return _gravar(p, destino)
