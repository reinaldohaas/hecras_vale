# -*- coding: utf-8 -*-
"""
Desenho: mapa, perfil e comparacao. Matplotlib puro, sem estado global.

Cada funcao recebe um Axes e desenha nele. A mesma rotina serve para a janela
do PySide6 e para exportar PNG em lote, sem duplicar codigo de plotagem.

RESTRICAO DE PRIMITIVAS -- ler antes de "melhorar" este arquivo. Nesta
instalacao (matplotlib 3.10.9 sobre numpy 2.5.1) varias chamadas comuns
derrubam o processo em codigo nativo: exit 127, sem traceback, levando a janela
junto. Foram medidas uma a uma, com dado sintetico, fora de qualquer
dependencia geoespacial:

    ax.axvline / ax.axhline      derruba
    ax.vlines (LineCollection)   derruba
    ax.annotate com arrowprops   derruba
    figure.tight_layout()        derruba no draw seguinte
    ax.legend(loc="best")        derruba

    ax.plot / ax.text / ax.grid              ok
    ax.imshow / legend com loc fixo          ok
    figure.subplots_adjust                   ok

Entao tudo aqui e feito com ax.plot: linha vertical e um plot de dois pontos, e
a barra de profundidade tambem. Nao e elegante; e o que nao quebra.
"""
import numpy as np

COR_OK = "#2e7d32"
COR_ATENCAO = "#ef6c00"
COR_CRITICA = "#c62828"
COR_INCERTO = "#6a1b9a"
CORES = {"OK": COR_OK, "ATENCAO": COR_ATENCAO,
         "CRITICA": COR_CRITICA, "INCERTO": COR_INCERTO}


def cor_status(st):
    return CORES.get(st, "#546e7a")


def _ajustar(ax, esq=0.13, dir_=0.98, baixo=0.14, cima=0.92):
    try:
        ax.figure.subplots_adjust(left=esq, right=dir_, bottom=baixo, top=cima)
    except Exception:                                        # noqa: BLE001
        pass


def _vertical(ax, x, y0, y1, **kw):
    """Linha vertical em coordenadas de dado (substitui axvline)."""
    return ax.plot([x, x], [y0, y1], **kw)


def _faixa_y(z, folga=0.08):
    z = np.asarray(z, float)
    f = z[np.isfinite(z)]
    if not f.size:
        return 0.0, 1.0
    lo, hi = float(f.min()), float(f.max())
    m = folga * ((hi - lo) or 1.0)
    return lo - m, hi + m


# ------------------------------------------------------------------- mapa
def mapa(ax, dem, eixo=None, secoes=(), selecionada=None, reducao=4):
    ax.clear()
    b = dem.bounds
    # masked em vez de NaN cru: o NoData fica transparente em vez de virar
    # cota falsa no extremo da escala de cores
    ax.imshow(np.ma.masked_invalid(dem.z[::reducao, ::reducao]),
              extent=(b.left, b.right, b.bottom, b.top), origin="upper",
              cmap="terrain", interpolation="nearest")
    if eixo is not None:
        for l in eixo.linhas:
            x, y = l.xy
            ax.plot(x, y, "-", color="#0d47a1", lw=1.4, zorder=3)
    for s in secoes:
        x, y = s.geom.xy
        st = s.qc.status if s.qc else None
        ax.plot(x, y, "-", lw=1.0, zorder=4,
                color=cor_status(st) if st else "#37474f", alpha=0.9)
    if selecionada is not None:
        x, y = selecionada.geom.xy
        ax.plot(x, y, "-", color="k", lw=3.0, zorder=6)
        ax.plot(x, y, "-", color="#ffeb3b", lw=1.8, zorder=7)
        i = selecionada.i_talvegue
        if i is not None and selecionada.xs is not None:
            ax.plot([selecionada.xs[i]], [selecionada.ys[i]], "o", ms=7,
                    mfc="#d50000", mec="k", zorder=8)
    ax.set_aspect("equal")
    ax.set_xlabel("E (m)")
    ax.set_ylabel("N (m)")
    _ajustar(ax)


def zoom_secao(ax, secao, folga=0.6):
    """Enquadra o mapa na secao selecionada."""
    x, y = secao.geom.xy
    dx = (max(x) - min(x)) or 100.0
    dy = (max(y) - min(y)) or 100.0
    m = folga * max(dx, dy)
    cx, cy = 0.5 * (max(x) + min(x)), 0.5 * (max(y) + min(y))
    ax.set_xlim(cx - dx / 2 - m, cx + dx / 2 + m)
    ax.set_ylim(cy - dy / 2 - m, cy + dy / 2 + m)


# ------------------------------------------------------------------ perfil
def perfil(ax, secao, mostrar_candidatos=True):
    ax.clear()
    if not secao.valida:
        ax.text(0.5, 0.5, "sem terreno valido nesta secao",
                ha="center", va="center", transform=ax.transAxes)
        _ajustar(ax)
        return
    sta, z = secao.sta, secao.z
    y0, y1 = _faixa_y(z)
    ax.plot(sta, z, "-", color="#37474f", lw=1.3, label="terreno")

    buraco = ~np.isfinite(z)
    if buraco.any():
        ax.plot(sta[buraco], np.full(int(buraco.sum()), y0), "|",
                color="#b71c1c", ms=8, label=f"NoData ({int(buraco.sum())})")

    t = secao.talvegue or {}
    i = t.get("i_talvegue")
    ia = t.get("i_min_abs")
    ip = t.get("i_principal")
    # os tres aparecem separados de proposito: o usuario precisa poder ver que
    # o minimo absoluto e o talvegue escolhido nao sao a mesma coisa, e
    # discordar da escolha do algoritmo
    if ia is not None and ia != i:
        ax.plot([sta[ia]], [z[ia]], "v", color=COR_INCERTO, ms=8,
                label=f"minimo absoluto ({z[ia]:.2f} m)")
    if ip is not None and ip not in (i, ia):
        ax.plot([sta[ip]], [z[ip]], "s", color="#00838f", ms=7,
                label="depressao principal")
    if mostrar_candidatos:
        c = [(x["sta"], x["z"]) for x in (t.get("candidatos") or [])[1:6]]
        if c:
            ax.plot([p[0] for p in c], [p[1] for p in c], ".",
                    color="#9e9e9e", ms=6, zorder=1, label="outros candidatos")
    if i is not None:
        _vertical(ax, sta[i], y0, y1, color="#d50000", lw=0.8, ls=":")
        ax.plot([sta[i]], [z[i]], "o", color="#d50000", ms=9,
                label=f"talvegue ({z[i]:.2f} m)")
    if secao.sta_eixo is not None:
        _vertical(ax, secao.sta_eixo, y0, y1, color="#0d47a1", lw=1.2,
                  ls="--", label="eixo do rio")

    prof = secao.profundidade_relativa
    if np.isfinite(prof) and i is not None and prof > 0:
        _vertical(ax, sta[i], z[i], z[i] + prof, color="#455a64", lw=1.2)
        ax.text(sta[i], z[i] + prof / 2, f"  {prof:.2f} m",
                va="center", fontsize=8, color="#455a64")

    st = secao.qc.status if secao.qc else ""
    ax.set_ylim(y0, y1)
    ax.set_title(f"{secao.rotulo}   [{st}]", color=cor_status(st), fontsize=10)
    ax.set_xlabel("distancia ao longo da secao (m)")
    ax.set_ylabel("elevacao (m)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
    _ajustar(ax)


def comparacao(ax, original, proposta):
    ax.clear()
    todos = []
    for s, cor, rot in ((original, "#78909c", "ORIGINAL"),
                        (proposta, "#1565c0", "PROPOSTA")):
        if s is None or not s.valida:
            continue
        todos.append(s.z)
        ax.plot(s.sta, s.z, "-", color=cor, lw=1.5,
                label=(f"{rot}  QC {s.qc.nota:.0f}  [{s.qc.status}]"
                       if s.qc else rot))
        i = s.i_talvegue
        if i is not None:
            ax.plot([s.sta[i]], [s.z[i]], "o", color=cor, ms=8, mec="k")
    if todos:
        ax.set_ylim(*_faixa_y(np.concatenate(todos)))
    ax.set_xlabel("distancia ao longo da secao (m)")
    ax.set_ylabel("elevacao (m)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("original x proposta", fontsize=10)
    _ajustar(ax)
