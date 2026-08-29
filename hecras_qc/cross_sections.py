# -*- coding: utf-8 -*-
"""
As secoes transversais: leitura, amostragem do perfil e metricas geometricas.

Uma Secao guarda a linha (no CRS metrico), o perfil amostrado do DEM e o
resultado da deteccao do talvegue. Ela nao decide se esta boa ou ruim -- isso e
do qc.py -- nem se corrige -- isso e do correction.py.

O perfil e sempre reamostrado do DEM com espacamento uniforme configuravel. Nao
se guarda cota "herdada" de arquivo nenhum: o objetivo do programa e conferir a
geometria CONTRA o terreno, e comparar contra uma cota que veio do mesmo lugar
que se quer auditar nao verifica nada.
"""
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString

# nomes que costumam trazer a River Station em arquivos de secao
COLUNAS_RS = ("rs", "river_sta", "riverstation", "river_station", "station",
              "xs_id", "id", "sta", "estaca")
COLUNAS_RIO = ("river", "rio", "reach", "name", "nome", "trecho")


class Secao:
    """Uma secao transversal com seu perfil extraido do terreno."""

    __slots__ = ("idx", "rio", "reach", "rs", "geom", "sta", "z", "xs", "ys",
                 "sta_eixo", "azimute", "talvegue", "qc", "origem", "atrib")

    def __init__(self, idx, geom, rio="", reach="", rs=None, atrib=None):
        self.idx = int(idx)
        self.geom = geom
        self.rio = rio or ""
        self.reach = reach or ""
        self.rs = rs
        self.atrib = dict(atrib or {})
        self.sta = self.z = self.xs = self.ys = None
        self.sta_eixo = None
        self.azimute = None
        self.talvegue = None
        self.qc = None
        self.origem = "original"

    # ------------------------------------------------------------- metricas
    @property
    def largura(self):
        return float(self.geom.length)

    @property
    def valida(self):
        return self.z is not None and np.isfinite(self.z).sum() >= 3

    @property
    def i_talvegue(self):
        return (self.talvegue or {}).get("i_talvegue")

    @property
    def z_talvegue(self):
        i = self.i_talvegue
        return float(self.z[i]) if i is not None and self.valida else float("nan")

    @property
    def posicao_relativa(self):
        """0 = comeco da secao, 1 = fim. O criterio central do QC."""
        i = self.i_talvegue
        if i is None or not self.valida:
            return float("nan")
        L = float(self.sta[-1] - self.sta[0]) or 1.0
        return float((self.sta[i] - self.sta[0]) / L)

    @property
    def dist_margem_esq(self):
        i = self.i_talvegue
        return float(self.sta[i] - self.sta[0]) if i is not None else float("nan")

    @property
    def dist_margem_dir(self):
        i = self.i_talvegue
        return float(self.sta[-1] - self.sta[i]) if i is not None else float("nan")

    @property
    def profundidade_relativa(self):
        return float((self.talvegue or {}).get("profundidade_relativa", float("nan")))

    @property
    def rotulo(self):
        base = f"RS {self.rs}" if self.rs is not None else f"#{self.idx}"
        return f"{self.rio} {self.reach} {base}".strip()

    # ------------------------------------------------------------ operacoes
    def extrair(self, dem, espacamento=2.0, eixo=None, proeminencia_min=0.5):
        """Amostra o DEM e detecta o talvegue."""
        from . import talweg
        self.sta, self.z, self.xs, self.ys = dem.perfil_linha(
            self.geom, espacamento, crs=None)
        self.sta_eixo = None
        if eixo is not None:
            p, dist, linha = eixo.ponto_no_eixo(self.geom)
            if dist <= max(2.0 * espacamento, 5.0):
                self.sta_eixo = float(self.geom.project(p))
            self.azimute = float(_angulo_entre(self.geom, eixo.direcao(p, linha)))
        self.talvegue = talweg.detectar(self.sta, self.z, self.sta_eixo,
                                        proeminencia_min=proeminencia_min)
        return self

    def perfil_df(self):
        import pandas as pd
        return pd.DataFrame({"River": self.rio, "River Station": self.rs,
                             "Station": self.sta, "Elevation": self.z})

    def copia_com_geometria(self, geom, origem):
        s = Secao(self.idx, geom, self.rio, self.reach, self.rs, self.atrib)
        s.origem = origem
        return s


def _angulo_entre(linha, direcao_eixo):
    """Angulo entre a secao e a direcao do eixo, em graus (0 a 90).

    Uma secao transversal deve estar proxima de 90 graus. Menos que isso e uma
    secao obliqua: ela mede uma largura maior que a real e, num 1D, superestima
    a area de escoamento.
    """
    (x0, y0), (x1, y1) = linha.coords[0], linha.coords[-1]
    vx, vy = x1 - x0, y1 - y0
    n = float(np.hypot(vx, vy)) or 1.0
    tx, ty = direcao_eixo
    cos = abs((vx * tx + vy * ty) / n)
    return float(np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))))


def carregar(caminho, crs_alvo, eixo=None):
    """Le as secoes e devolve uma lista de Secao ordenada rio abaixo."""
    gdf = gpd.read_file(caminho)
    if gdf.crs is None:
        raise ValueError(
            f"{caminho} nao tem CRS definido. Defina o CRS do arquivo antes "
            f"de usar -- adivinhar produziria geometria errada em silencio.")
    if crs_alvo is not None:
        gdf = gdf.to_crs(crs_alvo)

    col_rs = _achar(gdf.columns, COLUNAS_RS)
    col_rio = _achar(gdf.columns, COLUNAS_RIO)
    secoes = []
    for i, linha in enumerate(gdf.itertuples(index=False)):
        g = linha.geometry
        if g is None or g.is_empty:
            continue
        if isinstance(g, MultiLineString):
            g = max(g.geoms, key=lambda p: p.length)
        if not isinstance(g, LineString) or g.length <= 0:
            continue
        atrib = {c: getattr(linha, c, None) for c in gdf.columns
                 if c != "geometry"}
        rs = atrib.get(col_rs) if col_rs else None
        try:
            rs = float(rs) if rs is not None else None
        except (TypeError, ValueError):
            pass
        sec = Secao(i, g, rio=str(atrib.get(col_rio, "") or ""),
                    rs=rs, atrib=atrib)
        # trecho separado quando o arquivo traz a coluna: o agrupamento das
        # vizinhas no qc.avaliar_todas usa (rio, trecho)
        for c in ("reach", "trecho", "Reach"):
            if c in atrib and atrib[c] is not None:
                sec.reach = str(atrib[c])
                break
        secoes.append(sec)

    # ordem rio abaixo. Com RS numerica, ela manda (convencao do HEC-RAS: RS
    # decresce para jusante). Sem ela, a projecao no eixo -- os testes de salto
    # de talvegue e de largura comparam VIZINHAS, e sem ordem nao ha vizinha.
    if secoes and all(isinstance(s.rs, float) and np.isfinite(s.rs)
                      for s in secoes):
        secoes.sort(key=lambda s: -s.rs)
    elif eixo is not None:
        secoes.sort(key=lambda s: eixo.estacao(s.geom)[0])
    return secoes


def _achar(colunas, nomes):
    baixa = {str(c).lower(): c for c in colunas}
    for n in nomes:
        if n in baixa:
            return baixa[n]
    return None
