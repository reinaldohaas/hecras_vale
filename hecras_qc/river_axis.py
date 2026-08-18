# -*- coding: utf-8 -*-
"""
Eixo do rio: direcao local, cruzamento com a secao, estacao ao longo do rio.

O eixo e a referencia de tudo que este programa faz de inteligente. E ele que
diz onde o canal DEVERIA estar quando o DEM sozinho e ambiguo, e e a partir
dele que se gera uma secao perpendicular nova.

Cuidado com a direcao local: tomada entre dois vertices consecutivos ela oscila
com o serrilhado da digitalizacao, e as perpendiculares saem tortas. Aqui ela e
tomada numa janela de +-JANELA_DIR metros, que e o suficiente para ignorar o
serrilhado sem perder a curva do meandro.
"""
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge, nearest_points

JANELA_DIR = 50.0        # m, meia-janela da direcao local


class EixoRio:
    """Uma ou mais linhas de eixo, no CRS de trabalho."""

    def __init__(self, geometrias, crs):
        linhas = []
        for g in geometrias:
            if g is None or g.is_empty:
                continue
            if isinstance(g, MultiLineString):
                linhas.extend(list(g.geoms))
            elif isinstance(g, LineString):
                linhas.append(g)
        if not linhas:
            raise ValueError("nenhuma linha de eixo valida no arquivo")
        # linemerge junta os trechos que se tocam; o que sobrar solto continua
        # valendo como linha propria (afluentes, por exemplo)
        unido = linemerge(MultiLineString(linhas)) if len(linhas) > 1 else linhas[0]
        self.linhas = (list(unido.geoms) if isinstance(unido, MultiLineString)
                       else [unido])
        self.crs = crs

    # ------------------------------------------------------------------ IO
    @classmethod
    def ler(cls, caminho, crs_alvo):
        gdf = gpd.read_file(caminho)
        if gdf.crs is None:
            raise ValueError(
                f"{caminho} nao tem CRS definido. Defina o CRS do arquivo "
                f"antes de usar -- adivinhar aqui produziria geometria errada "
                f"em silencio.")
        if crs_alvo is not None:
            gdf = gdf.to_crs(crs_alvo)
        return cls(list(gdf.geometry), gdf.crs)

    # ------------------------------------------------------------ geometria
    def linha_mais_proxima(self, geom):
        return min(self.linhas, key=lambda l: l.distance(geom))

    def ponto_no_eixo(self, geom):
        """Onde a geometria encontra o eixo.

        Se cruzar, o cruzamento. Se nao cruzar (secao curta demais, ou eixo
        deslocado), o ponto do eixo mais proximo -- e quem chama recebe tambem
        a distancia, para saber se pode confiar.
        """
        linha = self.linha_mais_proxima(geom)
        inter = linha.intersection(geom)
        if not inter.is_empty:
            if inter.geom_type == "Point":
                return inter, 0.0, linha
            # varios cruzamentos: o eixo serpenteia dentro da secao. Fica o do
            # meio da secao, que e o que corresponde ao canal principal.
            pts = [g for g in getattr(inter, "geoms", []) if g.geom_type == "Point"]
            if pts:
                meio = geom.interpolate(0.5 * geom.length)
                return min(pts, key=lambda p: p.distance(meio)), 0.0, linha
        a, b = nearest_points(linha, geom)
        return a, float(a.distance(b)), linha

    def direcao(self, ponto, linha=None):
        """Vetor unitario da direcao local do eixo, suavizado."""
        linha = linha if linha is not None else self.linha_mais_proxima(ponto)
        s = linha.project(ponto)
        a = linha.interpolate(max(0.0, s - JANELA_DIR))
        b = linha.interpolate(min(linha.length, s + JANELA_DIR))
        tx, ty = b.x - a.x, b.y - a.y
        n = float(np.hypot(tx, ty)) or 1.0
        return tx / n, ty / n

    def azimute(self, ponto, linha=None):
        """Azimute da direcao local, em graus (0 = norte, horario)."""
        tx, ty = self.direcao(ponto, linha)
        return float(np.degrees(np.arctan2(tx, ty)) % 180.0)

    def perpendicular(self, ponto, meia_esq, meia_dir, linha=None):
        """LineString perpendicular ao eixo, passando por 'ponto'."""
        tx, ty = self.direcao(ponto, linha)
        rx, ry = ty, -tx                       # normal, apontando a direita
        return LineString([(ponto.x - meia_esq * rx, ponto.y - meia_esq * ry),
                           (ponto.x + meia_dir * rx, ponto.y + meia_dir * ry)])

    def estacao(self, geom):
        """Distancia ao longo do eixo -- serve para ordenar as secoes.

        Util quando o arquivo de secoes nao traz River Station: sem uma ordem
        ao longo do rio nao ha como comparar uma secao com as vizinhas, e os
        testes de salto de talvegue e de largura anormal dependem disso.
        """
        linha = self.linha_mais_proxima(geom)
        p = geom.interpolate(0.5 * geom.length) if geom.geom_type != "Point" else geom
        return float(linha.project(p)), linha

    @property
    def comprimento(self):
        return float(sum(l.length for l in self.linhas))
