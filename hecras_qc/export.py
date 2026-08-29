# -*- coding: utf-8 -*-
"""
Exportacao e registro de alteracoes.

Regra de seguranca do projeto: o arquivo original NUNCA e sobrescrito. Toda
saida vai para um arquivo novo, e toda substituicao aceita fica registrada num
log em texto, com o antes e o depois. Sem isso, uma sessao de correcao vira uma
caixa preta: passados alguns dias ninguem sabe quais secoes foram trocadas nem
por que.
"""
import csv
import datetime
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point


def _linha_atrib(s):
    i = s.i_talvegue
    return {
        "idx": s.idx,
        "river": s.rio,
        "reach": s.reach,
        "rs": s.rs,
        "origem": s.origem,
        "largura_m": round(s.largura, 2),
        "talvegue_pct": (round(100 * s.posicao_relativa, 2)
                         if np.isfinite(s.posicao_relativa) else None),
        "talvegue_sta": (round(float(s.sta[i]), 2) if i is not None else None),
        "talvegue_z": (round(s.z_talvegue, 3)
                       if np.isfinite(s.z_talvegue) else None),
        "dist_esq_m": (round(s.dist_margem_esq, 2)
                       if np.isfinite(s.dist_margem_esq) else None),
        "dist_dir_m": (round(s.dist_margem_dir, 2)
                       if np.isfinite(s.dist_margem_dir) else None),
        "prof_rel_m": (round(s.profundidade_relativa, 3)
                       if np.isfinite(s.profundidade_relativa) else None),
        "orientacao": (round(s.azimute, 2) if s.azimute is not None else None),
        "qc_status": s.qc.status if s.qc else "",
        "qc_nota": s.qc.nota if s.qc else None,
        "qc_motivos": s.qc.resumo if s.qc else "",
    }


def secoes_gdf(secoes, crs):
    if not secoes:
        return gpd.GeoDataFrame(geometry=[], crs=crs)
    return gpd.GeoDataFrame([_linha_atrib(s) for s in secoes],
                            geometry=[s.geom for s in secoes], crs=crs)


def talvegues_gdf(secoes, crs):
    """Camada de pontos com a posicao do talvegue de cada secao."""
    pts, atrib = [], []
    for s in secoes:
        i = s.i_talvegue
        if i is None or s.xs is None:
            continue
        pts.append(Point(float(s.xs[i]), float(s.ys[i])))
        atrib.append(_linha_atrib(s))
    return gpd.GeoDataFrame(atrib, geometry=pts, crs=crs)


def exportar_vetor(secoes, caminho, crs, incluir_talvegue=True):
    """Grava GeoJSON ou Shapefile, sem tocar no arquivo de entrada."""
    caminho = _nome_livre(caminho)
    gdf = secoes_gdf(secoes, crs)
    gdf.to_file(caminho)
    saidas = [caminho]
    if incluir_talvegue:
        base, ext = os.path.splitext(caminho)
        p = _nome_livre(f"{base}_talvegue{ext}")
        t = talvegues_gdf(secoes, crs)
        if len(t):
            t.to_file(p)
            saidas.append(p)
    return saidas


def exportar_csv_perfis(secoes, caminho):
    """CSV no formato que o HEC-RAS aceita importar.

    Colunas River / River Station / Station / Elevation, uma linha por ponto,
    secoes na ordem de rio abaixo. Pontos NoData sao OMITIDOS -- escrever cota
    inventada aqui seria contrabandear o problema para dentro do modelo.
    """
    caminho = _nome_livre(caminho)
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["River", "Reach", "River Station", "Station", "Elevation"])
        for s in secoes:
            if not s.valida:
                continue
            for sta, z in zip(s.sta, s.z):
                if np.isfinite(z):
                    w.writerow([s.rio, s.reach, s.rs,
                                f"{sta:.3f}", f"{z:.3f}"])
    return caminho


def exportar_tabela(secoes, caminho):
    caminho = _nome_livre(caminho)
    pd.DataFrame([_linha_atrib(s) for s in secoes]).to_csv(
        caminho, index=False, encoding="utf-8")
    return caminho


def registrar(caminho_log, secao_original, secao_nova, aceito=True):
    """Anota uma substituicao aceita (ou recusada) no log da sessao."""
    agora = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(caminho_log, "a", encoding="utf-8") as f:
        f.write(
            f"[{agora}] {'ACEITA' if aceito else 'RECUSADA'} "
            f"{secao_original.rotulo}\n"
            f"    antes : {secao_original.origem} | largura "
            f"{secao_original.largura:.1f} m | talvegue "
            f"{100*secao_original.posicao_relativa:.1f}% | QC "
            f"{secao_original.qc.nota if secao_original.qc else 0:.0f}\n"
            f"    depois: {secao_nova.origem} | largura "
            f"{secao_nova.largura:.1f} m | talvegue "
            f"{100*secao_nova.posicao_relativa:.1f}% | QC "
            f"{secao_nova.qc.nota if secao_nova.qc else 0:.0f}\n")
    return caminho_log


def _nome_livre(caminho):
    """Nunca sobrescreve: acrescenta _1, _2... se ja existir."""
    if not os.path.exists(caminho):
        return caminho
    base, ext = os.path.splitext(caminho)
    i = 1
    while os.path.exists(f"{base}_{i}{ext}"):
        i += 1
    return f"{base}_{i}{ext}"
