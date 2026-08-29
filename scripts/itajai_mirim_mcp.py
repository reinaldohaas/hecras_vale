"""Servidor MCP para construir a geometria 1D preliminar do Itajai-Mirim.

Entrada: MDT GeoTIFF e eixo do rio (GeoJSON/GPKG/SHP).
Saida: GeoPackage com River/Reach e XS Cut Lines, mais CSVs com os perfis.

O programa NAO inventa batimetria, margens ou cotas. Os perfis representam
somente a superficie observada no MDT e precisam de revisao hidraulica.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Annotated

import geopandas as gpd
import numpy as np
import rasterio
from mcp.server.fastmcp import FastMCP
from pyproj import CRS
from rasterio.warp import transform_geom
from shapely.geometry import LineString, MultiLineString, Point, mapping, shape
from shapely.ops import linemerge


mcp = FastMCP("itajai-mirim-hecras")


def _single_line(gdf: gpd.GeoDataFrame) -> LineString:
    merged = linemerge([g for g in gdf.geometry if g is not None and not g.is_empty])
    if isinstance(merged, MultiLineString):
        raise ValueError("O eixo tem trechos desconectados; corrija-o antes de gerar as secoes.")
    if not isinstance(merged, LineString):
        raise ValueError("A camada de eixo precisa conter uma LineString continua.")
    return merged


def _local_tangent(line: LineString, d: float, delta: float) -> tuple[float, float]:
    p0 = line.interpolate(max(0.0, d - delta))
    p1 = line.interpolate(min(line.length, d + delta))
    dx, dy = p1.x - p0.x, p1.y - p0.y
    norm = float(np.hypot(dx, dy))
    if norm == 0:
        raise ValueError(f"Nao foi possivel calcular a direcao do eixo em {d:.2f} m.")
    return dx / norm, dy / norm


def _sample_profile(
    src: rasterio.io.DatasetReader,
    xs: LineString,
    xs_crs: CRS,
    interval: float,
) -> list[tuple[float, float]]:
    distances = np.arange(0.0, xs.length, interval).tolist()
    if not distances or distances[-1] < xs.length:
        distances.append(xs.length)
    xy = [(xs.interpolate(d).x, xs.interpolate(d).y) for d in distances]

    # Rasterio amostra no CRS do raster.
    if xs_crs != CRS.from_user_input(src.crs):
        geom = transform_geom(xs_crs.to_string(), src.crs, mapping(xs))
        rline = shape(geom)
        ratio = rline.length / xs.length
        rxy = [(rline.interpolate(d * ratio).x, rline.interpolate(d * ratio).y) for d in distances]
    else:
        rxy = xy

    nodata = src.nodata
    values = [float(v[0]) for v in src.sample(rxy, indexes=1, masked=False)]
    profile: list[tuple[float, float]] = []
    for station, z in zip(distances, values):
        invalid = not np.isfinite(z) or (nodata is not None and np.isclose(z, nodata))
        if not invalid:
            profile.append((float(station), z))
    return profile


@mcp.tool()
def criar_geometria_itajai_mirim(
    dem_tif: Annotated[str, "Caminho absoluto do MDT GeoTIFF"],
    eixo_vetorial: Annotated[str, "Caminho absoluto do eixo do Itajai-Mirim"],
    pasta_saida: Annotated[str, "Pasta nova ou existente para os resultados"],
    espacamento_secoes_m: Annotated[float, "Distancia entre secoes"] = 250.0,
    meia_largura_secao_m: Annotated[float, "Extensao de cada lado do eixo"] = 750.0,
    intervalo_amostragem_m: Annotated[float, "Passo de amostragem do MDT"] = 5.0,
    eixo_esta_montante_para_jusante: Annotated[
        bool, "True quando o primeiro vertice do eixo esta a montante"
    ] = True,
) -> str:
    """Gera geometria 1D preliminar importavel no RAS Mapper.

    Cria um GeoPackage com as camadas River e XS_Cut_Lines e CSVs separados
    com station/elevation. A estacao fluvial cresce para montante. Nenhuma
    batimetria, bank station ou levee e inferida automaticamente.
    """
    if espacamento_secoes_m <= 0 or meia_largura_secao_m <= 0 or intervalo_amostragem_m <= 0:
        raise ValueError("Espacamentos e larguras devem ser positivos.")

    dem_path, eixo_path = Path(dem_tif), Path(eixo_vetorial)
    out = Path(pasta_saida)
    if not dem_path.is_file() or not eixo_path.exists():
        raise FileNotFoundError("Confira os caminhos do MDT e do eixo vetorial.")
    out.mkdir(parents=True, exist_ok=True)

    eixo_gdf = gpd.read_file(eixo_path)
    if eixo_gdf.crs is None:
        raise ValueError("O eixo nao possui CRS definido.")
    crs = CRS.from_user_input(eixo_gdf.crs)
    if not crs.is_projected:
        raise ValueError("Reprojete o eixo para um CRS projetado em metros (ex.: SIRGAS 2000 / UTM 22S).")
    unit = crs.axis_info[0].unit_name.lower() if crs.axis_info else ""
    if "metre" not in unit and "meter" not in unit:
        raise ValueError("O CRS do eixo precisa usar metros.")

    river = _single_line(eixo_gdf)
    if not eixo_esta_montante_para_jusante:
        river = LineString(list(river.coords)[::-1])

    # Evita secoes exatamente nas extremidades, onde a tangente e instavel.
    chainages = np.arange(espacamento_secoes_m, river.length, espacamento_secoes_m)
    if len(chainages) == 0:
        raise ValueError("O eixo e menor que o espacamento escolhido.")

    xs_rows: list[dict] = []
    profiles: list[dict] = []
    rejected: list[dict] = []
    tangent_delta = max(5.0, min(50.0, espacamento_secoes_m / 5.0))

    with rasterio.open(dem_path) as src:
        if src.crs is None:
            raise ValueError("O MDT nao possui CRS definido.")
        for idx, d in enumerate(chainages, start=1):
            center = river.interpolate(float(d))
            tx, ty = _local_tangent(river, float(d), tangent_delta)
            # Linha da esquerda para a direita, olhando para jusante.
            left = Point(center.x - ty * meia_largura_secao_m, center.y + tx * meia_largura_secao_m)
            right = Point(center.x + ty * meia_largura_secao_m, center.y - tx * meia_largura_secao_m)
            xs = LineString([left, right])
            profile = _sample_profile(src, xs, crs, intervalo_amostragem_m)
            expected = int(np.floor(xs.length / intervalo_amostragem_m)) + 1
            coverage = len(profile) / max(expected, 1)
            river_station = river.length - float(d)
            xs_id = f"IM_{idx:04d}"

            if coverage < 0.90:
                rejected.append({"xs_id": xs_id, "river_station_m": river_station,
                                 "coverage": coverage, "reason": "menos de 90% do perfil tem MDT valido"})
                continue

            xs_rows.append({
                "xs_id": xs_id,
                "river": "Itajai Mirim",
                "reach": "Principal",
                "river_sta": round(river_station, 3),
                "chainage": round(float(d), 3),
                "coverage": round(coverage, 4),
                "geometry": xs,
            })
            for station, elevation in profile:
                profiles.append({"xs_id": xs_id, "river_station_m": river_station,
                                 "station_m": station, "elevation_m": elevation})

    gpkg = out / "itajai_mirim_hecras.gpkg"
    if not xs_rows:
        raise ValueError(
            "Nenhuma secao atingiu 90% de cobertura valida no MDT; "
            "verifique extensao, CRS, NoData e largura das secoes."
        )
    gpd.GeoDataFrame([{"river": "Itajai Mirim", "reach": "Principal", "geometry": river}],
                     crs=crs).to_file(gpkg, layer="River", driver="GPKG")
    gpd.GeoDataFrame(xs_rows, crs=crs).to_file(gpkg, layer="XS_Cut_Lines", driver="GPKG")

    profiles_csv = out / "xs_station_elevation.csv"
    with profiles_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["xs_id", "river_station_m", "station_m", "elevation_m"])
        writer.writeheader()
        writer.writerows(profiles)

    rejected_csv = out / "xs_rejeitadas.csv"
    with rejected_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["xs_id", "river_station_m", "coverage", "reason"])
        writer.writeheader()
        writer.writerows(rejected)

    report = {
        "status": "ok",
        "gpkg": str(gpkg.resolve()),
        "profiles_csv": str(profiles_csv.resolve()),
        "rejected_csv": str(rejected_csv.resolve()),
        "accepted_cross_sections": len(xs_rows),
        "rejected_cross_sections": len(rejected),
        "crs": crs.to_string(),
        "warning": (
            "Geometria preliminar baseada apenas no relevo. Nao contem batimetria, "
            "bank stations, levees, ineffective flow areas, pontes nem Manning n. "
            "Revise no RAS Mapper antes de calcular."
        ),
    }
    (out / "relatorio_geometria.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return json.dumps(report, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    mcp.run(transport="stdio")
