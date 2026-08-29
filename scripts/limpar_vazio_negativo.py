# -*- coding: utf-8 -*-
"""Tira do terreno a celula envenenada por vazio negativo, e tapa com Copernicus.

    python scripts/limpar_vazio_negativo.py modelo/Terrain/MDT_SIGSC_30m.tif \
        --vrt modelo/Terrain/mirim30_sigsc_1m.vrt --nome mirim30

A REGRA: NO SIG-SC, VALOR NEGATIVO E VAZIO

  Ja se sabia de duas convencoes de vazio: `nodata=None` em 983 das 1019
  folhas, e `0.0` em 36. Ha uma terceira, e ela e pior: folhas que gravam o
  vazio como NUMERO NEGATIVO GRANDE. Em `MDT_SG-22-Z-B-IV-3-SO-D.tif` o vazio
  vale cerca de -644 m, numa faixa de -644,1 a -644,5, e a folha declara
  `nodata=None`.

  Nenhum `-srcnodata` resolve: ele aceita UM valor, e aqui o vazio e uma
  FAIXA. Por isso a limpeza e feita depois da reducao, e nao na leitura.

  O criterio e simples e nao custa nada de verdade neste vale: cota negativa e
  vazio. A foz esta no nivel do mar e o terreno nu nao tem cota abaixo de zero
  em lugar nenhum do dominio -- a unica cota negativa do modelo e -2,79 m, e e
  LEITO no Station-Elevation do HEC-RAS, que nao vem daqui.

POR QUE NAO BASTA APAGAR A CELULA QUE SAIU NEGATIVA

  Ao reduzir 1 m -> 30 m cada celula e a media de 900 pixels. Se 36% deles
  valem -644, a media sai negativa e a celula se denuncia -- foram 67 assim.
  Mas se apenas 5% valerem -644, a media de um terreno a 77 m sai por volta de
  41 m: POSITIVA, plausivel, e errada em 36 metros. Essa nao se denuncia.

  Entao o criterio nao pode ser o valor da MEDIA. E o MINIMO dos 900 pixels:
  se o minimo daquela celula for negativo, a media dela esta contaminada,
  qualquer que seja o valor que tenha saido. Dai a segunda passada com
  `-r min`, que custa outra leitura das folhas e e o preco de nao entregar
  cota plausivel e falsa.

O QUE ENTRA NO LUGAR

  A celula contaminada vira NoData -- nao se inventa cota para ela. Quem tapa
  o buraco e o Copernicus, entrando como SEGUNDO raster do terreno: o
  RasTerrain empilha na ordem recebida e o primeiro tem prioridade, entao o
  Copernicus so aparece onde o SIG-SC nao tem nada. E o mesmo arranjo que
  `vale/terreno.py:preparar()` ja usa, e mantem o erro do modelo de superficie
  (copa de mata sobre o leito) restrito aos vazios.
"""
import os
import sys

import numpy as np
import rasterio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vale.terreno import _gdal, _rodar, piramides, tamanho  # noqa: E402
from vale.config import WKT                                 # noqa: E402

VAZIO = -9999.0
PISO = 0.0       # cota negativa no SIG-SC e vazio (ver cabecalho)
COPERNICUS = "modelo/Terrain/so_mirim_Terreno.Terreno_Copernicus.tif"


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    tif = argv[0]
    vrt = _arg(argv, "--vrt")
    nome = _arg(argv, "--nome", "mirim30")
    cop = _arg(argv, "--copernicus", COPERNICUS)
    if not vrt or not os.path.exists(vrt):
        raise SystemExit("informe --vrt com o mosaico virtual de 1 m")
    pasta = os.path.dirname(tif) or "."

    with rasterio.open(tif) as s:
        b, tr = s.bounds, s.transform
        a = s.read(1)
        perfil = s.profile

    # ---- segunda passada: o MINIMO de cada celula, na MESMA grade
    piso_tif = os.path.join(pasta, f"{nome}_minimo.tif")
    print(f"segunda passada (-r min) sobre {os.path.basename(vrt)}")
    _rodar([_gdal("gdalwarp"),
            "-te", f"{b.left:.0f}", f"{b.bottom:.0f}",
            f"{b.right:.0f}", f"{b.top:.0f}",
            "-tr", f"{tr.a:g}", f"{-tr.e:g}", "-tap",
            "-r", "min", "-dstnodata", str(int(VAZIO)),
            "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-wm", "2048",
            "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
            "-co", "BIGTIFF=YES", "-overwrite", vrt, piso_tif], None, print)

    with rasterio.open(piso_tif) as s:
        m = s.read(1)
    if m.shape != a.shape:
        raise SystemExit(f"grades diferentes: {m.shape} x {a.shape}")

    tinha = a > VAZIO + 1
    ruim = tinha & (m > VAZIO + 1) & (m < PISO)
    neg = tinha & (a < PISO)
    print(f"\ncelulas com dado                 : {int(tinha.sum())}")
    print(f"   com minimo negativo             : {int(ruim.sum())}  "
          "<- media contaminada, viram NoData")
    print(f"      das quais a media ja era negativa: {int((ruim & neg).sum())}"
          "  <- as unicas que se denunciavam")
    print(f"      contaminadas com media POSITIVA : {int((ruim & ~neg).sum())}"
          "  <- as que passariam despercebidas")
    if ruim.any():
        err = a[ruim] - m[ruim]
        print(f"      amplitude interna delas: mediana {np.median(err):.1f} m,"
              f" maxima {err.max():.1f} m")

    a = a.copy()
    a[ruim] = VAZIO
    perfil.update(nodata=VAZIO, compress="deflate", tiled=True, bigtiff="YES")
    with rasterio.open(tif, "w", **perfil) as s:
        s.write(a, 1)
    piramides(tif, log=print)
    v = a > VAZIO + 1
    print(f"\nraster limpo: {tif}  ({tamanho(os.path.getsize(tif))})")
    print(f"   cota {a[v].min():.2f} a {a[v].max():.2f} m   "
          f"celulas negativas restantes: {int((v & (a < PISO)).sum())}")

    # ---- o .hdf: SIG-SC primeiro, Copernicus so para tapar buraco
    entradas = [tif]
    if cop and os.path.exists(cop):
        entradas.append(cop)
        print(f"   tapa-buraco: {os.path.basename(cop)} (entra em 2o lugar)")
    else:
        print(f"   AVISO: sem Copernicus em {cop} -- os buracos ficam abertos")

    from ras_commander import RasTerrain
    prj = os.path.join(pasta, f"{nome}.prj")
    with open(prj, "w", encoding="ascii") as f:
        f.write(WKT)
    destino = os.path.join(pasta, f"{nome}_Terreno.hdf")
    print(f"\nRasTerrain.create_terrain_hdf -> {os.path.basename(destino)}")
    RasTerrain.create_terrain_hdf(
        input_rasters=entradas, output_hdf=destino,
        projection_prj=prj, units="Meters", timeout_seconds=14400)
    try:
        os.remove(piso_tif)
    except OSError:
        pass
    print("\nCONFERENCIA")
    for p in (destino,):
        ok = os.path.exists(p)
        print(f"   {'OK   ' if ok else 'FALTA'} {os.path.basename(p)}"
              + (f"   {tamanho(os.path.getsize(p))}" if ok else ""))
    return destino


if __name__ == "__main__":
    main(sys.argv[1:])
