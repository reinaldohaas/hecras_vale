# -*- coding: utf-8 -*-
"""Cura o 0,0-como-cota do SIG-SC: vira nodata e e preenchido por vizinhos.

    python scripts/curar_zero_do_sigsc.py taha_ai_novo/Terrain/taha_ai_fundo_5m.tif
    python scripts/curar_zero_do_sigsc.py taha_ai_novo/Terrain/taha_ai_corredor_1m.tif

Sai um `<nome>_curado.tif` ao lado. A ENTRADA NAO E TOCADA.

O DEFEITO (visto pelo usuario no RAS Mapper, 26/08/2026)

  No SIG-SC a lamina d'agua e o colarinho de borda das folhas valem 0,0 --
  e 0,0 NAO e nodata declarado. O mosaico do fundo (e trechos do corredor)
  engoliu as bordas como cota valida: listras REBAIXADAS a zero, perfeitamente
  alinhadas a grade das folhas, cruzando o vale. No produto do RasTerrain
  elas aparecem como 'buracos entre os blocos de relevo' -- com dado, mas
  dado errado, o que e pior que buraco: nenhuma caca a nodata acha.

O QUE SE FAZ

  Por blocos (com halo de 64 px): todo 0,0 EXATO vira nodata e o
  `fillnodata` da rasterio interpola dos vizinhos validos. Na foz o 0,0
  do espelho do mar tambem e trocado por interpolacao das margens (~0-1 m);
  a hidraulica do estuario nunca veio do MDT -- vem da batimetria do 1D --
  e uma lamina de 0,5 m interpolada engana menos que uma cota de borda de
  folha atravessando o vale a zero.
"""
import os
import sys

import numpy as np


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import rasterio
    from rasterio.windows import Window
    from rasterio.fill import fillnodata

    entrada = argv[0]
    saida = os.path.splitext(entrada)[0] + "_curado.tif"
    src = rasterio.open(entrada)
    perfil = src.profile.copy()
    perfil.update(bigtiff="YES", compress="deflate", tiled=True)
    nod = src.nodata if src.nodata is not None else -9999.0
    perfil.update(nodata=nod)

    BLOCO, HALO = 2048, 64
    n_zero_tot = 0
    with rasterio.open(saida, "w", **perfil) as dst:
        for r0 in range(0, src.height, BLOCO):
            for c0 in range(0, src.width, BLOCO):
                h = min(BLOCO, src.height - r0)
                w = min(BLOCO, src.width - c0)
                rh0 = max(r0 - HALO, 0)
                ch0 = max(c0 - HALO, 0)
                rh1 = min(r0 + h + HALO, src.height)
                ch1 = min(c0 + w + HALO, src.width)
                jan = Window(ch0, rh0, ch1 - ch0, rh1 - rh0)
                z = src.read(1, window=jan)
                zeros = (z == 0.0)
                n_zero = int(zeros.sum())
                if n_zero:
                    n_zero_tot += int(zeros[(r0 - rh0):(r0 - rh0 + h),
                                            (c0 - ch0):(c0 - ch0 + w)].sum())
                    mask = (~zeros) & (z != nod)
                    z = np.where(zeros, nod, z)
                    z = fillnodata(z, mask=(z != nod).astype(np.uint8),
                                   max_search_distance=HALO,
                                   smoothing_iterations=0)
                rec = z[(r0 - rh0):(r0 - rh0 + h), (c0 - ch0):(c0 - ch0 + w)]
                dst.write(rec, 1, window=Window(c0, r0, w, h))
    print(f"{entrada}")
    print(f"   pixels 0,0 curados: {n_zero_tot}")
    print(f"   -> {saida}")


if __name__ == "__main__":
    main(sys.argv[1:])
