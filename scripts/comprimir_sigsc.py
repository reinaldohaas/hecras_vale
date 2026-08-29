# -*- coding: utf-8 -*-
"""Recomprime as folhas do SIG-SC sem perda (DEFLATE + preditor 3).

    python scripts/comprimir_sigsc.py            todas as nao comprimidas
    python scripts/comprimir_sigsc.py --limite 50   so as 50 primeiras

As folhas vem do portal sem compressao (float32 cru, ~126 MB); com
DEFLATE + preditor de gradiente caem para ~45 MB (35%), byte a byte
identicas na leitura. Escreve <folha>.tmp ao lado e so troca (replace)
depois de conferir uma amostra 500x500 — queda de energia no meio nao
corrompe nada. Ja-comprimidas sao puladas: rodar de novo e barato.
"""
import glob
import os
import sys

import numpy as np

PASTA = r'C:\Users\haas\Downloads\sigsc'


def main(argv):
    import rasterio
    limite = None
    if '--limite' in argv:
        limite = int(argv[argv.index('--limite') + 1])
    feitas = puladas = 0
    ganho = 0.0
    folhas = sorted(glob.glob(os.path.join(PASTA, 'MDT_*.tif')))
    for p in folhas:
        with rasterio.open(p) as s:
            if s.profile.get('compress') is not None:
                puladas += 1
                continue
            prof = s.profile
            prof.update(compress='deflate', predictor=3, tiled=True,
                        blockxsize=512, blockysize=512)
            dados = s.read()
        tmp = p + '.tmp'
        with rasterio.open(tmp, 'w', **prof) as o:
            o.write(dados)
        with rasterio.open(tmp) as s2:
            amostra = s2.read(1, out_shape=(500, 500))
        if not np.array_equal(amostra[:100],
                              dados[0][::dados.shape[1] // 500 or 1,
                                       ::dados.shape[2] // 500 or 1
                                       ][:100]):
            os.remove(tmp)
            print(f'  DIVERGIU (mantida original): {os.path.basename(p)}')
            continue
        antes = os.path.getsize(p)
        os.replace(tmp, p)
        ganho += (antes - os.path.getsize(p)) / 1e9
        feitas += 1
        if feitas % 50 == 0:
            print(f'  {feitas} comprimidas, {ganho:.1f} GB ganhos...',
                  flush=True)
        if limite and feitas >= limite:
            break
    print(f'{feitas} comprimidas, {puladas} ja estavam, '
          f'{ganho:.1f} GB liberados')


if __name__ == '__main__':
    main(sys.argv[1:])
