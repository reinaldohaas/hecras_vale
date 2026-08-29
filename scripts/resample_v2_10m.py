# -*- coding: utf-8 -*-
"""Reamostra o corredor v2 (1 m, corrigido) para 10 m -- pratico para
cortar secoes novas sem o peso do 1 m nem a costura de resolucao
contra o Copernico.

    python scripts/resample_v2_10m.py
"""
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RAIZ)
ORIGEM = os.path.join('taha_ai_novo', 'Terrain',
                      'taha_ai_corredor_1m_v2.tif')
SAIDA = os.path.join('taha_ai_novo', 'Terrain',
                     'taha_ai_corredor_10m_v2.tif')


def main():
    import rasterio
    from rasterio.enums import Resampling
    with rasterio.open(ORIGEM) as src:
        fator = 10
        h, w = src.height // fator, src.width // fator
        dados = src.read(1, out_shape=(h, w),
                         resampling=Resampling.average)
        t = src.transform * src.transform.scale(
            src.width / w, src.height / h)
        prof = src.profile.copy()
        prof.update(height=h, width=w, transform=t,
                    compress='deflate', predictor=3, tiled=True,
                    blockxsize=512, blockysize=512)
        with rasterio.open(SAIDA, 'w', **prof) as dst:
            dst.write(dados, 1)
    print(f'{SAIDA}: {w}x{h} px, '
          f'{os.path.getsize(SAIDA) / 1e6:.0f} MB')


if __name__ == '__main__':
    main()
