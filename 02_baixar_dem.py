"""
02 - Baixa o Modelo Digital de Elevacao (relevo) da Bacia do Itajai INTEIRA
via OpenTopography (Copernicus GLO-30, ~30 m).

Rode com:  python 02_baixar_dem.py

Saida: dem_bacia_itajai.tif  (GeoTIFF, EPSG:4326)

Obs.: use SUA chave de API do OpenTopography (gratuita em
https://opentopography.org). A area da bacia (~24.000 km2) esta dentro do
limite. O arquivo tem ~50-120 MB e pode levar alguns minutos.
"""
import os
import sys
import requests

# Chave de API do OpenTopography (troque pela sua se quiser)
API_KEY = os.environ.get("OPENTOPO_API_KEY", "cbe30884d238441ae080d36328459286")

# Mesma area do 01_baixar_rios.py (S, W, N, E) com pequena folga
SOUTH, WEST, NORTH, EAST = -27.80, -50.30, -26.35, -48.50

OUT = "dem_bacia_itajai.tif"
DEMTYPE = "COP30"   # Copernicus GLO-30 (~30 m)


def main():
    url = ("https://portal.opentopography.org/API/globaldem"
           f"?demtype={DEMTYPE}&south={SOUTH}&north={NORTH}"
           f"&west={WEST}&east={EAST}&outputFormat=GTiff&API_Key={API_KEY}")
    print(f"Baixando DEM {DEMTYPE} da bacia ({SOUTH},{WEST}) -> ({NORTH},{EAST})...")
    print("(pode levar alguns minutos)")
    try:
        r = requests.get(url, stream=True, timeout=600)
    except Exception as e:
        print("Erro de conexao:", e); sys.exit(1)

    if r.status_code != 200:
        print(f"Falha no download (status {r.status_code}):")
        print(r.text[:500])
        print("\nDica: verifique a chave de API e o limite de area.")
        sys.exit(1)

    total = 0
    with open(OUT, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 16):
            if chunk:
                f.write(chunk); total += len(chunk)
    mb = total / (1024 * 1024)
    print(f"OK: {OUT} ({mb:.1f} MB)")

    # resumo da extensao
    try:
        import rasterio
        with rasterio.open(OUT) as d:
            print(f"  CRS={d.crs}  shape={d.shape}  res={d.res}")
            print(f"  bounds={d.bounds}")
    except Exception:
        pass


if __name__ == "__main__":
    main()
