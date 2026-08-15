import requests
import os

def download_dem():
    api_key = "cbe30884d238441ae080d36328459286"
    print("Iniciando download do Modelo Digital de Elevação (Copernicus 30m) via OpenTopography...")
    
    # Bounding box para a região do Itajaí-Açu (Blumenau até Itajaí)
    south = -27.0
    north = -26.8
    west = -49.1
    east = -48.6
    
    url = f"https://portal.opentopography.org/API/globaldem?demtype=COP30&south={south}&north={north}&west={west}&east={east}&outputFormat=GTiff&API_Key={api_key}"
    
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        output_file = "dem_blumenau_itajai.tif"
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print(f"Download concluído com sucesso: {output_file}")
    else:
        print(f"Erro no download. Status code: {response.status_code}")
        print(response.text)

if __name__ == "__main__":
    download_dem()
