# RELATÓRIO TÉCNICO: IMPLEMENTAÇÃO DO MODELO HAND (HEIGHT ABOVE NEAREST DRAINAGE)
**Projeto**: Modelagem Hidrodinâmica do Vale do Itajaí (SC)  
**Objetivo do Documento**: Fornecer diagnóstico técnico completo e transparente do código desenvolvido para revisão por especialistas ou outra LLM.

---

## 1. O Que Foi Proposto
Desenvolver o modelo **HAND (Height Above Nearest Drainage)** para mapear a mancha 2D de inundação da Bacia do Rio Itajaí ($15.000\text{ km}^2$), onde a lâmina de inundação em cada célula $(x,y)$ no tempo $t$ é dada por:
$$\text{depth}(x,y,t) = \max\Big(0.0, \; H_{\text{rio}}(x_{\text{drenagem}}, t) - \text{HAND}(x,y)\Big)$$
sendo $\text{HAND}(x,y) = Z_{\text{DEM}}(x,y) - Z_{\text{DEM}}(x_{\text{drenagem}}, y_{\text{drenagem}})$.

---

## 2. Dados Disponíveis no Repositório (`C:\Users\haas\github\hecras_vale\`)

1. **Modelos Digitais de Elevação (DEM GeoTIFFs)**:
   - `dem_bacia_itajai.tif`: $5220 \times 6480$ pixels (resolução 30m, cobrindo toda a bacia: lon $[-50.30, -48.50]$, lat $[-27.80, -26.35]$).
   - `dem_blumenau_itajai.tif`: $720 \times 1800$ pixels (resolução 30m, cobrindo o Médio e Baixo Vale: lon $[-49.10, -48.60]$, lat $[-27.00, -26.80]$).
2. **Hidrografia Oficial**:
   - `app/itajai_ana_rios_alta_resolucao.geojson`: 1.127 trechos da base oficial ANA BHO 5k 1:5.000.
   - `app/itajai_real_dem_model.json`: Perfis altimétricos e coordenadas dos 10 rios principais.
3. **Modelos HEC-RAS Existentes**:
   - `Itajai_Bacia_Completa.p01.hdf`: Resultados não-permanentes (49 time steps, 140 seções transversais).
   - `Itajai_Bacia_Completa.rasmap`: Projeto do RAS Mapper.

---

## 3. O Que Foi Feito no Código

Foram criados os seguintes módulos em `itajai_flood_model/src/inundation/`:

### A. `spatiotemporal_hand.py` (`DynamicSpatiotemporalHAND`)
1. **Rasterização dos Rios**:
   - Transforma as coordenadas dos 10 rios em células da grade raster do DEM.
   - Associa a cada célula fluvial um ID de rio e uma fração de distância ao longo do canal ($s \in [0, 1]$).
2. **Cálculo da Distância e Célula Mais Próxima**:
   - Usa `scipy.ndimage.distance_transform_edt` sobre a máscara de rios para encontrar, para cada pixel $(r,c)$ do DEM, os índices da célula fluvial mais próxima `(nearest_r, nearest_c)`.
3. **Cálculo do HAND**:
   - $\text{HAND}(r,c) = \max(0.0, \; Z_{\text{DEM}}(r,c) - Z_{\text{DEM}}(\text{nearest\_r}, \text{nearest\_c}))$.
4. **Aplicação do Nível do Rio $H(x,t)$**:
   - Para cada hora $t$ ($0, 12, 18, 24, 30, 36, 48\text{h}$), extrai a cota do rio no trecho $H_{\text{river}}(\text{rio}, s, t)$.
   - Calcula $\text{depth}(r,c,t) = \max(0.0, \; H_{\text{river}}(\text{rio}, s, t) - \text{HAND}(r,c))$.
5. **Conectividade Hidráulica e Vetorização**:
   - Aplica rotulagem de componentes conectadas (`scipy.ndimage.label`) mantendo apenas feições ligadas ao rio.
   - Vetoriza em polígonos GeoJSON via `rasterio.features.shapes`.

---

## 4. Diagnóstico Técnico: Por Que o Método Atual Não Está Correto?

1. **Distância Euclidiana vs Direção de Fluxo Hidrológica Real (D8)**:
   - O código usou `distance_transform_edt` (distância euclidiana em linha reta) para encontrar o rio mais próximo.
   - **O HAND verdadeiro NÃO usa distância euclidiana geométrica**. O HAND exige a determinação da **direção de fluxo hidrológica (D8 / D-Infinity)** a partir de um DEM hidrologicamente condicionado (*Pit-Filling / Sink Removal* e *Stream Burning*). Cada pixel precisa drenar seguindo o caminho real do escoamento superficial (Flow Path) até encontrar o canal receptor correto, e não o ponto geometricamente mais próximo em linha reta.
2. **Atribuição e Propagação da Linha d'Água $H(x,t)$**:
   - A interpolação de cotas ao longo do talvegue foi simplificada e precisa estar rigidamente acoplada aos nós e seções transversais com perfil de remanso contínuo.
3. **Resolução e Vetorização**:
   - A vetorização direta de rasters em GeoJSON de múltiplos passos temporais gera polígonos complexos que necessitam de um pipeline otimizado (ex: raster dynamic canvas ou GeoJSON simplificado com topologia preservada).

---

## 5. O Que Pedir de Ajuda para o Especialista ou Outra LLM

Para implementar o HAND verdadeiro ou utilizar o RAS Mapper, sugerir a formulação das seguintes perguntas:

1. *"Como gerar o raster HAND hidrologicamente correto a partir de um GeoTIFF DEM 30m e de um GeoJSON de hidrografia em Python utilizando bibliotecas hidrológicas padrão (como `pyflwdir`, `pysheds`, `whitebox`, `richdem` ou `xarray`) com direção de fluxo D8?"*
2. *"Como interpolar a linha d'água dinâmica 1D/2D $H(x,t)$ dos rios para a grade HAND e gerar a mancha de inundação $\text{depth}(x,y,t) = \max(0, H_{\text{drain}}(x,y,t) - \text{HAND}(x,y))$ hora a hora de forma computacionalmente eficiente?"*
3. *"Alternativamente, como automatizar a exportação das manchas de inundação diretamente do HEC-RAS HDF5 (`Itajai_Bacia_Completa.p01.hdf` / `Itajai_Bacia_Completa.rasmap`) para GeoTIFF / GeoJSON?"*
