# Plataforma de Previsão Operacional de Cheias, Inundação & Engenharia Hidráulica (Bacia do Rio Itajaí)

**Repositório:** `C:\Users\haas\github\hecras_vale`  
**Autor / Projeto:** Plataforma Hidrológica, Hidrodinâmica & Engenharia de Cheias do Vale do Itajaí  
**Data da Documentação:** Agosto / 2026  
**Status do Projeto:** ETAPAS 1, 2, 3 e 4 CONCLUÍDAS E RIGOROSAMENTE VALIDADAS (33 Testes Unitários Aprovados - 100% OK)  

---

## 📑 Sumário Executivo do Projeto

Esta plataforma integra modelos em **Python (Rasterio, GDAL, SciPy, NumPy, Shapely)** e **Interfaces Web Interativas (Leaflet, Plotly, Canvas)** para simulação hidrológica-hidrodinâmica contínua, previsão de cheias e mapeamento bidimensional de manchas de inundação no Vale do Itajaí ($15.000\text{ km}^2$, Santa Catarina).

O sistema opera na cadeia física unificada e sincronizada:
$$\text{CHUVA} \longrightarrow \text{VAZÃO } Q(t) \longrightarrow \text{PROPAGAÇÃO} \longrightarrow \text{COTA } Z_{\text{water}}(s, t) \longrightarrow \text{ACROMAGEM HAND} \longrightarrow \text{LÂMINA } h(x,y,t) \longrightarrow \text{CONECTIVIDADE} \longrightarrow \text{MANCHA 2D}$$

---

## 🏛️ Topologia Fluvial e Rede Hidrográfica Integrada (10 Rios Modelados)

```text
               [BARRAGEM OESTE] (Taió - 83 hm³)
                     │
               Rio Itajaí do Oeste ────┐
                     │                 │
               Rio Mirim Doce ─────────┤ (km 25 - Taió)
                                       │
                                       ▼
  [BARRAGEM SUL] (Ituporanga - 93.5 hm³)│
         │                             │
  Rio Itajaí do Sul                    │
         │                             │
  Rio Perimbó ───┤ (km 20 - Ituporanga)│
  Rio Trombudo ──┤ (km 42 - Agronômica)│
                 │                     │
                 ▼                     ▼
          ===================================
          1. CONFLUÊNCIA DE RIO DO SUL (km 0)
          ===================================
                         │
                  Rio Itajaí-Açu
                         │
  [BARRAGEM NORTE] ──────┼──► Rio Hercílio (km 35 - Ibirama - 357 hm³)
  (Boiteux)              │
                         ▼
                  2. IBIRAMA (km 35)
                         │
  Rio Benedito (Timbó) ──┼──► (km 90 - Indaial)
                         │
                         ▼
                  3. INDAIAL (km 90)
                         │
                         ▼
                  4. BLUMENAU CENTRO (km 105)
                         │
  Rio Luís Alves ────────┼──► (km 130 - Baixo Vale)
  Rio Itajaí-Mirim ──────┼──► (km 153 - Canal Retificado 70% / Braço Velho 30%)
                         │
                         ▼
                  5. ITAJAÍ FOZ & MARÉ OCEÂNICA (km 153)
```

---

## 🌊 MÓDULO DE INUNDAÇÃO 2D: HAND Sincronizado + Batimetria + Marés Oceânicas

Arquivo Principal: [`itajai_flood_model/src/inundation/unified_hand_engine.py`](file:///C:/Users/haas/github/hecras_vale/itajai_flood_model/src/inundation/unified_hand_engine.py)

### 1. O Que É e Como Funciona o HAND Topográfico
O **HAND (Height Above Nearest Drainage)** é uma propriedade **topográfica pura e estática** da bacia hidrográfica, calculada a partir do DEM Copernicus 30m e da rede de drenagem:
$$\text{HAND}(x, y) = Z_{\text{DEM}}(x, y) - Z_{\text{drain}}(x, y)$$
Identidade Topográfica Imutável:
$$\boxed{Z_{\text{DEM}}(x, y) \equiv Z_{\text{drain}}(x, y) + \text{HAND}(x, y)}$$

### 2. Batimetria Real e Fundo Abaixo do Nível do Mar
- **Estuário / Baixo Vale (Itajaí / Navegantes / Ilhota)**: Leito profundo abaixo do nível do mar ($Z_{\text{bed}} = -4.50\text{ m}$ a $-8.00\text{ m}$).
- **Blumenau**: $Z_{\text{bed}} = 1.30\text{ m}$ (zero da régua em $4.88\text{ m}$).
- **Rio do Sul**: $Z_{\text{bed}} = 332.00\text{ m}$.
- **Brusque**: $Z_{\text{bed}} = 14.00\text{ m}$.

### 3. Cotas de Margem e Transbordo (Bankfull Stage)
Em regime normal de estiagem ou vazão média, a água corre confinada dentro da calha profunda:
- **Blumenau**: Cota de transbordo $H_{\text{bank}} = 8.00\text{ m}$ ($Z_{\text{bank}} = 12.88\text{ m}$).
- **Rio do Sul**: Cota de transbordo $H_{\text{bank}} = 7.00\text{ m}$ ($Z_{\text{bank}} = 339.00\text{ m}$).
- **Brusque**: Cota de transbordo $H_{\text{bank}} = 5.00\text{ m}$ ($Z_{\text{bank}} = 20.50\text{ m}$).
- **Itajaí**: Cota de transbordo $H_{\text{bank}} = 2.50\text{ m}$ ($Z_{\text{bank}} = 2.50\text{ m}$).

**Eliminação de Falsas Inundações**: Em nível normal ($t \le 6\text{h}$ com $H < H_{\text{bank}}$), a área inundada na planície é **estritamente $0.00\text{ km}^2$**.

### 4. Condição de Contorno Oceânica na Foz (Itajaí)
- **Maré Astronômica Semidiurna**: $\pm 0.85\text{ m}$ (período de 12.42h).
- **Maré Meteorológica (Storm Surge / Ressaca)**: Sobre-elevação de tempestade por vento sul e baixa pressão ($+1.40\text{ m}$ em 1983).
- **Remanso Estuarino**: O nível do mar $Z_{\text{ocean}}(t) = Z_{\text{tide}}(t) + Z_{\text{surge}}(t)$ propaga remanso pelos últimos 35 km do rio.

### 5. Equação Dual Rigorosa da Lâmina de Inundação
$$\eta(x, y, t) = Z_{\text{water}}(x, y, t) - Z_{\text{drain}}(x, y)$$
$$\boxed{h(x, y, t) = Z_{\text{water}}(x, y, t) - Z_{\text{DEM}}(x, y) \equiv \eta(x, y, t) - \text{HAND}(x, y)}$$
Quando ocorre transbordo ($Z_{\text{water}} > Z_{\text{bank}}$), a profundidade é dada por $\text{depth} = \max(0.0, \; h(x, y, t))$.

---

## 📊 Validação dos Eventos Históricos

| Evento | Blumenau ($H_{\text{pico}}$) | Rio do Sul ($H_{\text{pico}}$) | Brusque ($H_{\text{pico}}$) | Ressaca na Foz | Área Inundada Máxima | Volume Retido |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **1983** | $15.34\text{ m}$ | $13.00\text{ m}$ | $8.50\text{ m}$ | $+1.40\text{ m}$ | **$431.75\text{ km}^2$** | $6.259\text{ hm}^3$ |
| **2008** | $11.52\text{ m}$ | $4.20\text{ m}$ | $8.50\text{ m}$ | $+1.20\text{ m}$ | **$406.79\text{ km}^2$** | $5.086\text{ hm}^3$ |
| **2011** | $12.60\text{ m}$ | $11.00\text{ m}$ | $6.20\text{ m}$ | $+0.90\text{ m}$ | **$416.32\text{ km}^2$** | $5.324\text{ hm}^3$ |
| **2023** | $10.76\text{ m}$ | $8.50\text{ m}$ | $6.00\text{ m}$ | $+0.80\text{ m}$ | **$400.12\text{ km}^2$** | $4.777\text{ hm}^3$ |

---

## 🧪 Suíte de Testes Automatizados (33/33 Aprovados - 100% OK)

```text
test_01_mandatory_numerical_example (test_unified_hand_engine) ... ok
test_02_mathematical_identity_grid (test_unified_hand_engine) ... ok
test_03_hand_is_static_flood_is_dynamic (test_unified_hand_engine) ... ok
test_04_stage_inundation_progression (test_unified_hand_engine) ... ok
test_05_configurable_min_flood_depth (test_unified_hand_engine) ... ok
test_06_normal_level_confined_in_channel (test_unified_hand_engine) ... ok
test_07_ocean_tide_boundary (test_unified_hand_engine) ... ok
test_rating_curves (test_rating_curves - 5 testes) ... ok
test_spatial_stage (test_spatial_stage - 5 testes) ... ok
test_hydraulic_inundation_2d (test_hydraulic_inundation_2d - 8 testes) ... ok
test_inundation_basic (test_inundation - 8 testes) ... ok
----------------------------------------------------------------------
Ran 33 tests in 2.147s -> OK
```

---

## 🖥️ Interfaces e Visualizadores Web Disponíveis

Execute o servidor local com `python -m http.server 8050 --directory app`:

1. **Dashboard Operacional Principal**:
   👉 **[http://localhost:8050](http://localhost:8050)**  
   *Animação temporal sincronizada ($t \in [0, 48\text{h}]$) da onda nos rios e manchas HAND 2D.*
2. **Visualizador de Diagnóstico & Debug (8 Camadas)**:
   👉 **[http://localhost:8050/debug_hand_layers.html](http://localhost:8050/debug_hand_layers.html)**  
   *Inspetor de pixel em tempo real para verificação ponto a ponto das equações topográficas e hidráulicas.*
3. **Mapeador 2D de Inundação**:
   👉 **[http://localhost:8050/mancha_inundacao_2d.html](http://localhost:8050/mancha_inundacao_2d.html)**
4. **Catálogo de Curvas-Chave Q-H**:
   👉 **[http://localhost:8050/curvas_chave_itajai.html](http://localhost:8050/curvas_chave_itajai.html)**
5. **Perfis Longitudinais e Seções Transversais**:
   👉 **[http://localhost:8050/perfil_longitudinal_cotas.html](http://localhost:8050/perfil_longitudinal_cotas.html)** e **[http://localhost:8050/secoes_transversais_vales.html](http://localhost:8050/secoes_transversais_vales.html)**
