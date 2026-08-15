# Modelo de Propagação de Cheias da Bacia do Rio Itajaí
## Fase 1: Rio Itajaí-Mirim (Muskingum & SCS Unit Hydrograph)

Modelo hidrológico e hidrodinâmico modular em Python para simulação e previsão de propagação de ondas de cheia na Bacia do Rio Itajaí, iniciando pelo **Rio Itajaí-Mirim**.

---

## 1. Estrutura do Projeto

```text
itajai_flood_model/
│
├── data/
│   └── itajai_mirim/
│       ├── reaches.csv        # Discretização física dos trechos do rio
│       ├── stations.csv       # Estações de monitoramento e coordenadas
│       ├── rainfall.csv       # Séries de chuva de projeto (incremental e acumulada)
│       └── discharge.csv      # Série observada para calibração em Brusque
│
├── src/
│   ├── __init__.py
│   ├── unit_hydrograph.py     # Hidrograma Unitário Sintético / SCS NEH-4 Curvilíneo
│   ├── muskingum.py           # Propagador clássico de Muskingum com estabilidade
│   ├── muskingum_cunge.py     # Propagador físico de Muskingum-Cunge
│   ├── river.py               # Estrutura de dados RiverReach e RiverNetwork
│   ├── routing.py             # Controlador sequencial trecho a trecho (FloodRouter)
│   ├── calibration.py         # Métricas estatísticas (RMSE, NSE, Erro de Pico, Volume)
│   ├── visualization.py       # Gráficos de propagação, atenuação e espaço-tempo
│   └── mapping.py             # Mapa geográfico local dos trechos
│
├── examples/
│   └── itajai_mirim_demo.py   # Script de demonstração executável
│
├── tests/
│   └── test_muskingum.py      # Testes unitários de conservação de massa e estabilidade
│
├── requirements.txt
└── README.md
```

---

## 2. Fundamentos Matemáticos

### 2.1. Hidrograma Unitário do SCS (NRCS NEH-4)
A transformação chuva-vazão na sub-bacia é governada por:
- **Armazenamento Potencial**:
  $$S = \frac{25400}{CN} - 254$$
- **Chuva Efetiva**:
  $$P_e = \frac{(P - 0.2S)^2}{P + 0.8S} \quad (P > 0.2S)$$
- **Tempo de Pico**:
  $$t_p = \frac{\Delta t}{2} + 0.6\, t_c$$
- **Vazão de Pico Unitária**:
  $$q_p = \frac{2.08 \cdot A \cdot 1\text{ mm}}{t_p}$$
- **Convolução Linear Discreta**:
  $$Q(t) = \sum_{k=0}^{t} P_e(k) \cdot U(t - k) + Q_{\text{base}}$$

### 2.2. Método de Muskingum Clássico
A propagação da onda ao longo de cada trecho fluvial utiliza a equação de armazenamento:
$$S = K [X I + (1-X) Q]$$
e a equação da continuidade $\frac{dS}{dt} = I - Q$, resultando na forma discreta:
$$Q_{t+\Delta t} = C_0 I_{t+\Delta t} + C_1 I_t + C_2 Q_t$$
onde:
$$C_0 = \frac{\Delta t - 2KX}{2K(1-X) + \Delta t}, \quad C_1 = \frac{\Delta t + 2KX}{2K(1-X) + \Delta t}, \quad C_2 = \frac{2K(1-X) - \Delta t}{2K(1-X) + \Delta t}$$
com a garantia estrita de conservação de massa:
$$C_0 + C_1 + C_2 = 1.0$$

### 2.3. Formulação Física de Muskingum-Cunge
Para relacionar $K$ e $X$ com a geometria física do leito (declividade $S_0$, largura $B$ e rugosidade $n$ de Manning):
- **Celeridade da Onda**:
  $$c = \frac{5}{3} v = \frac{5}{3} \frac{Q_0}{B \cdot y_0}$$
- **Parâmetros Físicos**:
  $$K = \frac{\Delta x}{c}, \quad X = \frac{1}{2}\left(1 - \frac{Q_0}{B \cdot S_0 \cdot c \cdot \Delta x}\right)$$

---

## 3. Classificação dos Parâmetros

| Parâmetro | Fonte / Tipo | Classificação |
| :--- | :--- | :--- |
| Comprimento dos Trechos ($L$) | Medição georreferenciada da calha da ANA | **REAL** |
| Declividade média ($S_0$) | Modelo Digital de Elevação Copernicus DEM 30m | **REAL** |
| Coordenadas das Estações | Base ANA / HidroWeb e OpenStreetMap | **REAL** |
| Traçado do Canal Retificado | Geometria Oficial `canal_itajai_mirim.geojson` | **REAL** |
| Coeficiente $K$ (Tempo de Trânsito) | Estimativa preliminar baseada em velocidade média de $1.2\text{ m/s}$ | *DEMO / PLACEHOLDER* |
| Fator $X$ (Difusão) | Valor demonstrativo adotado ($0.20 \text{ a } 0.30$) | *DEMO / PLACEHOLDER* |
| Rugosidade $n$ de Manning | Valores tabelados para calhas naturais com meandros ($0.035 \text{ a } 0.045$) | *DEMO / PLACEHOLDER* |

---

## 4. Como Executar

Execute o script de demonstração:
```bash
python examples/itajai_mirim_demo.py
```

Para rodar a bateria de testes unitários:
```bash
python -m unittest tests/test_muskingum.py
```
