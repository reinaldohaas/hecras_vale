# Plataforma de Previsão Operacional de Cheias & Simulação Hidrodinâmica da Bacia do Rio Itajaí

**Repositório:** `C:\Users\haas\github\hecras_vale`  
**Autor / Projeto:** Plataforma Hidrológica e Hidrodinâmica do Vale do Itajaí / HEC-RAS Bridge  
**Data da Documentação:** Agosto / 2026  

---

## 📑 Sumário Executivo do Projeto

Este projeto consiste em uma plataforma integrada em **Python e Web/Leaflet/Plotly** para **modelagem hidrológica, previsão operacional em tempo real com faixas de incerteza, replay de cheias históricas (1983, 2008, 2011, 2023) e ponte de exportação para modelos hidrodinâmicos HEC-RAS 1D/2D**.

O sistema evoluiu da onda sintética inicial para um motor hidrológico operacional capaz de responder:
> **“Dada a chuva observada nas últimas 24h e a chuva prevista para as próximas horas, qual será a vazão e o tempo de chegada da crista em Rio do Sul, Ibirama, Indaial, Blumenau, Brusque e Itajaí?”**

---

## 🏛️ Topologia Fluvial e Rede Hidrográfica Integrada

A rede de drenagem modelada cobre integralmente as sub-bacias do Vale do Itajaí com **snapping geométrico exato ($0.000\text{ m}$)** nas confluências:

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
  Rio Itajaí-Mirim ──────┼──► (km 153 - Canal Retificado de Brusque)
                         │
                         ▼
                  5. ITAJAÍ FOZ & OCEANO (km 153)
```

### Detalhes das 10 Calhas Fluviais Modeladas:

1. **🌊 Rio Itajaí-Açu (Tronco Principal)**: $153.1\text{ km}$ desde a confluência de Rio do Sul até a Foz no Oceano Atlântico em Itajaí.
2. **💧 Rio Itajaí do Oeste**: $65.0\text{ km}$ ($2.480\text{ km}^2$), controlado pela Barragem Oeste (Taió).
3. **💧 Rio Mirim Doce**: $42.0\text{ km}$ ($640\text{ km}^2$), afluente da margem esquerda do Rio do Oeste.
4. **💧 Rio Itajaí do Sul**: $50.0\text{ km}$ ($1.250\text{ km}^2$), controlado pela Barragem Sul (Ituporanga).
5. **💧 Rio Perimbó**: $38.5\text{ km}$ ($510\text{ km}^2$), drena Petrolândia e deságua no Rio do Sul em Ituporanga.
6. **💧 Rio Trombudo**: $48.5\text{ km}$ ($720\text{ km}^2$), drena Agrolândia, Braço do Trombudo e Trombudo Central, desaguando em Agronômica.
7. **💧 Rio Hercílio / Norte**: $85.0\text{ km}$ ($3.450\text{ km}^2$), controlado pela Barragem Norte (José Boiteux) e deságua em Ibirama.
8. **💧 Rio Benedito**: $35.0\text{ km}$ ($1.540\text{ km}^2$), drena Timbó e deságua no Açú em Indaial.
9. **💧 Rio Itajaí-Mirim**: $110.0\text{ km}$ ($1.680\text{ km}^2$), drena Vidal Ramos, Botuverá, Brusque e escoa via Canal Retificado em Itajaí.
10. **💧 Rio Luís Alves**: $25.0\text{ km}$ ($580\text{ km}^2$), drena Luís Alves e deságua no Baixo Vale.

---

## 🔬 Fundamentação Matemática dos Modelos

### 1. Precipitação Efetiva (SCS-CN Dinâmico)
A chuva excedente $P_e(t)$ é calculada com abstração inicial $I_a = 0.2 S$:
$$S = \frac{25400}{CN} - 254$$
$$P_e(t) = \frac{(P(t) - 0.2 S)^2}{P(t) + 0.8 S} \quad \text{para } P(t) > 0.2 S$$

### 2. Condição de Umidade Antecedente (AMC I, II, III)
Baseada no acumulado dos últimos 5 dias ($P_5$):
- **AMC I ($P_5 < 35\text{ mm}$)**: $CN_{\text{I}} = \frac{CN_{\text{II}}}{2.281 - 0.0128 \cdot CN_{\text{II}}}$
- **AMC II ($35 \le P_5 \le 53\text{ mm}$)**: $CN_{\text{II}}$ (Condição média)
- **AMC III ($P_5 > 53\text{ mm}$)**: $CN_{\text{III}} = \frac{CN_{\text{II}}}{0.427 + 0.00573 \cdot CN_{\text{II}}}$

### 3. Hidrograma Unitário Curvilíneo do SCS
Pico unitário:
$$Q_p = \frac{2.08 \cdot A}{T_p}, \quad T_p = 0.5 \Delta t + 0.6 T_c$$
Convolução discreta:
$$Q(t) = Q_{\text{base}} + \sum_{k=0}^{t} P_e(k) \cdot U(t - k)$$

### 4. Propagação Fluvial Muskingum com Sub-stepping Estável
Para evitar coeficientes negativos ($C_0 < 0$) quando $2KX > \Delta t$, o trecho é subdividido em $n_{\text{sub}} = \lceil \frac{2KX}{\Delta t} \rceil$:
$$K_{\text{sub}} = \frac{K}{n_{\text{sub}}}$$
$$Q_{t+1} = C_0 I_{t+1} + C_1 I_t + C_2 Q_t$$
$$\sum C_i = 1.0, \quad C_0, C_1, C_2 \ge 0$$

### 5. Hidráulica de Operação de Barragens
Para as comportas de fundo com vazão ecológica/base $Q_{\text{base}}$:
$$Q_{\text{efluente}}(t) = Q_{\text{base}} + \left(\frac{N_{\text{abertas}}}{N_{\text{total}}}\right) \cdot \max(0, Q_{\text{afluente}}(t) - Q_{\text{base}})$$
Armazenamento no reservatório:
$$V(t+1) = V(t) + \max(0, Q_{\text{afluente}} - Q_{\text{efluente}}) \cdot \Delta t$$
Ao atingir $V_{\text{máx}} = \text{Capacidade (hm}^3)$, o excedente verte livremente.

---

## 📦 Estrutura de Arquivos do Projeto

```text
hecras_vale/
├── app/
│   ├── index.html                   # Interface Web com Modos Simulação, Previsão e Replay
│   ├── itajai_real_dem_model.json   # Perfis DEM 30m, 10 rios, coordenadas e distâncias
│   └── data.json                    # Base completa dos rios da bacia
│
├── itajai_flood_model/
│   ├── src/
│   │   ├── rainfall/
│   │   │   ├── provider.py          # Leitor de chuva CSV + Sintética + Qualidade
│   │   │   ├── spatial.py           # Interpolação espacial (Thiessen / IDW P_bacia(t))
│   │   │   ├── antecedent_moisture.py # P5 -> AMC I/II/III -> CN ajustado
│   │   │   └── forecast.py          # Timeline contínua [AGORA] + Cenários Low/Mean/High
│   │   ├── forecasting/
│   │   │   ├── engine.py            # Motor operacional com 10 rios e propagação em cascata
│   │   │   ├── assimilation.py      # Assimilação telemétrica com decaimento exponencial
│   │   │   └── alerts.py            # Sistema de alertas (NORMAL, ATENÇÃO, ALERTA, EMERGÊNCIA)
│   │   ├── historical_events/
│   │   │   ├── events_loader.py     # Séries horárias de 1983, 2008, 2011 e 2023
│   │   │   └── auto_calibration.py  # Otimização multievento e teste cego
│   │   └── export/
│   │       └── hecras_bridge.py     # Gerador de hidrogramas de contorno para HEC-RAS 1D/2D
│   ├── data/
│   │   └── thresholds.json          # Limiares de alerta e emergência por município
│   └── examples/
│       ├── forecast_itajai_mirim.py # Laboratório operacional no Rio Itajaí-Mirim
│       ├── replay_evento_historico.py # Replay dos 4 grandes eventos históricos
│       └── auto_calibration_demo.py # Calibração multievento
│
├── DOCUMENTACAO_TECNICA_COMPLETA.md # Esta documentação detalhada
└── README.md                        # Guia rápido de execução
```

---

## 🚀 Como Executar Localmente

### 1. Iniciar o Servidor Web da Aplicação
No terminal (PowerShell / Bash):
```powershell
python -m http.server 8050 --directory "app"
```
Acesse no navegador: **`http://localhost:8050`**

### 2. Executar os Scripts de Previsão e Calibração
```powershell
# Previsão operacional no Itajaí-Mirim:
python itajai_flood_model/examples/forecast_itajai_mirim.py

# Replay dos eventos históricos (1983, 2008, 2011, 2023):
python itajai_flood_model/examples/replay_evento_historico.py

# Calibração multievento e validação cega em 2023:
python itajai_flood_model/examples/auto_calibration_demo.py
```

---

## 📊 Validação da Propagação em Cascata

A verificação do motor hidrológico confirmou a cronologia física da onda de cheia:

$$\begin{aligned}
\text{1. Tributários Rápidos (Mirim Doce, Perimbó, Trombudo)} &\longrightarrow \text{Picos em } t = 7\text{h} \text{ a } 9\text{h} \\
\text{2. Confluência de Rio do Sul (Oeste + Sul)} &\longrightarrow \text{Pico em } t = 24\text{h} \quad (5.338\text{ m}^3/\text{s}) \\
\text{3. Ibirama (+ Rio Hercílio)} &\longrightarrow \text{Pico em } t = 30\text{h} \quad (7.809\text{ m}^3/\text{s}) \\
\text{4. Indaial (+ Rio Benedito)} &\longrightarrow \text{Pico em } t = 37\text{h} \quad (7.543\text{ m}^3/\text{s}) \\
\text{5. Blumenau Centro (Médio Vale)} &\longrightarrow \text{Pico em } t = 40\text{h} \quad (7.425\text{ m}^3/\text{s}) \\
\text{6. Itajaí Foz (+ Canal Mirim e Luís Alves)} &\longrightarrow \text{Pico em } t = 47\text{h} \quad (7.529\text{ m}^3/\text{s})
\end{aligned}$$
