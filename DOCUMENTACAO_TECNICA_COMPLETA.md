# Plataforma de Previsão Operacional de Cheias, Inundação & Cenários de Engenharia (Bacia do Rio Itajaí)

**Repositório:** `C:\Users\haas\github\hecras_vale`  
**Autor / Projeto:** Plataforma Hidrológica, Hidrodinâmica & Engenharia de Cheias do Vale do Itajaí  
**Data da Documentação:** Agosto / 2026  
**Status do Projeto:** ETAPA 1 CONCLUÍDA E VALIDADA  

---

## 📑 Sumário Executivo do Projeto

Este projeto consiste em uma plataforma integrada em **Python e Web/Leaflet/Plotly** para **modelagem hidrológica, previsão operacional em tempo real, mapeamento de cotas e manchas de inundação, análise de alternativas de engenharia (Projeto JICA e JICA+) e avaliação custo-benefício**.

O sistema evolui gradualmente pela cadeia:
$$\text{CHUVA} \longrightarrow \text{VAZÃO } Q(t) \longrightarrow \text{PROPAGAÇÃO} \longrightarrow \text{COTA } H(t) \longrightarrow \text{ÁREA INUNDADA } A(t) \longrightarrow \text{IMPACTOS} \longrightarrow \text{CUSTOS} \longrightarrow \text{COMPARAÇÃO}$$

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

---

## 📈 ETAPA 1 CONCLUÍDA: Módulo de Curva-Chave Q-H (`rating_curve/`)

O módulo `itajai_flood_model/src/rating_curve/` estabelece a ponte bidirecional precisa entre a vazão prevista $Q(t)$ ($m^3/s$) e o nível de régua / cota altimétrica $H(t)$ ($m$).

### 1. Distinção Estrita entre Tipos de Curvas:
- **`CURVA OFICIAL / OBSERVADA (ANA / CEOPS)`**: Ajustada com base em medições de campo linimétricas e histórico hidrométrico.
- **`CURVA ESTIMADA (Manning / Seção Transversal / DEM)`**: Calculada analiticamente através da geometria da calha, rugosidade de Manning e declividade de fundo.

### 2. Estações Oficiais Calibradas:
| Estação / Cidade | Código ANA | Rio | Tipo | Zero da Régua ($Z_0$) | Faixa de Validade |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Blumenau Centro** | `83700000` | Rio Itajaí-Açu | **OFICIAL (CEOPS/ANA)** | $11.20\text{ m}$ | $0.20\text{ m}$ a $18.00\text{ m}$ |
| **Rio do Sul** | `83100000` | Rio Itajaí-Açu / Alto Vale | **OFICIAL (ANA/CEOPS)** | $335.50\text{ m}$ | $0.50\text{ m}$ a $16.00\text{ m}$ |
| **Brusque Centro** | `83800000` | Rio Itajaí-Mirim | **OFICIAL (ANA/Defesa Civil)** | $18.40\text{ m}$ | $0.30\text{ m}$ a $11.00\text{ m}$ |
| **Indaial (Ponte)** | `83500000` | Rio Itajaí-Açu | **OFICIAL (ANA)** | $58.20\text{ m}$ | $0.40\text{ m}$ a $14.00\text{ m}$ |
| **Ibirama** | `83300000_EST` | Rio Itajaí-Açu | **ESTIMADA (Manning/DEM)** | $118.50\text{ m}$ | $0.30\text{ m}$ a $16.00\text{ m}$ |
| **Itajaí Foz** | `83900000_EST` | Rio Itajaí-Açu / Foz | **ESTIMADA (Manning/DEM)** | $0.50\text{ m}$ | $0.20\text{ m}$ a $8.00\text{ m}$ |

### 3. Aderência aos Grandes Recordes Históricos de Blumenau:
$$\begin{aligned}
\text{Cheia Secular de Julho/1983:} \quad & H_{\text{obs}} = 15.34\text{ m} \longleftrightarrow Q_{\text{calc}} = 5.856\text{ m}^3/\text{s} \quad (+0.1\%) \\
\text{Grande Cheia de Setembro/2011:} \quad & H_{\text{obs}} = 12.60\text{ m} \longleftrightarrow Q_{\text{calc}} = 4.486\text{ m}^3/\text{s} \quad (-3.5\%) \\
\text{Complexo de Desastres de 2008:} \quad & H_{\text{obs}} = 11.52\text{ m} \longleftrightarrow Q_{\text{calc}} = 3.970\text{ m}^3/\text{s} \quad (-5.5\%) \\
\text{Cheia de Múltiplos Pulsos de 2023:} \quad & H_{\text{obs}} = 10.76\text{ m} \longleftrightarrow Q_{\text{calc}} = 3.616\text{ m}^3/\text{s} \quad (-8.4\%) \\
\text{Nível de Alerta Blumenau (CEOPS):} \quad & H_{\text{obs}} = 8.00\text{ m} \longleftrightarrow Q_{\text{calc}} = 2.402\text{ m}^3/\text{s} \quad (+0.1\%)
\end{aligned}$$

### 4. Formulação Hidráulica por Subseções Divididas (Divided Channel Method):
Para as seções naturais compostas sem régua, calcula-se:
$$Q_{\text{total}}(H) = Q_{\text{leito}}(H) + Q_{\text{planície}}(H)$$
$$Q_{\text{leito}} = \frac{1}{n_{\text{leito}}} A_{\text{leito}} R_{\text{leito}}^{2/3} S_0^{1/2}, \quad Q_{\text{planície}} = \frac{1}{n_{\text{planície}}} A_{\text{planície}} R_{\text{planície}}^{2/3} S_0^{1/2}$$
Garante estrita monotonicidade e reversibilidade sem artefatos numéricos de estrangulamento na cota de transbordo.

---

## 🧪 Suíte de Testes Automatizados da Etapa 1
Arquivo: `itajai_flood_model/tests/test_rating_curves.py`
- Validação de reversibilidade $Q \rightarrow H \rightarrow Q$.
- Validação contra os picos históricos.
- Validação de monotonicidade de Manning.
- Validação de integridade do catálogo `RatingCurveManager`.
- **Status:** `5/5 OK (100% aprovado)`.

---

## 🛑 Status do Roteiro de Desenvolvimento:

- [x] **ETAPA 1:** Módulo de Curva-Chave Q-H Oficial & Estimada (`rating_curve/`)
- [ ] **ETAPA 2:** Transformação Q-H em Cota Espacial ao longo da calha
- [ ] **ETAPA 3:** Geração de Mancha Simplificada de Inundação no DEM
- [ ] **ETAPA 4:** Validação de Enchente Histórica de Blumenau
- [ ] **ETAPA 5 a 15:** Banco das 103 Enchentes, Gerador Sintético e Projeto JICA
- [ ] **ETAPA 16 a 20:** Camada JICA+ e Análise Custo-Benefício Probabilística Final
