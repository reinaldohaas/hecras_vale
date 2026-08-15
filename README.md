# Modelagem Hidrodinâmica & Previsão Operacional de Cheias — Vale do Itajaí

Este repositório contém a **Plataforma de Previsão Operacional de Cheias, Replay Histórico e Simulação Hidrodinâmica da Bacia do Rio Itajaí**, integrando modelos em Python e aplicação web com Leaflet e Plotly.

Consulte o documento completo em:
📖 **[`DOCUMENTACAO_TECNICA_COMPLETA.md`](DOCUMENTACAO_TECNICA_COMPLETA.md)**

---

## 🌊 Destaques da Plataforma

1. **Previsão Operacional com Faixas de Incerteza**:
   - Linha do tempo contínua separando o passado observado do futuro previsto com cenários de incerteza ($Q_{\text{low}}, Q_{\text{mean}}, Q_{\text{high}}$).
   - Assimilação de telemetria em tempo real com decaimento exponencial suave.
2. **Replay dos Grandes Eventos Históricos**:
   - Calibração e reprodução das cheias de **1983** (15.34m Blumenau), **2008** (Baixo Vale), **2011** e **2023** (múltiplos pulsos e barragens).
3. **Topologia Fluvial Completa de 10 Rios**:
   - **Rio Itajaí-Açu**, **Rio do Oeste**, **Rio Mirim Doce**, **Rio do Sul**, **Rio Perimbó**, **Rio Trombudo**, **Rio Hercílio / Norte**, **Rio Benedito**, **Rio Itajaí-Mirim** e **Rio Luís Alves**.
   - Conexão geométrica com **snapping exato ($0.0\text{ m}$)** nas confluências.
4. **Controle Dinâmico de Barragens**:
   - Operação de comportas das barragens **Oeste (Taió)**, **Sul (Ituporanga)** e **Norte (José Boiteux)** com estrita independência física a montante.
5. **Ponte de Exportação HEC-RAS 1D/2D**:
   - Geração automática de hidrogramas de contorno para importação em modelos não permanentes `.u01` e HDF/DSS.

---

## 🚀 Como Rodar Localmente

```powershell
# 1. Iniciar a Interface Web Interativa (Porta 8050)
python -m http.server 8050 --directory "app"

# Acesse no navegador:
# http://localhost:8050

# 2. Executar Demonstrações em Python
python itajai_flood_model/examples/forecast_itajai_mirim.py
python itajai_flood_model/examples/replay_evento_historico.py
python itajai_flood_model/examples/auto_calibration_demo.py
```
