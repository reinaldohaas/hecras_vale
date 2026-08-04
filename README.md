# Modelagem Hidráulica do Vale do Itajaí no HEC-RAS (com 3 Barragens)

Este repositório contém os scripts de automação em Python e a geometria hidráulica multi-trecho para a **Bacia do Rio Itajaí-Açu**, incluindo as 3 barragens de contenção de cheias da Defesa Civil de Santa Catarina.

---

## 🌊 Estrutura do Modelo

A rede hidráulica abrange os principais afluentes e junções:

- **Rio Itajaí do Sul** (com a **Barragem Sul** em Ituporanga - 5 comportas)
- **Rio Itajaí do Oeste** (com a **Barragem Oeste** em Taió - 7 comportas)
- **Rio Itajaí do Norte / Hercílio** (com a **Barragem Norte** em José Boiteux - 2 comportas)
- **Rio Itajaí-Açu** (Trecho principal passando por Rio do Sul, Indaial, Blumenau e Gaspar)
- **Rio Benedito** e **Rio Itajaí-Mirim**
- **Foz no Oceano Atlântico** (Itajaí / Navegantes)

---

## 📁 Arquivos do Projeto

- `create_full_geometry.py`: Script Python que constrói do zero a geometria multi-trecho com as 3 barragens (`Itajai_Bacia_Completa.prj` e `.g01`).
- `plot_basin_and_results.py`: Script didático que gera figuras ilustrativas sobre o comportamento hidráulico e a contenção de cheias.
- `run_hecras.py`: Controlador em Python que conecta ao HEC-RAS via COM Interface (`RAS701.HECRASController`) para rodar simulações em segundo plano.
- `fetch_full_basin_geojson.py`: Utilitário para baixar vetores da bacia via OpenStreetMap.
- `figuras/`: Pasta com os gráficos e diagramas gerados.

---

## 🎨 Figuras Didáticas Geradas

### 1. Rede Hidrográfica e Localização das Barragens
![Rede Hidrográfica](figuras/figura_1_rede_de_rios_e_barragens.png)

### 2. Perfil Longitudinal de Elevação da Bacia
![Perfil Longitudinal](figuras/figura_2_perfil_longitudinal_cotas.png)

### 3. Efeito Prático do Amortecimento pelas Barragens
![Amortecimento de Cheia](figuras/figura_3_hidrogramas_e_amortecimento.png)

---

## 🚀 Como Executar Localmente

### 1. Gerar as Figuras Didáticas
```bash
python plot_basin_and_results.py
```

### 2. Gerar o Projeto HEC-RAS da Bacia Completa
```bash
python create_full_geometry.py
```

### 3. Rodar a Simulação no HEC-RAS via Python
> Use o **Python 3.10** (que tem o `pywin32` instalado). O Python 3.14 padrão não tem `win32com`.
```bash
py -3.10 run_hecras.py
```
Saída esperada (com as 3 barragens laminando os picos; confluência no Itajaí-Açu):
```
=== Simulacao concluida com SUCESSO ===
Trecho                Q max (m3/s)   Nivel max (m)
Trecho_Norte                2000.0          153.46
Trecho_Oeste                1500.0          133.12
Trecho_Sul                  1200.0          142.65
Trecho_Principal            4257.7           57.52   ← laminado (sem barragens seria 4480.9)
```

---

## ✅ Estado do modelo

**O modelo completo roda de ponta a ponta no HEC-RAS 7.0.1**, headless via COM: rede 1D (3 afluentes → confluência → Itajaí-Açu → mar) **com as 3 barragens** (inline structures com vertedouro + comportas), hidrogramas de cheia a montante e profundidade normal a jusante. Resultados em `Itajai_Bacia_Completa.p01.hdf` (erro de balanço de volume ~0,0015%).

- `INCLUDE_DAMS=True` / `DAM_GATES=True` em `create_full_geometry.py` liga as barragens; as comportas usam abertura fixa (`GATE_OPEN`), o que reduz o pico no Itajaí-Açu de 4480,9 para ~4257,7 m³/s (laminação).

### Detalhes técnicos que destravaram o modelo
Todos os formatos foram validados contra os **projetos-exemplo oficiais** do HEC-RAS (repos `neeraip/hecras-example-models` e `gpt-cmdr/ras-commander`):

1. **Séries em colunas fixas de 8 caracteres, 10 valores por linha** — vale para os hidrogramas E para o perfil do vertedouro (`#Inline Weir SE`). Era a causa raiz do *"missing data"* / do diálogo de load. Valores em linha única eram mal-lidos.
2. `Boundary Location` com 6 campos e padding fixo (Rio=16, Reach=16, RS=8).
3. Profundidade normal a jusante = `Friction Slope=<decl>` (um valor).
4. Condição inicial unsteady = `Use Restart= 0` + `Initial Flow Loc=Rio,Reach,RS_montante,Q`.
5. Junção liga por **nome**: `Up River,Reach=` / `Dn River,Reach=` / `Junc L&A=`.
6. Comportas exigem controle no `.u01`: bloco `Gate Name=` + `Gate Openings=` (série de aberturas).
7. Abrir o projeto via COM com caminho **Windows** (`\`); barras `/` corrompem o caminho do DSS.
