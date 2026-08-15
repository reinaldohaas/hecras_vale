# Documentação Técnica e Configuração HEC-RAS: Bacia do Itajaí

**Data de Geração:** 2026-08-07 22:04:55  
**Orquestrador Mestre:** `run_hecras.py`  
**Fonte de Dados Hidrográficos:** Agência Nacional de Águas e Saneamento Básico (ANA BHO 2017)  
**Modelo Digital de Elevação:** Copernicus DEM GLO-30 (~30m)  

---

## 1. Resumo Executivo da Bacia

A modelagem hidrodinâmica 1D/2D para a **Bacia do Itajaí** foi estruturada automaticamente com base em dados vetoriais oficiais online. A calha principal e os afluentes conectados foram isolados topologicamente a partir da exutória principal no oceano.

- **Nome da Bacia / Sistema:** Itajaí
- **Vazão Sintética Simulação HEC-RAS:** 1785.0 m³/s
- **Resolução do Relevo DEM:** 30 metros (Copernicus GLO-30)
- **Formato de Saída do Projeto HEC-RAS:** `.prj`, `.g01`, `.u01`, `.p01`

---

## 2. Estruturas Hidráulicas e Sistema de Contenção de Cheias

| Nome da Estrutura | Rio de Localização | Capacidade de Acumulação | Comportas | Função Hidráulica |
|---|---|---|---|---|
| **Barragem Sul (Ituporanga)** | Rio Itajaí do Sul | 110.0 milhões m³ | 5 comportas | Retém os picos de cheia do Rio Itajaí do Sul antes da junção em Rio do Sul |
| **Barragem Oeste (Taió)** | Rio Itajaí do Oeste | 110.0 milhões m³ | 7 comportas | Retém os picos de cheia do Alto Vale vindos de Taió e Rio do Oeste |
| **Barragem Norte (José Boiteux)** | Rio Hercílio / Itajaí do Norte | 357.0 milhões m³ | 2 comportas | Maior reservatório de contenção do estado de SC; protege Ibirama, Apiúna, Indaial e Blumenau |


---

## 3. Arquivos de Geometria e Projeto Gerados

Os seguintes arquivos do modelo HEC-RAS foram criados em `C:\Users\haas\github\hecras_vale`:

1. **[`Itajaí_Bacia_Real.prj`](file:///C:\Users\haas\github\hecras_vale/Itajaí_Bacia_Real.prj):** Arquivo de projeto mestre do HEC-RAS.
2. **[`Itajaí_Bacia_Real.g01`](file:///C:\Users\haas\github\hecras_vale/Itajaí_Bacia_Real.g01):** Geometria 1D contendo a amostragem de elevação do DEM e estacas de margem.
3. **[`Itajaí_Bacia_Real.u01`](file:///C:\Users\haas\github\hecras_vale/Itajaí_Bacia_Real.u01):** Condições de contorno de hidrograma de cheia não-permanente.
4. **[`Itajaí_Bacia_Real.p01`](file:///C:\Users\haas\github\hecras_vale/Itajaí_Bacia_Real.p01):** Plano de execução e parâmetros do solver hidrodinâmico.

---

## 4. Gráficos e Resultados de Inundação

Os gráficos gerados pela simulação foram salvos no diretório `figuras/`:

- **Rede Hidrográfica e Barragens:** [`figuras/figura_1_rede_de_rios_e_barragens.png`](file:///C:\Users\haas\github\hecras_vale/figuras/figura_1_rede_de_rios_e_barragens.png)
- **Perfil Longitudinal do Rio:** [`figuras/figura_2_perfil_longitudinal.png`](file:///C:\Users\haas\github\hecras_vale/figuras/figura_2_perfil_longitudinal.png)
- **Hidrogramas de Vazão:** [`figuras/figura_3_hidrogramas_cheia.png`](file:///C:\Users\haas\github\hecras_vale/figuras/figura_3_hidrogramas_cheia.png)
