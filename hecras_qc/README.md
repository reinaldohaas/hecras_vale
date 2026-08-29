# hecras_qc — controle de qualidade geométrico de seções transversais

Detecta seções transversais geometricamente ruins **a partir do terreno** e
propõe uma seção melhor. Não é um visualizador de perfis: é uma ferramenta de
QC e correção assistida.

O caso que motivou o programa: no modelo do Itajaí-Mirim, a seção
`R1 RS 104738.3` tinha `Bank Sta 511,54 / 616,95` num corte que termina em
620 m — três metros de planície à direita do rio. O HEC-RAS só revelou isso
abortando a simulação. Rodando o QC sobre o terreno, essa mesma seção aparece
como **CRÍTICA** antes de qualquer simulação.

## Instalação

```bash
pip install -r requirements.txt
```

Todas as dependências já estão no ambiente miniforge deste repositório.

## Execução

Interface gráfica:

```bash
python -m hecras_qc.main
```

Lote (sem tela, para milhares de seções):

```bash
python -m hecras_qc.main --lote --dem terreno.tif --eixo rio.geojson --secoes xs.geojson --saida qc.csv
```

Se as seções estiverem num modelo HEC-RAS em vez de shapefile, converta antes:

```bash
python -m hecras_qc.ras_geometry modelo.g01 EPSG:31982 Itajai_Mirim
```

Isso grava `modelo_eixo.geojson` e `modelo_secoes.geojson`.

## O algoritmo do talvegue

O ponto delicado. **O menor ponto do DEM não é necessariamente o talvegue** —
um buraco de uma célula na margem vence o canal real se o critério for só
"cota mínima". A decisão usa três informações:

1. **Proeminência** (`scipy.signal.find_peaks` no perfil invertido). É quanto
   a depressão afunda em relação ao entorno: separa vale de cova. O buraco de
   uma célula tem centímetros; o canal tem metros.
2. **Proximidade ao eixo do rio** — a única informação que não vem do DEM.
   Mesmo com terreno ruim, o canal está perto do eixo.
3. **Penalidade de borda**, para não premiar um mínimo que é só a seção
   cortada ao meio.

A suavização (Savitzky-Golay) serve **apenas para detectar**; as cotas
exportadas continuam sendo as amostradas do DEM. Onde há NoData, fica NaN — não
se inventa terreno, e o CSV de saída omite esses pontos em vez de preenchê-los.

O programa mostra sempre os **três** pontos separados — mínimo absoluto,
depressão principal e talvegue escolhido — para você poder discordar.

## Os testes

| | teste | reprova quando |
|---|---|---|
| **A** | talvegue na extremidade | posição relativa fora de 20–80% (ATENÇÃO até 10–90%, CRÍTICA além) |
| **B** | canal não identificado | profundidade relativa (média das extremidades − talvegue) abaixo do limiar |
| **C** | salto de talvegue | cota fora da **tendência** local do trecho |
| **D** | largura anormal | largura fora de uma razão contra a mediana das vizinhas |
| **E** | orientação | ângulo com o eixo do rio longe de 90° |

Todos os limiares são configuráveis na interface e no `qc.Limiares`.

Duas decisões do teste C que valem registro, porque a versão ingênua não
funciona:

- comparar com a **mediana** das vizinhas só serve em rio plano. No Mirim, com
  150 m entre seções e 1,6% de declividade, a queda legítima entre vizinhas já
  é 2,4 m, e o teste acusava quase tudo. Agora ajusta-se uma reta pelas
  vizinhas e mede-se o resíduo;
- só o limiar absoluto ainda marcava **21%** das seções — ruído de amostragem
  sobre DEM de 30 m apresentado como anomalia. O teste é de *ponto fora da
  curva*, então o limite também considera a dispersão local (MAD). Com isso,
  de 223 seções do Mirim: **209 OK, 7 atenção, 1 incerto, 6 críticas**.

## Correção: geométrica, nunca por cota

Quando uma seção é ruim, o programa **não altera as elevações**. Isso apagaria
o problema do relatório e o manteria no modelo. Duas estratégias, nesta ordem:

1. **Estender** para os dois lados (10/20/30/50/100 m, configurável) e
   reamostrar o terreno. Preserva posição e orientação originais, que muitas
   vezes foram escolhidas por quem conhece o rio.
2. **Perpendicular ao eixo** — se estender não resolve, o problema não é
   comprimento, é posição ou ângulo. Acha-se o cruzamento com o eixo, toma-se a
   direção local e gera-se uma seção nova centrada nele.

A proposta vem com o QC dela e uma pontuação de 0 a 100, lado a lado com a
original. Quem aceita é você.

## Segurança

- o arquivo de entrada **nunca** é sobrescrito;
- toda saída vai para arquivo novo (com sufixo `_1`, `_2`… se já existir);
- cada substituição aceita ou recusada fica registrada em
  `<secoes>_qc_alteracoes.log`, com o antes e o depois;
- `Ctrl+Z` desfaz.

## Não mascarar problemas

O objetivo **não** é deixar tudo verde. Quando o terreno não permite concluir,
o resultado é `INCERTO — revisão manual`, que é diferente de OK e diferente de
CRÍTICA. Nenhuma correção artificial é feita para melhorar a pontuação.

## Módulos

```
dem.py             leitura do GeoTIFF e amostragem bilinear com NoData
river_axis.py      eixo do rio: direção local, cruzamento, perpendicular
cross_sections.py  leitura das seções e extração do perfil
talweg.py          detecção do talvegue provável
qc.py              os cinco testes, classificação e pontuação
correction.py      estender / recriar perpendicular
export.py          GeoJSON, Shapefile, CSV e log de alterações
plotting.py        mapa, perfil e comparação
gui.py             interface PySide6
main.py            ponto de entrada (GUI e lote)
ras_geometry.py    ponte: lê .g01 do HEC-RAS e exporta eixo/seções
```

`qc.py` e `correction.py` são funções puras — a interface e o modo lote chamam
exatamente as mesmas, de propósito: se divergissem, não daria para confiar em
nenhuma das duas.

## Nota sobre o ambiente

Neste ambiente (matplotlib 3.10.9 sobre numpy 2.5.1) várias chamadas do
matplotlib derrubam o processo em código nativo — `axvline`, `vlines`,
`annotate` com seta, `tight_layout()` e `legend(loc="best")`. Foram medidas uma
a uma com dado sintético. O `plotting.py` usa só primitivas que sobrevivem, e o
cabeçalho dele traz a lista. Se você mudar de máquina ou de versão, isso deixa
de ser necessário — mas não custa nada manter.
