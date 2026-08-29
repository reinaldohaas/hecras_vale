# Validação da baseline — QC da geometria HEC-RAS do Itajaí-Mirim

Reprodução independente dos números de `doc/QC_so_mirim_handoff.md` a partir
dos arquivos atuais do repositório. **Nenhuma correção foi feita, nenhum script
de correção foi executado, nenhum arquivo de geometria foi alterado.**

Métrica de referência: `scripts/qc_geometria.py`. O contador do
`Validate Geometry` do RAS Mapper **não foi usado** — o handoff documenta que o
mesmo arquivo, byte a byte, produziu 240 e 454.

---

## A. Ambiente utilizado

```
python        3.12.4   C:\Users\haas\github\hecras_vale\.venv\Scripts\python.exe
plataforma    Windows-11-10.0.26300-SP0
numpy 2.5.2 · rasterio 1.5.1 · shapely 2.1.2 · h5py 3.16.0 · pyproj 3.7.2
matplotlib 3.11.1 · ras-commander 0.99.1

osgeo (GDAL)  ausente
gdalbuildvrt  ausente
```

Confirmado: **não há GDAL nesta máquina**. Toda a amostragem do MDT é feita com
`rasterio` + `numpy`, despachando cada ponto para a folha que o contém
(`scripts/mdt_sigsc.py`). Nenhum mosaico foi montado.

---

## B. Arquivos e SHA1

Todos os SHA1 do handoff **conferem**:

| arquivo | sha1 (10) | MB | handoff |
| --- | --- | --- | --- |
| `modelo/so_mirim.g01` | `cac7971762` | 2,77 | confere |
| `modelo/so_mirim.g03` | `aeba7a8403` | 2,77 | confere |
| `modelo/so_mirim.g04` | `c08482779d` | 2,77 | confere — **não usado**, conforme instrução |
| `modelo/so_mirim.g05` | `4115903669` | 2,77 | confere |
| `modelo/so_mirim.g06` | `c0196b254e` | 2,77 | confere |
| `modelo/mirim_mdt/mirim_mdt.g01` | `ea99ffac74` | 2,72 | confere |

Arquivos originais, inalterados nesta etapa:

```
modelo/so_mirim.prj       1cfee3dada
modelo/so_mirim.g01       cac7971762      <- o original
modelo/so_mirim.g01.hdf   32030431a8
modelo/so_mirim.p01       8f8fc57630
modelo/so_mirim.u01       78fe31bfc9
```

A cadeia `g01 → g03 → g05 → g06 → mirim_mdt` está íntegra e reproduzível.

---

## C. MDT utilizado

**SIG-SC, `C:\Users\haas\Downloads\sigsc\MDT_SG-22-*.tif`.** Verificado folha a
folha, nas 1019:

```
CRS          EPSG:31982    -- 1019 de 1019, o mesmo do modelo, sem reprojecao
resolucao    1,0 x 1,0 m   -- 1019 de 1019
dtype        float32       -- 1019 de 1019
nodata       None em 983 folhas · 0.0 em 36 folhas
extensao     573642 6900301 .. 766209 7068825   (193 x 169 km)
```

Domínio das cutlines: `666491 6980160 .. 730731 7024001` — 64,2 × 43,8 km.
**122 folhas** intersectam, listadas em `modelo/sigsc_tiles_so_mirim.txt`.

Cobertura sobre as 1418 cross sections, medida ponto a ponto nas mesmas estacas
do HEC-RAS:

```
fracao de pontos com dado:  mediana 1,000    minima 0,463
secoes com cobertura < 90%:  7
```

O Copernicus **não** foi usado como referência nesta etapa.

---

## D. Tratamento do 0.0

Amostrando 11 folhas do domínio (subamostragem de 1/10 em cada eixo):

```
folha                            nodata    %zeros     %<0       max
MDT_SG-22-Z-B-IV-3-SE-A.tif      None      93,51%    0,00%    443,26
MDT_SG-22-Z-B-IV-4-SO-A.tif      None      93,97%    0,00%     61,96
MDT_SG-22-Z-B-V-3-SO-A.tif       None      94,75%    0,00%     87,35
MDT_SG-22-Z-D-I-1-NE-A.tif       None      93,50%    0,00%    163,74
MDT_SG-22-Z-D-I-2-NE-A.tif       None      94,28%    0,00%    309,76
MDT_SG-22-Z-D-I-2-SE-A.tif       None      94,43%    0,00%    881,94
MDT_SG-22-Z-D-I-3-NE-A.tif       None      93,49%    0,00%    692,50
MDT_SG-22-Z-D-II-1-NE-E.tif      None      95,06%    0,00%     97,42
MDT_SG-22-Z-D-II-1-SE-E.tif      None      95,05%    0,00%    150,05
MDT_SG-22-Z-D-II-2-NO-E.tif      None      95,37%    0,00%    219,79
MDT_SG-22-Z-D-II-4-NO-A.tif      None      95,36%    0,00%     50,19

zeros: mediana 94,43%   min 93,49%   max 95,37%
```

**`0.0` é tratado como vazio** (`NaN`), em todas as folhas, declarem nodata ou
não. Justificativa quantitativa: 94% dos pixels são exatamente `0.00` e
**nenhum pixel é negativo**. Se `0.0` fosse cota válida, o terreno teria um
platô perfeitamente plano ao nível zero cobrindo 94% da área e nada abaixo
dele, o que não é topografia. O custo é perder o nível do mar exato, irrelevante
num rio cujo leito modelado está entre −2,79 m e 211 m.

Implementação: `scripts/mdt_sigsc.py`, `MosaicoSigsc.cota()` — aplica também
`nodata` declarado e descarta `< -1000`.

---

## E. Metodologia

Amostragem do MDT: para cada cross section, a cutline é percorrida por
comprimento de arco normalizado **nas mesmas estacas do HEC-RAS**, de modo que
cada ponto do perfil tenha o seu par no terreno:

```
f_i = (station_i - station_0) / (station_N - station_0)      f in [0,1]
P_i = cutline.interpolate(f_i, normalized=True)
z_MDT_i = MDT(P_i)
```

Eixo do rio: lido do **próprio `.gNN`** (`Reach XY=`, 10 664 pontos), não de
arquivo auxiliar.

Tangente local do rio: **janela adaptativa**, nunca fixa. A meia-janela é
`clip(2 × largura_do_canal, 20, 150)` metros, e a tangente é a corda entre
`s − janela` e `s + janela` sobre o eixo. Motivo: num meandro fechado uma janela
longa mede a corda, não a tangente.

Critérios, com o método matemático de cada um, em `scripts/qc_geometria.py`:

| categoria | método |
| --- | --- |
| Station × cutline | `abs(L_cutline − (station[-1] − station[0]))`, com `L_cutline` = soma dos segmentos da polilinha |
| interseções múltiplas | `cutline.intersection(eixo_do_proprio_reach)`; conta pontos do resultado (`Point` = 1, `MultiPoint` = n) |
| orientação angular | `ang = min(d, 180−d)` com `d = |azimute_cutline − azimute_tangente| mod 180`; ideal 90°; mede-se `|90 − ang|` |
| overlap com vizinha | `cutline_i.intersects(cutline_{i±1})` |
| dog-leg | maior deflexão entre segmentos consecutivos da própria cutline |
| talvegue na extremidade | `(station[argmin z] − station[0]) / (station[-1] − station[0])` |
| spikes fora do canal | `|z − mediana_movel(z, janela=7)| > 3 m`, apenas em estacas fora de `[lb, rb]` |
| overbank × MDT | mediana de `z_HEC − z_MDT` nas estacas fora de `[lb, rb]` |
| degraus de fundo | `|min(z)_i − min(z)_{i+1}|` entre seções consecutivas, ordenadas por RS decrescente |

Limiares: `station_length` 0,05 / 0,50 m · ângulo 30° / 50° de desvio ·
overbank 1,0 / 3,0 m · dog-leg 30° · spike 3 m. WARNING no primeiro, CRITICAL
no segundo.

---

## F. Números reproduzidos

### Status, pelas 1418 seções

```
                n       OK   WARNING   CRITICAL
g01          1418       52       259       1107
g03          1418      574       678        166
g05          1418      574       678        166
g06          1418      574       678        166
mirim_mdt    1418      574       678        166
```

`g01` e `g06` **conferem exatamente** com o handoff (52/259/1107 e
574/678/166).

### Categorias

```
                            g01    g03    g05    g06   mirim_mdt
overbank fora do MDT       1297     65     65     65      65
spikes fora do canal        805    571    571    571     571
angulo com a tangente       253    253    253    253     253
cruza o eixo 2+ vezes        67     67     67     67      67
cruza cutline vizinha        66     66     66     66      66
cobertura do MDT              7      7      7      7       7
```

Todas conferem com o handoff.

### HEC-RAS menos MDT SIG-SC 1 m (p10 / mediana / p90 / máx|·|)

```
        overbank                            canal                      acima do MDT
g01     +0,40  +7,15  +17,74   41,62        -2,81  +1,50  +9,64  56,83      92%
g03     -0,00  -0,00   +0,00   16,63        -2,82  -0,00  +0,00  56,83      22%
g06     -0,00  -0,00   +0,00   16,63        -2,81  -0,00  +0,00  56,83      22%
```

O achado central do handoff está reproduzido: no `g01` o overbank está **+7,15 m
acima** do terreno real na mediana, em 92% das seções, e o canal **+1,50 m
acima** — não é batimetria.

### Station × comprimento da cutline

```
|L_cutline - (station[-1] - station[0])|
   mediana 0,00335   p90 0,00771   MAXIMO 0,0156 m
   acima de 0,05 m: 0      acima de 0,50 m: 0
```

**Categoria vazia.** As três maiores: RS 135238.76 (0,0156 m), RS 40602.9
(0,0154 m), RS 80967.87 (0,0134 m). Confere com o handoff ("máx 0,015 m").

### Degraus de fundo

```
                mediana    p90     max    >1 m   >2 m
g01              0,040    0,38    1,33      10      0
g03/g05/g06      0,100    0,55    2,77      44      5
```

Confere. Os piores pares no `g06`:

```
RS  127448.69 -> 127354.94   137,49 -> 134,72   -2,77 m
RS   94739.50 ->  94601.23    60,01 ->  57,44   -2,57 m
RS   51017.91 ->  50773.21    11,50 ->   9,17   -2,33 m
RS  123528.74 -> 123503.74   124,08 -> 121,95   -2,13 m
RS   94601.23 ->  94462.96    57,44 ->  59,50   +2,06 m
```

---

## G. Diferenças em relação ao handoff

Duas, ambas de redação, nenhuma de resultado:

**G.1 — cobertura mínima do MDT: 46,3%, e não 68%.**

**G.2 — o efeito atribuído ao `g06` já está inteiro no `g03`.** O handoff
apresenta a coluna como "g06"; a medição mostra `g03`, `g05`, `g06` e
`mirim_mdt` com status e categorias **idênticos**.

---

## H. Explicação de cada diferença

**G.1.** O número do handoff (68%) veio de uma amostra de **60 seções**
espaçadas uniformemente ao longo do rio, usada para estimar o custo antes de
rodar as 1418. A validação mediu todas: a mínima real é **0,463**, e há
**7 seções abaixo de 90%**. O valor do handoff não estava errado para a amostra
— estava errado como afirmação sobre o conjunto. Matematicamente: o mínimo de
uma amostra de 60 é um estimador enviesado para cima do mínimo de 1418, e
nenhuma das 7 seções problemáticas caiu na amostra. **O algoritmo não foi
ajustado**; o texto é que precisa dizer "mínima 46%".

**G.2.** Os três critérios que separam `g03` de `g06` não entram em nenhuma
categoria do auditor:

- `g03 → g05`: move a `Bank Sta` de **13** seções. Isso desloca a fronteira
  entre canal e overbank, mas nas 13 o overbank continua colado no MDT (o
  recorte já o fizera), então nenhuma muda de status. Como 13/1418 = 0,9% e
  nenhuma cruza limiar, o total não se move.
- `g05 → g06`: mexe apenas em `XS HTab Starting El`, que é parâmetro de
  **cálculo hidráulico**, não geometria. O auditor não o lê.

Ou seja: a diferença entre `g03` e `g06` é real e necessária para o HEC-RAS
rodar (o handoff documenta 45 s → 8 h), mas é **invisível para o QC
geométrico**. Isso não é falha do auditor: são camadas diferentes do modelo.

---

## I. Seções problemáticas

Lista completa por categoria em **`doc/QC_so_mirim_secoes_problematicas.csv`**.
Tabela por River Station, com as 40 colunas, em `modelo/so_mirim_qc.csv`.

### I.1 Interseções múltiplas — 67 (4,7%)

```
RS 116184.01   5 cruzamentos   angulo 34,4   largura 135,9 m
RS 127991.53   4               angulo  1,0   largura 120,0 m
RS  95517.90   4               angulo  1,8   largura 142,9 m
RS  91690.34   4               angulo  2,4   largura 194,6 m
RS 127941.63   3               angulo  5,7   largura 120,0 m
RS 104878.43   3               angulo 27,3   largura 177,0 m
RS 101698.58   3               angulo 26,3   largura 132,6 m
RS  92348.02   3               angulo 26,0   largura 259,2 m
RS  89956.85   3               angulo 33,1   largura 350,2 m
RS  80967.87   3               angulo 20,7   largura 142,9 m
```

Todos no mesmo River/Reach — não são braço secundário nem ilha. Note a
correlação com a categoria I.3: os ângulos são baixos.

### I.2 Overlap com cutline vizinha — 66 (4,7%)

```
RS 14619.52   1 vizinha    largura 1550,4 m   dist ant 917,1   prox 480,5
RS 15536.60   2            largura 1354,2 m   dist ant 673,4   prox 917,1
RS 20057.47   1            largura 1351,5 m   dist ant 301,7   prox 451,7
RS 19605.79   2            largura 1272,1 m   dist ant 451,7   prox 451,7
RS 16210.04   2            largura 1211,8 m   dist ant 653,8   prox 673,4
RS 14138.98   1            largura 1196,1 m   dist ant 480,5   prox 480,5
RS 19154.12   2            largura 1192,8 m   dist ant 451,7   prox 451,8
RS 18702.28   2            largura 1170,8 m   dist ant 451,8   prox 451,8
RS 20359.14   1            largura 1152,6 m   dist ant 301,7   prox 301,7
RS 16863.82   1            largura 1148,7 m   dist ant 653,8   prox 653,8
```

Concentrados no baixo vale (RS 14–20 km) e todos em seções muito largas —
1148 a 1550 m contra a mediana de 132 m do modelo.

### I.3 Orientação angular — 253 (17,8%)

```
RS 102974.17   angulo  0,0 graus   janela da tangente 119,4 m   canal 59,7 m
RS   1695.94   angulo  0,2         janela 150,0 m               canal 80,0 m
RS   6717.02   angulo  0,2         janela 150,0 m               canal 85,7 m
RS  76244.85   angulo  0,3         janela 150,0 m               canal 76,1 m
RS 103021.32   angulo  0,7         janela 125,7 m               canal 62,8 m
RS 114674.12   angulo  0,7         janela 110,5 m               canal 55,3 m
RS 114636.78   angulo  0,8         janela 104,2 m               canal 52,1 m
RS 127991.53   angulo  1,0         janela  67,9 m               canal 34,0 m
RS  53597.06   angulo  1,0         janela 149,4 m               canal 74,7 m
RS   1221.01   angulo  1,1         janela 150,0 m               canal 83,4 m
```

Distribuição do ângulo nas 1418: p05 20,1° · p25 69,9° · **mediana 82,3°** ·
p75 87,1° · p95 89,6°. A massa está perto de 90°, como deve ser; o problema é a
cauda. **42 seções abaixo de 10°** — praticamente paralelas ao fluxo. Das 253,
**188 cruzam o eixo exatamente uma vez**, ou seja não é artefato do cálculo da
tangente sobre cruzamentos múltiplos.

### I.4 Talvegue próximo da extremidade — 6 (0,4%)

Critério estrito (fração < 5% ou > 95% da largura):

```
RS  75255.62   talvegue em  4,9%   estaca   7,1 de 143,2   cota  33,57
RS  74583.12   talvegue em  4,9%   estaca   7,1 de 143,4   cota  33,00
RS  63522.38   talvegue em  3,0%   estaca  12,3 de 417,4   cota  21,03
RS     75.00   talvegue em  0,8%   estaca   5,1 de 650,4   cota  -2,79
RS 118529.66   talvegue em  0,0%   estaca   0,0 de 158,2   cota 110,00
RS  73981.17   talvegue em  0,0%   estaca   0,0 de 150,2   cota  33,00
```

Com o critério mais frouxo de "fora do terço central" seriam 352 (25%), mas
esse número **não indica erro**: a seção é do vale, e não há razão física para o
rio estar no meio dela. Só as 6 acima merecem exame.

### I.5 Cobertura do MDT abaixo de 90% — 7 (0,5%)

Ver `doc/QC_so_mirim_secoes_problematicas.csv`. Mínima 46,3%.

### I.6 Dog-leg — 0

Todas as 1418 cutlines têm exatamente 2 pontos; a deflexão interna é 0 por
construção. Categoria vazia, e permanecerá vazia enquanto nenhuma cutline
receber vértices intermediários.

### I.7 Station × cutline — 0

Já reportada em F.

---

## J. Análise detalhada da RS 127448.69

Figura: `modelo/figuras/validacao_RS127448.svg`.

### J.1 As três seções

```
papel         RS            largura     lb      rb   talvegue g01   talvegue g06
anterior      127542.44      120,0    30,00   69,47      143,52         137,63
ALVO          127448.69      120,0    45,00   82,11      143,51         137,49
posterior     127354.94      120,0    45,00   82,11      143,50         134,72
posterior+1   127261.18      120,0    45,00   85,26      143,49         134,51
```

Diferença de elevação entre alvo e posterior, em 93,8 m de eixo:

```
g01:  143,51 -> 143,50   =  -0,01 m
g06:  137,49 -> 134,72   =  -2,77 m
```

Localização do talvegue: alvo na estaca 48,0 de 120,0 (**40%** da seção),
posterior na estaca 51,0 de 120,0 (**42%**). **As duas dentro das margens** —
não é caso de talvegue na extremidade.

### J.2 O teste que separa artefato de feição real

Amostrei o MDT **ao longo do eixo do rio**, 300 pontos entre as duas seções
mais 150 m de folga de cada lado, cobertura 100%:

```
cota no eixo:  min 134,83   max 137,73   amplitude 2,90 m
ENTRE as duas secoes:  137,51 -> 136,18
   queda 1,33 m em 93,8 m  =  1,42%
   maior degrau interno:  0,11 m em 3,3 m
```

Interpretação, com os dois lados:

- **Mudança topográfica real: sim, parcialmente.** O terreno cai mesmo, 1,33 m
  em 93,8 m — 1,42%, plausível a 137 m de altitude. E cai **suavemente**: o
  maior degrau interno ao longo do eixo é 0,11 m.
- **Erro de geometria: sim, o resto.** O degrau registrado é 2,77 m, o real ao
  longo do eixo é 1,33 m. **1,44 m — 52% do degrau — não está no terreno entre
  as duas seções.** Vem de as duas cutlines pegarem mínimos que não estão sobre
  a mesma linha longitudinal: no alvo o mínimo (137,49) coincide com o eixo
  (137,51, diferença 0,02 m), mas na posterior o mínimo (134,72) está **1,46 m
  abaixo** do eixo no cruzamento dela (136,18). A seção posterior alcança algo
  mais fundo fora do eixo.

Ou seja: **o degrau é real na ordem de 1,3 m e artificial na ordem de 1,4 m.**

### J.3 O que isto NÃO decide

Nenhuma das três alternativas do handoff foi escolhida. A medição acima é o
insumo para essa escolha, não a escolha. Em particular, ela **não** autoriza
alisar o leito: metade do degrau é terreno de verdade.

---

## Conclusão da etapa

A baseline é reprodutível. Todos os números de status, categoria, distribuição
HEC × MDT, `station_length_error` e degraus de fundo conferem com o handoff. As
duas diferenças encontradas (cobertura mínima 46,3% e o efeito estar já no
`g03`) são de redação do handoff, foram explicadas matematicamente, e **nenhum
algoritmo foi ajustado para fazer número coincidir**.

As três categorias espaciais continuam abertas e intocadas: 253 de orientação,
67 de interseção múltipla, 66 de overlap.

**Arquivos gerados nesta etapa** — apenas relatório e listas, nenhuma geometria:

```
doc/QC_so_mirim_validation.md                 este documento
doc/QC_so_mirim_secoes_problematicas.csv      RS por categoria
modelo/figuras/validacao_RS127448.svg         figura da J.2
modelo/so_mirim_qc.csv, .geojson              regravados pelo qc_geometria.py
```
