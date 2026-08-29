# `g07` — análise de conflito geométrico e dos degraus de fundo

Análise das causas dos problemas restantes. **Nenhuma seção foi alterada,
nenhum perfil foi suavizado, nenhuma cota, Manning ou HTab foi tocado. Não há
`g08`.**

Base: `modelo/so_mirim.g07`, `sha1 c42a6f744c`.
Métrica: MDT SIG-SC 1 m, a mesma da baseline validada.

---

## 0. Rodada de integridade — e o que ela achou

O `g07` foi executado numa cópia isolada
(`%TEMP%\...\scratchpad\run_g07`), sem tocar em `modelo\`.

**O HEC-RAS recusou a geometria.** Não computou nada; escreveu o motivo em
`so_mirim.p01.data_errors.txt`, 1565 linhas, com apenas dois tipos:

```
521  - Right bank station not in station elevation data.
519  - Left bank station not in station elevation data.
```

São 521 das 541 seções recortadas. A causa é minha, e é a **terceira regra do
formato `.g01`** da mesma família das duas já documentadas no handoff:

> **A estaca de margem tem de ser uma das estacas do perfil.** Não basta estar
> dentro da faixa: o valor de `Bank Sta=` precisa coincidir com um valor
> presente no bloco `#Sta/Elev`.

Verificado nas três geometrias:

```
geometria         secoes   lb fora da grade   rb fora   pior erro
so_mirim.g01        1418            0             0      0,0000 m
so_mirim.g06        1418            0             0      0,0000 m
so_mirim.g07        1418          519           521      2,3100 m
```

Exemplo, RS 134746.0: as estacas novas são `0,00 · 1,67 · 3,34 … 128,71`
(passo 1,67 m), e eu gravei `lb = 43,27`, que cai entre `43,46` e a anterior —
erro de 0,19 m. Ao recortar, gerei as estacas com `linspace(0, L, n)` e pus as
margens em `L/2 ± largura_do_canal/2`, valores que quase nunca caem na grade.

**Consequência para esta etapa:** a integridade do `g07` está reprovada. As
análises abaixo continuam válidas — elas medem a geometria espacial, que é
independente desse defeito de formato — mas o arquivo não roda, e o teste
hidráulico fica pendente até isso ser corrigido.

O `.g01.hdf` foi recompilado normalmente (779 KB), então o arquivo **abre**; é
a checagem de dados que barra o compute.

---

## 1. Análise de conflito dos overlaps

106 seções envolvidas, formando **63 pares** de cutlines que se cruzam. Método
em `scripts/analisar_overlaps.py`.

### Critério de classificação, na ordem em que é testado

| classe | regra |
| --- | --- |
| **A** erro inequívoco de geometria | desvio da normal > 30° em alguma das duas, **ou** centro da seção a mais de 25% da largura de distância do eixo |
| **B** consequência inevitável da curvatura | as duas praticamente normais (desvio ≤ 15°), as duas centradas, e **meia-largura > raio de curvatura local** |
| **C** seção excessivamente larga | não é B, e a largura de uma delas passa de 3× a mediana do modelo (131,9 m) |
| **D** problema de espaçamento | não é A/B/C, e a distância longitudinal é menor que a largura média do canal das duas |
| **E** ambíguo | o que sobra |

A ordem importa: A antes de B porque seção torta explica o cruzamento sozinha;
B antes de C porque curvatura é razão física e largura é escolha.

Raio local: `R = ds/dθ`, com θ medido entre duas cordas de meia janela
adaptativa (`clip(2 × largura_do_canal, 20, 150)` m).

### Resultado

```
A  erro inequivoco de geometria           16   25,4%
B  consequencia inevitavel da curvatura   23   36,5%
C  secao excessivamente larga             15   23,8%
D  problema de espacamento                 5    7,9%
E  ambiguo                                 4    6,3%
                                          --
                                          63
```

Assinaturas por classe, que confirmam que a separação é real:

```
      desvio max      meia-largura / R      intersecao no canal
      (mediana)          (mediana)
A       55,4 graus          0,81                 9 de 16
B        4,3 graus          1,84                 7 de 23
C        2,8 graus          0,75                 1 de 15
D        8,6 graus          0,85                 0 de  5
E        5,1 graus          0,79                 1 de  4
```

**A** tem desvio mediano de 55° — são seções tortas, o cruzamento é sintoma.
**B** tem desvio de 4,3° e meia-largura **1,84 vezes o raio local**: são seções
corretas cruzando por geometria pura. As demais têm desvio baixo mas
meia-largura menor que o raio, ou seja a curva não explica.

### Exemplos, os 6 mais largos de cada classe

**A — erro inequívoco (16)**

```
RS i        RS j          dx   desvI  desvJ   largI  largJ      R  lado
14138.98    13658.43     481     0,0   55,8    1196    842   1713  esq
13658.43    13397.10     261    55,8    0,0     842    855   1018  esq
40968.47    40602.90     366     0,0    1,5     727    574    176  dir
2164.73     1695.94      469    55,2    0,9     708    740   1012  esq
7505.47     7111.25      394    32,2    0,0     623    622    168  dir
3952.87     3696.55      256    58,0   44,8     543    490    541  esq
```

Note que em cinco dos seis o par tem **uma** seção torta e a outra perfeita
(0,0°): são as 25 que sobraram da categoria de ângulo, arrastando a vizinha.

**B — curvatura (23)**

```
RS i        RS j          dx   desvI  desvJ   largI  largJ      R  lado
15536.60    14619.52     917     6,5    2,8    1354   1550    157  esq
20057.47    19605.79     452     9,8    0,0    1351   1272    208  esq
16210.04    15536.60     673     8,5    6,5    1212   1354     93  dir
19154.12    18702.28     452     6,6    0,0    1193   1171    160  dir
20359.14    20057.47     302     0,0    9,8    1153   1351    162  dir
16863.82    16210.04     654    12,3    8,5    1149   1212     80  esq
```

Todas no baixo vale, RS 14–21 km, com larguras de 1149 a 1550 m e raios de 80 a
208 m. Meia-largura de ~600 m contra raio de ~150 m: duas seções normais a um
arco desse raio **têm** de se cruzar. Não há orientação que evite.

**C — largura demais (15)**

```
RS i        RS j          dx   desvI  desvJ   largI  largJ      R  lado
14619.52    14138.98     481     2,8    0,0    1550   1196   3805  dir
18702.28    18250.43     452     0,0   11,2    1171   1149   1693  dir
36089.93    35964.93     125     0,0    0,0    1074   1070    565  esq
33105.05    32605.05     500     2,7    0,0    1073   1089    558  dir
35616.26    35392.59     224     0,0    0,0    1060   1053    938  dir
35168.92    34945.25     224     0,0   22,2    1047   1040    155  esq
```

Seções perfeitamente normais (0,0°) e raio grande — 3805 m no primeiro caso.
Aqui a curva não tem culpa nenhuma: o cruzamento vem de a seção ter 1550 m
numa planície onde as vizinhas estão a 481 m.

**D — espaçamento (5)** e **E — ambíguo (4)**

```
D:  103115.63 / 103068.48   dx  47 m   largura ~150 m   R  98
    129022.95 / 128997.95   dx  25 m   largura ~148 m   R  76
    141297.34 / 141272.34   dx  25 m   largura  120 m   R  63
E:   56099.99 /  55953.81   dx 146 m   largura 371/301  R 214
     99417.97 /  99257.80   dx 160 m   largura 211/189  R 162
```

Tabela completa: **`modelo/overlaps_g07.csv`**, 21 colunas por par.

### Leitura

Os overlaps não são um problema, são **três**:

- **B (23 pares, 37%)** é geometria pura e não tem conserto por orientação. Ou
  se reduz a largura nesses trechos, ou se aceita.
- **C (15 pares, 24%)** é largura excessiva com curva folgada — tem conserto,
  mas mexer em largura é decisão sua e você já rejeitou regra fixa de largura.
- **A (16 pares, 25%)** vem das 25 seções que ainda estão tortas. É o resíduo
  do problema que o `g07` resolveu em 90%.
- **D + E (9 pares, 14%)** é cauda.

---

## 2. Análise dos degraus de fundo

81 degraus acima de 1 m, em 1417 pares. Método em
`scripts/analisar_degraus.py`.

### O teste que decide

Compara o degrau **registrado nas seções** com a queda do terreno **medida ao
longo do eixo** entre elas:

```
d_secoes = |min(z_i) - min(z_j)|                    o que o modelo diz
d_eixo   = |MDT(eixo em s_i) - MDT(eixo em s_j)|    o que o terreno diz
explica  = d_eixo / d_secoes
```

Registra-se também o **maior salto interno do MDT ao longo do eixo** entre as
duas: degrau real costuma aparecer ali como salto localizado, não como rampa.

| classe | regra |
| --- | --- |
| **C** erro de geometria | o mínimo de alguma das duas cai **fora das margens**, ou a seção está a mais de 25% da largura do eixo — antes de discutir o degrau, a seção está errada |
| **A** provavelmente real | `explica ≥ 70%` |
| **B** provavelmente artefato | `explica ≤ 30%` |
| **D** inconclusivo | entre 30% e 70% |

### Resultado

```
A  provavelmente real                     28   34,6%
B  provavelmente artefato de amostragem   23   28,4%
C  provavelmente erro de geometria        18   22,2%
D  inconclusivo                           12   14,8%
                                          --
                                          81
```

As assinaturas separam limpo, o que é a melhor evidência de que o critério
mede o que diz medir:

```
      n    degrau mediana   explicado pelo eixo   salto interno do MDT   minimo fora do canal
A    28        1,31 m            99%                    0,21 m                  0 de 28
B    23        1,42 m             2%                    0,03 m                  0 de 23
C    18        1,49 m            14%                    0,03 m                 15 de 18
D    12        1,37 m            52%                    0,11 m                  0 de 12
```

Em **A** o terreno explica 99% do degrau e há salto localizado de 0,21 m no
MDT. Em **B** explica 2% e o terreno é liso (0,03 m). Em **C**, 15 de 18 têm o
mínimo fora do canal declarado. A magnitude do degrau é **igual nas quatro
classes** (1,31 a 1,49 m) — ou seja, o tamanho do degrau **não** diz nada sobre
a natureza dele. Só a comparação com o eixo diz.

### Exemplos

**A — provavelmente real (28)**

```
RS i        RS j         dx   degrau   d_eixo  explica  saltoMDT  minI%  minJ%
51017.91    50773.21    245    -2,37     2,59    109%     0,43     52%    48%
132580.53   132555.53    25    -1,99     1,40     70%     0,14     57%    42%
140922.34   140872.34    50    -1,97     1,75     89%     0,29     58%    58%
86358.99    86081.03    278    -1,92     2,32    121%     0,32     53%    38%
137623.17   137593.17    30    -1,89     1,67     88%     0,22     45%    45%
46230.36    46134.69     96    -1,80     1,64     91%     0,44     46%    56%
35168.92    34945.25    224    -1,78     2,24    126%     0,96     49%    53%
```

Em vários o eixo cai **mais** que as seções (109%, 121%, 126%) — a seção está
até conservadora. Os mínimos ficam perto do meio da seção (38–58%).

**B — provavelmente artefato (23)**

```
RS i        RS j         dx   degrau   d_eixo  explica  saltoMDT  minI%  minJ%
3036.88     2633.52     403    -2,59     0,02      1%     0,02     55%    28%
3440.24     3036.88     403    +2,50     0,02      1%     0,05     47%    55%
12459.53    12044.61    415    -1,63     0,02      1%     0,02     51%    52%
13135.77    12874.44    261    -1,57     0,02      1%     0,04     49%    53%
12874.44    12459.53    415    +1,55     0,04      3%     0,03     53%    51%
81638.39    81496.43    142    -1,52     0,00      0%     0,00     42%    27%
```

O terreno ao longo do eixo é **plano** (0,00 a 0,04 m) e as seções registram
2,59 m. Note os pares alternados — `−2,59` seguido de `+2,50` — que é a
assinatura de serrilha: uma seção pega um ponto fundo, a seguinte não, a
terceira pega de novo.

**C — provavelmente erro de geometria (18)**

```
RS i        RS j         dx   degrau   d_eixo  explica  saltoMDT  minI%  minJ%
1221.01     746.07      475    -2,75     0,03      1%     0,02     85%    26%
1695.94     1221.01     475    +2,66     0,00      0%     0,01     28%    85%
2164.73     1695.94     469    -2,65     0,02      1%     0,02     33%    28%
3696.55     3440.24     256    -2,45     0,01      0%     0,02     10%    47%
94877.77    94739.50    138    -2,00     0,65     33%     0,39     58%    73%
```

Concentrados na foz (RS 0,7–4 km), com o mínimo em 85% ou 10% da seção — fora
do canal. São as mesmas seções que aparecem na classe A dos overlaps: seção mal
posta produz degrau e cruzamento ao mesmo tempo.

**D — inconclusivo (12)**

```
RS i        RS j         dx   degrau   d_eixo  explica  saltoMDT
127448.69   127354.94    94    -2,65     1,34     51%     0,10
123528.74   123503.74    25    -2,13     0,98     46%     0,12
88219.76    88064.53    155    -1,86     0,78     42%     0,16
122093.18   122071.27    22    -1,65     0,88     53%     0,10
```

**A RS 127448.69 cai aqui, com 51%** — consistente com a validação anterior,
que mediu 1,33 m de terreno real contra 2,77 m registrados (48%). Os dois
métodos, feitos separadamente, chegam ao mesmo lugar: metade real, metade
artefato.

Tabela completa: **`modelo/degraus_so_mirim.csv`**, 19 colunas por degrau.

---

## 3. Conclusões

1. **O `g07` não roda**, por um defeito de formato meu — estaca de margem fora
   da grade de estacas em 521 seções. Terceira regra da mesma família das duas
   já documentadas. O teste hidráulico fica pendente.

2. **Os 63 overlaps são três problemas distintos.** Só 25% são erro de
   geometria (resíduo das 25 seções ainda tortas). **37% são consequência
   inevitável da curvatura** e não têm conserto por orientação: meia-largura
   1,84 vezes o raio local. 24% são largura excessiva com curva folgada.

3. **Os 81 degraus não são um fenômeno só.** 35% são terreno real (o eixo
   explica 99%), 28% são artefato de amostragem (o eixo explica 2% e é plano),
   22% são seção mal posta, 15% ficam indefinidos. **O tamanho do degrau não
   distingue** — a mediana é 1,3 a 1,5 m nas quatro classes.

4. Isso confirma, com número, que **alisar o fundo seria errado**: 28 dos 81
   degraus são topografia de verdade, e 18 são seção fora do lugar — que
   alisamento nenhum conserta.

---

## Arquivos gerados

```
doc/QC_so_mirim_g07_analysis.md    este documento
modelo/overlaps_g07.csv            63 pares, 21 colunas, com a classe
modelo/degraus_so_mirim.csv        81 degraus, 19 colunas, com a classe
scripts/analisar_overlaps.py       metodo dos overlaps
scripts/analisar_degraus.py        metodo dos degraus
```

Nenhuma geometria foi alterada. `so_mirim.g01` continua em `cac7971762`.
