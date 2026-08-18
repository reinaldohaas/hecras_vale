# vale — modelo do Vale do Itajaí a partir do MDT de 1 m do SIG-SC

Programa que constrói o modelo hidrodinâmico 1D do Vale do Itajaí no HEC-RAS,
**um passo por vez**, com todas as decisões expostas e nenhuma tomada às
escondidas.

Usa as ferramentas do próprio HEC-RAS (via `ras-commander`) onde elas existem,
e código próprio só onde não existem.

## Como rodar

```bash
python -m vale
```

Lista os dez passos e marca o que já foi feito. Depois:

```bash
python -m vale 1
```
```bash
python -m vale 1-4
```
```bash
python -m vale tudo --sim
```

Sem `--sim` o programa pergunta antes de cada passo (`S` roda, `n` pula,
`q` para). Cada passo grava o resultado em `modelo/estado.pkl` e o próximo lê de
lá — dá para parar, conferir no RAS Mapper, mudar uma opção e retomar do meio.

Opções na linha de comando, `chave=valor`:

```bash
python -m vale 2-6 selecao=todos corredor_m=500 escavar=false
```

```bash
python -m vale opcoes
```

mostra as 50 opções com os valores atuais. Todas estão em `vale/config.py`, com
o motivo de cada padrão escrito ao lado.

## Os passos

| | passo | o que faz |
|---|---|---|
| 1 | `rios` | catálogo da BHO 2017 da ANA, numerado, e a seleção |
| 2 | `eixos` | eixos dos rios montados das linhas da ANA |
| 3 | `terreno` | mosaico do SIG-SC: 1 m no corredor, 5 m no fundo, e o `.hdf` |
| 4 | `secoes` | corta as seções transversais no terreno |
| 5 | `perfil` | condiciona o perfil longitudinal |
| 6 | `calha` | escava a calha e o pilot channel |
| 7 | `escrever` | grava `.g01`, `.u01`, `.p01`, `.prj`, `.rasmap` |
| 8 | `corrigir` | ferramentas de correção do HEC-RAS e auditoria |
| 9 | `rodar` | computa e lê o log do solver |
| 10 | `visual` | página interativa da cheia |

## Seleção de rios

```bash
python -m vale.rios                    # os 12 atuais (padrão)
python -m vale.rios --sel todos        # os 36 acima de 100 km²
python -m vale.rios --sel 1,2,3,8,10   # por número
python -m vale.rios --sel 1-6,krauel   # faixa e nome, misturados
python -m vale.rios --area 50          # outro limiar
```

O catálogo agrupa por nome **normalizado** (sem acento, minúsculo). Isso corrige
um erro que passava calado: a ANA grava `Rio Itajaí do Oeste` (3.007 km², 68 km,
jusante) e `Rio Itajai do Oeste` (1.103 km², 56 km, montante) — o mesmo rio, com
e sem acento. O casamento por nome acentuado truncava o Oeste em 68 km e
descartava 56 km de cabeceira. A área na foz continuava certa, a vazão entrava
normal, e nada acusava a falta. O mesmo vale para `Ribeirão Dollmann`/`Dollman`.

## Terreno: por que duas resoluções

O SIG-SC são 995 tiles de MDT a 1 m, 118 GB. Nada disso cabe num terreno único
do HEC-RAS, e a maior parte é encosta que a cheia nunca alcança. Então:

- **corredor** a 1 m numa faixa em torno dos eixos (`corredor_m`, padrão
  1.000 m — 1.020 km de eixo dão ~6,7 GB);
- **fundo** a 5 m no resto, para a seção que passar do corredor encontrar
  terreno em vez de NoData.

O `RasTerrain` empilha na ordem recebida, e o primeiro tem prioridade: o
corredor vem antes do fundo.

**MDT não é MDS, e isso muda premissa.** O Copernicus GLO-30 usado antes é
modelo de *superfície*: inclui mata, ponte e a lâmina d'água, gravada como um
plano na cota do espelho. Daí vinham degraus de 12 m nas seções, um corcovo de
9 m no Itajaí do Sul (uma soleira que era copa de mata) e a impossibilidade de
escavar sem contar a profundidade duas vezes. Com MDT o leito submerso continua
ausente do dado, mas o que aparece é terreno — e a batimetria sintética passa a
ser a melhor aproximação disponível em vez de uma duplicação.

## Ferramentas do HEC-RAS usadas

| ferramenta | para quê |
|---|---|
| `RasTerrain.create_terrain_hdf` | terreno do RAS a partir dos GeoTIFF |
| `GeomCrossSection.build_cross_section` | montar as seções (insere ponto na margem) |
| `GeomHtabUtils.calculate_optimal_xs_htab` | tabela hidráulica por seção |
| `RasFixit.fix_bank_stations` | margem que não casa com a tabela |
| `RasFixit.fix_ineffective_flow` | escoamento inefetivo nas seções largas |
| `RasFixit.fix_htab_starting_elevations` | tabela começando na cota errada |
| `RasCheck.run_all` | checagem do RAS antes de computar |
| `RasCmdr.compute_plan` + `HdfResultsPlan` | rodar e ler o log sem a GUI |

Cada uma pode ser desligada (`usar_fixit=false` etc.), e o que faltar na versão
instalada é reportado sem derrubar o resto.

## Duas coisas que o programa não faz sozinho

O log do solver **não** está em `computeMsgs.txt` — o `Compute` via COM não
escreve esse arquivo, só a GUI escreve. O passo 9 lê de dentro do `.p01.hdf`.

E o passo 9 **isola** o projeto numa pasta própria antes de computar, sem levar
`.p01.hdf`, `.u01.hdf` nem `.g01.hdf`. Sem isso: `compute_plan('01')` resolve o
plano dentro da pasta e pode computar outro projeto devolvendo `SUCCESS`; e um
`.u01.hdf` velho faz o solver ler os contornos antigos, produzindo um resumo
idêntico ao da rodada anterior — o que leva à conclusão de que a correção não
fez efeito, quando o que houve foi ler o log de antes.

## Visualização

O passo 10 grava um HTML único, sem servidor e sem biblioteca externa: abre com
duplo clique, roda offline, pode ser mandado por e-mail. Mapa com as seções
coloridas pela **profundidade** (não pela cota — mapa de cota mostra o relevo, o
rio desce 300 m e a cheia desaparece dentro dessa variação), perfil
longitudinal com a envoltória de máxima, seção transversal clicada, hidrograma
no ponto, e linha do tempo com play.

## Estado

Programa escrito e com a sintaxe conferida. **Não foi executado ainda** — o
passo 3 lê 118 GB e é a decisão do usuário quando gastar esse tempo.
