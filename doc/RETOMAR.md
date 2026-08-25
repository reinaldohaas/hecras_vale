# Prompt para a próxima sessão

Copie tudo abaixo da linha.

---

Continuo o modelo HEC-RAS 1D não permanente da bacia do Itajaí (SC), em
`C:\Users\haas\github\hecras_vale`, branch `vale-perfil-isotonico`. Leia
primeiro este arquivo e depois `git log --oneline -15`, que os commits explicam
o porquê de cada decisão.

**Objetivo final:** rede com os seis rios, as três barragens, calibrada para a
enchente de 1983 — evento em que a Barragem Norte AINDA NÃO EXISTIA (só Sul e
Oeste). Depois: projeto JICA e um sistema de visualização e decisão para
minimizar cheias.

## O que já funciona, medido

- **Verificação desta retomada:** os scripts citados existem e compilam. Branch
  atual: `vale-perfil-isotonico`; último commit: `fc1bde9 feat(pipeline): a
  batimetria entra no UM comando, atras do porteiro de eixo`. Atenção: existe
  `modelo/itajai_acu/itajai_acu.g02` no disco, mas ele NÃO deve ser tomado como
  válido: o porteiro atual ainda mede eixo alto no Açu
  `(11,25–31,80 km, rebaixamento até 117,4 m)`. Se rodar o pipeline atual, esse
  g02 velho deve ser removido até o eixo ser refeito.

- **`medir()` do `rodar_rios.py` consertado.** Lê do `.bco01` (o HEC-RAS não
  gera `computeMsgs.txt` aqui): erro de volume no bloco "Total Volume
  Accounting", iteração no `<i>` que abre cada linha `<i> <Reach> ...`. E o pool
  virou de PROCESSOS: em threads o estado global do `ras_commander` fazia os
  projetos se atropelarem (o init do norte lia o u01 do acu). Os cinco rios
  isolados rodam limpos em 3,1 min.

- **Rios isolados (medido, tabela de veredito):**
  - Mirim `iter 6  0,02%`, Oeste `iter 4  0,02%`, Sul `iter 3  0,02%` — convergem.
  - Norte `iter 40  0,01%` — fecha o volume mas itera no teto (37 seções que a
    batimetria não pôde rebaixar deixam degraus locais).
  - Açu `iter 40  252%` — **falha**, ver "eixo alto" abaixo.

- **UM comando faz tudo, do MDT ao g02 — usar SEMPRE ele, nunca os scripts
  soltos** (regra do usuário: fazer um software e usá-lo, não SER o software):

      .venv/Scripts/python.exe scripts/construir_rio.py --todos

  O `construir_rio.py` agora inclui o passo 5: preenche a batimetria do
  legado, passa o PORTEIRO DE EIXO ALTO, aplica o g02 e reaponta o projeto.
  O porteiro barra rio com trecho contíguo ≥ 1 km pedindo rebaixar > 25 m
  (critério medido: sadios têm blips de 0,0–0,2 km; quebrados têm 20,6 e
  24,6 km), gera a figura e remove g02 velho — duas rodadas dão o mesmo
  estado. Verificado: Mirim sai `g02 aplicado` (77 casados, 117
  contradeclives→0); Benedito sai `REFAZER EIXO (0–25 km, reb até 284 m)`
  sem g02. O passo 8 (religar terreno) foi corrigido para reapontar ao g02
  quando existe — antes reescrevia `Geom File=g01` calado.

- **A rede MONTA e VALIDA (mas o solver não fecha ainda).**

      .venv/Scripts/python.exe scripts/montar_rede.py   # costura + aperta ate 0 Fatal

  Topologia do legado, não inventada: Rio_do_Sul (Sul+Oeste→Açu R1), Ibirama
  (Norte→R2), Itajaí (Mirim→R3→mar). O Açu é o tronco, partido nos pontos de RS
  dos reaches do legado. Só entram junções cujo afluente tem g02, então Indaial
  (Benedito) fica inativa e o Açu passa reto — sem cair na junção 1-entrada-
  1-saída que o HEC-RAS recusa. O laço descarta as seções que o próprio
  validador rejeita por cruzarem o reach vizinho na confluência (22 seções:
  Oeste 7, Sul 7, Açu 6, Mirim 1, Norte 1) e chega a **0 mensagens, 0 Fatal**.

- **Reservatórios (barragens) roteando hidrograma**, com dados reais de
  `dados_estruturas/barragens_itajai.json` (capacidade, vazão do extravasor):

      .venv/Scripts/python.exe scripts/reservatorio.py --barragem barragem_sul --rio Itajai_Sul

  Detenção por balanço de massa; corte ótimo calculado do próprio hidrograma.
  Pico: Sul −40%, Oeste −37%, Norte −75%. É a alavanca do sistema de decisão.

## A descoberta que mudou tudo: o LEGADO E FICCAO nas cabeceiras

Nas cabeceiras do Açu (RS ~141–164 km) e do Benedito (RS 18–44 km) — e no
montante do Mirim (RS 103–114 km) — o "fundo levantado" do legado é uma RETA
DESENHADA: declive exatamente 8,00 m/km com resíduo rms de 1–2 milímetros por
dezenas de seções, a mesma constante nos dois rios ("rede real ANA + relevo
DEM", diz o título do próprio legado). Rio real dá rms de 1–15 m (controles).
O diagnóstico anterior ("eixo pela encosta") estava ERRADO e está corrigido em
`diagnostico_eixo_alto.py`: buscado o caminho conectado de menor cota num
corredor de ±1500 m, não há vale naquela cota — o eixo está certo, o legado é
que mente ali. Ancorar nele pediria cavar até 284 m e foi o que instabilizou o
Açu ("Solution Solver Failed" às 03:14) e, por ele, a rede.

O conserto está NO PIPELINE, medido nos seis rios:

- `batimetria_do_legado.py`: detector de ficção (reta local com rms < 5 cm E
  declive ≥ 1 m/km — canal dragado é reta PLANA e é real) + filtro de
  rebaixamento implausível (> 25 m numa calha de ~11) + casamento por RS
  (folga 800 m; XY só de sanidade 1200 m).
- `batimetria.py aplicar`: interpolação absoluta entre âncoras (zera
  contradeclives: Mirim 117→0), e rebaixamento ZERADO por intervalo em
  aglomerados de pontos sem âncora ≥ 3 km, com rampa de 1 km — sem ponte por
  baixo do terreno. Medido: Açu ajuste máx 117→24 m, contradeclives 125→17;
  Norte 62→0; Sul 50→0; Oeste 59→0; Benedito máx 5,5 m.
- `rio_do_relevo.py`: onde o MDT não cobre (< 50% do centro) ou onde a zona é
  de maré (talvegue < 2 m e legado com fundo > 2 m abaixo), a seção levantada
  do legado entra INTEIRA (cutline, largura e margens dela — o recorte no grid
  de ±400 m cortava o canal do porto de 2.116 m pela metade); o pente de
  lâmina entre levantadas cai; ponta degenerada (restinga de 76 m recebendo a
  maré) cai. Foz do Açu: 8 últimas seções levantadas, fundo −10,1→−10,8,
  4,2–4,4 km de largura.

## Estado 25/08 (tarde): Açu aprovável, contrato calibrado pela referência

Rodada do usuário 14:27 com o gerador redesenhado (500 m adaptativo, subida
20 m): **222 seções, 0 pontas n'água, GRAVES 7 (só vãos), bank lines 0/0,
41 dobras de edge**. Referência (legado 1983, mesma régua): 348 dobras, 556
bank×eixo, 129 Fatal, 12 GRAVES — e roda. O aviso do Mapper sobre edge lines
é da superfície de MAPA, não do solver; a referência o exibe 8× pior.

Decisões do usuário: portões DUROS = GRAVES 0 + bank 0/0; dobras/Fatal são
MACIOS (aceitos ≤ referência/5 = 70/26, relatados no veredito). Os 7 vãos
eram podas reabrindo espaçamento depois da seleção — o gerador agora TAPA
vãos > 1400 m reinstalando a melhor candidata medida. Seguir para a REDE.

## O que fazer agora, nesta ordem

1. **Rodar o pipeline inteiro** (`construir_rio.py --todos`) — o usuário roda
   no terminal; o programa imprime o aceite (validador 0/0, banks 0, adotadas
   na foz, tabela com "g02"). Depois `montar_rede.py` — com o Benedito agora
   COM g02, a junção Indaial reativa sozinha (Açu em 4 reaches, como o
   legado).
2. **Rodar a rede** (`rodar_rios.py itajai_rede`) e ler o veredito pelo MCP
   (`get_plan_results_summary`: foi "Unsteady Went Unstable" que desmascarou o
   run abortado). Vigiar a condição inicial (a rede partia com 185 hm³ e
   drenava o excesso).
## Estado dos scripts novos desta rodada

- `construir_rede.py` — costura a rede, parte o tronco, escreve as junções.
- `projeto_rede.py` — prj/p01/u01 da rede; `montar(geom, q_override=)` reusável.
- `montar_rede.py` — laço costura→valida→descarta seção que cruza vizinho.
- `reservatorio.py` — roteia hidrograma pelas barragens (dados reais).
- `diagnostico_eixo_alto.py` — mede e figura o eixo alto (Açu R1, Benedito).

## Regras que valem, do usuário

- **NÃO invente cotas. NÃO escave artificialmente o canal. NÃO crie pilot
  channel trapezoidal. NÃO force as extremidades das seções para cima.**
- Não desloque nem estenda seções automaticamente só porque o talvegue do DEM
  está perto da extremidade.
- **Não altere `so_mirim.g01`, `so_mirim.g01.hdf` nem qualquer arquivo
  original** em `legado/`. Ler pode; escrever não.
- No MDT SIG-SC, **0.0 é vazio**, e há folhas que gravam vazio como número
  negativo grande. Não use GDAL para montar mosaicos com `-srcnodata 0` sobre
  fonte que já declara nodata.
- **Tudo por script**, nunca conserto na mão, e um script que sirva a qualquer
  rio — não uma correção para um caso.
- **Todo erro detectado vira figura**, de preferência gerada pelo programa.
- Rodar com `py -3.10` ou `.venv/Scripts/python.exe`, caminho Windows.
- O HEC-RAS lê número com a **cultura do Windows**: com vírgula decimal no
  sistema, 227.19 vira 22719 e nenhum modelo abre.

## MCP

`.mcp.json` configura o servidor `hecras` (ras-commander-mcp) apontando para
`.venv/Scripts/ras-commander-mcp.exe`, com `mcp<2` e o h5py do venv — as duas
condições sem as quais ele crasha. **Use-o** para ler resultado de plano,
mensagens de cômputo e estrutura de HDF (foi assim que se viu que a rede ia
instável: `get_plan_results_summary` deu "Unsteady Went Unstable"), em vez de
abrir HDF na mão com `h5py`.

## Como me cobrar

Não aceite "está pronto" sem número. Peça a medida antes e depois, e diga se a
causa foi diagnosticada ou só contornada. Eu já errei nesta bacia afirmando que
um servidor subia sem ter falado o protocolo com ele, dando batimetria como
solução de convergência quando a causa era a condição inicial, e dizendo que a
rede "fechava em 13%" quando na verdade o solver tinha ido instável e abortado.
