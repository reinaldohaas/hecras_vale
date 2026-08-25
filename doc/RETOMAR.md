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

- **Batimetria do legado aplicada em Mirim, Oeste, Sul, Norte, Açu.** Fluxo:

      .venv/Scripts/python.exe scripts/batimetria_do_legado.py doc/batimetria_<rio>.csv --rio <Rio>
      .venv/Scripts/python.exe scripts/batimetria.py aplicar modelo/<rio>/<rio>.g01 --pontos doc/batimetria_<rio>.csv --saida g02
      .venv/Scripts/python.exe scripts/projeto_rio_avulso.py modelo/<rio>/<rio>.g02 --rio-fonte <Rio>

  Oeste/Sul/Norte: contradeclives → 0, saudável. Açu e Benedito: ver abaixo.

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

## O bloqueio de verdade: EIXO ALTO no Açu R1 e no Benedito

`scripts/diagnostico_eixo_alto.py` mede e figura (`doc/figuras/eixo_alto.png`):
o `eixos_do_relevo.geojson` É a linha esquemática do legado. No vale largo ela
cai sobre o canal e a batimetria é um rebaixamento sadio (~5–10 m). No vale
encaixado de **cabeceira** ela corre pela ENCOSTA, e a seção que o
`rio_do_relevo` corta do MDT pega o talvegue dezenas a centenas de metros ACIMA
do canal real:

- **Benedito**, montante toda (RS 22–44 km): rebaixamento mediana 104 m. Sem
  g02, fora da rede (por isso Indaial inativa).
- **Açu R1** (cabeceira, RS ~143 000–174 000, Rio do Sul→Ibirama): rebaixamento
  mediana 55 m, até 118 m. Ancorar ali cava um cânion, o Açu vai instável
  sozinho, e a rede vai junto — **"Solution Solver Failed" às 03:14 do 1º dia**,
  erro de volume 12,96% de um run ABORTADO (só 4 passos de saída; a solução
  explode desde o 1º passo, sempre no Açu R1/R2). O Açu de JUSANTE (R2/R3, RS <
  140 000) tem batimetria boa (5–10 m).

**Não é problema de batimetria — é de eixo.** Filtrar ponto a ponto não resolve
(tentei; a interpolação entre o que sobra piora o degrau). O conserto é REFAZER
O EIXO de cabeceira do Açu e do Benedito seguindo o talvegue do MDT (acumulação
de fluxo), recortar as seções nesse eixo, e só então ancorar no levantamento —
que a jusante já bate.

## O que fazer agora, nesta ordem

1. **Refazer o eixo de cabeceira do Açu (R1) e do Benedito pelo talvegue do
   MDT** (acumulação de fluxo do relevo), não pela linha esquemática do legado.
   Recortar as seções, reaplicar batimetria (a jusante já bate). É o que
   destrava o Açu — e, com ele, a rede inteira.
2. **Rodar a rede e fechar o balanço.** Com o Açu estável, `montar_rede.py` +
   `rodar_rios.py itajai_rede`. Vigiar a condição inicial: a rede começou com
   185 hm³ nos canais (plausível para calha funda de 174 km) e drenava o
   excesso — se sobrar erro depois de estável, é warmup/restart.
3. **Benedito na rede** depois do eixo: ele reativa a junção Indaial (Açu passa
   a ter 4 reaches, como no legado).
4. **Barragens na rede + calibração de 1983 sem a Norte.** As barragens entram
   como hidrograma roteado nas cabeceiras (`reservatorio.py` → `montar()` de
   `projeto_rede.py` aceita `q_override`). Cuidado: o Canal Retificado precisa
   de nome de rio próprio (`Canal_Retif`), senão dá 1.401 erros de bankline.
5. **Sistema de decisão / visualização.** Cenários (natural / 1983 Sul+Oeste /
   atual três) → rodar a rede → pico de cota nos pontos críticos
   (`barragens_itajai.json`: Rio do Sul, Blumenau, Itajaí, com cotas de
   alerta/emergência). Comparar RELATIVO (redução de pico) é robusto ao erro; o
   ABSOLUTO precisa do zero da régua de cada ponto (dado externo, não no repo).
6. **Projeto JICA:** dados externos (regras de operação, estruturas) não estão
   no repo — pedir/localizar antes.

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
