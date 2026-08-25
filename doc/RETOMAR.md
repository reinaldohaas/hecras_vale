# Prompt para a próxima sessão

Copie tudo abaixo da linha.

---

Continuo o modelo HEC-RAS 1D não permanente da bacia do Itajaí (SC), em
`C:\Users\haas\github\hecras_vale`, branch `vale-perfil-isotonico`. Último
commit `2797832`. Leia primeiro este arquivo e depois `git log --oneline -12`,
que os commits explicam o porquê de cada decisão.

**Objetivo final:** rede com os seis rios, as três barragens, calibrada para a
enchente de 1983 — evento em que a Barragem Norte AINDA NÃO EXISTIA (só Sul e
Oeste).

## O que já funciona, medido

Um comando faz geometria, projeto, validação sem solver, pedido de batimetria,
terreno da bacia e conferência da edge line:

    .venv/Scripts/python.exe scripts/construir_rio.py --todos

Os seis rios estão com **0 mensagens e 0 Fatal** no validador do próprio
HEC-RAS, e **0 cruzamentos** na edge line lida do HDF do RAS.

O **Itajaí-Mirim roda**: 40 iterações em 6.892 de 6.900 passos virou **máximo
6**, erro de volume **92,38% → 0,01771%**, tempo de horas sem terminar para
**1,1 min**. Geometria em uso: `modelo/itajai_mirim/itajai_mirim.g02`.

Três defeitos meus em série causavam aquilo, cada um escondendo o seguinte:

1. `Initial RS=` não é chave do HEC-RAS. A certa é `Initial Flow Loc=`, com a
   RS em campo de 8 sem decimais. Ele ignorava a condição inicial CALADO e o
   rio de 114 km partia com vazão ~zero.
2. Profundidade normal na foz, com leito a −9,81 m — contorno de rio de
   montanha num estuário. Agora usa a **maré** copiada de
   `legado/Itajai_Rede_1983.u01` (192 h, −0,20 a +0,80 m), e só nos rios cujo
   fundo da última seção fica abaixo de 2 m: Açu (0,00) e Mirim (−9,81). Norte,
   Sul, Oeste e Benedito terminam a 51–334 m e ficam com profundidade normal —
   o contorno certo deles é a junção.
3. `Geom File=g01` fixo no `.p01`: o plano rodava a geometria SEM batimetria.

**A batimetria estava no repositório.** `legado/Itajai_Rede_1983.g01` tem 1.240
seções levantadas dos seis rios, com calha de 7,5 a 10,9 m de mediana, contra
0,02 m na tirada do MDT — que vê a lâmina, não o fundo. O eixo é o mesmo: a
distância entre `eixos_do_relevo.geojson` e o `Reach XY` do legado é ZERO em
cinco dos seis. Fluxo:

    .venv/Scripts/python.exe scripts/batimetria_do_legado.py doc/batimetria_<rio>.csv --rio <Rio>
    .venv/Scripts/python.exe scripts/batimetria.py aplicar modelo/<rio>/<rio>.g01 --pontos doc/batimetria_<rio>.csv --saida g02
    .venv/Scripts/python.exe scripts/projeto_rio_avulso.py modelo/<rio>/<rio>.g02 --rio-fonte <Rio>

No Mirim: 77 de 88 pontos casados, distância mediana 73 m, rebaixamento
mediano 9,02 m, **117 contradeclives → 0**, declividade máxima 11,60% → 1,33%.

## O que fazer agora, nesta ordem

1. **Consertar `scripts/rodar_rios.py`**: a função `medir()` procura
   `<nome>.p01.computeMsgs.txt`, que o HEC-RAS não gera aqui. Os números estão
   no `.bco01` — erro de volume em "Total Volume Accounting", iterações nas
   linhas `^ *(\d+) <Rio>`. Sem isso a tabela sai toda com `-`.
2. **Aplicar batimetria e rodar os outros cinco**, como o Mirim. Rodar em
   paralelo: `scripts/rodar_rios.py --workers 3 --cores 2`.
3. **Juntar em rede com junções.** O legado mostra como: o Açu tem 4 reaches e
   os afluentes terminam em junção, sem contorno próprio. Cuidado com dois
   defeitos já medidos: junção com 1 entrada e 1 saída o HEC-RAS recusa, e dar
   ao Canal Retificado o mesmo NOME DE RIO do Mirim gera 1.401 erros de
   bankline (ele precisa de nome próprio, `Canal_Retif`).
4. **Barragens** e a calibração de 1983 sem a Norte.

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
condições sem as quais ele crasha. Testado: `initialize OK`, 6 ferramentas.
**Use-o** para ler resultado de plano, mensagens de cômputo e estrutura de HDF,
em vez de abrir HDF na mão com `h5py`.

## Como me cobrar

Não aceite "está pronto" sem número. Peça a medida antes e depois, e diga se a
causa foi diagnosticada ou só contornada. Eu já errei nesta bacia afirmando que
um servidor subia sem ter falado o protocolo com ele, e dando batimetria como
solução de convergência quando a causa era a condição inicial.
