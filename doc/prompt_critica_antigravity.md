# Prompt de crítica para o Antigravity — copie tudo abaixo da linha

---

Seu modelo `modelos/_anti/mirim` (gerado por `a_scripts/gerar_mirim.py`) roda
muito bem — 192 h completas, erro de volume 0,0015%, convergência na 1ª
iteração, estável até com dt de 60 s. Crédito dado. Mas ele foi **medido
contra o terreno e contra o dado real**, e tem defeitos que o desqualificam
para o objetivo do projeto (reproduzir a enchente de 1983). Você deve
corrigi-los SEM perder o que fez de bom.

## Regras que valem (do dono do projeto, inegociáveis)

- NÃO invente cotas. NÃO escave canal artificial. NÃO crie pilot channel
  trapezoidal. NÃO force as extremidades das seções para cima.
- Tudo por script, geral (que sirva a qualquer rio), nunca conserto de caso.
- Não altere nada em `legado/` (ler pode, escrever não).
- Arquivos HEC-RAS em CRLF; atenção à cultura do Windows (vírgula decimal).
- No MDT SIG-SC, 0.0 é vazio (água); há folhas com vazio como negativo grande.
- Aceite é NÚMERO medido no HDF que o HEC-RAS constrói — não screenshot,
  não "parece certo".

## Defeitos medidos no seu modelo (números, não opinião)

1. **NÃO HÁ PLANÍCIE DE INUNDAÇÃO.** Suas seções têm 140–280 m de largura
   total (fórmula `70 + frac*70` por margem) num vale que tem QUILÔMETROS de
   várzea real no baixo Mirim — e você ainda LEVANTA as bordas com
   `np.maximum(z, z_lob)`, `z_lob ≥ talvegue+2,5` e paredes de `4,50 m` no
   canal. Isso apaga o armazenamento e o espraiamento — exatamente o que a
   enchente de 1983 fez. Um modelo de cheia sem várzea não modela cheia:
   ele canaliza o pico e entrega níveis errados com convergência bonita.
2. **O LEITO É INVENTADO.** Oito pontos de controle escritos no código
   (275,0 → −2,68) + PCHIP + monotonia forçada, e calha-parábola
   (`z = alvo + (margem−alvo)·dist²`). Nada disso é medido. O repositório TEM
   batimetria real: `legado/Itajai_Rede_1983.g01`, 258 seções levantadas do
   Mirim (cuidado: o trecho RS ~103–114 km do legado é reta desenhada de
   8,00 m/km com resíduo de milímetros — ficção; use o detector em
   `scripts/batimetria_do_legado.py` ou reimplemente o critério: reta local
   com rms < 5 cm E declive ≥ 1 m/km não é levantamento).
3. **A HIDROLOGIA É INVENTADA.** Seu hidrograma de montante tem pico de
   109,64 m³/s; o REAL de 1983 (em `legado/Itajai_Rede_1983.u01`, chave
   `Flow Hydrograph` do `Itajai_Mirim`) tem pico de **1.671 m³/s** — 15×
   maior. Sua "maré" é uma senoide; a real (192 h, −0,20 a +0,80 m) está no
   mesmo arquivo (`Stage Hydrograph` no `Itajai_Acu` RS 75). Sua lateral de
   532 m³/s é chute. Com a cheia real, suas seções de 280 m transbordam no
   primeiro dia — e não há várzea para receber.
4. **`Initial RS=` É CHAVE MORTA.** O HEC-RAS a ignora em silêncio. A chave
   certa é `Initial Flow Loc=`, com a RS em campo de 8 sem decimais (está
   provado no próprio repo; seu modelo só partiu porque o dt fino engoliu a
   partida a frio).
5. **Manning inventado por zona** (0,025/0,045) sem dado de rugosidade.
   Não há levantamento de rugosidade nesta bacia: declare constantes
   justificadas ou use as do projeto (0,032 calha / 0,055 planície), mas não
   apresente número inventado como calibrado.
6. **Não há portões de aceite.** Seu fluxo não mede nada depois de gerar.
   Use os que já existem no repo (ou reimplemente equivalentes):
   `scripts/qc_perfis.py` (ponta n'água, vãos, margens),
   `scripts/ler_erros_geometria.py` (validador real do HEC-RAS, sem solver),
   `scripts/conferir_edge_lines.py` (edge + bank lines lidas do `.g01.hdf`).

## O que preservar do seu trabalho (não regrida)

- A unificação do eixo com o Canal Retificado (já adotada pelo projeto).
- Bank Sta por construção (lb/centro/rb injetados na lista de estacas).
- As chaves de estabilidade do p01 (WFStab/SFStab, DZMax Abort).
- O empacotamento autocontido (um script → projeto completo).

## Implementação obrigatória

1. **Planície do MDT**: cada seção se estende até conter a cheia — a ponta
   sobe ao menos 8–20 m acima do talvegue OU alcança encosta real, medida no
   MDT SIG-SC 1 m (vazio = água; não interpole vazio como terreno). Remova
   TODOS os `np.maximum` de borda e as paredes de 4,50 m — dique só onde o
   MDT mostrar dique.
2. **Leito ancorado em medida**: batimetria do legado onde crível (com o
   detector de ficção), MDT onde não há levantamento. Zero número de leito
   escrito em código.
3. **Hidrologia real de 1983**: hidrograma e maré copiados de
   `legado/Itajai_Rede_1983.u01`. Zero série sintética.
4. **`Initial Flow Loc=`** no lugar de `Initial RS=`.
5. Rode os portões e o solver e entregue os números.

## Aceite numérico obrigatório

- `qc_perfis.py`: GRAVES 0 (nenhuma ponta a < 1 m do talvegue; sem vão >
  1600 m).
- `ler_erros_geometria.py`: 0 Fatal.
- `conferir_edge_lines.py` no `.g01.hdf`: TOTAL 0.
- Solver com a cheia REAL de 1983: 192 h completas, erro de volume ≤ 0,02%,
  iteração máxima ≪ 40.
- **Prova da várzea**: no pico, a largura molhada no baixo Mirim tem de ser
  MUITO maior que a calha (mostre top width no pico por trecho); um modelo
  que passa o pico de 1.671 m³/s com 200 m de largura molhada está errado
  por construção.
- Relatório final: antes/depois de cada número, qual HDF foi medido, e
  arquivos/commits alterados.
