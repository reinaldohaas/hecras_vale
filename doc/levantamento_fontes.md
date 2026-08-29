# Levantamento de fontes — batimetria, curvas-chave e 1983

*Noite de 26→27/08/2026. O que existe, onde está, e o que já foi trazido
para o repositório.*

## JICA (baixado — `doc/fontes/jica/`)

**Estudo Preparatório para o Projeto de Prevenção e Mitigação de
Desastres na Bacia do Rio Itajaí** (JICA/Nippon Koei, nov/2011) — PDFs
públicos do [openjicareport](https://openjicareport.jica.go.jp), séries
12043618, 12043683 e 12043691.

O que interessa ao modelo:

- **143 seções transversais levantadas em campo (jun–ago/2010)** no canal
  do Itajaí — perfil transversal E longitudinal. A batimetria REAL que o
  modelo não tem. Localizar o anexo com as pranchas (Volume IV "DATA BOOK
  CD" tem os dados brutos; os PDFs Vol. II/III trazem as figuras).
- **Anexo A (Hidrologia)** = `12043691_01.pdf`: enchente de 1983 (§6.3.1),
  curvas H-Q e H-V (Tab. 7.4.6), vazões máximas nas barragens
  (Tabs. 4.5.1/4.6.1), calibração chuva-vazão (Muskingum-Cunge).
- **Fichas das barragens** (§7.5.5–7.5.6, Figs. 7.5.5–7.5.9):
  - *Barragem Oeste* (Taió, 1973): vertedouro na cota **360,0 m**,
    7 comportas nos condutos, capacidade **163 m³/s**, reservatório
    ~83 hm³; TR50 verte +0,9 m.
  - *Barragem Sul* (Ituporanga, 1976): vertedouro **399,0 m**, condutos
    **194 m³/s**, ~93 hm³; altura 43,5 m → base ≈ 355,5 m (bate com o
    talvegue do modelo em RS 32 km do Itajaí do Sul: 354,9 m).
  - *Barragem Norte* (José Boiteux): vertedouro 295,0 m, condutos
    347 m³/s — **concluída só em 1992, fora de julho/1983**.
- Estudo anterior JICA/DNOS (1986–1990): propôs retificação e alargamento
  do canal do Açu (contexto do [[canal-retificado-ausente]]).

## Projeto CRISE (a buscar na FURB)

Pós-enchentes de 1983/84, FURB + CELESC: origem do sistema de alerta,
monitoramento de níveis, modelos de previsão e mapas de risco. Os
relatórios têm os dados hidrológicos de 1983/84 analisados na época.
Sucedido pelo Projeto Itajaí e pelo CEOPS. **Onde pedir**: biblioteca da
FURB / CEOPS (ceops.furb.br).

## CEOPS/FURB

- Artigo **"Curva-chave de Blumenau" (2003)** — estimativa da curva com
  vazões de Indaial para níveis acima de 4,5 m (download do portal CEOPS
  falhou por certificado; baixar manualmente em
  ceops.furb.br → Publicações → Artigos).
- **Cotas de enchente de Blumenau** (ABRH 2013) — baixado:
  `doc/fontes/ceops/cotas_enchente_blumenau_2013.pdf`.
- Registros contínuos de nível em Blumenau desde 1939 (observador) e
  telemetria desde 1984.

## Curvas-chave próprias (feito — `doc/curvas_chave/`)

Ajustadas por potência `Q = a(h−h0)^b` com pares cota×vazão diários da
ANA 1980–1985 (`scripts/ajustar_curva_chave.py`), 11 estações, erro
mediano 1–21%. Blumenau: `Q = 90,6(h+0,51)^1,44`, válida até os 15,19 m
de julho/1983. Figura: `doc/figuras/curvas_chave.png`.

## ANA — amanhã, com o token do usuário

`HidroSeriePerfilTransversal` (perfis medidos nas seções das estações) e
`HidroSerieResumoDescarga` (medições de descarga: largura, área, nível —
a matéria-prima da curva-chave) via HidroWebService OAuth. Estações
prioritárias: Benedito (83660000, 83664000, 83677000/1, 83680000),
Blumenau (83800002), Indaial (83690000), Rio do Sul (83300200).

## Outras pistas

- Dissertação Speckhann (UFSC/LabHidro) — mapeamento de inundação
  Blumenau; Tachini (FURB) — danos de enchente.
- FBDS/RapidEye 5 m: hidrografia vetorial usada nas larguras
  (`doc/fbds/`, [[larguras-do-sigsc]]).
- Base documental do Comitê do Itajaí em aguas.sc.gov.br
  ("base-documental-rio-itajai") — inclui o estudo preparatório em PT.
