# `g07` — recorte do terreno nos vãos longos: resultado medido

Opção **B** do plano de correção geométrica, limiar de vão **300 m**, aplicada
sobre o `g06`. Nenhum arquivo anterior foi alterado.

Métrica: `scripts/qc_geometria.py` contra o MDT SIG-SC de 1 m, a mesma da
baseline validada em `doc/QC_so_mirim_validation.md`. O contador do
`Validate Geometry` não foi usado.

---

## O que foi feito

541 seções **interpoladas em vão maior que 300 m** deixaram de ser
interpoladas e passaram a ser recortadas do terreno:

- centro no ponto do eixo, na estaca da própria seção;
- direção perpendicular à tangente local, janela adaptativa
  `clip(2 × largura_do_canal, 20, 150)` m;
- **largura preservada** — a que a seção já tinha, sem regra nova;
- estacas = comprimento acumulado da polilinha, logo
  `L_cutline ≡ station[-1] − station[0]` por construção;
- cotas amostradas do MDT 1 m; 140 pontos sem dado herdaram a cota da mesma
  fração da largura;
- **largura de canal preservada**, centrada no cruzamento com o eixo;
- **valores de Manning preservados**; só as quebras andaram, com `n` em 3 casas;
- HTab 2 cm acima do novo talvegue, incremento e contagem preservados.

Verificado: RS, comprimentos de trecho e contagem de pontos **inalterados**;
541 cutlines movidas (mediana 35,9 m, p90 127,8 m, máx 806,8 m), e as listas de
"cutline movida" e "perfil alterado" coincidem exatamente.

---

## Resultado

```
                                        g06        g07
status OK                               574        758
status WARNING                          678        628
status CRITICAL                         166         32

angulo fora da perpendicular            253         25     <- alvo
cruza o eixo 2+ vezes                    67         21
cruza cutline vizinha                    66        106     <- PIOROU
UNIAO das tres                          278        118

angulo: mediana (90 = ideal)           82,3       88,7
angulo: p05                            20,1       73,4
secoes com angulo < 10 graus             42         12
cobertura MDT minima                  0,463      0,685
|L_cut - estacas| maximo (m)         0,0156     0,0156

degrau de fundo: mediana (m)          0,100      0,120
degrau de fundo: maximo (m)            2,77       2,75
degraus acima de 1 m                     44         81     <- PIOROU
degraus acima de 2 m                      5         11     <- PIOROU
```

---

## Leitura

**Funcionou no que era o alvo.** A categoria de ângulo, que subsumia as outras
duas, caiu de 253 para 25 — 90% de redução. A mediana do ângulo foi de 82,3°
para 88,7°, e o p05 de 20,1° para 73,4°: a cauda de cutlines quase paralelas ao
fluxo praticamente desapareceu. Interseções múltiplas caíram de 67 para 21. A
união das três categorias foi de 278 para 118, e o CRITICAL de 166 para 32.

**Piorou em dois lugares, e os dois são explicáveis.**

*Overlap com a vizinha: 66 → 106.* Consequência direta de acertar o ângulo. Uma
seção perpendicular ao eixo, num meandro fechado, converge com a vizinha pelo
lado de dentro da curva. Antes elas estavam tortas e não se encontravam; agora
estão certas e se cruzam. É a mesma tensão entre largura e raio de curvatura já
medida antes: com meia-largura mediana de 66 m e raios de meandro menores que
isso, o cruzamento é geométrico, não corrigível pela orientação.

*Degraus de fundo: 44 → 81 acima de 1 m, 5 → 11 acima de 2 m.* Cada seção
recortada pega o mínimo do MDT na sua própria linha, que agora está em lugar
diferente. É o mesmo mecanismo de amostragem independente que já produziu o
degrau da RS 127448.69 — validado ali como 1,33 m de terreno real e 1,44 m de
artefato. O máximo não subiu (2,77 → 2,75 m), mas a quantidade dobrou.

---

## O que isto não resolve, e o que não foi testado

- As 106 de overlap continuam abertas, e agora são a maior categoria. Não há
  correção de orientação que as elimine.
- Os degraus de fundo pioraram; nenhum foi alisado, conforme a regra.
- **A hidráulica não foi testada.** O `g07` não foi rodado no HEC-RAS.

---

## Como testar

Projeto pronto em `modelo\mirim_g07\`:

```
mirim_g07.prj      Proj Title=mirim_g07
mirim_g07.g01      Geom Title=mirim_g07 (perfil do MDT SIG-SC 1 m)
                              + margens + htab + vaos longos do terreno
mirim_g07.p01      o mesmo plano (Mixed Flow, dt 5SEC, 192 h)
mirim_g07.u01      o mesmo fluxo
mirim_g07.rasmap   + Terrain\ (Copernicus, para fundo de tela)
```

Sem resultados: nada foi computado. Abrir e rodar, ou pedir que se rode.

**Advertência de tela**, a mesma do `mirim_mdt`: a camada de terreno é o
Copernicus, ~7 m acima do MDT de 1 m de onde vieram os perfis. As seções vão
parecer enterradas. Isso é o esperado.

---

## Cadeia completa, com SHA1

```
so_mirim.g01   cac7971762   original, nunca tocado
so_mirim.g03   aeba7a8403   perfil recortado no MDT 1 m
so_mirim.g05   4115903669   g03 + 13 margens (Manning junto)
so_mirim.g06   c0196b254e   g05 + HTab reancorado
so_mirim.g07   c42a6f744c   g06 + 541 secoes recortadas em vao > 300 m
                            (copiada em modelo\mirim_g07\mirim_g07.g01)
so_mirim.g04   c08482779d   nao usar
```
