# -*- coding: utf-8 -*-
"""
Configuracao de uma execucao. Tudo que se decide fica aqui, num objeto so.

O programa nao esconde escolha nenhuma dentro do codigo. Cada constante que
muda o resultado esta neste arquivo, com o valor padrao e o motivo dele, e
pode ser trocada pela linha de comando ou por um arquivo .json. Quem roda
decide; o programa executa e mostra o que fez.

Os valores padrao nao sao chutes: sao o que sobrou de uma reconstrucao inteira
do modelo do Vale, e cada comentario abaixo registra o que acontece quando o
valor esta errado.
"""
import json
import os
from dataclasses import dataclass, asdict, field, fields

EPSG = 31982                # SIRGAS 2000 / UTM 22S -- o CRS do SIG-SC

# WKT NO DIALETO ESRI, que e o que o HEC-RAS pede: o proprio dialogo diz
# "the ESRI Projection file (*.prj)". Duas versoes anteriores falharam --
# a escrita a mao (sem nos AUTHORITY) e a do GDAL (com eles). Nas duas o
# RAS Mapper carregava o arquivo e mostrava a caixa "Definition" VAZIA, e
# depois estourava com "Referencia de objeto nao definida".
#
# O ESRI e outro dialeto, nao um WKT melhor ou pior: "GCS_SIRGAS_2000" e
# "D_SIRGAS_2000" com sublinhado, sem AUTHORITY, numeros com casa decimal
# explicita. O GDAL gera com to_wkt("WKT1_ESRI").
#
# E O ARQUIVO NAO PODE SE CHAMAR "<projeto>.prj": esse nome ja e o do
# arquivo de PROJETO do HEC-RAS. Ver NOME_SRS abaixo.
WKT = ('PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",'
       'DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,'
       '298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",'
       '0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
       'PARAMETER["False_Easting",500000.0],'
       'PARAMETER["False_Northing",10000000.0],'
       'PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",'
       '0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')

# A EXTENSAO TEM DE SER .prj. Fugir do choque de nomes trocando a extensao
# para ".projection" resolve a colisao e cria outra falha: o RAS Mapper abre
# a pagina de projecao e estoura em
# MapperOptionProjection.DisplaySRSFileText -> "Referencia de objeto nao
# definida", porque nao reconhece o arquivo como SRS. O mesmo ja estava
# escrito em terreno.py sobre o RasProcess ("passar um .projection faz o
# RasProcess falhar"), e a licao nao tinha sido aplicada aqui.
#
# A saida e manter .prj e trocar o NOME: o conflito e com "<projeto>.prj",
# nao com a extensao. Um nome fixo, do sistema de coordenadas, tambem deixa
# claro que o arquivo nao pertence a um projeto so.
NOME_SRS = "SIRGAS2000_UTM22S.prj"

RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"


@dataclass
class Opcoes:
    # ----------------------------------------------------------- entrada
    projeto: str = "vale"
    sigsc: str = r"C:\Users\haas\Downloads\sigsc"
    bho: str = "rios_itajai.geojson"
    saida: str = "modelo"          # pasta de trabalho

    # ------------------------------------------------------------- rios
    area_minima: float = 100.0     # km2 para entrar no catalogo
    selecao: str = "atuais"        # 'todos' | 'atuais' | '1,3,5' | '1-6,krauel'

    # ---------------------------------------------------------- terreno
    # FONTE do terreno. O que muda nao e so a resolucao: muda o TIPO de modelo,
    # e com ele uma premissa do resto do programa.
    #
    #   'copernicus'  GLO-30, um GeoTIFF de 30 m que ja esta no disco. Segundos.
    #                 E modelo de SUPERFICIE: inclui mata, ponte e a LAMINA
    #                 D'AGUA, gravada como um plano na cota do espelho.
    #   'sigsc'       995 tiles de MDT a 1 m, 118 GB. Minutos a horas -- o custo
    #                 e LEITURA, nao reamostragem: o corredor toca centenas de
    #                 tiles e o fundo sobre a bacia toca quase todos.
    #   'misto'       corredor a 1 m do SIG-SC, fundo do Copernicus. Corta a
    #                 passada mais cara (o fundo sobre a bacia inteira) e
    #                 mantem 1 m onde a cheia e decidida.
    #
    # ATENCAO, e o programa cuida disso sozinho: com fonte de SUPERFICIE a
    # escavacao da calha e desligada. O dado ja contem a lamina, entao escavar
    # conta a profundidade DUAS VEZES -- o modelo parte seco e o assentamento
    # tem de encher centenas de hm3 de canal inventado. Ver vale/calha.py.
    fonte: str = "copernicus"
    # VRT em vez de mosaico fisico. Os tiles ja estao a 1 m e no CRS do
    # modelo: nao ha o que reamostrar. O mosaico fisico escrevia 21,9
    # bilhoes de pixels para guardar 1,67 bilhao -- a caixa envolvente do
    # corredor cobre a bacia toda. Duas horas davam 10%.
    # Resolucao do SIG-SC no modelo. A 1 m o dado nao cabe em lugar nenhum e
    # obriga leitura por janela; a 10 m a bacia inteira sao ~105 milhoes de
    # pixels (420 MB) que cabem em memoria de uma vez. E continua sendo MDT --
    # o ganho sobre o Copernicus nao esta na resolucao, esta em NAO ter 25 m de
    # copa de mata sobre o leito.
    res_sigsc: float = 10.0
    # VRT (sem copiar nada) so faz sentido na resolucao NATIVA. Reamostrando,
    # o mosaico fisico e melhor: 420 MB lidos uma vez contra reamostrar a cada
    # acesso.
    vrt: bool = True
    # Gerar o .hdf de terreno do HEC-RAS? So o RAS MAPPER precisa dele: as
    # secoes sao cortadas do VRT e o solver nao o le. O RasTerrain materializa
    # tudo num GeoTIFF proprio -- sobre 10.494 km2 a 1 m sao ~42 GB e uma hora.
    # Desligue quando o objetivo for so rodar o modelo.
    terreno_hdf: bool = True
    # Tempo limite do RasTerrain, em segundos. O padrao DELE e 600, e sobre
    # 10.494 km2 a 1 m isso estoura: 55 minutos e 41 GB escritos sem produzir
    # o arquivo. Materializar tudo leva horas -- e por isso que o padrao aqui
    # e nao gerar o .hdf quando a fonte e o SIG-SC inteiro.
    terreno_timeout: int = 7200
    copernicus: str = "Terrain/Terreno_Copernicus.tif"

    # Corredor em 1 m em torno dos eixos, fundo mais grosseiro no resto (so
    # vale para 'sigsc' e 'misto'). 1.020 km de eixo com 1.000 m de
    # meia-largura dao ~6,7 GB.
    corredor_m: float = 1000.0
    res_corredor: float = 1.0
    res_fundo: float = 5.0
    fundo: str = "bacia"           # 'bacia' | 'mosaico' | 'nenhum'

    # -------------------------------------------------------------- 2D
    # O modelo bidimensional nao usa NADA do bloco de secoes abaixo: nao ha
    # secao, nao ha entalhe piloto, nao ha vazao inicial. Ver vale/malha.py.
    #
    # ATENCAO ao corredor: `corredor_m` acima tem de ser >= `buffer_2d`, senao
    # a borda da area 2D cai fora do terreno de 1 m e a malha assenta sobre o
    # fundo grosseiro justamente onde a agua espalha. Com buffer de 1.500 m,
    # corredor_m=1000 nao serve.
    buffer_2d: float = 1500.0      # meia-largura da area 2D, a partir do eixo
    perimetro_max_pontos: int = 1200

    # Celula grossa sobre terreno fino e o ponto do 2D: as tabelas de
    # sub-grade guardam a curva cota-volume do MDT de 1 m DENTRO de cada
    # celula de 100 m. Quem representa a calha e a tabela, nao a celula.
    celula: float = 100.0
    refino_2d: float = 25.0        # celula na faixa da calha; 0 desliga
    refino_largura: float = 300.0  # meia-largura da faixa refinada
    face_minima: float = 0.05      # min_face_length_ratio do gerador
    n_2d: float = 0.06             # Manning da planicie

    # Contornos. Todos sao arcos DO PERIMETRO -- o HEC-RAS nao aceita contorno
    # por dentro da area. A tampa inteira de um buffer de 1.500 m e um
    # semicirculo de 4,7 km, e lancar a vazao nos 4,7 km poe agua na meia
    # encosta dos dois lados; por isso o arco e recortado pelo meio.
    bc_largura: float = 1000.0     # arco de montante e da foz
    n_laterais: int = 8            # por lado; 0 concentra tudo na cabeceira
    bc_lateral_largura: float = 300.0
    bc_minima: float = 150.0       # arco menor que isto e descartado
    decl_bc: float = 0.0010        # Flow Hydrograph Slope dos contornos

    # DWE (onda difusiva) na primeira rodada, de proposito: SWE-ELM e mais
    # fisico e muito menos tolerante, e se o modelo nao fecha volume em DWE
    # nao ha por que atribuir a diferenca a fisica.
    equacao_2d: str = "DWE"        # 'DWE' | 'SWE-ELM'
    intervalo_2d: str = "1MIN"
    courant_2d: float = 1.0        # 0 usa passo fixo
    ztol_2d: float = 0.01
    max_iter_2d: int = 20
    ic_horas: float = 0.0          # aquecimento; 2D comeca SECO

    # ----------------------------------------------------------- secoes
    # O Acu cai 195 m em 13 km na garganta do Salto Pilao. A 1 km de
    # espacamento sao 8 m de queda ENTRE SECOES VIZINHAS e o solver falha no
    # primeiro passo; o criterio dx <~ 0,15*D/S da ~75 m ali. Por isso o
    # espacamento sai da declividade local, e nao de um numero fixo.
    #
    # E DEPOIS O CRITERIO FICOU SO NO COMENTARIO. Estas quatro linhas
    # interpolavam entre 1.000 m e 150 m conforme a declividade, e 150 m e o
    # DOBRO dos 75 m que o proprio comentario acima calcula. No Benedito o
    # resultado foi 134 dos 147 trechos (91%) fora do criterio, com dx exigido
    # mediano de 18 m contra os 150 m usados -- oito vezes grosso demais, ao
    # longo de um rio que nunca terminou de rodar. Agora `samuels` faz valer o
    # que estava escrito: a interpolacao continua sendo o alvo, e o criterio e
    # um TETO por cima dela.
    espacamento: float = 1000.0
    espacamento_min: float = 150.0
    decl_plano: float = 0.0010
    decl_ingreme: float = 0.0060

    # Samuels (1989), dx <= k*D/S. D e a profundidade caracteristica de calha
    # cheia; k=0,15 e o valor classico. O PISO existe porque o criterio nao
    # tem fundo: a 5% ele pede 4,5 m, o que daria mais de 10 mil secoes no
    # Benedito. O piso e so um limite de custo.
    #
    # E NAO E CAUSA DE INSTABILIDADE, ao contrario do que este comentario
    # afirmava. O Benedito tem 56% dos vaos abaixo do piso e, sozinho, completa
    # as 192 h com 0,024% de erro -- com as MESMAS 819 secoes que ele tem
    # dentro da rede, onde a simulacao cai. Mesmo espacamento, dois desfechos.
    #
    # NASCE DESLIGADO, por decisao do dono do projeto, e a medicao sustenta.
    # Ligado, ele densifica o CORTE do terreno: no Mirim produziu 1.553 secoes
    # cortadas com mediana de 58 m -- 79% dos vaos abaixo dos 150 m pedidos e
    # 16% mais finos que o pixel de 30 m do Copernicus, ou seja amostrando
    # mais fino que o dado. Desligado, saem 432 cortadas com minimo de 150 m.
    #
    # E a densidade que ele impunha NAO e exigida pela numerica. Medido no
    # so_mirim.g08, com o dt de 15 s do plano e as velocidades da propria
    # rodada (mediana 0,31 a 0,40 m/s no canal, p90 0,69 a 0,92):
    #
    #     Courant = V*dt/dx = 0,10 a 0,29        -- com folga para 1
    #     distancia percorrida por passo: 4,5 a 13,5 m, contra dx de 47 m
    #
    # Ou seja: o espacamento denso nunca veio de exigencia de passo de tempo.
    # Veio deste criterio, e o preco dele e amostrar terreno inexistente.
    samuels: bool = False
    samuels_k: float = 0.15
    samuels_D: float = 1.5         # m; profundidade de calha cheia tipica
    samuels_leopold: bool = True   # D = kh*A^eh por posicao, e nao fixo

    # AMOSTRAR o terreno e RESOLVER o escoamento sao coisas diferentes. O corte
    # do terreno usa `espacamento`/`espacamento_min` (geometria: de quanto em
    # quanto o vale muda); a densidade que o solver exige e obtida DEPOIS,
    # interpolando entre secoes reais com o interpolador da propria biblioteca
    # (GeomCrossSection.interpolate_station_elevation). Cortar do terreno na
    # densidade numerica dava 1.553 secoes no Mirim -- e a declividade que
    # pedia isso vinha do dossel do Copernicus, nao do rio.
    interpolar: bool = True
    interp_max: int = 40           # teto de intermediarias por par
    # Salto maximo de espacamento entre vaos vizinhos. Samuels nao limita a
    # planicie (S->0 => dx->infinito), e o Cedros saiu com vaos de 944 m ao
    # lado de 25 m -- razao de 14x. Foi no trecho plano e grosseiro que o
    # par Benedito+Cedros instabilizou, com 2 cm de lamina.
    razao_dx: float = 2.0

    # Janela em que a DIRECAO do eixo e medida, e com ela a perpendicular
    # da secao. Era `espacamento_min * 1.7` = 255 m, herdado de quando as
    # secoes eram cortadas a 150 m; hoje a mediana e 25 m e 255 m e a corda
    # de um meandro inteiro, nao a tangente. Medido no Benedito com 255 m:
    # mediana de 13,4 graus fora da perpendicular, p90 de 35,8 e MAXIMO DE
    # 90 -- secao paralela ao rio. As bank lines e edge lines saem da secao,
    # entao saem tortas junto, e foi assim que isto apareceu: olhando o
    # RAS Mapper, nao o log.
    #
    # Medido por janela (mediana / p90 / max / fora de 30 graus):
    #     255 m   11,8   32,2   86,9   12%
    #     120 m    5,6   18,2   58,4    2%
    #      60 m    2,5   10,4   27,6    0%
    #      30 m    1,1    4,8   16,6    0%
    #
    # 60 m e o maior que zera o descontrole. Menor que isso comeca a pegar
    # o serrilhado da digitalizacao da ANA, que e por que a suavizacao
    # existe -- sem ela o RAS acusa "edge lines have self intersections".
    janela_direcao: float = 60.0

    # ------------------------------------------------- calha no terreno
    # A calha passa a ser escrita NO TERRENO por
    # RasTerrainModWriter.add_channel_modification, e nao desenhada secao a
    # secao. Ver vale/canal.py. Como a funcao recebe UM depth e UMA width por
    # polilinha, o rio e partido em segmentos de profundidade quase constante.
    calha_no_terreno: bool = True
    tolerancia_prof: float = 0.5    # m; variacao aceita dentro de um segmento
    passo_polilinha: float = 50.0   # m entre pontos da polilinha do canal
    talude_canal: float = 3.0       # H:V
    prof_minima_canal: float = 0.5  # m
    espacamento_piso: float = 25.0  # m; abaixo disto o modelo 1D nao vale
    n_pontos: int = 280            # limite do HEC-RAS e 450
    # Espacamento ALVO entre pontos da secao. Numero fixo de pontos numa
    # secao estreita da espacamento menor que o pixel do terreno (0,62 m
    # numa secao de 175 m contra 30 m do Copernicus): sao pontos
    # inventados, e foi um deles que colidiu com a estaca da margem e
    # fez o HEC-RAS recusar a geometria por ponto duplicado.
    espacamento_pontos: float = 5.0
    n_pontos_min: int = 40
    # Numa curva de raio R as perpendiculares vizinhas se encontram a R do eixo
    # do lado concavo. Passar disso e o que cruza as cutlines -- 24% dos pares
    # na primeira tentativa -- e a mancha do RAS Mapper perde sentido.
    folga_curva: float = 0.70

    # O MESMO CRITERIO DE CURVA, AGORA SOBRE A GEOMETRIA FINAL. O limitador
    # acima age dentro do cortar(), sobre as secoes CORTADAS e com as larguras
    # de ANTES do recorte pela cheia. Depois dele entram a densificacao (que
    # insere secoes com outro espacamento) e o recorte (que muda a largura
    # toda) -- entao a condicao de nao se cruzarem nunca chegava a valer para
    # o que o HEC-RAS recebe.
    #
    # Medido no Mirim: meia-largura mediana de 66 m em meandros cujo raio e
    # menor que isso. As vizinhas se cruzam por dentro da curva e a edge line
    # que liga as pontas da laco, com o RAS avisando "The generated edge lines
    # have self intersections".
    #
    # E NASCE DESLIGADO, porque foi medido e PIOROU. Com ele, as cutlines que
    # se cruzam caem de 65 para 34 -- e o Validate Geometry sobe de 240 para
    # 995, com o defeito mudando de lugar: sai das edge lines e vai para as
    # proprias secoes, que ficam curtas e desencontradas (largura minima de
    # 120 m para 69 m).
    #
    # A licao, que ja custou tres tentativas neste dia -- fixar o lado da
    # estaca 0, apertar `folga_curva` para 0,45, e este limitador -- e que
    # melhorar uma secao ISOLADA piora as linhas que passam por ela. Bank line
    # e edge line ligam ponta a ponta entre vizinhas: o que conta e a
    # concordancia com a vizinha, nao a qualidade de cada uma.
    #
    # Fica aqui, desligado, com o numero ao lado: quem tentar de novo comeca
    # sabendo o que ja foi medido.
    curva_pos: bool = False
    curva_piso: float = 45.0       # m; meia-largura minima ao apertar

    # PUXAR O EIXO PARA O TALVEGUE antes de cortar. O tracado da BHO 2017 da
    # ANA e esquematico: no Mirim o talvegue lido do terreno fica a 16 m dele
    # na mediana, 42 m no p90 e 296 m no pior caso, contra 26 m de meia-calha.
    # Em 28% das secoes o rio real esta fora da calha declarada e em 12% as
    # duas margens caem do mesmo lado do eixo -- e e por isso que as bank
    # lines geradas pelo HEC-RAS cruzam o eixo 320 vezes.
    #
    # E NASCE DESLIGADO POR UM MOTIVO MEDIDO: no Copernicus nao ha talvegue
    # lateral para seguir. Medindo 1.070 perfis de +-80 m no Mirim, o quanto o
    # terreno DESCE do eixo ate o minimo da janela:
    #
    #     mediana 0,00 m    p75 0,84 m    p90 3,98 m
    #     65% dos perfis descem menos de 10 cm
    #
    # O eixo da ANA ja esta no fundo em dois tercos do rio. O GLO-30 e modelo
    # de SUPERFICIE e traz a LAMINA D'AGUA, que e plana: dentro do rio a cota
    # nao varia, o minimo lateral e um empate e o `argmin` devolve uma posicao
    # arbitraria. Uma penalidade de 0,001 m por metro -- 8 cm ao longo da
    # janela inteira -- ja anula a correcao, o que mostra o tamanho do sinal.
    #
    # Ou seja: os "16 m de afastamento entre eixo e talvegue" que eu havia
    # medido nao sao o eixo fora do rio, sao empates num plano de agua. Ligar
    # isto no Copernicus ajusta RUIDO, e deixa o eixo mais quebrado (giro
    # maximo entre vertices de 117 para 167 graus).
    #
    # COM O MDT DO SIG-SC A 1 m (fonte=sigsc), que e solo exposto, o talvegue
    # existe de verdade e vale ligar. Meca antes de trocar o padrao.
    #
    # As pontas NAO se movem: sao a conexao com a juncao, onde o snapping e
    # exato. `eixo_taper` e o comprimento em que a correcao afunila ate zero.
    eixo_talvegue: bool = False
    eixo_passo: float = 50.0       # m; de quanto em quanto se procura o fundo
    eixo_janela: float = 80.0      # m; meia-janela de busca lateral
    eixo_res: float = 5.0          # m; passo da amostragem dentro da janela
    eixo_desloc_max: float = 120.0  # m; teto do deslocamento
    eixo_alisar: float = 250.0     # m; janela da media movel do deslocamento
    eixo_taper: float = 300.0      # m; afunilamento junto as pontas
    eixo_penalidade: float = 0.02  # m de cota por m de afastamento do centro
    # RECORTE PELA COTA DE CHEIA. A largura saia so do porte do rio
    # (180*sqrt(A/100), com piso de 500 m de meia-largura), e o piso e que
    # mandava: 129 das 148 secoes do Benedito estavam nele. Medido, a cheia de
    # pico molhava 13% da largura na mediana e 6% na cabeceira -- 62 m de agua
    # numa secao de 966 m. O resto e encosta de vale, e ela ESTRAGA: a secao
    # vira bacia fechada (40 m de profundidade mediana) e a conducao entre
    # vizinhas chegou a variar 2.809 vezes, com a vazao invertendo de sinal no
    # primeiro passo do solver.
    recortar_secao: bool = True
    folga_secao: float = 3.0       # m acima da cota de cheia de projeto
    margem_secao: float = 1.5      # multiplica a meia-largura necessaria
    meia_largura_min: float = 60.0  # piso, para nao estrangular a planicie
    razao_lados: float = 2.5       # um lado no maximo 2,5x o outro
    minimo_lado: float = 120.0

    # ------------------------------------------------------------ calha
    # Com MDT (solo exposto) a lamina d'agua NAO esta no dado: o leito real
    # esta abaixo do que o MDT mostra no espelho d'agua. Escavar passa a ser
    # correto -- diferente do MDS de 30 m, onde o dado ja continha a superficie
    # da agua e escavar contava a mesma profundidade duas vezes.
    escavar: bool = True
    canal_kh: float = 0.277        # profundidade = kh * A^eh  (Leopold)
    canal_eh: float = 0.35
    # LARGURA CALIBRADA CONTRA OBSERVACAO DE CAMPO, nao contra Leopold cru.
    # Com kw = 5,0 (o valor original, ajustado para rios norte-americanos de
    # planicie) o modelo escavava 62 m na cabeceira do Benedito e 93 m na foz.
    # Quem conhece o rio observa 20 a 60 m, mesmo nas partes planas. Com
    # kw = 2,5 sai 31 m na cabeceira e 47 m na foz -- dentro da faixa -- e o
    # Itajai-Acu fica com 103 m em Blumenau e 117 m na foz, que tambem fecha.
    #
    # Ajustar uma potencia direta aos dois extremos observados daria expoente
    # 1,08 -- largura crescendo mais que a area, o que nao e geometria
    # hidraulica (os expoentes de Leopold ficam entre 0,4 e 0,5). Ou seja: os
    # 20 a 60 m sao variacao LOCAL, de meandro e garganta, e o que se calibra e
    # o coeficiente, mantendo o expoente.
    #
    # MEDIDO: isto NAO muda a hidraulica. A lamina de base e governada pelo
    # pilot channel (0,41 m na cabeceira com kw de 5,0 ou 2,0), e a de pico e a
    # mesma 2,69 m na foz em todos os casos, porque no pico a agua ocupa a
    # secao inteira de qualquer jeito. O que muda e a FIDELIDADE: terreno
    # preservado passa de 42% para 55%, e a escavacao mediana de 11,0 para
    # 8,1 m. E correcao de representacao, nao de estabilidade.
    canal_kw: float = 2.5          # largura = kw * A^ew
    canal_ew: float = 0.40
    # ENTALHE DIMENSIONADO PELA VAZAO DE BASE, e nao por constante. Com 25 m
    # fixos, a lamina de base na cabeceira do Benedito dava 6 cm: 0,58 m3/s
    # espalhados por 25 m a 5% de declividade. E 25 m nao entalha nada ali --
    # a calha de margens plenas por Leopold, com os 75 km2 da cabeceira, tem
    # 28 m. O entalhe tinha a largura da calha inteira.
    #
    # Medido na cabeceira, com a vazao de base de 0,58 m3/s:
    #     25 m -> 0,058 m     10 m -> 0,103 m      3 m -> 0,235 m
    # Estreitar de 25 para 3 m rende o mesmo que multiplicar a vazao por
    # QUINZE -- e multiplicar a vazao por quinze significa 115 L/s/km2, que
    # nao e escoamento de base, e cheia.
    pilot_largura: float = 25.0     # TETO da largura do entalhe
    pilot_largura_min: float = 3.0  # piso
    pilot_prof_alvo: float = 0.35   # profundidade que a vazao de base deve ter
    pilot_prof: float = 1.5

    # ----------------------------------------------------------- perfil
    decl_minima: float = 1e-4
    decl_maxima: float = 0.008     # piso do teto POR RIO
    decl_teto: float = 0.05        # limite de validade de Jarrett
    escavacao_maxima: float = 12.0
    # Secao cujo TERRENO tem menos que isto de desnivel esta inteira no fundo
    # plano do vale e nao contem cheia nenhuma. Alargar so essas: alargar todas
    # cruza as cutlines uma com as outras (o Acu vai de 0 para 213 cruzamentos
    # a 2x), e o RAS ja avisa de auto-interseccao na largura atual.
    desnivel_minimo: float = 3.0   # m
    fator_alargar: float = 2.0     # quantas vezes, so nas secoes sem desnivel
    # Acima desta cota a foz NAO esta no mar, e o contorno de jusante tem de ser
    # profundidade normal em vez de mare (que vai de -0,3 a +0,9 m). Rodando um
    # rio isolado a foz dele vira a saida do modelo: no Benedito sozinho a mare
    # foi imposta a uma secao com fundo em 50 m, o HEC-RAS recusou os dados
    # antes de computar e a rodada anunciou "NENHUM PROBLEMA DETECTADO".
    cota_mare: float = 3.0         # m
    # O aparo da cabeceira nao pode engolir a primeira confluencia: sem esse
    # limite o Itajai_Norte perdeu a cabeceira acima da foz do Iraputa, a
    # juncao ficou com um trecho entrando e um saindo, e o HEC-RAS recusou a
    # geometria antes de computar.
    corte_max_fracao: float = 0.35
    secoes_acima_juncao: int = 4

    # ------------------------------------------------------------ fluxo
    evento: str = ""               # '' (sintetico) | '1983' | '2008' | ...
    barragens: bool = True
    horas: int = 192
    fracao_cabeceira: float = 0.05  # o resto da area entra como lateral
    # Escoamento de base como fracao do pico. 2% dao 7,7 L/s/km2 na
    # cabeceira do Benedito, que e o valor plausivel para a regiao; 30%
    # dariam 115 L/s/km2, que e cheia. Nao e por aqui que se resolve
    # lamina rasa -- ver pilot_largura.
    base_frac: float = 0.02

    # ------------------------------------------------------------ plano
    # LPI ligado: a rede tem trechos de serra com 6 a 10% de declividade, onde
    # o escoamento e torrencial por fisica. Sem ele o solver oscila desde o
    # aquecimento. ZTol de 2 cm porque exigir os 6 mm do padrao num modelo
    # sobre DEM e pedir precisao abaixo da que o dado tem.
    lpi: bool = True
    ztol: float = 0.02
    max_iter: int = 40
    # INTERVALO DE CALCULO. Medido no Benedito com a geometria recortada:
    # com 1MIN o numero de Courant fica em 1,07 de mediana e 3,29 no
    # maximo -- 95 das 148 secoes acima de 1. Com 15SEC cai para 0,27 de
    # mediana e 0,82 no maximo, NENHUMA acima de 1. O solver e implicito e
    # tolera Courant acima de 1, mas nao em lamina rasa sobre leito
    # inclinado, que e exatamente o regime aqui: a agua esvaziava o
    # modelo na primeira hora e a vazao invertia de sinal.
    intervalo: str = "15SEC"

    # ------------------------------------------------------- ferramentas
    # Ferramentas do proprio HEC-RAS, via ras-commander. Ligadas por padrao:
    # sao o comportamento de referencia, e desligar e que precisa de motivo.
    usar_build_xs: bool = True     # GeomCrossSection.build_cross_section
    usar_htab: bool = True         # GeomHtabUtils.calculate_optimal_xs_htab
    # DESLIGADO desde 18/08/2026, e o motivo e escala. O que sobrou aqui e o
    # fix_htab_starting_elevations, e ele nao escala: com 2.077 secoes levava
    # segundos; com 11.251 passou de 44 min, e com 11.684 travou a rodada de
    # novo -- Python de thread unica, um nucleo de 20, sem uma linha de log no
    # meio. Cresce MUITO acima do proporcional. E o passo 7 ja calcula a tabela
    # hidraulica otima de cada secao pelo GeomHtabUtils, entao esta passagem
    # corrige o que ja estava certo. Religue com usar_fixit=true se editar a
    # geometria por fora do passo 7.
    usar_fixit: bool = False       # RasFixit.fix_htab_starting_elevations
    # fix_bank_stations SEPARADO e DESLIGADO. Medido em 18/08/2026: 17,5 min e
    # 26,1 GB de .bak para nao mudar NENHUMA das 2.077 secoes -- toda linha do
    # log dizia "original: 280, interpolated: 0". Nao podia mudar: quem escreve
    # a geometria e o build_cross_section, e ele ja INSERE as estacas das
    # margens na tabela, que e exatamente a invariante que esta ferramenta
    # confere. Sao 45% do tempo da rodada inteira gastos confirmando o que o
    # passo 7 garantiu. Ligue (corrigir_margens=true) se editar a geometria por
    # fora do builder, que e quando a invariante deixa de ser garantida.
    corrigir_margens: bool = False  # RasFixit.fix_bank_stations
    usar_check: bool = True        # RasCheck.run_all
    # DESLIGADO, e tem flag PROPRIA: usar_fixit=false NAO o desliga, o que ja
    # custou uma rodada. Medido no Benedito: a geometria saiu com ZERO areas
    # inefetivas depois de rodar -- gastou o tempo para nao mudar nada. Area
    # inefetiva e remedio para secao larga que atravessa meandro; criar uma so
    # porque o solver nao converge mascara a causa em vez de trata-la.
    usar_ineffective: bool = False  # RasFixit.fix_ineffective_flow

    # --------------------------------------------------------- execucao
    ras_exe: str = RAS_EXE
    confirmar: bool = True         # perguntar antes de cada passo

    # ------------------------------------------------------------ util
    def superficie(self):
        """A fonte inclui a lamina d'agua? (MDS em vez de MDT)"""
        return self.fonte in ("copernicus",)

    def coerir(self, log=print):
        """Ajusta o que nao pode ficar inconsistente, e AVISA.

        Escolha de fonte e decisao de escavar estao acopladas por fisica, nao
        por gosto: nao da para deixar as duas soltas e esperar que quem roda se
        lembre. O que o programa muda por conta propria, ele diz.
        """
        if self.superficie() and self.escavar:
            self.escavar = False
            log("   AVISO: fonte 'copernicus' e modelo de SUPERFICIE -- inclui "
                "a lamina d'agua.")
            log("          escavar DESLIGADO: escavar em cima dela conta a "
                "profundidade duas vezes,")
            log("          o modelo parte seco e o assentamento tem de encher "
                "canal inventado.")
            log("          Para forcar: escavar=true")
        if self.fonte == "misto" and not self.escavar:
            log("   nota: fonte 'misto' tem corredor de MDT; escavar=true "
                "faz sentido no corredor.")
        return self

    def caminho(self, *p):
        return os.path.join(self.saida, *p)

    def dict(self):
        return asdict(self)

    def gravar(self, caminho=None):
        caminho = caminho or self.caminho("opcoes.json")
        os.makedirs(os.path.dirname(caminho) or ".", exist_ok=True)
        with open(caminho, "w", encoding="utf-8") as f:
            json.dump(self.dict(), f, indent=2, ensure_ascii=False)
        return caminho

    @classmethod
    def ler(cls, caminho):
        with open(caminho, encoding="utf-8") as f:
            d = json.load(f)
        validos = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in validos})

    def aplicar(self, pares):
        """Aplica 'chave=valor' vindos da linha de comando, com o tipo certo."""
        tipos = {f.name: f.type for f in fields(self)}
        for par in pares or []:
            if "=" not in par:
                raise ValueError(f"esperado chave=valor, veio {par!r}")
            k, v = par.split("=", 1)
            k = k.strip()
            if k not in tipos:
                raise ValueError(
                    f"opcao desconhecida: {k!r}. "
                    f"Veja 'python -m vale opcoes' para a lista.")
            atual = getattr(self, k)
            if isinstance(atual, bool):
                v = str(v).strip().lower() in ("1", "true", "sim", "s", "yes")
            elif isinstance(atual, int):
                v = int(v)
            elif isinstance(atual, float):
                v = float(v)
            setattr(self, k, v)
        return self
