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

WKT = ('PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",'
       'DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",'
       'SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],'
       'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
       'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],'
       'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
       'PARAMETER["false_northing",10000000],UNIT["metre",1]]')

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

    # ----------------------------------------------------------- secoes
    # O Acu cai 195 m em 13 km na garganta do Salto Pilao. A 1 km de
    # espacamento sao 8 m de queda ENTRE SECOES VIZINHAS e o solver falha no
    # primeiro passo; o criterio dx <~ 0,15*D/S da ~75 m ali. Por isso o
    # espacamento sai da declividade local, e nao de um numero fixo.
    espacamento: float = 1000.0
    espacamento_min: float = 150.0
    decl_plano: float = 0.0010
    decl_ingreme: float = 0.0060
    n_pontos: int = 280            # limite do HEC-RAS e 450
    # Numa curva de raio R as perpendiculares vizinhas se encontram a R do eixo
    # do lado concavo. Passar disso e o que cruza as cutlines -- 24% dos pares
    # na primeira tentativa -- e a mancha do RAS Mapper perde sentido.
    folga_curva: float = 0.70
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
    canal_kw: float = 5.0          # largura = kw * A^ew
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
    intervalo: str = "1MIN"

    # ------------------------------------------------------- ferramentas
    # Ferramentas do proprio HEC-RAS, via ras-commander. Ligadas por padrao:
    # sao o comportamento de referencia, e desligar e que precisa de motivo.
    usar_build_xs: bool = True     # GeomCrossSection.build_cross_section
    usar_htab: bool = True         # GeomHtabUtils.calculate_optimal_xs_htab
    usar_fixit: bool = True        # RasFixit.fix_*
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
    usar_ineffective: bool = True  # RasFixit.fix_ineffective_flow

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
