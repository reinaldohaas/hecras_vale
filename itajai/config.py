# -*- coding: utf-8 -*-
"""
Configuracao compartilhada do modelo.

Um lugar so para o que mais de um modulo precisa saber. O gerador antigo tinha
essas constantes espalhadas e duplicadas -- o WKT da projecao aparecia em tres
arquivos, e bastava um ficar para tras numa edicao para o RAS Mapper abrir
desalinhado.
"""

# Nome dos projetos HEC-RAS gerados por esta reescrita. Os projetos antigos
# chamavam-se Itajai_* e estao em legado/; o prefixo diferente deixa claro, no
# disco e no RAS Mapper, o que veio de qual geracao do modelo.
PROJETO = "Tajai"

EPSG = 31982                   # SIRGAS 2000 / UTM 22S

WKT = ('PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",'
       'DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",'
       'SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],'
       'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
       'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],'
       'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
       'PARAMETER["false_northing",10000000],UNIT["metre",1]]')


def nome_projeto(evento=None):
    """Tajai para o cenario sintetico, Tajai_1983 para um evento real."""
    return f"{PROJETO}_{evento}" if evento else PROJETO
