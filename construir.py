# -*- coding: utf-8 -*-
"""
Constroi o modelo da bacia do Itajai, do relevo ao projeto HEC-RAS.

Reescrita limpa. O gerador anterior tinha 1.376 linhas com quinze correcoes
empilhadas que interagiam entre si -- o resultado oscilava entre 1 e 48 dos 192
passos conforme se mexia num parametro, sem direcao. Aqui cada etapa e um
modulo pequeno, verificavel sozinho:

    itajai/terreno.py     relevo Copernicus: UTM, amostragem, terreno .hdf
    itajai/tracado.py     eixo dos rios pelo relevo (priority-flood + D8)
    itajai/topologia.py   quem desagua em quem, e area de drenagem (ANA)

O que veio do trabalho anterior, porque esta VALIDADO:
  - series e #Sta/Elev em colunas de 8 caracteres, 10 por linha;
  - Boundary Location em 6 campos com padding;
  - juncao por nome, com Junc L&A medido da geometria;
  - Viewing Rectangle real (sem ele o RAS Mapper abre vazio);
  - rugosidade de Jarrett (1984) nas gargantas;
  - espacamento adaptado a declividade.

O que foi deixado para tras: a mistura de fontes de terreno, e as chaves
vestigiais (LATERAIS, AREA_CABECEIRA_ACU, CORTAR_NO_REENCONTRO, canais
enxertados) que ja nao faziam nada ou faziam duas coisas ao mesmo tempo.

Uso:  python construir.py            traca eixos e monta a rede
      python construir.py --terreno  refaz tambem o terreno do RAS Mapper
"""
import sys
import time

import numpy as np
from shapely.geometry import Point

from itajai import terreno, tracado, topologia

WKT = ('PROJCS["SIRGAS 2000 / UTM zone 22S",GEOGCS["SIRGAS 2000",'
       'DATUM["Sistema_de_Referencia_Geocentrico_para_las_Americas_2000",'
       'SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],'
       'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
       'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],'
       'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
       'PARAMETER["false_northing",10000000],UNIT["metre",1]]')


def main():
    t0 = time.time()
    print("=" * 70)
    print("BACIA DO ITAJAI  |  construcao a partir do relevo")
    print("=" * 70)

    print("\n[1] topologia (ANA BHO 2017)")
    rede = topologia.carregar()
    receptor, filhos = topologia.arvore(rede)
    for k, v in sorted(rede.items(), key=lambda x: -x[1]["area"]):
        alvo = receptor.get(k)
        print(f"    {v['nome']:<14} {v['area']:8.1f} km2"
              + (f"  -> {rede[alvo]['nome']}" if alvo else "   (calha principal)"))

    print("\n[2] relevo Copernicus -> UTM 22S")
    caminho = terreno.preparar_utm()
    print(f"    {caminho}")

    print("\n[3] eixo dos rios pelo relevo (priority-flood + D8)")
    tr = tracado.Tracador()
    for k, v in sorted(rede.items(), key=lambda x: -x[1]["area"]):
        ln = tr.eixo(v["cabeceira"], v["foz"])
        if ln is None:
            print(f"    {v['nome']:<14} FALHOU")
            continue
        v["linha"] = ln
        d = np.array([ln.distance(Point(p))
                      for p in list(v["linha_ana"].coords)[::20]])
        print(f"    {v['nome']:<14} {ln.length/1000:7.1f} km "
              f"(ANA {v['linha_ana'].length/1000:6.1f})   "
              f"afastamento mediano {np.median(d):5.0f} m")

    if "--terreno" in sys.argv:
        print("\n[4] terreno do RAS Mapper (RasProcess CreateTerrain)")
        print(f"    {terreno.preparar_hdf(WKT)}")

    print(f"\nconcluido em {time.time()-t0:.0f} s")
    return rede


if __name__ == "__main__":
    main()
