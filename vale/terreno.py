# -*- coding: utf-8 -*-
"""
Terreno: mosaico do SIG-SC e o .hdf que o HEC-RAS le.

O SIG-SC entrega 995 tiles de MDT a 1 m, 118 GB. Nada disso cabe num terreno
unico do HEC-RAS, e a maior parte e encosta que a cheia nunca alcanca. A saida
e um mosaico de duas resolucoes:

    CORREDOR   1 m numa faixa em torno dos eixos, onde a lamina decide a mancha
    FUNDO      5 m no resto da bacia, para a secao que passar do corredor ainda
               encontrar terreno em vez de NoData

MDT E DIFERENTE DE MDS, e isso muda premissa. O Copernicus GLO-30 usado antes e
modelo de SUPERFICIE: inclui mata, ponte e -- o pior -- a LAMINA D'AGUA, que
ele grava como um plano na cota do espelho. Dai vinham os degraus de 12 m nas
secoes, o corcovo de 9 m no Itajai do Sul (uma soleira que era copa de mata) e
a impossibilidade de escavar a calha sem contar a profundidade duas vezes. O
MDT do SIG-SC e solo exposto: o leito submerso continua ausente, mas o que
aparece e terreno de verdade.

Nada e reprojetado: o SIG-SC ja vem em EPSG:31982, o mesmo CRS do modelo.
Reprojetar 118 GB para nada seria a operacao mais cara do programa.
"""
import glob
import os

import numpy as np

SIMPLIFICA_CORTE = 5.0   # m; tolerancia da linha de corte do mosaico


def tiles(pasta):
    """Todos os tiles do SIG-SC, com a extensao de cada um lida do .tfw.

    O .tfw basta para saber onde o tile esta -- seis numeros num arquivo de
    texto. Abrir os 995 GeoTIFF so para descobrir a extensao custaria minutos
    e nao acrescenta nada.
    """
    saida = []
    for tif in sorted(glob.glob(os.path.join(pasta, "*.tif"))):
        tfw = tif[:-4] + ".tfw"
        if not os.path.exists(tfw):
            continue
        v = [float(x) for x in open(tfw).read().split()]
        px, py = v[0], abs(v[3])
        x0, y1 = v[4], v[5]
        saida.append({"tif": tif, "px": px, "py": py, "x0": x0, "y1": y1})
    return saida


def _dimensoes(tif):
    import rasterio
    with rasterio.open(tif) as ds:
        return ds.width, ds.height


def extensao(pasta, amostra=8):
    """Extensao do mosaico inteiro, medindo o tamanho de poucos tiles.

    Os tiles do SIG-SC sao todos do mesmo tamanho; medir alguns e suficiente,
    e evita abrir 995 arquivos.
    """
    t = tiles(pasta)
    if not t:
        raise FileNotFoundError(f"nenhum tile .tif com .tfw em {pasta}")
    w, h = _dimensoes(t[0]["tif"])
    for d in t[:amostra]:
        w2, h2 = _dimensoes(d["tif"])
        w, h = max(w, w2), max(h, h2)
    x0 = min(d["x0"] for d in t)
    x1 = max(d["x0"] + w * d["px"] for d in t)
    y1 = max(d["y1"] for d in t)
    y0 = min(d["y1"] - h * d["py"] for d in t)
    return (x0, y0, x1, y1), len(t), (w, h)


def tiles_que_tocam(pasta, geom, folga=0.0):
    """So os tiles que interceptam a geometria -- o resto nem e aberto."""
    from shapely.geometry import box
    alvo = geom.buffer(folga) if folga else geom
    t = tiles(pasta)
    w, h = _dimensoes(t[0]["tif"])
    dentro = []
    for d in t:
        cx = box(d["x0"], d["y1"] - h * d["py"],
                 d["x0"] + w * d["px"], d["y1"])
        if cx.intersects(alvo):
            d["caixa"] = cx
            dentro.append(d)
    return dentro


def corredor(eixos, meia_largura):
    """Poligono do corredor em torno dos eixos."""
    from shapely.ops import unary_union
    return unary_union([d["linha"].buffer(meia_largura) for d in eixos])


GDAL_BIN = os.path.join(os.path.dirname(os.__file__), "..", "Library", "bin")


def _gdal(nome):
    """Caminho do executavel do GDAL.

    Usa os EXECUTAVEIS, nao os bindings 'osgeo'. O rasterio traz o GDAL
    compilado mas nao o modulo Python, e instalar so para isto mexeria no
    ambiente de quem roda. Os .exe ja estao la, e ainda imprimem progresso.
    """
    for base in (GDAL_BIN, os.path.join(os.path.dirname(__file__), "..")):
        c = os.path.normpath(os.path.join(base, nome + ".exe"))
        if os.path.exists(c):
            return c
    return nome            # no PATH, se estiver


def _rodar(cmd, prog=None, log=print):
    """Roda o comando e converte o progresso do GDAL em tempo restante.

    O GDAL escreve '0...10...20...' sem quebra de linha; le-se caractere a
    caractere e cada '10' vira 10% na barra.
    """
    import subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    buf, erros = "", []
    for ch in iter(lambda: p.stdout.read(1), ""):
        buf += ch
        if prog is not None and buf.endswith("."):
            n = buf.rstrip(".").rsplit(".", 1)[-1]
            if n.isdigit() and 0 <= int(n) <= 100:
                prog.fracao(int(n) / 100.0)
        if ch == "\n":
            if "ERROR" in buf or "FAILURE" in buf:
                erros.append(buf.strip())
            buf = ""
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"{os.path.basename(cmd[0])} falhou "
                           f"({p.returncode}): {' | '.join(erros[-3:])}")
    return erros


def vrt(pasta_sigsc, geom, destino, log=print):
    """Mosaico VIRTUAL dos tiles que tocam 'geom'. Segundos, sem copiar nada.

    Por que isto substitui o mosaico fisico: os tiles do SIG-SC ja estao a 1 m
    e ja em EPSG:31982 -- nao ha o que reamostrar nem reprojetar. Um .vrt e um
    XML que o GDAL trata como raster unico, lendo direto dos originais.

    O mosaico fisico era pior por um motivo que so aparece medindo: com
    -crop_to_cutline o GDAL cria o raster do tamanho da CAIXA ENVOLVENTE do
    corte. O corredor e uma cobra fina de 1.670 km2, mas a caixa dele cobre
    163 x 135 km -- 21,9 BILHOES de pixels para guardar 1,67 bilhao, 92% de
    vazio, comprimido a DEFLATE. Duas horas produziram 1 GB e 10% do trabalho.

    -srcnodata 0 porque no SIG-SC zero e vazio (ver mosaico()).
    """
    dentro = tiles_que_tocam(pasta_sigsc, geom)
    if not dentro:
        raise ValueError("nenhum tile do SIG-SC toca a area pedida")
    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    lista = destino + ".tiles.txt"
    with open(lista, "w", encoding="utf-8") as f:
        f.write(chr(10).join(d["tif"] for d in dentro))
    _rodar([_gdal("gdalbuildvrt"), "-srcnodata", "0", "-vrtnodata", "-9999",
            "-input_file_list", lista, destino], None, log)
    b = sum(os.path.getsize(d["tif"]) for d in dentro)
    log(f"      VRT sobre {len(dentro)} tiles ({tamanho(b)} de origem), "
        f"{tamanho(os.path.getsize(destino))} de indice")
    return destino


def mosaico(pasta_sigsc, geom, destino, resolucao, log=print, prog=True):
    """Recorta e reamostra os tiles que tocam 'geom' num GeoTIFF unico.

    Via gdalwarp, e nao por leitura no Python: a entrada pode ter centenas de
    arquivos e dezenas de GB, e o warp resolve recorte, reamostragem e mosaico
    numa passada sem trazer nada para a memoria. Leitura por janela do rasterio
    sobre estes GeoTIFF derruba o interpretador em codigo nativo.

    ZERO E VAZIO no SIG-SC. Os tiles nao declaram nodata (o .tif.xml diz
    -3.4e38, mas essa marca nao esta no GeoTIFF) e os vazios saem como 0,00.
    Em tiles cuja altitude vai de 446 a 937 m o valor 0 ocupa 1,5 a 1,9% dos
    pixels e e CEM VEZES mais frequente que qualquer outro: e preenchimento,
    nao cota. Sem -srcnodata 0 entrariam buracos no nivel do mar dentro de um
    rio a 380 m de altitude.
    """
    import geopandas as gpd
    from .config import EPSG
    from .progresso import Progresso

    dentro = tiles_que_tocam(pasta_sigsc, geom)
    if not dentro:
        raise ValueError("nenhum tile do SIG-SC toca a area pedida")
    log(f"      {len(dentro)} tiles, reamostrando para {resolucao:g} m")

    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    corte = destino + ".corte.geojson"
    lista = destino + ".tiles.txt"
    # SIMPLIFICAR O CORTE. O gdalwarp rasteriza a mascara da cutline bloco a
    # bloco; com o corredor bruto -- uniao de 12 buffers, 36.878 vertices --
    # isso domina o custo e nao a leitura: 10% em 5 minutos num raster de 1,67
    # bilhao de pixels. A 5 m de tolerancia sobram 24,8% dos vertices e a area
    # muda 1,96 km2 em 1.670 (0,1%), muito abaixo da escala do corredor de
    # 1.000 m de meia-largura.
    simples = geom.simplify(SIMPLIFICA_CORTE)
    gpd.GeoDataFrame(geometry=[simples], crs=EPSG).to_file(corte,
                                                           driver="GeoJSON")
    with open(lista, "w", encoding="utf-8") as f:
        f.write("\n".join(d["tif"] for d in dentro))

    p = (Progresso(1000, f"mosaico {os.path.basename(destino)}", log=None)
         if prog else None)
    cmd = [_gdal("gdalwarp"),
           "-tr", str(resolucao), str(resolucao),
           # media ao reduzir: 'near' num MDT de 1 m levado a 5 m sorteia UM
           # pixel entre 25, e o ruido do laser vira o valor da celula
           "-r", "average" if resolucao > 1.0 else "near",
           "-cutline", corte, "-crop_to_cutline",
           "-srcnodata", "0", "-dstnodata", "-9999",
           "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-wm", "1024",
           "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "ZLEVEL=6",
           "-co", "BIGTIFF=YES", "-co", "NUM_THREADS=ALL_CPUS",
           "-overwrite", "--optfile", lista, destino]
    # --optfile poe a lista de entradas num arquivo: 468 caminhos estouram o
    # limite de tamanho da linha de comando do Windows
    cmd = [c for c in cmd if c]
    _rodar(cmd, p, log)
    if p:
        p.fim(tamanho(os.path.getsize(destino)))
    for f in (corte, lista):
        try:
            os.remove(f)
        except OSError:
            pass
    return destino


def piramides(caminho, niveis=(2, 4, 8, 16, 32, 64), log=print):
    """Piramides internas. Sem elas o RAS Mapper trava ao dar zoom out."""
    from .progresso import Progresso
    p = Progresso(1000, f"piramides {os.path.basename(caminho)}", log=None)
    _rodar([_gdal("gdaladdo"), "-r", "average", caminho]
           + [str(n) for n in niveis], p, log)
    p.fim()
    return caminho


def tamanho(b):
    from .progresso import tamanho as _t
    return _t(b)


def preparar_copernicus(op, log=print):
    """Usa o GeoTIFF de 30 m que ja esta no disco. Nada e reamostrado.

    E a fonte RAPIDA: um arquivo, ja em EPSG:31982, com piramides. Custa
    segundos contra as horas do SIG-SC. O preco esta na resolucao (30 m) e no
    TIPO -- modelo de superficie, com a lamina d'agua dentro do dado. Ver
    Opcoes.coerir(), que desliga a escavacao por causa disso.
    """
    if not os.path.exists(op.copernicus):
        raise FileNotFoundError(
            f"nao achei {op.copernicus}. Aponte com copernicus=CAMINHO.tif "
            f"ou use fonte=sigsc.")
    import rasterio
    with rasterio.open(op.copernicus) as ds:
        log(f"   Copernicus: {ds.width}x{ds.height} px, "
            f"{ds.res[0]:g} m, {ds.crs}, {tamanho(os.path.getsize(op.copernicus))}")
        piram = bool(ds.overviews(1))
    if not piram:
        try:
            piramides(op.copernicus, log=log)
        except Exception as e:                                # noqa: BLE001
            log(f"      piramides nao geradas ({e}); o RAS Mapper fica lento "
                f"no zoom out, mas o modelo roda")
    return [op.copernicus]


def estimativa(op, eixos, log=print):
    """Quanto vai custar o passo do terreno, ANTES de comecar.

    Conta tiles e bytes que serao LIDOS, que e onde o tempo esta. A conversao
    para tempo usa uma taxa de leitura conservadora e e apresentada como ordem
    de grandeza -- serve para decidir entre as fontes, nao para cronometrar.
    """
    from shapely.ops import unary_union
    TAXA_MB_S = 120.0            # leitura+descompressao, disco local

    if op.fonte == "copernicus":
        b = os.path.getsize(op.copernicus) if os.path.exists(op.copernicus) else 0
        log(f"   fonte 'copernicus': 1 arquivo, {tamanho(b)}  ->  segundos")
        return {"tiles": 1, "bytes": b, "seg": b / 1e6 / TAXA_MB_S}

    total_t, total_b = 0, 0
    faixa = corredor(eixos, op.corredor_m)
    areas = [("corredor", faixa, op.res_corredor)]
    if op.fonte == "sigsc" and op.fundo != "nenhum":
        if op.fundo == "bacia":
            areas.append(("fundo",
                          unary_union([d["linha"].buffer(8000.0) for d in eixos]),
                          op.res_fundo))
        else:
            from shapely.geometry import box
            (x0, y0, x1, y1), _, _ = extensao(op.sigsc)
            areas.append(("fundo", box(x0, y0, x1, y1), op.res_fundo))
    for nome, geom, res in areas:
        dentro = tiles_que_tocam(op.sigsc, geom)
        b = sum(os.path.getsize(d["tif"]) for d in dentro)
        total_t += len(dentro)
        total_b += b
        log(f"   {nome:<10} {geom.area/1e6:8.0f} km2   {len(dentro):>4} tiles   "
            f"{tamanho(b):>9}   ~{hms_(b / 1e6 / TAXA_MB_S)}")
    if op.fonte == "misto":
        log(f"   fundo      Copernicus de 30 m (nao le tile do SIG-SC)")
    seg = total_b / 1e6 / TAXA_MB_S
    if getattr(op, "vrt", False):
        # com VRT nao ha leitura antecipada: os tiles sao lidos sob demanda,
        # ponto a ponto, quando as secoes forem cortadas
        log(f"   TOTAL      {total_t} tiles indexados por VRT -- segundos "
            f"(sem VRT seriam ~{hms_(seg)} de leitura)")
    else:
        log(f"   TOTAL      {total_t} tiles, {tamanho(total_b)}, "
            f"ordem de ~{hms_(seg)} de leitura")
    return {"tiles": total_t, "bytes": total_b, "seg": seg}


def hms_(s):
    from .progresso import hms
    return hms(s)


def preparar(op, eixos, log=print):
    """Monta os GeoTIFF do terreno. Devolve a lista, do fino para o grosso.

    A ordem importa: o RasTerrain empilha e o PRIMEIRO tem prioridade onde
    houver sobreposicao.
    """
    from shapely.ops import unary_union
    if op.fonte == "copernicus":
        return preparar_copernicus(op, log)
    saidas = []
    faixa = corredor(eixos, op.corredor_m)

    res = float(getattr(op, "res_sigsc", 1.0))
    area_bacia = unary_union([faixa,
                              unary_union([d["linha"].buffer(8000.0)
                                           for d in eixos])])
    if res > 1.5:
        # MOSAICO FISICO quando ha reamostragem. A 10 m a bacia inteira cabe
        # num arquivo pequeno que o amostrador carrega em memoria de uma vez:
        # some a leitura por janela, e com ela a divisao recursiva que a caixa
        # envolvente do rio (40.575 x 100.177 px, 30 GiB) obrigou a existir.
        # A 1 m isso seria impossivel -- dai o VRT continuar sendo o caminho
        # para a resolucao nativa.
        fino = op.caminho("Terrain", f"{op.projeto}_sigsc_{res:g}m.tif")
        log(f"   mosaico do SIG-SC a {res:g} m sobre {area_bacia.area/1e6:.0f} km2")
        mosaico(op.sigsc, area_bacia, fino, res, log)
        piramides(fino, log=log)
        saidas.append(fino)
        if os.path.exists(op.copernicus):
            log("   Copernicus ao final, so para preencher vazio do MDT")
            saidas += preparar_copernicus(op, log)
        return saidas

    if op.vrt:
        # Um VRT so, sobre tudo que a bacia toca. Nao ha corredor e fundo
        # separados: o VRT le os tiles no 1 m nativo em qualquer lugar que eles
        # existam, e custa segundos em vez de horas.
        area = area_bacia
        v = op.caminho("Terrain", f"{op.projeto}_sigsc.vrt")
        log(f"   VRT do SIG-SC sobre {area.area/1e6:.0f} km2 (1 m nativo)")
        vrt(op.sigsc, area, v, log)
        saidas.append(v)
        if os.path.exists(op.copernicus):
            log("   Copernicus ao final, so para preencher vazio do MDT")
            saidas += preparar_copernicus(op, log)
        return saidas

    fino = op.caminho("Terrain", f"{op.projeto}_corredor_{op.res_corredor:g}m.tif")
    log(f"   corredor de {op.corredor_m:g} m a {op.res_corredor:g} m "
        f"({faixa.area/1e6:.0f} km2)")
    mosaico(op.sigsc, faixa, fino, op.res_corredor, log)
    piramides(fino, log=log)
    saidas.append(fino)

    if op.fonte == "misto":
        # o fundo vem do Copernicus: corta a passada mais cara, que e ler tile
        # do SIG-SC sobre a bacia inteira so para produzir 5 m
        saidas += preparar_copernicus(op, log)
    elif op.fundo != "nenhum":
        if op.fundo == "bacia":
            # a bacia, aproximada pela envoltoria dos eixos com folga larga.
            # Nao e o divisor de aguas real; e o suficiente para a secao nunca
            # cair em NoData, que e para o que o fundo existe.
            area = unary_union([d["linha"].buffer(8000.0) for d in eixos])
        else:
            from shapely.geometry import box
            (x0, y0, x1, y1), _, _ = extensao(op.sigsc)
            area = box(x0, y0, x1, y1)
        grosso = op.caminho("Terrain", f"{op.projeto}_fundo_{op.res_fundo:g}m.tif")
        log(f"   fundo '{op.fundo}' a {op.res_fundo:g} m ({area.area/1e6:.0f} km2)")
        mosaico(op.sigsc, area, grosso, op.res_fundo, log)
        piramides(grosso, log=log)
        saidas.append(grosso)

    # COPERNICUS POR ULTIMO, como tapa-buraco. O SIG-SC tem 1,5 a 1,9% de
    # vazios por tile, e agora eles sao NoData de verdade em vez de zeros
    # falsos -- o que deixaria furos na secao. Como a lista e consultada em
    # ordem (o primeiro que tiver dado vence, aqui e no RasTerrain), por na
    # ultima posicao faz o Copernicus preencher SO onde o MDT nao tem nada.
    # O erro do MDS (25 m de copa de mata sobre o leito) fica restrito aos
    # vazios, em vez de valer para o modelo inteiro.
    if op.fonte == "sigsc" and os.path.exists(op.copernicus):
        log("   Copernicus ao final, so para preencher vazio do MDT")
        saidas += preparar_copernicus(op, log)
    return saidas


def hdf(op, geotiffs, log=print):
    """Terreno do HEC-RAS, pelo RasProcess (ferramenta do proprio RAS).

    O RasTerrain empilha os GeoTIFF na ordem recebida: o PRIMEIRO tem
    prioridade onde houver sobreposicao. Por isso o corredor de 1 m vem antes
    do fundo de 5 m.

    Exige um .prj ESRI de verdade. Passar um .projection faz o RasProcess
    chegar a PROGRESS=100 e falhar com "Referencia de objeto nao definida" --
    erro do .NET que nao diz nada sobre a causa.
    """
    from ras_commander import RasTerrain
    from .config import WKT

    pasta = op.caminho("Terrain")
    os.makedirs(pasta, exist_ok=True)     # com fonte=copernicus a pasta pode
    prj = os.path.join(pasta, f"{op.projeto}.prj")   # nem existir ainda
    with open(prj, "w", encoding="ascii") as f:
        f.write(WKT)
    destino = os.path.join(pasta, f"{op.projeto}_Terreno.hdf")
    log(f"   RasTerrain.create_terrain_hdf -> {os.path.basename(destino)}")
    # UNITS="Meters", explicito. O padrao da funcao e "Feet", e num modelo em
    # metros isso daria um terreno 3,28 vezes errado -- sem erro nenhum, so
    # cotas absurdas. Os nomes dos argumentos tambem sao estes, e nao os que
    # pareciam obvios (input_files/output_path/projection_file).
    RasTerrain.create_terrain_hdf(
        input_rasters=list(geotiffs), output_hdf=destino,
        projection_prj=prj, units="Meters",
        timeout_seconds=int(getattr(op, "terreno_timeout", 7200)))
    return destino


class Amostrador:
    """Le cota do terreno, com o primeiro raster tendo prioridade.

    NAO usa rasterio.sample(). Sobre estes GeoTIFF o sample() do rasterio 1.5
    derruba o interpretador em codigo NATIVO -- sem excecao, sem traceback, o
    processo simplesmente termina (exit 127). Leitura por janela do rasterio
    faz o mesmo. Foi assim que o passo das secoes morreu calado.

    Em vez disso: banda inteira em memoria quando cabe (indexacao por
    transformada, que e so aritmetica), e leitura por janela pelo GDAL quando
    nao cabe -- o GDAL nao tem esse defeito.
    """

    LIMITE_MEM = 3_000_000_000        # bytes; acima disto, janela pelo GDAL

    def __init__(self, caminhos, log=print):
        import rasterio
        self.fontes = []
        for c in caminhos:
            ds = rasterio.open(c)
            n = ds.width * ds.height * 4
            f = {"ds": ds, "tr": ~ds.transform, "w": ds.width, "h": ds.height,
                 "nodata": ds.nodata, "banda": None, "caminho": c}
            if n <= self.LIMITE_MEM:
                f["banda"] = ds.read(1)
                log(f"      {os.path.basename(c)}: banda em memoria "
                    f"({tamanho(n)})")
            else:
                from osgeo import gdal
                gdal.UseExceptions()
                f["gdal"] = gdal.Open(c)
                log(f"      {os.path.basename(c)}: {tamanho(n)}, leitura por "
                    f"janela pelo GDAL")
            self.fontes.append(f)

    MAX_JANELA = 16_000_000        # pixels por leitura (~64 MB em float32)

    def _ler(self, f, xs, ys):
        col, lin = f["tr"] * (xs, ys)
        c = np.floor(np.asarray(col, float)).astype(np.int64)
        l = np.floor(np.asarray(lin, float)).astype(np.int64)
        dentro = (c >= 0) & (c < f["w"]) & (l >= 0) & (l < f["h"])
        out = np.full(len(c), np.nan)
        if not dentro.any():
            return out
        if f["banda"] is not None:
            out[dentro] = f["banda"][l[dentro], c[dentro]]
        else:
            idx = np.flatnonzero(dentro)
            self._janela(f, c, l, idx, out)
        nod = f["nodata"]
        if nod is not None:
            out[out == nod] = np.nan
        out[out < -1000.0] = np.nan
        return out

    def _janela(self, f, c, l, idx, out, prof=0):
        """Le uma janela do raster, DIVIDINDO quando a caixa e grande demais.

        A caixa envolvente de um corte transversal e estreita -- os pontos
        estao numa reta de algumas centenas de metros. Mas o mesmo amostrador
        e chamado para percorrer o EIXO INTEIRO do rio, e ai a caixa e a do
        rio: 40.575 x 100.177 pixels, 30 GiB, para colher alguns milhares de
        valores. Foi assim que o passo das secoes morreu.

        Dividir recursivamente pelo lado maior resolve os dois padroes sem
        precisar saber qual deles esta em uso.
        """
        if len(idx) == 0:
            return
        cd, ld = c[idx], l[idx]
        c0, c1 = int(cd.min()), int(cd.max())
        l0, l1 = int(ld.min()), int(ld.max())
        w, h = c1 - c0 + 1, l1 - l0 + 1
        if w * h > self.MAX_JANELA and len(idx) > 1 and prof < 40:
            eixo = cd if w >= h else ld
            meio = (eixo.min() + eixo.max()) // 2
            esq = idx[eixo <= meio]
            dir_ = idx[eixo > meio]
            if len(esq) == 0 or len(dir_) == 0:      # tudo num lado: parte no
                meio_i = len(idx) // 2               # numero de pontos
                esq, dir_ = idx[:meio_i], idx[meio_i:]
            self._janela(f, c, l, esq, out, prof + 1)
            self._janela(f, c, l, dir_, out, prof + 1)
            return
        b = f["gdal"].GetRasterBand(1)
        arr = b.ReadAsArray(c0, l0, w, h)
        # float32, nao float64: converter dobra a memoria da janela sem ganho
        out[idx] = np.asarray(arr, dtype=np.float32)[ld - l0, cd - c0]

    def cota(self, xs, ys):
        xs = np.asarray(xs, float).ravel()
        ys = np.asarray(ys, float).ravel()
        out = np.full(xs.shape, np.nan)
        for f in self.fontes:                 # o primeiro que tiver dado vence
            falta = ~np.isfinite(out)
            if not falta.any():
                break
            out[falta] = self._ler(f, xs[falta], ys[falta])
        return out

    def fechar(self):
        for f in self.fontes:
            f["ds"].close()
