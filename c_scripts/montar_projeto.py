# -*- coding: utf-8 -*-
"""Monta um projeto HEC-RAS avulso a partir de uma geometria, e confere.

    python scripts/montar_projeto.py modelo/so_mirim.g08 --nome mirim_g08 \
        --desc "Itajai-Mirim -- MDT 1 m + vaos longos + margens na grade"

    python scripts/montar_projeto.py modelo/so_mirim.g08 --nome g08_teste \
        --sem-renomear        # mantem os nomes de arquivo do projeto de origem

Copia plano, fluxo, terreno e (se houver) os resultados de uma rodada, e grava
um projeto que ABRE. Escrito depois de quatro tentativas falhas de abrir um
projeto montado a mao, cada uma com uma causa diferente. As quatro estao
tratadas aqui:

  1. TERMINACAO DE LINHA. O HEC-RAS espera CRLF. Gravar LF produz arquivo com
     conteudo correto que ele nao le: os campos Plan, Geometry e Unsteady Flow
     abrem VAZIOS e ele reclama de "files not found". Tudo passa por
     `ras_io.escrever`, e no fim a pasta e auditada.

  2. O NOME DO TERRENO NAO SE RENOMEIA. O `<terreno>.hdf` guarda o nome do
     raster como GRUPO INTERNO (`/Terrain/<nome>.Terreno_Copernicus/0..5`).
     Renomear o arquivo por fora nao renomeia o grupo, e o RAS Mapper procura
     o raster pelo nome que esta la dentro. O terreno vem com o nome de origem.

  3. O FORMATO DO .prj E O DA VERSAO INSTALADA. Copiar o formato de outro
     projeto qualquer nao serve: ha dialetos (`X Axis Title(PF)/(XS)` contra
     `(PR)/(CS)`, com ou sem `RASMap Filename=`). O .prj do projeto de ORIGEM
     e usado como molde, trocando so titulo e descricao.

  4. TITULO CURTO. `Geom Title`, `Plan Title`, `Short Identifier` e
     `Flow Title` recebem o nome do projeto e nada mais. Titulos longos
     nascem sozinhos quando cada etapa da cadeia acrescenta um sufixo -- os
     daqui chegaram a 97 caracteres. A linhagem vai para a descricao do .prj.

Alem disso, remove do .rasmap as camadas cujo arquivo nao existe (menos os
`.hdf`, que o proprio RAS gera ao abrir), e confere no fim que todo caminho
resolve.

`--terreno <hdf>` embarca um terreno ADICIONAL e o poe em PRIMEIRO lugar na
lista `<Terrains>`. A ordem nao e cosmetica: o RAS Mapper desenha o perfil do
terreno no editor de secao a partir do primeiro terreno da lista, entao quem
vem primeiro e o que se ve por tras da Station-Elevation. O terreno de origem
continua embarcado logo abaixo, e pode ser ligado no mesmo painel.
"""
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever, conferir_crlf                    # noqa: E402

GERADOS = (".p01.hdf", ".u01.hdf", ".g01.hdf", ".O01", ".O02", ".r01", ".x01",
           ".bco01", ".ic.o01", ".dss", ".b01")


def _titulo(caminho, chave, valor):
    t = open(caminho, encoding="latin-1", errors="replace").read()
    novo = re.sub(r"(?m)^" + chave + r"=.*$", chave + "=" + valor, t, count=1)
    if novo != t:
        escrever(caminho, novo)
        return True
    return False


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    geom = argv[0]
    raiz = os.path.dirname(geom) or "."
    origem = os.path.basename(geom).split(".")[0]          # so_mirim
    nome = argv[argv.index("--nome") + 1] if "--nome" in argv else origem + "_novo"
    desc = argv[argv.index("--desc") + 1] if "--desc" in argv else nome
    manter = "--sem-renomear" in argv
    res = argv[argv.index("--resultados") + 1] if "--resultados" in argv else None
    extra = argv[argv.index("--terreno") + 1] if "--terreno" in argv else None
    base = origem if manter else nome
    destino = os.path.join(raiz, nome)

    print(f"geometria : {geom}")
    print(f"projeto   : {destino}   (arquivos como '{base}.*')")
    if os.path.isdir(destino):
        shutil.rmtree(destino)
    os.makedirs(os.path.join(destino, "Terrain"))

    # ---- geometria, plano e fluxo
    shutil.copy2(geom, os.path.join(destino, base + ".g01"))
    for ext in (".p01", ".u01"):
        o = os.path.join(raiz, origem + ext)
        if os.path.exists(o):
            shutil.copy2(o, os.path.join(destino, base + ext))

    # ---- resultados de uma rodada, se houver
    n_res = 0
    if res and os.path.isdir(res):
        for f in os.listdir(res):
            if not f.startswith(origem):
                continue
            if not any(f.endswith(e) for e in GERADOS):
                continue
            shutil.copy2(os.path.join(res, f),
                         os.path.join(destino, base + f[len(origem):]))
            n_res += 1
    print(f"resultados copiados: {n_res}")

    # ---- terreno: NOME DE ORIGEM, sempre (armadilha 2)
    T = os.path.join(raiz, "Terrain")
    n_t = 0
    if os.path.isdir(T):
        for f in os.listdir(T):
            if f.startswith(origem + "_Terreno"):
                shutil.copy2(os.path.join(T, f),
                             os.path.join(destino, "Terrain", f))
                n_t += 1
    # ---- terreno adicional, COM O NOME QUE ELE JA TEM (armadilha 2 de novo:
    # o .hdf guarda o nome do raster como grupo interno, e renomear por fora
    # deixa o RAS Mapper procurando um raster que nao existe)
    extra_nome = None
    if extra:
        po = os.path.dirname(extra) or "."
        fo = os.path.basename(extra)
        extra_nome = fo[:-4] if fo.lower().endswith(".hdf") else fo
        n_e = 0
        for f in os.listdir(po):
            if f.startswith(extra_nome):
                shutil.copy2(os.path.join(po, f),
                             os.path.join(destino, "Terrain", f))
                n_e += 1
        if not n_e:
            raise SystemExit(f"nao achei arquivo algum de '{extra_nome}' em {po}")
        print(f"terreno extra: {n_e} arquivo(s) de '{extra_nome}' -- ira em 1o lugar")

    for f in ("SIRGAS2000_UTM22S.prj",):
        o = os.path.join(raiz, f)
        if os.path.exists(o):
            shutil.copy2(o, os.path.join(destino, f))
    print(f"terreno: {n_t} arquivo(s), com o nome de origem '{origem}_Terreno.*'")

    # ---- .prj: molde do projeto de origem (armadilha 3)
    molde = os.path.join(raiz, origem + ".prj")
    if not os.path.exists(molde):
        raise SystemExit(f"falta o molde {molde}")
    t = open(molde, encoding="latin-1", errors="replace").read()
    t = re.sub(r"(?m)^Proj Title=.*$", "Proj Title=" + base, t, count=1)
    m = re.search(r"BEGIN DESCRIPTION:\r?\n(.*?)\r?\nEND DESCRIPTION:", t, re.S)
    if m:
        t = t.replace(m.group(1), desc, 1)
    escrever(os.path.join(destino, base + ".prj"), t)

    # ---- rasmap
    rm = os.path.join(raiz, origem + ".rasmap")
    if os.path.exists(rm):
        t = open(rm, encoding="latin-1", errors="replace").read()
        if base != origem:
            t = re.sub(r"\b" + re.escape(origem) + r"(?!_Terreno)", base, t)
        escrever(os.path.join(destino, base + ".rasmap"), t)

    # ---- registra o terreno extra em PRIMEIRO lugar em <Terrains>
    p = os.path.join(destino, base + ".rasmap")
    if extra_nome and os.path.exists(p):
        t = open(p, encoding="latin-1", errors="replace").read()
        # NAO DECLARAR DUAS VEZES. Depois de uma rodada o proprio RAS grava o
        # terreno no .rasmap do projeto de origem; inserir de novo produz duas
        # camadas com o mesmo nome, e o passo "Computing Stored Results Maps"
        # morre com System.ArgumentException "Ja foi adicionado um item com a
        # mesma chave" em StoreAllMapsCommand -- exit code -532462766.
        ja = re.search(r'<Layer\s+Name="%s"[^>]*Type="TerrainLayer"'
                       % re.escape(extra_nome), t) is not None
        bloco = "\n".join((
            r'    <Layer Name="%s" Type="TerrainLayer" '
            r'Filename=".\Terrain\%s.hdf">' % (extra_nome, extra_nome),
            "      <ResampleMethod>near</ResampleMethod>",
            '      <Surface On="True" />',
            "    </Layer>"))
        # a tag pode vir com atributos: depois de uma rodada o proprio RAS
        # reescreve o .rasmap e ela vira `<Terrains Checked="True" ...>`.
        # Casar a string nua deixava de achar o bloco.
        m = re.search(r"<Terrains\b[^>]*?/>|<Terrains\b[^>]*>", t)
        if not m:
            raise SystemExit("o .rasmap nao tem bloco <Terrains>")
        if ja:
            print(f"rasmap: '{extra_nome}' ja estava declarado -- nao repito")
        elif m.group(0).endswith("/>"):
            t = t[:m.start()] + "<Terrains>\n" + bloco + "\n  </Terrains>" \
                + t[m.end():]
        else:
            t = t[:m.end()] + "\n" + bloco + t[m.end():]
        escrever(p, t)
        if not ja:
            print("rasmap: '" + extra_nome + "' registrado como primeiro terreno")

    # ---- titulos curtos (armadilha 4)
    for ext, chave in ((".g01", "Geom Title"), (".p01", "Plan Title"),
                       (".p01", "Short Identifier"), (".u01", "Flow Title")):
        p = os.path.join(destino, base + ext)
        if os.path.exists(p):
            _titulo(p, chave, base)

    # ---- limpa camadas do rasmap sem arquivo
    p = os.path.join(destino, base + ".rasmap")
    if os.path.exists(p):
        try:
            tree = ET.parse(p); root = tree.getroot(); n = 0
            for pai in root.iter():
                for filho in list(pai):
                    v = filho.get("Filename")
                    if not v or v.endswith(".hdf"):
                        continue
                    q = os.path.normpath(os.path.join(destino,
                                                      v.lstrip("." + os.sep)))
                    if not os.path.exists(q):
                        pai.remove(filho); n += 1
            if n:
                escrever(p, ET.tostring(root, encoding="unicode"))
                print(f"rasmap: {n} camada(s) sem arquivo removida(s)")
        except ET.ParseError as e:
            print(f"rasmap: XML nao pode ser lido ({e}) -- deixado como esta")

    # ---------------------------------------------------------- conferencia
    print()
    print("CONFERENCIA")
    conferir_crlf(destino, corrigir=True)
    faltam = []
    for ext in (".prj", ".g01", ".p01", ".u01"):
        if not os.path.exists(os.path.join(destino, base + ext)):
            faltam.append(base + ext)
    p = os.path.join(destino, base + ".rasmap")
    if os.path.exists(p):
        try:
            ET.parse(p); print("   rasmap: XML valido")
        except ET.ParseError as e:
            faltam.append(f"rasmap invalido: {e}")
        t = open(p, encoding="latin-1", errors="replace").read()
        for v in sorted(set(re.findall(r'Filename="([^"]+)"', t))):
            q = os.path.normpath(os.path.join(destino, v.lstrip("." + os.sep)))
            if not os.path.exists(q) and not v.endswith(".hdf"):
                faltam.append("rasmap aponta para " + v)
    for ext, chave in ((".prj", "Proj Title"), (".g01", "Geom Title"),
                       (".p01", "Plan Title"), (".u01", "Flow Title")):
        q = os.path.join(destino, base + ext)
        if not os.path.exists(q):
            continue
        t = open(q, encoding="latin-1", errors="replace").read()
        mm = re.search(r"(?m)^" + chave + r"=(.*)$", t)
        v = mm.group(1).strip() if mm else ""
        print(f"   {chave:<18} {v}")
        if len(v) > 32:
            faltam.append(f"{chave} com {len(v)} caracteres")
    print("   FALHAS:", faltam if faltam else "nenhuma")
    return destino


if __name__ == "__main__":
    main(sys.argv[1:])
