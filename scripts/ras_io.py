# -*- coding: utf-8 -*-
"""Escrita e conferencia de arquivos de texto do HEC-RAS.

Existe por um defeito que custou varias tentativas de abrir um projeto: o
HEC-RAS grava e ESPERA CRLF nos arquivos de texto (.prj, .pNN, .uNN, .gNN,
.rasmap, .vrt). Gravar com LF produz um arquivo cujo conteudo esta correto --
caminhos certos, titulos certos -- e que o HEC-RAS simplesmente nao le: abre o
projeto com os campos Plan, Geometry e Unsteady Flow VAZIOS e reclama de
"files not found".

O modo de errar e sutil e se repete: le-se o arquivo com `.replace("\\r","")`
para facilitar o parsing por linha, e grava-se com `newline=""`, que preserva a
string como esta. O conteudo fica com LF.

Use `escrever(caminho, texto)` para qualquer arquivo de texto do HEC-RAS, e
`conferir_crlf(pasta)` para auditar uma pasta inteira antes de entregar.
"""
import os

EXT_TEXTO = (".prj", ".rasmap", ".vrt", ".b01", ".bco01")
EXT_NUM = ("p", "u", "g", "f", "s")        # .p01, .u01, .g01, .f01, .s01...


def e_texto_ras(nome):
    n = nome.lower()
    if n.endswith(EXT_TEXTO):
        return True
    raiz, ext = os.path.splitext(n)
    return (len(ext) == 4 and ext[1] in EXT_NUM and ext[2:].isdigit())


def escrever(caminho, texto, encoding="latin-1"):
    """Grava texto do HEC-RAS SEMPRE com CRLF, venha ele como vier."""
    t = texto.replace("\r\n", "\n").replace("\r", "\n")
    with open(caminho, "wb") as f:
        f.write(t.replace("\n", "\r\n").encode(encoding, errors="replace"))
    return caminho


def terminacoes(caminho):
    b = open(caminho, "rb").read()
    crlf = b.count(b"\r\n")
    return crlf, b.count(b"\n") - crlf


def conferir_crlf(pasta, corrigir=False, log=print):
    """Audita (e opcionalmente conserta) as terminacoes de uma pasta."""
    maus = []
    for raiz, _, arquivos in os.walk(pasta):
        for f in arquivos:
            if not e_texto_ras(f):
                continue
            p = os.path.join(raiz, f)
            crlf, lf = terminacoes(p)
            if lf == 0:
                continue
            maus.append((os.path.relpath(p, pasta), crlf, lf))
            if corrigir:
                b = open(p, "rb").read()
                open(p, "wb").write(
                    b.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    if maus:
        log(f"   {len(maus)} arquivo(s) com LF"
            f"{' -- corrigidos' if corrigir else ''}:")
        for n, c, l in maus:
            log(f"      {n:<44} CRLF {c:6d}  LF {l:6d}")
    else:
        log("   terminacoes: todas CRLF")
    return maus
