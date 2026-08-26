# -*- coding: utf-8 -*-
"""Restaura secoes NOMEADAS a partir de uma geometria de referencia.

    python scripts/restaurar_secoes.py taha_ai_novo/taha_ai.g01 \
        --de taha_ai_novo/taha_ai.g01.antes_do_reparo \
        --alvo "Itajai_Mirim,R1,75.00" --saida g18

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

POR QUE EXISTE

  O reparo automatico apara perfil que cruza eixo alheio -- mas ha lugares
  onde a LARGURA e decisao de modelagem, nao defeito: na foz do Itajai
  Mirim, dentro de Itajai, a planicie dos dois rios e uma so, e encolher o
  perfil ali corta a area de extravasamento da cidade (veto do usuario,
  26/08/2026). Este script devolve o bloco ORIGINAL da secao, inteiro:
  estacas, cotas, margens, Manning, HTab e cutline.

  O preco e explicito: a secao restaurada volta a disparar as mensagens do
  Validate Geometry que motivaram a apara ("XS must intersect exactly one
  Reach", quando a cutline cruza o eixo do rio vizinho). E um custo aceito,
  registrado, e que nao impede o solver: erro de VALIDACAO de mapa nao e
  erro de DADOS.

O QUE MUDA, E O QUE NAO

  So os blocos nomeados. O `Type RM Length L Ch R` fica com os COMPRIMENTOS
  ATUAIS (se uma vizinha foi removida depois do backup, o comprimento atual
  e o correto); todo o resto do bloco vem da referencia.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever            # noqa: E402


def _alvos(argv):
    out = []
    for i, a in enumerate(argv):
        if a == "--alvo":
            rio, rch, rs = argv[i + 1].split(",")
            out.append((rio.strip(), rch.strip(), float(rs)))
    return out


def blocos(caminho):
    """[(ch, rs, i_ini, i_fim, linhas)] das secoes do arquivo."""
    linhas = open(caminho, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out, ch = [], None
    ini = None
    for i, l in enumerate(linhas):
        eh_reach = l.startswith("River Reach=")
        eh_sec = l.startswith("Type RM Length L Ch R")
        if (eh_reach or eh_sec) and ini is not None:
            out[-1][3] = i
            ini = None
        if eh_reach:
            p = l.split("=", 1)[1].split(",")
            ch = (p[0].strip(), p[1].strip())
        elif eh_sec:
            rs = float(l.split("=", 1)[1].split(",")[1])
            ini = i
            out.append([ch, rs, i, None])
    if ini is not None:
        out[-1][3] = len(linhas)
    return linhas, out


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ref = argv[argv.index("--de") + 1]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g18"
    alvos = _alvos(argv)
    if not alvos:
        raise SystemExit("preciso de ao menos um --alvo Rio,Reach,RS")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    lin_e, sec_e = blocos(entrada)
    lin_r, sec_r = blocos(ref)
    print(f"entrada   : {entrada}   (intocada)")
    print(f"referencia: {ref}")
    print(f"saida     : {novo}\n")

    trocas = {}
    for rio, rch, rs in alvos:
        de = [s for s in sec_e if s[0] == (rio, rch)
              and abs(s[1] - rs) < 0.05]
        dr = [s for s in sec_r if s[0] == (rio, rch)
              and abs(s[1] - rs) < 0.05]
        if len(de) != 1 or len(dr) != 1:
            raise SystemExit(f"alvo {rio},{rch},{rs}: {len(de)} na entrada, "
                             f"{len(dr)} na referencia -- recusado")
        bloco_ref = list(lin_r[dr[0][2]:dr[0][3]])
        # comprimentos ATUAIS no cabecalho (vizinhas podem ter mudado)
        cab_atual = lin_e[de[0][2]].split("=", 1)[1].split(",")
        cab_ref = bloco_ref[0].split("=", 1)[1].split(",")
        bloco_ref[0] = ("Type RM Length L Ch R = %s,%s,%s,%s,%s"
                        % (cab_ref[0].strip(), cab_ref[1].strip(),
                           cab_atual[2].strip(), cab_atual[3].strip(),
                           cab_atual[4].strip()))
        trocas[de[0][2]] = (de[0][3], bloco_ref)
        print(f"   restaura {rio:14s} {rch:3s} RS {rs:10.2f}   "
              f"({dr[0][3]-dr[0][2]} linhas no lugar de "
              f"{de[0][3]-de[0][2]})")

    saida, i = [], 0
    while i < len(lin_e):
        if i in trocas:
            fim, bloco = trocas[i]
            saida += bloco
            i = fim
            continue
        saida.append(lin_e[i])
        i += 1
    escrever(novo, "\n".join(saida))

    # -------------------------------------------------------- conferencia
    from qc_secoes import ler_secoes
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    R2 = ler_secoes(ref)
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    print(f"   secoes: {len(A2)} -> {len(B2)}   (nao pode mudar)")
    Rm = {(d['rio'], d['reach'], round(d['rs'], 2)): d for d in R2}
    for rio, rch, rs in alvos:
        b = next(d for d in B2 if d['rio'] == rio and d['reach'] == rch
                 and abs(d['rs'] - rs) < 0.05)
        r = Rm[(rio, rch, round(rs, 2))]
        larg_b = b['sta'][-1] - b['sta'][0]
        larg_r = r['sta'][-1] - r['sta'][0]
        ok = "ok" if abs(larg_b - larg_r) < 0.01 else "DIFERE"
        print(f"   {rio} RS {rs:.2f}: largura {larg_b:.1f} m "
              f"(referencia {larg_r:.1f})  {ok}")


if __name__ == "__main__":
    main(sys.argv[1:])
