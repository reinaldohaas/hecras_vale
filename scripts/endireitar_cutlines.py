# -*- coding: utf-8 -*-
"""Poe as cutlines gravadas ao contrario no sentido esquerda->direita.

    python scripts/endireitar_cutlines.py modelo/mirim_t30/mirim_t30.g02 --saida g03

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  O HEC-RAS exige que a linha de corte seja percorrida da margem ESQUERDA
  para a DIREITA, olhando para JUSANTE. Medido no `g01`, 45 das 1.418 secoes
  (3,2%) estao gravadas ao contrario -- o produto vetorial entre a tangente do
  rio e a direcao da cutline sai com o sinal trocado.

  A consequencia se ve no RAS Mapper: a bank line liga o ponto `Bank Sta`
  esquerdo de uma secao ao da seguinte, e onde a orientacao inverte esse ponto
  cai do OUTRO lado do rio. A linha atravessa o canal e volta -- 94 travessias
  ao longo dos 141 km. E o zigue-zague vermelho da imagem.

  E NAO E SO DESENHO. Com a cutline invertida, o `#Mann` da margem esquerda e
  aplicado a planicie da DIREITA, e vice-versa. Neste modelo as duas planicies
  usam n = 0,162 e o efeito hidraulico e nulo -- mas a geometria fica dizendo
  o contrario do que representa, e qualquer valor distinto por margem, agora
  ou depois, sairia trocado.

A INVERSAO E EXATA, E NAO APROXIMADA

  Nenhuma cota muda e nenhum ponto entra ou sai. A secao inteira e espelhada:

      estaca      s  ->  L - s     (e a lista volta a ficar crescente)
      cota        z  ->  z         (so muda de ordem)
      Bank Sta   (lb, rb) -> (L - rb, L - lb)
      #Mann      as faixas invertem de ordem, e o inicio de cada uma vira
                 L menos o FIM da faixa correspondente
      cutline    troca os dois extremos

  O talvegue, a largura do canal e a area da secao sao os mesmos por
  construcao -- a conferencia no fim mede isso e falha se nao for.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from ras_io import escrever            # noqa: E402


def _col(v, larg, dec):
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += ("%*.*f" % (larg, dec, x))
        if (i + 1) % 10 == 0:
            saida.append(linha)
            linha = ""
    if linha:
        saida.append(linha)
    return saida


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def invertidas(S, eixo):
    """Indices das secoes cuja cutline nao vai de esquerda para direita."""
    from shapely.geometry import Point
    fora = []
    for i, d in enumerate(S):
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        v = B - A
        u = v / max(float(np.hypot(*v)), 1e-9)
        M = 0.5 * (A + B)
        s = eixo.project(Point(*M))
        q0 = np.array(eixo.interpolate(max(s - 25.0, 0.0)).coords[0])
        q1 = np.array(eixo.interpolate(min(s + 25.0, eixo.length)).coords[0])
        t = q1 - q0
        t /= max(float(np.hypot(*t)), 1e-9)
        if t[0] * u[1] - t[1] * u[0] > 0:      # a esquerda virou direita
            fora.append(i)
    return fora


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g03"
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    eixo = list(ler_eixos(entrada).values())[0]
    # A LISTA VEM DO PROPRIO HEC-RAS quando ha `--erros`.
    # O detector daqui (`invertidas`) usa a tangente do eixo numa janela de
    # +-25 m e erra onde a cutline e quase paralela ao fluxo: apontou 42
    # secoes onde o Validate Geometry aponta 16, e o excedente PIOROU as bank
    # lines. Preferir o veredito do proprio motor elimina o palpite.
    erros = _arg(argv, "--erros")
    if erros:
        import csv
        import re as _re
        # DUAS ARMADILHAS NO CASAMENTO, ambas medidas:
        #  - o RAS EXIBE o RS arredondado (114599.4 para 114599.36), entao
        #    igualdade exata perde 6 das 16;
        #  - o RS SE REPETE entre reaches (R2 vai de 20057 a 1221 e o canal de
        #    7500 a 0), entao casar so por RS pega a secao errada.
        # Casa-se por (rio, reach) e pela tolerancia do proprio arredondamento.
        cab = open(entrada, encoding="latin-1", errors="replace").read() \
            .replace("\r", "").split("\n")
        ch, mapa = None, []
        for l in cab:
            if l.startswith("River Reach="):
                p = l.split("=", 1)[1].split(",")
                ch = (p[0].strip(), p[1].strip())
            elif l.startswith("Type RM Length L Ch R"):
                mapa.append(ch)
        marcados, perdidos = [], []
        for r in csv.DictReader(open(erros, encoding="utf-8"), delimiter=";"):
            if "reversed" not in r["mensagem"]:
                continue
            m = _re.match(r"\s*([^,]+),\s*(\S+)\s*\(([\d.]+)\)",
                          r["onde"].strip())
            if not m:
                continue
            rio_e, rch_e, rs_e = m.group(1).strip(), m.group(2), float(m.group(3))
            # O RAS TRUNCA, nao arredonda: 111072.16 aparece como "111072.1",
            # a 0,06 do valor real -- fora de uma tolerancia de meia casa.
            # A janela cobre os dois comportamentos: [-meia casa, +uma casa).
            # Sem ambiguidade, porque a secao mais proxima esta a 185 m.
            dec = len(m.group(3).split(".")[-1]) if "." in m.group(3) else 0
            passo = 10 ** (-dec)
            cand = [i for i, d in enumerate(S)
                    if mapa[i] == (rio_e, rch_e)
                    and -0.5 * passo - 1e-9 <= float(d["rs"]) - rs_e
                    < passo + 1e-9]
            if len(cand) == 1:
                marcados.append(cand[0])
            else:
                perdidos.append((rio_e, rch_e, rs_e, len(cand)))
        alvo = set(marcados)
        print(f"lista de erros: {erros}")
        print(f"   marcadas pelo HEC-RAS: {len(marcados)+len(perdidos)}   "
              f"casadas: {len(alvo)}")
        for x in perdidos:
            print(f"   NAO CASOU: {x[0]} {x[1]} RS {x[2]}  "
                  f"({x[3]} candidatas)")
    else:
        alvo = set(invertidas(S, eixo))
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"secoes : {len(S)}   gravadas ao contrario: {len(alvo)} "
          f"({100*len(alvo)/len(S):.1f}%)\n")
    if not alvo:
        raise SystemExit("nada a endireitar")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    comp = {}                       # comprimento da secao, por indice
    for i, d in enumerate(S):
        comp[i] = float(np.asarray(d["sta"], float)[-1])

    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        if i_sec in alvo:
            d = S[i_sec]
            L = comp[i_sec]
            if l.startswith("XS GIS Cut Line"):
                A = np.asarray(d["cut"][0], float)
                B = np.asarray(d["cut"][-1], float)
                saida.append("XS GIS Cut Line= 2")
                saida.append("".join("%16.2f" % x for x in
                                     (B[0], B[1], A[0], A[1])))
                j += 1
                while j < len(linhas) and linhas[j].strip() and \
                        linhas[j][:1] in " -0123456789":
                    j += 1
                continue
            if l.startswith("#Sta/Elev"):
                st = L - np.asarray(d["sta"], float)[::-1]
                z = np.asarray(d["z"], float)[::-1]
                st[0] = 0.0                     # tira o -0.00 do espelhamento
                v = []
                for a, b in zip(st, z):
                    v += [a, b]
                saida.append("#Sta/Elev= %d " % len(st))
                saida += _col(v, 8, 2)
                cnt = int(l.split("=")[1])
                j += 1
                lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(x), 8)
                                  if x[c:c + 8].strip()])
                    j += 1
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%s,%s"
                             % (_fmt(L - float(d["rb"])),
                                _fmt(L - float(d["lb"]))))
                j += 1
                continue
            if l.startswith("#Mann="):
                cnt = int(l.split("=")[1].split(",")[0])
                bruto, k2 = [], j + 1
                while k2 < len(linhas) and len(bruto) < 3 * cnt:
                    x = linhas[k2]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                              if x[c:c + 8].strip()]
                    k2 += 1
                val = [float(x) for x in bruto[:3 * cnt]]
                ini = [val[3 * k] for k in range(cnt)]
                nn = [val[3 * k + 1] for k in range(cnt)]
                fl = [val[3 * k + 2] for k in range(cnt)]
                # cada faixa vai de ini[k] ate ini[k+1] (a ultima ate L);
                # espelhada, ela passa a comecar em L - fim
                fim = ini[1:] + [L]
                novo_ini = [L - f for f in fim][::-1]
                novo_ini[0] = 0.0
                nn = nn[::-1]
                fl = fl[::-1]
                saida.append(l)
                lin, corpo = "", []
                seq = []
                for k in range(cnt):
                    seq += [novo_ini[k], nn[k], fl[k]]
                for t, x in enumerate(seq):
                    lin += ("%8.2f" % x if t % 3 == 0 else
                            "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                    if (t + 1) % 9 == 0:
                        corpo.append(lin)
                        lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                j = k2
                continue
        saida.append(l)
        j += 1

    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    from shapely.geometry import Point
    print("CONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print(f"   secoes                 : {len(A2)} -> {len(B2)}")
    za = np.array([float(x["z"].min()) for x in A2])
    zb = np.array([float(x["z"].min()) for x in B2])
    print(f"   talvegue mudou em      : max {np.abs(zb - za).max():.6f} m "
          "(tem de ser zero)")
    ca = np.array([float(x["rb"] - x["lb"]) for x in A2])
    cb = np.array([float(x["rb"] - x["lb"]) for x in B2])
    print(f"   largura do canal mudou : max {np.abs(cb - ca).max():.6f} m "
          "(tem de ser zero)")
    la = np.array([float(x["sta"][-1] - x["sta"][0]) for x in A2])
    lb_ = np.array([float(x["sta"][-1] - x["sta"][0]) for x in B2])
    print(f"   largura da secao mudou : max {np.abs(lb_ - la).max():.6f} m "
          "(tem de ser zero)")
    aa = np.array([float(np.trapezoid(x["z"], x["sta"])) for x in A2])
    ab = np.array([float(np.trapezoid(x["z"], x["sta"])) for x in B2])
    print(f"   area sob o perfil mudou: max {np.abs(ab - aa).max():.6f} m2 "
          "(tem de ser zero)")
    # A conferencia de "quantas continuam ao contrario" e de travessia da bank
    # line so vale com um eixo; com quatro reaches ela mediria a coisa errada.
    # Quem responde isso agora e o proprio HEC-RAS, por `ler_erros_geometria`.
    print("   rode `ler_erros_geometria.py` sobre o .hdf para o veredito")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
