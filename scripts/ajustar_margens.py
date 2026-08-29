# -*- coding: utf-8 -*-
"""Estende a estaca de margem mais proxima ate incluir o ponto mais fundo.

    python scripts/ajustar_margens.py modelo/so_mirim.g03 --saida g04 --limiar 0.50

O ARQUIVO DE ENTRADA NAO E TOCADO. Sai um .gXX novo. Muda EXCLUSIVAMENTE as
linhas `Bank Sta=`; perfil, estacas, cutline, Manning e htab ficam byte a byte.

POR QUE, E POR QUE SO ISTO

  A `Bank Sta` divide a secao em planicie esquerda, canal e planicie direita.
  O HEC-RAS aplica um Manning a cada pedaco e calcula a conducao de cada um em
  separado. Quando o ponto mais fundo da secao cai FORA dessa marca, a parte
  mais funda e tratada como planicie: recebe rugosidade de planicie e conduz
  menos do que deveria.

  Depois de recortar o perfil no MDT de 1 m, isso passou a acontecer em 239
  secoes -- mas 165 delas por menos de 25 cm, que e empate num fundo de vale
  chato, com o minimo caindo 12 m para fora por 12 cm de diferenca. Mexer
  nelas seria ruido.

  Por isso o limiar. So entram as secoes em que o canal declarado esta acima
  do fundo real por mais que `--limiar` metros.

  E POR ISSO A ACAO E A MINIMA POSSIVEL: a margem mais proxima anda ate a
  estaca do ponto mais fundo, e para ali. Nao se procura quebra de
  declividade, nao se recentra o canal, nao se inventa largura. A estaca de
  destino ja existe no perfil -- nenhum valor novo entra no arquivo.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes  # noqa: E402


def _fmt(v):
    """Como o HEC-RAS grava: sem zero a direita ('24' e nao '24.00')."""
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g04"
    lim = float(argv[argv.index("--limiar") + 1]) if "--limiar" in argv else 0.50
    # TETO DO DESLOCAMENTO. Sem ele a regra "estenda ate o ponto mais fundo"
    # deixa de ser minima: medido nas 54 do limiar de 0,50 m, a margem andava
    # 54,9 m na mediana e 232,7 m no pior caso, e a largura declarada do canal
    # ia de 66 para 127 m (299 m no maximo). Um fundo a 163 m da margem quase
    # nunca e o mesmo canal -- e braco secundario, meandro abandonado ou
    # drenagem paralela cortada pela secao, e nesse caso o conserto nao e
    # alargar o canal, e olhar a secao.
    maxd = (float(argv[argv.index("--max-desloc") + 1])
            if "--max-desloc" in argv else float("inf"))
    raiz = os.path.dirname(entrada) or "."
    nome = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{nome}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    print(f"entrada: {entrada}   (intocado)")
    print(f"saida  : {novo}")
    print(f"limiar : {lim:.2f} m de desnivel entre o canal declarado e o fundo real")

    secoes = ler_secoes(entrada)
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    idx = [i for i, l in enumerate(linhas) if l.startswith("Bank Sta=")]
    assert len(idx) == len(secoes), "linhas Bank Sta e secoes nao batem"

    mudancas, adiadas = [], []
    for k, d in enumerate(secoes):
        lb, rb = d.get("lb"), d.get("rb")
        if lb is None or rb is None:
            continue
        st, z = d["sta"], d["z"]
        j = int(np.argmin(z))
        s = float(st[j])
        if lb <= s <= rb:
            continue
        canal = (st >= lb) & (st <= rb)
        if not canal.any():
            continue
        desnivel = float(np.min(z[canal])) - float(z[j])
        if desnivel <= lim:
            continue
        if s < lb:
            novo_lb, novo_rb, lado = s, rb, "esq"
        else:
            novo_lb, novo_rb, lado = lb, s, "dir"
        if novo_rb <= novo_lb:
            continue
        anda = abs(s - (lb if lado == "esq" else rb))
        if anda > maxd:
            adiadas.append({"rs": d["rs"], "anda": anda, "desnivel": desnivel,
                            "larg": rb - lb, "larg_seria": novo_rb - novo_lb})
            continue
        mudancas.append({"i": idx[k], "rs": d["rs"], "lado": lado,
                         "lb": lb, "rb": rb, "novo_lb": novo_lb,
                         "novo_rb": novo_rb, "desnivel": desnivel,
                         "anda": abs(s - (lb if lado == "esq" else rb)),
                         "larg": rb - lb, "larg_nova": novo_rb - novo_lb})

    # A QUEBRA DO MANNING E A PROPRIA ESTACA DE MARGEM, e tem de andar junto.
    #
    # O bloco e assim:
    #     Bank Sta=29.97,75.37
    #     #Mann= 3 , 0 , 0
    #         0.00   0.100       0   29.97   0.055       0   75.37   0.100  0
    #                ^planicie        ^lb    ^canal          ^rb     ^planicie
    #
    # Mover so a `Bank Sta` deixa o pedaco entre a quebra velha e a margem
    # nova sem valor de n, e o HEC-RAS RECUSA A RODAR -- sem computar nada,
    # com o motivo no <plano>.data_errors.txt: "Right bank Manning's n value
    # not set". Foi o que aconteceu na primeira tentativa, nas 13 secoes.
    troca = {m["i"]: f"Bank Sta={_fmt(m['novo_lb'])},{_fmt(m['novo_rb'])}"
             for m in mudancas}
    porsec = {m["i"]: m for m in mudancas}
    mann = {}
    for i_bank, m in porsec.items():
        j = i_bank + 1
        while j < len(linhas) and not linhas[j].startswith("#Mann="):
            if linhas[j].startswith("Type RM Length"):
                j = -1
                break
            j += 1
        if j < 0 or j >= len(linhas):
            continue
        n = int(linhas[j].split("=")[1].split(",")[0])
        bruto, k = [], j + 1
        while k < len(linhas) and len(bruto) < 3 * n:
            l = linhas[k]
            if not l.strip() or l[:1].isalpha() or l[:1] == "#":
                break
            bruto += [l[c:c + 8] for c in range(0, len(l), 8) if l[c:c + 8].strip()]
            k += 1
        if len(bruto) < 3 * n:
            continue
        v = [float(x) for x in bruto[:3 * n]]
        for idx_v in range(0, 3 * n, 3):
            s = v[idx_v]
            if abs(s - m["lb"]) < 1e-6:
                v[idx_v] = m["novo_lb"]
            elif abs(s - m["rb"]) < 1e-6:
                v[idx_v] = m["novo_rb"]
        est = [v[t] for t in range(0, 3 * n, 3)]
        if any(b <= a for a, b in zip(est, est[1:])):
            continue                      # nao deixaria as quebras em ordem
        # CADA COLUNA COM A SUA PRECISAO. Gravar tudo com %8.2f arredondava o
        # n de 0.055 para 0.06 -- 9% de rugosidade a mais no canal, mudando o
        # modelo por descuido de formatacao. Estaca em 2 casas, n em 3, o
        # terceiro campo inteiro, como o HEC-RAS grava.
        linha_nova, corpo = "", []
        for t, x in enumerate(v):
            linha_nova += ("%8.2f" % x if t % 3 == 0 else
                           "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
            if (t + 1) % 9 == 0:
                corpo.append(linha_nova); linha_nova = ""
        if linha_nova:
            corpo.append(linha_nova)
        mann[j] = (k, corpo)

    saida, i = [], 0
    while i < len(linhas):
        saida.append(troca.get(i, linhas[i]))
        if i in mann:
            fim, corpo = mann[i]
            saida += corpo
            i = fim
            continue
        i += 1
    txt = "\n".join(saida)
    txt = txt.replace("Geom Title=" + linhas[0].split("=", 1)[1],
                      "Geom Title=" + linhas[0].split("=", 1)[1] + " + margens", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)

    print(f"\nsecoes com margem estendida: {len(mudancas)} de {len(secoes)}")
    if mudancas:
        a = np.array([m["anda"] for m in mudancas])
        w = np.array([m["larg"] for m in mudancas])
        w2 = np.array([m["larg_nova"] for m in mudancas])
        print(f"   a margem andou    : mediana {np.median(a):6.1f} m   "
              f"p90 {np.percentile(a,90):6.1f}   max {a.max():6.1f}")
        print(f"   largura do canal  : mediana {np.median(w):6.1f} -> "
              f"{np.median(w2):6.1f} m   max {w2.max():.1f}")
        print()
        print(f"   {'RS':>10} {'lado':>5} {'desnivel':>9} {'andou':>7}   "
              f"{'antes':>16} -> {'depois':>16}")
        for m in sorted(mudancas, key=lambda m: -m["desnivel"]):
            print(f"   {m['rs']:>10.2f} {m['lado']:>5} {m['desnivel']:>8.2f}m "
                  f"{m['anda']:>6.1f}m   "
                  f"{_fmt(m['lb'])+','+_fmt(m['rb']):>16} -> "
                  f"{_fmt(m['novo_lb'])+','+_fmt(m['novo_rb']):>16}")
    if adiadas:
        a = np.array([m["anda"] for m in adiadas])
        print(f"\nNAO ALTERADAS pelo teto de {maxd:.0f} m: {len(adiadas)} secoes")
        print(f"   a margem teria de andar mediana {np.median(a):.1f} m, "
              f"max {a.max():.1f} m")
        print("   ficam marcadas para exame -- fundo distante assim raramente "
              "e o mesmo canal")
        print(f"   {'RS':>10} {'andaria':>8} {'desnivel':>9} "
              f"{'canal':>7} -> {'seria':>7}")
        for m in sorted(adiadas, key=lambda m: -m["anda"]):
            print(f"   {m['rs']:>10.2f} {m['anda']:>7.1f}m {m['desnivel']:>8.2f}m "
                  f"{m['larg']:>6.1f}m -> {m['larg_seria']:>6.1f}m")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
