# -*- coding: utf-8 -*-
"""Troca interpolacao por terreno onde o vao entre secoes REAIS e longo demais.

    python scripts/recortar_vaos_longos.py modelo/so_mirim.g06 --saida g07 --vao 300

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DIAGNOSTICO QUE JUSTIFICA (validado em doc/QC_so_mirim_validation.md)

  As tres categorias espaciais abertas sao uma so: 97% das secoes com
  intersecao multipla e 62% das com overlap tambem estao marcadas por angulo.
  Uniao: 278 secoes, 19,6%.

  E 94% delas sao INTERPOLADAS. A taxa de reprovacao das interpoladas e 23,6%
  contra 5,7% das cortadas do terreno -- quatro vezes maior.

  A causa nao e o meandro: a divergencia angular entre as duas cortadas que
  cercam cada interpolada e 35,7 graus de mediana TANTO nas marcadas quanto
  nas nao marcadas. O que separa e o VAO entre as secoes reais:

      vao entre cortadas    n     reprovadas    taxa
      150 - 300 m         557          22        4%
      300 - 500 m         164          57       35%
      500 - 800 m         266         128       48%
      800 -1200 m         111          53       48%
      CORTADAS            318          18        6%

  Degrau limpo em ~300 m. Das 278 marcadas, 238 (86%) sao interpoladas em vao
  maior que 300 m. Interpolar a direcao de uma cutline atraves de um vao longo
  nao funciona: o rio vira demais no meio e a media das direcoes das pontas
  nao acompanha.

O QUE SE FAZ, E O QUE NAO SE INVENTA

  Para cada interpolada num vao > `--vao`, a secao e RECORTADA do terreno:

    centro     ponto do eixo na estaca `s` da propria secao
    direcao    perpendicular a tangente local do eixo, janela ADAPTATIVA
               (2x a largura do canal, entre 20 e 150 m) -- nunca fixa
    largura    A QUE A SECAO JA TINHA. Nao ha regra de largura nova, nem
               teto de 350 m, nem envelope inventado
    estacas    comprimento acumulado da propria polilinha, entao
               L_cutline == station[-1] - station[0] por construcao
    cotas      amostradas do MDT SIG-SC 1 m; onde falta dado, mantem-se a
               cota que a secao tinha naquela fracao da largura
    margens    a largura de canal QUE ELA JA TINHA, centrada no cruzamento
               com o eixo
    Manning    os MESMOS tres valores de n; so as quebras andam para as
               novas lb/rb, com a precisao por coluna (n em 3 casas)
    HTab       inicio 2 cm acima do novo talvegue; incremento e contagem
               mantidos

  NAO se mexe em RS, em comprimento de trecho, nem em nenhuma secao cortada
  do terreno, nem em interpolada de vao curto.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                    # noqa: E402
from qc_geometria import ler_eixos, tangente_local  # noqa: E402
from mdt_sigsc import MosaicoSigsc                  # noqa: E402

FOLGA_HTAB = 0.02


def _col(v, larguras):
    """Coluna fixa de 8 chars, 10 por linha."""
    saida, linha = [], ""
    for i, x in enumerate(v):
        linha += "%8.2f" % x
        if (i + 1) % 10 == 0:
            saida.append(linha); linha = ""
    if linha:
        saida.append(linha)
    return saida


def _fmt_bank(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString, Point
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g07"
    VAO = float(argv[argv.index("--vao") + 1]) if "--vao" in argv else 300.0
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    ordem = sorted(range(len(S)), key=lambda i: -S[i]["rs"])
    eixo = list(ler_eixos(entrada).values())[0]
    mdt = MosaicoSigsc(tiles=open(os.path.join(
        raiz, f"sigsc_tiles_{base}.txt")).read().split("\n"))

    # quais sao REAIS: as que o estado do pipeline marca como nao interpoladas
    import pickle
    est = pickle.load(open(os.path.join(raiz, f"estado_{base}.pkl"), "rb"))
    chave = next(iter(est["xs_pronto"]))
    real = {round(float(x["rs"]), 2) for x in est["xs_pronto"][chave]
            if not x.get("interpolada")}

    rs = np.array([S[i]["rs"] for i in ordem])
    ereal = np.array([round(r, 2) in real for r in rs])
    idxr = np.flatnonzero(ereal)

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"vao    : {VAO:.0f} m")
    print(f"secoes : {len(S)}   reais {ereal.sum()}   interpoladas {(~ereal).sum()}")

    alvo = []
    for k in range(len(ordem)):
        if ereal[k]:
            continue
        a = idxr[idxr < k]; b = idxr[idxr > k]
        if not len(a) or not len(b):
            continue
        if (rs[a[-1]] - rs[b[0]]) > VAO:
            alvo.append(k)
    print(f"a recortar: {len(alvo)} interpoladas em vao > {VAO:.0f} m")

    # ---------------------------------------------------------------- recorte
    novas = {}
    sem_mdt = 0
    for k in alvo:
        i = ordem[k]
        d = S[i]
        st, z = d["sta"], d["z"]
        L = float(st[-1] - st[0])
        n = len(st)
        larg_canal = float(d["rb"] - d["lb"])

        # centro: o ponto do eixo na estaca desta secao
        g = LineString(d["cut"]).intersection(eixo)
        if g.is_empty:
            continue
        p = g if g.geom_type == "Point" else list(g.geoms)[0]
        s_ = float(eixo.project(p))
        C = np.array(eixo.interpolate(s_).coords[0])

        # direcao: perpendicular a tangente local, janela adaptativa
        jan = float(np.clip(2.0 * larg_canal, 20.0, 150.0))
        t = tangente_local(eixo, s_, jan)
        t = t / max(float(np.hypot(*t)), 1e-9)
        u = np.array([t[1], -t[0]])          # normal, para a direita

        # a largura E A QUE ELA JA TINHA; o centro do eixo fica no meio
        A = C - 0.5 * L * u
        B = C + 0.5 * L * u
        # estacas: comprimento acumulado -> L_cut == station range, exato
        ns = np.linspace(0.0, L, n)
        P = [A + (x / L) * (B - A) for x in ns]
        zm = mdt.cota([q[0] for q in P], [q[1] for q in P])
        falta = ~np.isfinite(zm)
        if falta.any():
            sem_mdt += int(falta.sum())
            # mantem a cota que a secao tinha na MESMA fracao da largura
            zm[falta] = np.interp(ns[falta], st - st[0], z)
        nz = zm

        # margens: a mesma largura de canal, centrada no eixo
        c = 0.5 * L
        nlb = max(0.0, c - 0.5 * larg_canal)
        nrb = min(L, c + 0.5 * larg_canal)
        novas[i] = {"sta": ns, "z": nz, "lb": nlb, "rb": nrb,
                    "cut": np.array([A, B]), "n": n}

    print(f"pontos sem MDT, herdados da secao anterior: {sem_mdt}")

    # ------------------------------------------------------------- reescrita
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    i_sec = -1
    saida, j = [], 0
    n_ok = 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        if i_sec in novas:
            nv = novas[i_sec]
            if l.startswith("XS GIS Cut Line"):
                saida.append("XS GIS Cut Line= 2 ")
                saida.append("".join("%16.4f" % x for x in
                                     (nv["cut"][0][0], nv["cut"][0][1],
                                      nv["cut"][1][0], nv["cut"][1][1])))
                j += 1
                while j < len(linhas) and linhas[j][:1] in " -0123456789" \
                        and linhas[j].strip():
                    j += 1
                continue
            if l.startswith("#Sta/Elev"):
                saida.append("#Sta/Elev= %d " % nv["n"])
                v = []
                for a, b in zip(nv["sta"], nv["z"]):
                    v.append(a); v.append(b)
                saida += _col(v, None)
                cnt = int(l.split("=")[1])
                j += 1; lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(x), 8)
                                  if x[c:c + 8].strip()])
                    j += 1
                n_ok += 1
                continue
            if l.startswith("Bank Sta="):
                saida.append("Bank Sta=%s,%s" % (_fmt_bank(nv["lb"]),
                                                 _fmt_bank(nv["rb"])))
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
                vv = [float(x) for x in bruto[:3 * cnt]]
                # os MESMOS n; so as quebras andam
                if cnt >= 3:
                    vv[0] = 0.0
                    vv[3] = nv["lb"]
                    vv[6] = nv["rb"]
                saida.append(l)
                lin, corpo = "", []
                for t_, x in enumerate(vv):
                    lin += ("%8.2f" % x if t_ % 3 == 0 else
                            "%8.3f" % x if t_ % 3 == 1 else "%8d" % int(x))
                    if (t_ + 1) % 9 == 0:
                        corpo.append(lin); lin = ""
                if lin:
                    corpo.append(lin)
                saida += corpo
                j = k2
                continue
            if l.startswith("XS HTab Starting El and Incr"):
                p_ = [x.strip() for x in l.split("=", 1)[1].split(",")]
                inc, cnt = float(p_[1]), int(p_[2])
                el = float(np.min(nv["z"])) + FOLGA_HTAB
                saida.append(f"XS HTab Starting El and Incr={el:.2f},"
                             f"{inc:.3f}, {cnt} ")
                j += 1
                continue
        saida.append(l)
        j += 1

    txt = "\n".join(saida)
    t0 = linhas[0].split("=", 1)[1]
    txt = txt.replace("Geom Title=" + t0,
                      "Geom Title=" + t0 + " + vaos longos do terreno", 1)
    open(novo, "w", encoding="latin-1", newline="\r\n").write(txt)
    print(f"secoes recortadas: {n_ok}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
