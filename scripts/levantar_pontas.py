# -*- coding: utf-8 -*-
"""Estende a secao que termina dentro d'agua ate achar barranco no MDT.

    python scripts/levantar_pontas.py taha_ai_novo/taha_ai.g04 \
        --dem taha_ai_novo/Terrain/taha_ai_fundo_5m.tif --saida g05

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Secao cuja PONTA fica a menos de 1 m do talvegue termina dentro d'agua:
  nao ha barranco para segurar a cheia, o HEC-RAS extrapola a tabela
  hidraulica na primeira lamina e a agua "vaza" pela borda. O qc_perfis
  acusa como GRAVE ("ponta n'agua"). Medido no g04 do taha_ai_novo: 16
  secoes, todas com barranco real a 5-165 m da ponta no MDT.

O QUE SE FAZ

  So nas secoes com ponta abaixo de `--minimo` m sobre o talvegue: a
  cutline e prolongada NA PROPRIA DIRECAO, amostrando o MDT a cada
  `--passo` m, ate a cota alcancar talvegue + `--folga` m (ou o teto de
  `--teto` m). Os pontos amostrados ENTRAM no perfil -- e terreno real,
  nao rampa inventada. Se o prolongamento cruzaria algum eixo (o proprio
  ou o de outro reach), ele para antes, com recuo; se nem assim a ponta
  sair do grave, a secao fica como esta e entra no relatorio.

  Estender pela ESQUERDA desloca o zero das estacas: Bank Sta e as
  quebras do #Mann andam junto, como no corrigir_cutlines (mesmo
  mecanismo, sinal oposto). Nenhuma cota existente muda.

QUANDO ESTENDER E PROIBIDO, RECUA-SE

  Ha pontas que terminam na agua DA OUTRA PASSADA do proprio rio (meandro
  vizinho) ou do rio receptor: estender para la cruzaria um eixo -- erro
  fatal do validador -- e representaria a mesma agua duas vezes. Nesses
  casos a ponta RECUA para a estaca mais alta que ja existe no proprio
  perfil, fora do canal: o barranco que o perfil ja tem passa a ser a
  borda. So se recua se esse barranco tirar a secao do grave; se nem ele
  existir (foz colada no rio receptor), a secao fica e sai no relatorio.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from qc_geometria import ler_eixos                     # noqa: E402
from corrigir_cutlines import (mapa_reaches, travessias, _col, _fmt,
                               _arg, TOL)              # noqa: E402
from ras_io import escrever                            # noqa: E402


class DemJanela:
    """Amostra por janela de 1 pixel; `src.sample` e `src.index` derrubam o
    processo quando o numpy.linalg da rasterio nao acha o BLAS -- a conta
    afim aqui e feita a mao, sem algebra linear."""

    def __init__(self, caminho):
        import rasterio
        self.src = rasterio.open(caminho)
        self.T = self.src.transform
        self.nod = self.src.nodata

    def cota(self, xy):
        from rasterio.windows import Window
        x, y = xy
        c = int((x - self.T.c) / self.T.a)
        r = int((y - self.T.f) / self.T.e)
        if not (0 <= r < self.src.height and 0 <= c < self.src.width):
            return np.nan
        v = float(self.src.read(1, window=Window(c, r, 1, 1))[0, 0])
        if self.nod is not None and v == self.nod:
            return np.nan
        return v


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    dem_p = _arg(argv, "--dem", None)
    ext = _arg(argv, "--saida", "g05")
    minimo = _arg(argv, "--minimo", 1.0, float)
    folga = _arg(argv, "--folga", 4.0, float)
    passo = _arg(argv, "--passo", 5.0, float)
    teto = _arg(argv, "--teto", 400.0, float)
    if not dem_p:
        raise SystemExit("preciso de --dem <mdt.tif>")
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    dem = DemJanela(dem_p)
    S = ler_secoes(entrada)
    eixos = ler_eixos(entrada)
    mapa = mapa_reaches(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    novos, falhas = {}, []
    for i, d in enumerate(S):
        z = np.asarray(d["z"], float)
        st = np.asarray(d["sta"], float)
        zt = float(z.min())
        lados = []
        if z[0] - zt < minimo:
            lados.append("esq")
        if z[-1] - zt < minimo:
            lados.append("dir")
        if not lados:
            continue
        A = np.asarray(d["cut"][0], float)
        B = np.asarray(d["cut"][-1], float)
        L = float(np.hypot(*(B - A)))
        u = (B - A) / max(L, 1e-9)
        ch = mapa[i]
        lb, rb = float(d["lb"]), float(d["rb"])

        def estender(lado):
            """Lista [(dist, cota)] ate o barranco, ou None."""
            P0, dirv = (A, -u) if lado == "esq" else (B, u)
            alvo = zt + folga
            pontos, dd = [], passo
            while dd <= teto:
                P = P0 + dd * dirv
                seg_ini = P0 + 0.01 * dirv
                if any(travessias(seg_ini, dirv, 0.0, dd, e)
                       for e in eixos.values()):
                    break               # cruzaria um eixo: erro fatal
                zp = dem.cota(P)
                if np.isnan(zp):
                    break
                if zp <= zt + 0.05:
                    break               # terreno abaixo do talvegue na frente
                pontos.append((dd, zp))
                if zp >= alvo:
                    return pontos
                dd += passo
            if pontos:
                melhor = max(p[1] for p in pontos)
                if melhor - zt >= minimo:
                    k = int(np.argmax([p[1] for p in pontos]))
                    return pontos[:k + 1]     # sai do grave, fica no aviso
            return None

        def recuar(lado):
            """Estaca do barranco mais alto ja presente no perfil, ou None."""
            if lado == "esq":
                m = st < lb - 0.01
            else:
                m = st > rb + 0.01
            if not m.any():
                return None
            k = int(np.argmax(np.where(m, z, -np.inf)))
            if z[k] - zt < minimo:
                return None
            return float(st[k])

        modos = {}
        for lado in lados:
            p = estender(lado)
            if p is not None:
                modos[lado] = ("ext", p)
                continue
            s = recuar(lado)
            if s is not None:
                modos[lado] = ("rec", s)
            else:
                falhas.append((ch, d["rs"], lado,
                               "sem barranco: nem na frente (MDT/eixo) "
                               "nem no proprio perfil"))
        if not modos:
            continue

        ns, nz = list(st), list(z)
        desl_e = 0.0        # positivo = recuo, negativo = extensao
        corte_d = 0.0
        rot = []
        if "esq" in modos:
            modo, val = modos["esq"]
            if modo == "ext":
                ns = [-p[0] for p in reversed(val)] + ns
                nz = [p[1] for p in reversed(val)] + nz
                desl_e = -val[-1][0]
                rot.append(f"+{val[-1][0]:.0f} m esq")
            else:
                keep = [(a, b) for a, b in zip(ns, nz) if a >= val - 0.01]
                ns = [a for a, _ in keep]
                nz = [b for _, b in keep]
                desl_e = val
                rot.append(f"-{val:.0f} m esq (recuo)")
        if "dir" in modos:
            modo, val = modos["dir"]
            if modo == "ext":
                ns = ns + [L + p[0] for p in val]
                nz = nz + [p[1] for p in val]
                corte_d = -val[-1][0]
                rot.append(f"+{val[-1][0]:.0f} m dir")
            else:
                keep = [(a, b) for a, b in zip(ns, nz) if a <= val + 0.01]
                ns = [a for a, _ in keep]
                nz = [b for _, b in keep]
                corte_d = L - val
                rot.append(f"-{L - val:.0f} m dir (recuo)")
        ns = np.asarray(ns, float) - desl_e    # zero no primeiro ponto
        nz = np.asarray(nz, float)
        novos[i] = {"sta": ns, "z": nz,
                    "lb": lb - desl_e, "rb": rb - desl_e,
                    "desl": desl_e,
                    "cut": (A + desl_e * u, B - corte_d * u)}
        print(f"   {ch[0]:14s} {ch[1]:3s} RS {d['rs']:9.2f}  "
              f"{'  '.join(rot):26s} ({len(ns)} pontos)")

    for ch, rs, lado, m in falhas:
        print(f"   FICOU: {ch[0]} {ch[1]} RS {rs:.2f} ({lado}): {m}")
    if not novos:
        print("nenhuma secao precisou de extensao")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
            if l.startswith("XS GIS Cut Line"):
                saida.append("XS GIS Cut Line= 2")
                saida.append("".join("%16.2f" % x for x in
                                     (nv["cut"][0][0], nv["cut"][0][1],
                                      nv["cut"][1][0], nv["cut"][1][1])))
                j += 1
                while j < len(linhas) and linhas[j].strip() and \
                        linhas[j][:1] in " -0123456789":
                    j += 1
                continue
            if l.startswith("#Sta/Elev"):
                v = []
                for a, b in zip(nv["sta"], nv["z"]):
                    v += [a, b]
                saida.append("#Sta/Elev= %d " % len(nv["sta"]))
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
                             % (_fmt(nv["lb"]), _fmt(nv["rb"])))
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
                topo = float(nv["sta"][-1])
                for t in range(0, 3 * cnt, 3):
                    val[t] = min(max(val[t] - nv["desl"], 0.0), topo)
                # o HEC-RAS exige n na PRIMEIRA estaca: apos estender pela
                # esquerda a primeira quebra iria parar em ext_e, deixando o
                # trecho novo sem n ("A horizontal Manning's n value needs
                # to be specified on first station")
                val[0] = 0.0
                saida.append(l)
                lin, corpo = "", []
                for t, x in enumerate(val):
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

    # -------------------------------------------------------- conferencia
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print(f"   secoes: {len(A2)} -> {len(B2)}")
    tal = max(abs(float(a['z'].min()) - float(b['z'].min()))
              for a, b in zip(A2, B2))
    print(f"   talvegue mudou no maximo {tal:.6f}  (tem de ser zero)")
    graves = sum(1 for b in B2
                 if min(b['z'][0], b['z'][-1]) - b['z'].min() < minimo)
    print(f"   pontas n'agua: {graves}   "
          f"(era {sum(1 for a in A2 if min(a['z'][0], a['z'][-1]) - a['z'].min() < minimo)})")
    ruim = sum(1 for b in B2 if len(b['sta']) > 500)
    print(f"   secoes com mais de 500 pontos: {ruim}   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
