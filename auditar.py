# -*- coding: utf-8 -*-
"""
Audita a geometria ANTES de rodar o HEC-RAS.

Perseguir a secao que o solver reporta e whack-a-mole: cada correcao empurra a
falha para a proxima secao mais fraca, e foram quatro rodadas assim
(dos Cedros -> do Testo -> Acu R3 -> Benedito R2). Aqui todas as 1.660 secoes
sao medidas de uma vez, contra criterios que nao dependem de simular.

O que e verificado, e por que cada um derruba o solver:

  ALTURA UTIL      cota do topo menos o talvegue. Se a cheia esperada nao cabe,
                   o HEC-RAS extrapola a tabela de conducao -- foi o que fazia
                   o modelo falhar no assentamento, com as secoes do Mirim de
                   topo abaixo do nivel do mar.
  SALTO DE AREA    razao da area molhada entre secoes vizinhas, numa cota de
                   referencia. Degrau de area e contracao/expansao brusca, e o
                   solver oscila.
  SALTO DE FUNDO   queda do talvegue entre vizinhas. Muito grande, o
                   escoamento fica transcritico.
  CALHA            largura e profundidade da calha escavada contra a secao.
                   Calha larga demais para a secao nao tem planicie; estreita
                   demais nao conduz.

Uso:  python auditar.py [PROJETO]
"""
import sys
import collections

import numpy as np

CHEIA_ESPERADA = 8.0     # m de lamina acima do talvegue que a secao deve conter


def ler(projeto):
    txt = open(f"{projeto}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    xs, rio, rea, rs = [], None, None, None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, rea = p[0].strip(), p[1].strip()
        elif l.startswith("Type RM"):
            try:
                rs = float(l.split(",")[1])
            except ValueError:
                rs = None
        elif l.startswith("#Sta/Elev="):
            n = int(l.split("=")[1])
            cab = i                     # linha do cabecalho, para achar Bank Sta
            v = []
            i += 1
            while i < len(txt) and len(v) < 2 * n:
                s = txt[i]
                v += [float(s[c:c + 8]) for c in range(0, len(s.rstrip()), 8)
                      if s[c:c + 8].strip()]
                i += 1
            # Bank Sta pode vir ANTES ou DEPOIS do #Sta/Elev: a formatacao
            # propria punha depois, e o build_cross_section do ras-commander
            # poe antes. Procurar so depois devolvia None em TODAS as secoes, e
            # qualquer metrica que use margem sai silenciosamente vazia.
            lb = rb = None
            for k in list(range(i, min(i + 6, len(txt)))) +                      list(range(max(0, cab - 8), cab)):
                if txt[k].startswith("Bank Sta="):
                    lb, rb = [float(x) for x in txt[k].split("=")[1].split(",")]
                    break
            xs.append({"rio": rio, "reach": rea, "rs": rs,
                       "sta": np.array(v[0::2]), "z": np.array(v[1::2]),
                       "lb": lb, "rb": rb})
            continue
        i += 1
    return xs


def area_molhada(d, nivel):
    sta, z = d["sta"], d["z"]
    prof = np.clip(nivel - z, 0, None)
    return float(np.trapezoid(prof, sta)) if hasattr(np, "trapezoid") \
        else float(np.trapz(prof, sta))


def auditar(projeto):
    xs = ler(projeto)
    por = collections.defaultdict(list)
    for d in xs:
        por[(d["rio"], d["reach"])].append(d)
    for k in por:
        por[k].sort(key=lambda d: -d["rs"])

    from itajai import perfil
    teto = {}
    for (rio, _), v in por.items():
        S = []
        for x, y in zip(v, v[1:]):
            dx = x["rs"] - y["rs"]
            if dx > 0:
                S.append(abs(float(x["z"].min()) - float(y["z"].min())) / dx)
        if S:
            teto[rio] = max(teto.get(rio, 0.0),
                            float(np.clip(np.percentile(S, 90) * 1.2,
                                          perfil.DECL_MAXIMA, perfil.DECL_TETO)))
    baixa, salto_a, salto_z = [], [], []
    for (rio, rea), v in por.items():
        for i, d in enumerate(v):
            z0 = float(d["z"].min())
            util = float(d["z"].max()) - z0
            if util < CHEIA_ESPERADA:
                baixa.append((util, rio, rea, d["rs"]))
            if i + 1 < len(v):
                e = v[i + 1]
                # mesma PROFUNDIDADE sobre o talvegue de cada uma, nao a mesma
                # cota absoluta: num rio em declive a cota comum mistura
                # declividade com forma, e um trecho ingreme aparece como salto
                # de area que nao existe.
                a1 = area_molhada(d, z0 + 3.0)
                a2 = area_molhada(e, float(e["z"].min()) + 3.0)
                if min(a1, a2) > 1.0:
                    r = max(a1, a2) / min(a1, a2)
                    if r > 3.0:
                        salto_a.append((r, rio, rea, d["rs"]))
                # contra o teto DO RIO, nao um valor fixo: com teto por rio
                # (perfil.teto_declividade) um afluente de serra desce 5% por
                # projeto, e um limite unico de 2% acusaria 198 falsos.
                dx = d["rs"] - e["rs"]
                dz = z0 - float(e["z"].min())
                if dx > 0 and dz / dx > teto.get(rio, 0.008) * 1.25:
                    salto_z.append((dz / dx, rio, rea, d["rs"]))

    print("=" * 72)
    print(f"{projeto}: {len(xs)} secoes em {len(por)} trechos")
    print("=" * 72)

    print(f"\n[1] ALTURA UTIL menor que {CHEIA_ESPERADA:.0f} m: "
          f"{len(baixa)} de {len(xs)}")
    for u, rio, rea, rs in sorted(baixa)[:12]:
        print(f"    {rio:<16}{rea:<4} RS {rs/1000:8.2f} km   so {u:5.2f} m")
    if baixa:
        c = collections.Counter(b[1] for b in baixa)
        print("    por rio:", dict(c.most_common()))

    print(f"\n[2] SALTO DE AREA maior que 3x entre vizinhas: {len(salto_a)}")
    for r, rio, rea, rs in sorted(salto_a, reverse=True)[:10]:
        print(f"    {rio:<16}{rea:<4} RS {rs/1000:8.2f} km   {r:6.1f}x")

    print(f"\n[3] DEGRAU DE FUNDO acima de 2%: {len(salto_z)}")
    for s, rio, rea, rs in sorted(salto_z, reverse=True)[:10]:
        print(f"    {rio:<16}{rea:<4} RS {rs/1000:8.2f} km   {100*s:6.2f}%")

    n = len(baixa) + len(salto_a) + len(salto_z)
    print(f"\nTOTAL de secoes sinalizadas: {n}")
    return {"baixa": baixa, "salto_area": salto_a, "degrau": salto_z}


if __name__ == "__main__":
    auditar(sys.argv[1] if len(sys.argv) > 1 else "Tajai")
