# -*- coding: utf-8 -*-
"""
Motor hidraulico quasi-permanente da bacia do Itajai.

MOTIVO
------
O modelo unsteady do HEC-RAS resolve Saint-Venant completo por esquema
implicito acoplado no tempo. Nesta bacia ele fica numa borda tao fina que muda
de comportamento com 3 secoes de diferenca na geometria: o mesmo modelo rodou
99 e 29 dos 192 passos com alteracoes minimas. A causa nao e o solver -- e a
bacia. O Alto Vale e plano (0,13 m/km) e drena pela garganta do Salto Pilao,
que cai 195 m em 13 km. Transicao plano<->ingreme com Froude perto de 1 e o
pior caso possivel para aquele esquema.

Aqui a mesma fisica e resolvida em duas etapas desacopladas, do jeito que o
proprio HEC-RAS faz no modo quasi-unsteady (o que ele usa para sedimentos):

  1. VAZAO: roteada por Muskingum-Cunge ao longo das secoes, com os
     parametros tirados da hidraulica local.
  2. COTA: a cada instante, um remanso PERMANENTE por passo padrao com a
     vazao ja roteada, de jusante para montante.

O que se ganha: nao existe divergencia. O passo padrao e uma raiz de funcao
secao a secao, resolvida por bisseccao -- sempre converge; e onde o escoamento
fica supercritico usa profundidade critica em vez de estourar. E faz REMANSO de
verdade, que e o mecanismo que alaga Rio do Sul e Ituporanga (barramento da
garganta subindo o vale plano) e que o HAND nao reproduz.

O que se perde: propagacao transiente dentro do passo horario. Para mancha de
inundacao isso e irrelevante perto da incerteza do terreno e da chuva.

A geometria vem do .g01 ja validado -- as mesmas secoes cortadas do SIG-SC de
1 m, com calha escavada e Manning de Jarrett nas gargantas.

Uso:   python motor_hidraulico.py Itajai_Rede_1983
Saida: <PROJETO>_motor.npz   (cota e vazao por secao e por instante)
"""
import os
import sys

import numpy as np

G = 9.80665
NIVEIS = 80          # niveis da tabela de propriedades por secao


# ================================================================== GEOMETRIA
def ler_geometria(projeto):
    """Secoes e juncoes do .g01."""
    txt = open(f"{projeto}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    secoes, juncoes = [], []
    rio = reach = rs = None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, reach = p[0].strip(), p[1].strip()
        elif l.startswith("Junct Name="):
            j = {"nome": l.split("=", 1)[1].strip(), "up": [], "dn": None}
            k = i + 1
            while k < len(txt) and txt[k].strip():
                if txt[k].startswith("Up River,Reach="):
                    p = txt[k].split("=", 1)[1].split(",")
                    j["up"].append((p[0].strip(), p[1].strip()))
                elif txt[k].startswith("Dn River,Reach="):
                    p = txt[k].split("=", 1)[1].split(",")
                    j["dn"] = (p[0].strip(), p[1].strip())
                k += 1
            juncoes.append(j)
            i = k
            continue
        elif l.startswith("Type RM"):
            try:
                rs = float(l.split(",")[1])
            except ValueError:
                rs = None
        elif l.startswith("#Sta/Elev="):
            n = int(l.split("=")[1])
            v = []
            i += 1
            while i < len(txt) and len(v) < 2 * n:
                s = txt[i]
                v += [float(s[c:c + 8]) for c in range(0, len(s.rstrip()), 8)
                      if s[c:c + 8].strip()]
                i += 1
            sta, z = np.array(v[0::2]), np.array(v[1::2])
            n_esq = n_ch = n_dir = 0.035
            lb = rb = None
            for k in range(i, min(i + 8, len(txt))):
                if txt[k].startswith("#Mann="):
                    w = txt[k + 1]
                    val = [w[c:c + 8].strip()
                           for c in range(0, len(w.rstrip()), 8)]
                    val = [x for x in val if x]
                    if len(val) >= 9:
                        n_esq, n_ch, n_dir = (float(val[1]), float(val[4]),
                                              float(val[7]))
                elif txt[k].startswith("Bank Sta="):
                    lb, rb = [float(x) for x in txt[k].split("=")[1].split(",")]
                    break
            if lb is None:
                lb, rb = sta[0], sta[-1]
            secoes.append({"rio": rio, "reach": reach, "rs": rs,
                           "sta": sta, "z": z, "lb": lb, "rb": rb,
                           "n_esq": n_esq, "n_ch": n_ch, "n_dir": n_dir})
            continue
        i += 1
    return secoes, juncoes


def ler_fluxo(projeto):
    """Contornos do .u01: hidrogramas, laterais uniformes e mare."""
    txt = open(f"{projeto}.u01", encoding="utf-8", errors="ignore").read().splitlines()

    def serie(i0, n):
        v, j = [], i0
        while len(v) < n and j < len(txt):
            s = txt[j]
            v += [float(s[c:c + 8]) for c in range(0, len(s.rstrip()), 8)
                  if s[c:c + 8].strip()]
            j += 1
        return np.array(v[:n]), j

    loc, itens = None, []
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("Boundary Location="):
            p = l.split("=", 1)[1].split(",")
            loc = (p[0].strip(), p[1].strip(),
                   float(p[2]) if p[2].strip() else None,
                   float(p[3]) if len(p) > 3 and p[3].strip() else None)
        elif l.startswith("Flow Hydrograph="):
            v, i = serie(i + 1, int(l.split("=")[1]))
            itens.append({"tipo": "vazao", "loc": loc, "v": v})
            continue
        elif l.startswith("Stage Hydrograph="):
            v, i = serie(i + 1, int(l.split("=")[1]))
            itens.append({"tipo": "cota", "loc": loc, "v": v})
            continue
        elif l.startswith("Uniform Lateral Inflow Hydrograph="):
            v, i = serie(i + 1, int(l.split("=")[1]))
            itens.append({"tipo": "uniforme", "loc": loc, "v": v})
            continue
        elif l.startswith("Lateral Inflow Hydrograph="):
            v, i = serie(i + 1, int(l.split("=")[1]))
            itens.append({"tipo": "lateral", "loc": loc, "v": v})
            continue
        i += 1
    return itens


# ================================================================ HIDRAULICA
def _fatia(sta, z, zw, a, b):
    """Area, perimetro e largura molhados entre as estacas a e b."""
    m = (sta >= a) & (sta <= b)
    if m.sum() < 2:
        return 0.0, 0.0, 0.0
    s, e = sta[m], z[m]
    d = zw - e
    A = P = T = 0.0
    for k in range(len(s) - 1):
        d1, d2, dx = d[k], d[k + 1], s[k + 1] - s[k]
        if d1 <= 0 and d2 <= 0:
            continue
        if d1 > 0 and d2 > 0:                      # painel todo molhado
            A += 0.5 * (d1 + d2) * dx
            P += float(np.hypot(dx, e[k + 1] - e[k]))
            T += dx
        else:                                      # painel parcialmente seco
            f = d1 / (d1 - d2) if (d1 - d2) != 0 else 0.0
            xw = dx * (f if d1 > 0 else 1.0 - f)
            prof = d1 if d1 > 0 else d2
            A += 0.5 * prof * xw
            P += float(np.hypot(xw, prof))
            T += xw
    return A, P, T


def _hidraulica(sec, zw):
    """(A, T, K) na cota zw, com conducao COMPOSTA das tres zonas."""
    sta, z = sec["sta"], sec["z"]
    A_t = T_t = K_t = 0.0
    for a, b, n in ((sta[0], sec["lb"], sec["n_esq"]),
                    (sec["lb"], sec["rb"], sec["n_ch"]),
                    (sec["rb"], sta[-1], sec["n_dir"])):
        if b <= a:
            continue
        A, P, T = _fatia(sta, z, zw, a, b)
        if A <= 0 or P <= 0:
            continue
        K_t += (1.0 / max(n, 1e-3)) * A * (A / P) ** (2.0 / 3.0)
        A_t += A
        T_t += T
    return A_t, T_t, K_t


def tabela(sec):
    """Tabela (cota -> A, T, K) montada UMA vez por secao.

    Resolver a hidraulica painel a painel dentro de cada bisseccao sairia caro:
    sao ~1.200 secoes x 192 instantes x dezenas de iteracoes. Com a tabela,
    todo o resto vira interpolacao. E o que o HEC-RAS faz em "Preprocessing
    Geometry".
    """
    z0, z1 = sec["z"].min(), sec["z"].max()
    lev = z0 + np.linspace(1e-3, max(z1 - z0, 2.0) * 1.05, NIVEIS)
    A = np.empty(NIVEIS)
    T = np.empty(NIVEIS)
    K = np.empty(NIVEIS)
    for i, zw in enumerate(lev):
        A[i], T[i], K[i] = _hidraulica(sec, zw)
    sec["_lev"], sec["_A"], sec["_T"], sec["_K"] = lev, A, T, K
    return sec


def hid(sec, zw):
    """(A, T, K) por interpolacao. Acima do topo tabelado, alarga em prisma."""
    lev = sec["_lev"]
    if zw <= lev[0]:
        return 0.0, 0.0, 0.0
    if zw >= lev[-1]:
        T = max(sec["_T"][-1], 1.0)
        A = sec["_A"][-1] + T * (zw - lev[-1])
        K = sec["_K"][-1] * (A / max(sec["_A"][-1], 1e-6)) ** (5.0 / 3.0)
        return A, T, K
    return (float(np.interp(zw, lev, sec["_A"])),
            float(np.interp(zw, lev, sec["_T"])),
            float(np.interp(zw, lev, sec["_K"])))


def z_normal(sec, Q, S):
    """Cota normal: inverte a tabela de conducao em K = Q/sqrt(S)."""
    alvo = Q / np.sqrt(max(S, 1e-5))
    K, lev = sec["_K"], sec["_lev"]
    if alvo <= K[0]:
        return float(lev[0])
    if alvo >= K[-1]:
        A = sec["_A"][-1] * (alvo / max(K[-1], 1e-9)) ** 0.6
        return float(lev[-1] + (A - sec["_A"][-1]) / max(sec["_T"][-1], 1.0))
    return float(np.interp(alvo, K, lev))


def z_critica(sec, Q):
    """Menor cota tabelada com Fr <= 1."""
    A, T, lev = sec["_A"], sec["_T"], sec["_lev"]
    with np.errstate(all="ignore"):
        fr = Q * Q * T / (G * np.maximum(A, 1e-6) ** 3)
    ok = np.flatnonzero(fr <= 1.0)
    return float(lev[ok[0]]) if len(ok) else float(lev[-1])


def energia_montante(sec_dn, sec_up, Q_dn, Q_up, z_dn, dx):
    """Cota a montante pela equacao da energia (passo padrao).

        z_up + V_up^2/2g = z_dn + V_dn^2/2g + hf + he

    Bisseccao no residuo, que e monotono no ramo subcritico. Se nem na cota
    critica a energia fecha, o controle ali E o critico: devolve ele. E assim
    que a garganta do Salto Pilao deixa de derrubar o calculo.
    """
    A_dn, _, K_dn = hid(sec_dn, z_dn)
    if A_dn <= 0 or K_dn <= 0:
        return z_dn + 0.01, True
    V_dn = Q_dn / A_dn
    E_dn = z_dn + V_dn * V_dn / (2 * G)
    Sf_dn = (Q_dn / K_dn) ** 2
    z_cr = z_critica(sec_up, Q_up)

    def resid(z):
        A, _, K = hid(sec_up, z)
        if A <= 0 or K <= 0:
            return 1e6
        V = Q_up / A
        hf = 0.5 * ((Q_up / K) ** 2 + Sf_dn) * dx
        he = 0.3 * abs(V * V - V_dn * V_dn) / (2 * G)     # expansao/contracao
        return (z + V * V / (2 * G)) - (E_dn + hf + he)

    lo = max(z_cr, float(sec_up["_lev"][0]))
    if resid(lo) > 0:
        return z_cr, True
    hi = lo + 50.0
    if resid(hi) < 0:
        return hi, True
    for _ in range(45):
        m = 0.5 * (lo + hi)
        if resid(m) < 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi), False


# ================================================================ ROTEAMENTO
def muskingum_cunge(Q_up, Q_up_ant, Q_dn_ant, dx, dt, sec, So):
    """Um passo de Muskingum-Cunge com parametros da hidraulica local.

    A celeridade cinematica e estimada em 5/3 da velocidade (canal largo), e
    dai saem K = dx/c e X = 0,5(1 - Q/(B So c dx)). Os coeficientes sao
    limitados a nao-negativos e renormalizados: e isso que impede oscilacao
    mesmo num trecho quase plano como o Itajai do Oeste (0,13 m/km), onde o
    esquema implicito secava a calha e divergia.
    """
    Qr = max(0.5 * (Q_up + Q_dn_ant), 1.0)
    A, T, _ = hid(sec, z_normal(sec, Qr, So))
    if A <= 0 or T <= 0:
        return Q_up
    c = max(5.0 / 3.0 * (Qr / A), 0.05)
    K = dx / c
    X = 0.5 * (1.0 - Qr / (max(T, 1.0) * max(So, 1e-5) * c * dx))
    X = min(max(X, 0.0), 0.5)
    den = 2.0 * K * (1.0 - X) + dt
    C = [(dt + 2.0 * K * X) / den,
         (dt - 2.0 * K * X) / den,
         (2.0 * K * (1.0 - X) - dt) / den]
    C = [max(c_, 0.0) for c_ in C]
    s = sum(C)
    if s <= 0:
        return Q_up
    C = [c_ / s for c_ in C]
    return C[0] * Q_up + C[1] * Q_up_ant + C[2] * Q_dn_ant


# ====================================================================== REDE
def montar(secoes, juncoes):
    """Trechos ordenados e ordem topologica de montante para jusante."""
    tre = {}
    for s in secoes:
        tre.setdefault((s["rio"], s["reach"]), []).append(s)
    for k in tre:
        tre[k].sort(key=lambda d: -d["rs"])          # montante -> jusante
    alim = {k: [] for k in tre}                       # quem entra em cada trecho
    for j in juncoes:
        if j["dn"] in alim:
            alim[j["dn"]] = [u for u in j["up"] if u in tre]
    ordem, visto = [], set()

    def visita(k):
        if k in visto:
            return
        visto.add(k)
        for u in alim.get(k, ()):
            visita(u)
        ordem.append(k)

    for k in tre:
        visita(k)
    return tre, alim, ordem


def declividade(xs):
    """Declividade do leito entre secoes vizinhas, limitada a faixa util."""
    z = np.array([s["z"].min() for s in xs])
    rs = np.array([s["rs"] for s in xs])
    if len(xs) < 2:
        return np.array([1e-3])
    dx = np.maximum(-np.diff(rs), 1.0)
    S = np.abs(np.diff(z)) / dx
    return np.clip(np.append(S, S[-1]), 1e-5, 0.05)


# ================================================================= SIMULACAO
def simular(projeto, verbose=True):
    secoes, juncoes = ler_geometria(projeto)
    itens = ler_fluxo(projeto)
    tre, alim, ordem = montar(secoes, juncoes)
    if verbose:
        print(f"{projeto}: {len(secoes)} secoes, {len(tre)} trechos, "
              f"{len(juncoes)} juncoes")
    for s in secoes:
        tabela(s)

    NT = max((len(it["v"]) for it in itens), default=0)
    if NT == 0:
        raise SystemExit("nenhuma serie no .u01")
    dt = 3600.0

    # quem esta a JUSANTE de cada trecho (para o remanso subir a rede)
    jus = {}
    for j in juncoes:
        for u in j["up"]:
            if u in tre and j["dn"] in tre:
                jus[u] = j["dn"]
    saida = [k for k in tre if k not in jus]

    # --- contornos, casados com o trecho pela estaca
    entrada = {k: np.zeros(NT) for k in tre}     # vazao imposta na cabeceira
    unif = {k: np.zeros(NT) for k in tre}        # lateral uniforme no trecho
    mare = None
    for it in itens:
        if it["loc"] is None:
            continue
        rio, reach, rs = it["loc"][0], it["loc"][1], it["loc"][2]
        k = (rio, reach)
        if k not in tre:
            continue
        v = np.resize(it["v"], NT)
        if it["tipo"] == "cota":
            mare = v
        elif it["tipo"] == "uniforme":
            unif[k] = unif[k] + v
        elif it["tipo"] in ("vazao", "lateral"):
            rs_topo = tre[k][0]["rs"]
            if rs is not None and abs(rs - rs_topo) < 1.0:
                entrada[k] = entrada[k] + v       # cabeceira
            else:
                unif[k] = unif[k] + v             # aporte no meio do trecho
    if mare is None:
        mare = np.full(NT, 0.3)

    decl = {k: declividade(xs) for k, xs in tre.items()}
    dxs = {k: np.maximum(-np.diff([s["rs"] for s in xs]), 1.0)
           for k, xs in tre.items()}

    # --- malha de ROTEAMENTO, separada da malha de secoes
    # Muskingum-Cunge so e valido com dx ~ c*dt. Aplicado no espacamento das
    # secoes (150 m nas gargantas) com dt de 1 h, dt fica muito maior que
    # K = dx/c: os coeficientes degeneram para C1 ~ C2 ~ 0,5, ou seja meio
    # passo de atraso POR SECAO. Com 1.202 secoes isso virava centenas de
    # horas de lag e o pico nao chegava na foz dentro da simulacao.
    # Com c da ordem de 1-3 m/s e dt de 1 h, o passo certo e alguns km.
    DX_ROTEIO = 5000.0
    nos = {}
    for k, xs in tre.items():
        idx, ultimo = [0], xs[0]["rs"]
        for i in range(1, len(xs)):
            if ultimo - xs[i]["rs"] >= DX_ROTEIO:
                idx.append(i)
                ultimo = xs[i]["rs"]
        if idx[-1] != len(xs) - 1:
            idx.append(len(xs) - 1)
        nos[k] = np.array(idx)
    Q = {k: np.zeros((NT, len(xs))) for k, xs in tre.items()}
    Z = {k: np.zeros((NT, len(xs))) for k, xs in tre.items()}

    n_crit = 0
    for t in range(NT):
        # ---------- 1. VAZAO: montante -> jusante, na malha de roteamento ----
        for k in ordem:
            xs, S, idx = tre[k], decl[k], nos[k]
            ta = t - 1 if t else t
            q0 = float(entrada[k][t])
            for u in alim.get(k, ()):
                q0 += float(Q[u][t, -1])
            rs_no = np.array([xs[i]["rs"] for i in idx])
            L = max(rs_no[0] - rs_no[-1], 1.0)
            q_no = np.empty(len(idx))
            q_no[0] = Q[k][t, idx[0]] = max(q0, 1.0)
            for m in range(1, len(idx)):
                i, ip = idx[m], idx[m - 1]
                dxr = max(rs_no[m - 1] - rs_no[m], 1.0)
                q = muskingum_cunge(q_no[m - 1], Q[k][ta, ip], Q[k][ta, i],
                                    dxr, dt, xs[i], float(S[i]))
                # lateral uniforme rateada pelo comprimento do sub-trecho
                q_no[m] = max(q + unif[k][t] * dxr / L, 1.0)
            # interpola a vazao dos nos para TODAS as secoes (a malha fina
            # existe para o remanso, nao para o roteamento)
            rs_tudo = np.array([s["rs"] for s in xs])
            Q[k][t, :] = np.interp(-rs_tudo, -rs_no, q_no)

        # ---------- 2. COTA: jusante -> montante ----------
        for k in reversed(ordem):
            xs, dx = tre[k], dxs[k]
            if k in jus:
                kd = jus[k]
                z_dn = Z[kd][t, 0]                # cota no topo do trecho de jusante
            else:
                z_dn = float(mare[t])
            # a secao mais a jusante nao pode ficar abaixo do proprio leito
            z_dn = max(z_dn, xs[-1]["z"].min() + 0.05)
            # se a cota imposta e menor que a normal, quem manda e a normal
            z_nor = z_normal(xs[-1], Q[k][t, -1], float(decl[k][-1]))
            Z[k][t, -1] = max(z_dn, z_nor)
            for i in range(len(xs) - 2, -1, -1):
                z, crit = energia_montante(xs[i + 1], xs[i],
                                           Q[k][t, i + 1], Q[k][t, i],
                                           Z[k][t, i + 1], float(dx[i]))
                Z[k][t, i] = z
                n_crit += int(crit)
        if verbose and (t % 24 == 0 or t == NT - 1):
            zz = np.concatenate([Z[k][t] for k in tre])
            qq = np.concatenate([Q[k][t] for k in tre])
            print(f"  t={t:4d} h   Q max {qq.max():8.0f} m3/s   "
                  f"cota max {zz.max():7.2f} m")

    # --- empacota na mesma ordem das secoes do .g01
    chaves, riv, rch, rs_l = [], [], [], []
    for k in sorted(tre):
        for i, s in enumerate(tre[k]):
            chaves.append((k, i))
            riv.append(s["rio"]); rch.append(s["reach"]); rs_l.append(s["rs"])
    WS = np.column_stack([Z[k][:, i] for k, i in chaves])
    QQ = np.column_stack([Q[k][:, i] for k, i in chaves])
    if verbose:
        print(f"  controle critico em {n_crit} avaliacoes "
              f"({100*n_crit/max(WS.size,1):.1f}% das secoes-instante)")
    return {"ws": WS, "q": QQ, "river": np.array(riv), "reach": np.array(rch),
            "rs": np.array(rs_l, dtype=float)}


def main():
    projeto = sys.argv[1] if len(sys.argv) > 1 else "Itajai_Rede_1983"
    if not os.path.exists(f"{projeto}.g01"):
        raise SystemExit(f"{projeto}.g01 nao encontrado")
    r = simular(projeto)
    np.savez_compressed(f"{projeto}_motor.npz", **r)
    print(f"\n[OK] {projeto}_motor.npz  "
          f"({r['ws'].shape[0]} instantes x {r['ws'].shape[1]} secoes)")


if __name__ == "__main__":
    main()
