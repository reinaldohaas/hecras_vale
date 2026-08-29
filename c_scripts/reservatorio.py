# -*- coding: utf-8 -*-
"""Roteia um hidrograma por uma barragem de contencao (detencao com alvo).

    python scripts/reservatorio.py --listar
    python scripts/reservatorio.py --barragem barragem_sul --rio Itajai_Sul

O QUE E, E O QUE NAO E

  Nao e uma inline structure do HEC-RAS: nao ha cota de crista no datum das
  secoes, nem curva cota-volume medida -- so o que `dados_estruturas/
  barragens_itajai.json` traz: CAPACIDADE de acumulacao e VAZAO MAXIMA do
  extravasor. Com isso da para um modelo de DETENCAO por balanco de massa, que
  e o que importa para cheia: quanto do pico o reservatorio segura.

  Os parametros fisicos (capacidade, vazao maxima) sao reais. A REGRA DE
  OPERACAO -- quanto soltar a cada hora -- e a variavel de decisao, e por isso
  fica exposta:

    alvo Q_alvo   solta ate Q_alvo enquanto couber; o excedente fica retido.
                  Cheio o reservatorio, solta o que entra (limitado ao
                  extravasor). Vazio, nunca solta mais que a entrada.

  ALVO OTIMO. Existe um Q_alvo que da a MAIOR reducao de pico possivel com
  aquela capacidade: o nivel Qc tal que o volume do hidrograma ACIMA de Qc
  iguala a capacidade. Abaixo dele o reservatorio transbordaria; acima,
  sobraria espaco. `alvo_otimo` calcula esse Qc do proprio hidrograma -- nao e
  numero escolhido a dedo, e o limite fisico do que a barragem consegue.

  Nao substitui a operacao real de 1983, que nao esta no repositorio. E o
  teto do que cada barragem pode fazer, e a base para o sistema de decisao
  comparar estrategias.
"""
import argparse
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))
from projeto_rio_avulso import hidrograma    # noqa: E402

DADOS = "dados_estruturas/barragens_itajai.json"
HIDRO = "legado/Itajai_Rede_1983.u01"


def carregar():
    d = json.load(open(DADOS, encoding="utf-8"))
    return {b["id"]: b for b in d["barragens"]}, d.get(
        "pontos_criticos_inundacao", [])


def alvo_otimo(Q_in, dt_s, capacidade_m3):
    """Qc tal que o volume acima de Qc == capacidade (maior corte possivel).

    Se nem retendo tudo acima do menor pico a capacidade e atingida, devolve o
    Qc que enche exatamente o reservatorio -- o corte maximo viavel.
    """
    Q = np.asarray(Q_in, float)
    lo, hi = 0.0, float(Q.max())
    for _ in range(60):
        qc = (lo + hi) / 2
        vol = np.sum(np.clip(Q - qc, 0, None)) * dt_s
        if vol > capacidade_m3:
            lo = qc            # precisa cortar menos (soltar mais)
        else:
            hi = qc
    return hi


def rotear(Q_in, dt_s, capacidade_m3, q_max_m3s, q_alvo_m3s):
    """Balanco de massa horario. Devolve (Q_out, S) do mesmo tamanho de Q_in.

        S_{t+1} = S_t + (Q_in - Q_out) * dt,   0 <= S <= capacidade
        Q_out   = min(q_alvo, extravasor), mas sobe se o reservatorio encheria,
                  e nunca passa do extravasor nem do que ha para soltar.
    """
    Q = np.asarray(Q_in, float)
    out = np.zeros_like(Q)
    Sserie = np.zeros_like(Q)
    S = 0.0
    for t in range(len(Q)):
        alvo = min(q_alvo_m3s, q_max_m3s)
        # se soltar `alvo` estouraria a capacidade, solta o suficiente p/ nao
        # passar de capacidade (limitado ao extravasor)
        excesso = S + (Q[t] - alvo) * dt_s
        if excesso > capacidade_m3:
            o = Q[t] + (excesso - capacidade_m3) / dt_s
        else:
            o = alvo
        o = min(o, q_max_m3s)                 # nao passa do extravasor
        o = min(o, Q[t] + S / dt_s)           # nao solta agua que nao tem
        o = max(o, 0.0)
        S = S + (Q[t] - o) * dt_s
        S = min(max(S, 0.0), capacidade_m3)
        out[t] = o
        Sserie[t] = S
    return out, Sserie


def resumo(nome, Q_in, Q_out, dt_s, cap):
    pin, pout = Q_in.max(), Q_out.max()
    tpico_in = int(np.argmax(Q_in))
    tpico_out = int(np.argmax(Q_out))
    ret = np.max(np.cumsum((Q_in - Q_out) * dt_s))
    print(f"{nome}")
    print(f"   pico entrada {pin:7.0f} m3/s  ->  saida {pout:7.0f} m3/s   "
          f"reducao {100*(1-pout/pin):4.1f}%")
    print(f"   pico atrasa  {tpico_out - tpico_in:+d} h")
    print(f"   volume retido no maximo {ret/1e6:6.1f} hm3  "
          f"(capacidade {cap/1e6:.0f} hm3, uso {100*ret/cap:4.0f}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--barragem", default=None)
    ap.add_argument("--rio", default=None, help="reach do legado p/ hidrograma")
    ap.add_argument("--alvo", type=float, default=None,
                    help="vazao alvo m3/s; padrao = alvo otimo do hidrograma")
    ap.add_argument("--listar", action="store_true")
    a = ap.parse_args()

    barr, criticos = carregar()
    if a.listar or not a.barragem:
        print("barragens em", DADOS)
        for bid, b in barr.items():
            print(f"   {bid:16} {b['nome']:28} cap {b['capacidade_acumulacao_m3']/1e6:.0f} hm3"
                  f"   extravasor {b['vazao_maxima_extravasor_m3s']:.0f} m3/s")
        print("\npontos criticos:")
        for p in criticos:
            print(f"   {p['nome']:12} alerta {p['cota_alerta_m']} m  "
                  f"emergencia {p['cota_emergencia_m']} m")
        if not a.barragem:
            return
    b = barr[a.barragem]
    rio = a.rio
    if rio is None:
        raise SystemExit("informe --rio (reach do legado, ex.: Itajai_Sul)")
    Q = hidrograma(HIDRO, rio)
    if Q is None:
        raise SystemExit(f"sem hidrograma de {rio}")
    dt = 3600.0
    cap = b["capacidade_acumulacao_m3"]
    qmax = b["vazao_maxima_extravasor_m3s"]
    qalvo = a.alvo if a.alvo is not None else alvo_otimo(Q, dt, cap)
    Qout, _ = rotear(Q, dt, cap, qmax, qalvo)
    print(f"\n{b['nome']}  (rio {rio})")
    print(f"   alvo de liberacao {qalvo:.0f} m3/s"
          + ("  (otimo do hidrograma)" if a.alvo is None else ""))
    resumo(b["nome"], Q, Qout, dt, cap)
    return Qout


if __name__ == "__main__":
    main()
