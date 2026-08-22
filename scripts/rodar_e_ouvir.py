# -*- coding: utf-8 -*-
"""Roda um projeto no HEC-RAS e guarda o que ELE diz, sem interpretar.

    python scripts/rodar_e_ouvir.py modelo/mirim_t30/mirim_t30.prj

Existe para parar de discutir o CONTADOR de erros do Validate Geometry e ler o
que o solver realmente reporta. O projeto e rodado ONDE ESTA -- ja e uma copia
isolada, montada por `montar_projeto.py`, e nenhum arquivo original e tocado.

Grava ao lado do projeto:
  <nome>_mensagens.txt   as mensagens de computo do proprio HEC-RAS
  <nome>_resumo.txt      o que se conseguiu extrair do .pNN.hdf

`clear_geompre=True` porque a superficie de interpolacao ARMAZENADA na
geometria e de outra rodada; sem limpar, o RAS reprova secao por secao com
"Stored Interpolation Surface does not contain XS'(s) at:" -- erro que fala de
cache velho e nao de geometria.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vale.terreno import HECRAS_DIR      # noqa: E402

RAS_EXE = os.path.join(HECRAS_DIR, "Ras.exe")


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    prj = os.path.abspath(argv[0])
    plano = argv[argv.index("--plano") + 1] if "--plano" in argv else "01"
    if not os.path.exists(prj):
        raise SystemExit(f"nao achei {prj}")
    raiz = os.path.dirname(prj)
    base = os.path.splitext(os.path.basename(prj))[0]

    from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan
    print(f"projeto : {prj}")
    print(f"plano   : {plano}")
    print(f"Ras.exe : {RAS_EXE}   existe={os.path.exists(RAS_EXE)}")
    t0 = time.time()
    p = init_ras_project(prj, RAS_EXE)
    ok = RasCmdr.compute_plan(plano, ras_object=p, force_rerun=True,
                              clear_geompre=True)
    dt = time.time() - t0
    print(f"\ncompute_plan devolveu {ok!r}  em {dt/60:.1f} min")

    hdf = os.path.join(raiz, f"{base}.p{plano}.hdf")
    saida = os.path.join(raiz, f"{base}_mensagens.txt")
    try:
        msgs = str(HdfResultsPlan.get_compute_messages(hdf))
    except Exception as e:                                   # noqa: BLE001
        msgs = f"(sem mensagens: {e})"
    with open(saida, "w", encoding="utf-8", errors="replace") as f:
        f.write(msgs)
    print(f"mensagens ({len(msgs)} caracteres) -> {saida}")

    # o .computeMsgs.txt fica ao lado quando o RAS o escreve
    for ext in (".computeMsgs.txt", f".p{plano}.computeMsgs.txt", ".bco01"):
        q = os.path.join(raiz, base + ext)
        if os.path.exists(q):
            print(f"   tambem existe: {os.path.basename(q)}  "
                  f"({os.path.getsize(q)} bytes)")
    return ok


if __name__ == "__main__":
    main(sys.argv[1:])
