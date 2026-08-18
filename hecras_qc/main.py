# -*- coding: utf-8 -*-
"""
Ponto de entrada.

Sem argumentos, abre a interface. Com --lote, roda a mesma analise sem GUI e
grava a tabela -- util para conferir milhares de secoes de uma vez, e para
rodar num servidor sem tela. Os dois caminhos chamam exatamente as mesmas
funcoes de qc/correction: se divergissem, nao daria para confiar em nenhum.

    python -m hecras_qc.main
    python -m hecras_qc.main --lote --dem t.tif --eixo rio.geojson \\
                             --secoes xs.geojson --saida qc.csv
"""
import argparse
import sys


def lote(args):
    import numpy as np
    from . import correction, cross_sections, export, qc
    from .dem import DEM
    from .river_axis import EixoRio

    dem = DEM(args.dem)
    print("DEM:", dem.resumo())
    eixo = EixoRio.ler(args.eixo, dem.crs_metrico) if args.eixo else None
    secoes = cross_sections.carregar(args.secoes, dem.crs_metrico, eixo)
    print(f"{len(secoes)} secoes")

    lim = qc.Limiares(espacamento=args.espacamento,
                      proeminencia_min=args.proeminencia)
    for s in secoes:
        s.extrair(dem, lim.espacamento, eixo, lim.proeminencia_min)
    qc.avaliar_todas(secoes, lim)

    c = qc.contagem(secoes)
    print(f"OK {c[qc.OK]} | atencao {c[qc.ATENCAO]} | "
          f"incerto {c[qc.INCERTO]} | critica {c[qc.CRITICA]}")

    ruins = [s for s in secoes if s.qc.status in (qc.CRITICA, qc.INCERTO)]
    print(f"\n{'RS':>12} {'largura':>8} {'talv%':>7} {'prof':>7} "
          f"{'orient':>7} {'QC':>5}  motivo")
    for s in ruins[:40]:
        p = 100 * s.posicao_relativa
        print(f"{s.rs if s.rs is not None else s.idx:>12} "
              f"{s.largura:>8.0f} {p:>7.1f} {s.profundidade_relativa:>7.2f} "
              f"{(s.azimute if s.azimute is not None else float('nan')):>7.1f} "
              f"{s.qc.nota:>5.0f}  {s.qc.resumo[:60]}")
    if len(ruins) > 40:
        print(f"... e mais {len(ruins)-40}")

    if args.corrigir and eixo is not None:
        print("\npropondo correcao para as problematicas...")
        melhorou = 0
        for s in ruins:
            nova, _ = correction.propor(s, dem, eixo, lim)
            if nova is not None:
                i = secoes.index(s)
                secoes[i] = nova
                melhorou += 1
        qc.avaliar_todas(secoes, lim)
        c = qc.contagem(secoes)
        print(f"{melhorou} secoes com proposta melhor.  depois: "
              f"OK {c[qc.OK]} | atencao {c[qc.ATENCAO]} | "
              f"incerto {c[qc.INCERTO]} | critica {c[qc.CRITICA]}")

    if args.saida:
        print("tabela:", export.exportar_tabela(secoes, args.saida))
    if args.saida_vetor:
        print("vetor:", export.exportar_vetor(secoes, args.saida_vetor,
                                              dem.crs_metrico))
    if args.saida_csv:
        print("perfis:", export.exportar_csv_perfis(secoes, args.saida_csv))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        description="QC geometrico de secoes transversais para HEC-RAS")
    p.add_argument("--lote", action="store_true", help="rodar sem interface")
    p.add_argument("--dem")
    p.add_argument("--eixo")
    p.add_argument("--secoes")
    p.add_argument("--saida", help="CSV da tabela de QC")
    p.add_argument("--saida-vetor", help="GeoJSON/SHP das secoes com QC")
    p.add_argument("--saida-csv", help="CSV dos perfis no formato HEC-RAS")
    p.add_argument("--espacamento", type=float, default=2.0)
    p.add_argument("--proeminencia", type=float, default=0.5)
    p.add_argument("--corrigir", action="store_true",
                   help="propor correcao geometrica para as problematicas")
    a = p.parse_args(argv)

    if a.lote:
        if not (a.dem and a.secoes):
            p.error("--lote exige --dem e --secoes")
        return lote(a)
    from .gui import rodar
    return rodar()


if __name__ == "__main__":
    sys.exit(main())
