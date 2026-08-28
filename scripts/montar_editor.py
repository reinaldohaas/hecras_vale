# -*- coding: utf-8 -*-
"""Monta a pagina do editor de linhas a partir do template + linhas + fundo.

    python scripts/montar_editor.py
        (le doc/linhas_hand/linhas_editor.json, doc/figuras/fundo_editor.jpg
         e doc/editor/editor_template.html; grava doc/editor/editor_edges.html)

Depois de montar, o Claude publica no artifact (mesmo link).
"""
import base64
import json
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    tpl = open(os.path.join(RAIZ, 'doc/editor/editor_template.html'),
               encoding='utf-8').read()
    img = base64.b64encode(
        open(os.path.join(RAIZ, 'doc/figuras/fundo_editor.jpg'),
             'rb').read()).decode()
    tpl = tpl.replace('@@IMG@@', 'data:image/jpeg;base64,' + img)
    dados = json.load(open(os.path.join(RAIZ,
                                        'doc/linhas_hand/linhas_editor.json')))
    if not any(l['nome'].startswith('Divisor') for l in dados['linhas']):
        print('aviso: sem linha do Divisor no JSON')
    saida = tpl.replace('@@DADOS@@', json.dumps(dados)) \
               .replace('@@TPL@@', base64.b64encode(
                   tpl.encode('utf-8')).decode())
    destino = os.path.join(RAIZ, 'doc/editor/editor_edges.html')
    open(destino, 'w', encoding='utf-8').write(saida)
    print(f'{destino}: {len(saida)/1e6:.2f} MB, '
          f"{len(dados['linhas'])} linhas")


if __name__ == '__main__':
    main()
