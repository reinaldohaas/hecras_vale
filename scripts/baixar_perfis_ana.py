# -*- coding: utf-8 -*-
"""Baixa PERFIS TRANSVERSAIS e RESUMOS DE DESCARGA do HidroWebService.

    python scripts/baixar_perfis_ana.py [--saida doc/perfis_ana]

CREDENCIAIS (cadastro do usuario na ANA): arquivo doc/ana_credenciais.txt
com DUAS linhas -- Identificador na 1a, Senha na 2a. O arquivo esta no
.gitignore; nunca sobe para o repositorio.

FLUXO (API oficial https://www.ana.gov.br/hidrowebservice):
  1. GET /EstacoesTelemetricas/OAUth/v1  (headers Identificador, Senha)
     -> token Bearer (expira; renovado a cada rodada)
  2. GET /EstacoesTelemetricas/HidroSeriePerfilTransversal/v1
  3. GET /EstacoesTelemetricas/HidroSerieResumoDescarga/v1
     por estacao, gravando o JSON bruto + um CSV achatado por tipo.

Estacoes: as 11 reguas da campanha + Benedito extras + barragens.
"""
import csv
import json
import os
import sys
import urllib.request

BASE = 'https://www.ana.gov.br/hidrowebservice'
ESTACOES = ['83800002', '83690000', '83300200', '83440000', '83345000',
            '83250000', '83105000', '83050000', '83900000', '83660000',
            '83664000', '83677000', '83680000', '83675000', '83880000',
            '83030000', '83029900', '83145000', '83145100']


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def http_json(url, headers):
    req = urllib.request.Request(url, headers=headers)
    import ssl
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE          # certificado gov...
    with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
        return json.loads(r.read().decode('utf-8', 'replace'))


def main(argv):
    pasta = _arg(argv, '--saida', 'doc/perfis_ana')
    os.makedirs(pasta, exist_ok=True)
    cred = 'doc/ana_credenciais.txt'
    if not os.path.exists(cred):
        raise SystemExit(
            f'falta {cred}: 1a linha = Identificador, 2a = Senha '
            '(do cadastro hidrowebservice da ANA)')
    ident, senha = [l.strip() for l in
                    open(cred, encoding='utf-8').read().split('\n')
                    if l.strip()][:2]

    tok = http_json(BASE + '/EstacoesTelemetricas/OAUth/v1',
                    {'Identificador': ident, 'Senha': senha})
    token = (tok.get('items') or {}).get('tokenautenticacao') \
        or tok.get('tokenautenticacao')
    if not token:
        raise SystemExit(f'autenticacao falhou: {tok}')
    print('token OK')
    H = {'Authorization': f'Bearer {token}'}

    SERVICOS = [
        ('HidroSeriePerfilTransversal', 'perfil'),
        ('HidroSerieResumoDescarga', 'resumo_descarga'),
    ]
    for cod in ESTACOES:
        for servico, rot in SERVICOS:
            url = (f'{BASE}/EstacoesTelemetricas/{servico}/v1'
                   f'?C%C3%B3digo%20da%20Esta%C3%A7%C3%A3o={cod}')
            try:
                d = http_json(url, H)
            except Exception as e:
                print(f'   {cod} {rot}: FALHOU ({e})')
                continue
            itens = d.get('items') or []
            arq = os.path.join(pasta, f'{cod}_{rot}.json')
            json.dump(d, open(arq, 'w', encoding='utf-8'),
                      ensure_ascii=False)
            if itens and isinstance(itens, list):
                chaves = sorted({k for it in itens for k in it})
                with open(os.path.join(pasta, f'{cod}_{rot}.csv'),
                          'w', newline='', encoding='utf-8') as f:
                    w = csv.DictWriter(f, chaves, delimiter=';')
                    w.writeheader()
                    w.writerows(itens)
            print(f'   {cod} {rot}: {len(itens) if isinstance(itens, list) else "?"} itens')


if __name__ == '__main__':
    main(sys.argv[1:])
