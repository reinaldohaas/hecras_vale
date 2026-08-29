# -*- coding: utf-8 -*-
"""
Status de andamento com tempo restante.

Um passo que le 118 GB nao pode ficar calado. Sem estimativa nao ha como
distinguir "esta demorando" de "travou", e a unica saida e matar o processo --
perdendo tudo que ja foi lido.

A estimativa e por taxa observada, nao por media desde o inicio: a leitura dos
tiles do SIG-SC nao e uniforme (tile na borda da bacia entra quase vazio, tile
no meio entra inteiro), e a media global fica presa no comeco. A janela movel
acompanha.

O que se ve:

    [4] secoes   372/1415   26%   decorrido 0:41   restante ~1:58   9,1/s
        Itajai_Mirim RS 106.460

Numa unica linha, reescrita no lugar quando o terminal permite; em log de
arquivo, uma linha a cada intervalo, para nao encher o arquivo de milhares de
linhas quase iguais.
"""
import sys
import time


def hms(seg):
    """Duracao legivel. Sem casas decimais: ninguem decide nada com elas."""
    if seg is None or seg != seg or seg < 0:
        return "--:--"
    seg = int(seg)
    h, r = divmod(seg, 3600)
    m, s = divmod(r, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


def tamanho(bytes_):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if abs(bytes_) < 1024.0 or u == "TB":
            return f"{bytes_:.1f} {u}" if u != "B" else f"{bytes_:.0f} B"
        bytes_ /= 1024.0
    return ""


class Progresso:
    """Contador com tempo restante, para laco de N itens conhecidos.

    total=None serve para trabalho de tamanho desconhecido: mostra o que ja
    passou e a taxa, sem inventar um restante que seria mentira.
    """

    def __init__(self, total, rotulo="", cada=1.0, log=None, janela=40):
        self.total = total
        self.rotulo = rotulo
        self.cada = cada
        self.log = log
        self.n = 0
        self.t0 = time.time()
        self.ultimo = 0.0
        self.marcas = [(self.t0, 0)]
        self.janela = janela
        self.tty = getattr(sys.stdout, "isatty", lambda: False)()

    # ------------------------------------------------------------- interno
    def _taxa(self):
        """Itens por segundo na janela movel."""
        if len(self.marcas) < 2:
            return 0.0
        t_a, n_a = self.marcas[0]
        t_b, n_b = self.marcas[-1]
        dt = t_b - t_a
        return (n_b - n_a) / dt if dt > 0 else 0.0

    def restante(self):
        if not self.total:
            return None
        taxa = self._taxa()
        if taxa <= 0:
            return None
        return max(self.total - self.n, 0) / taxa

    def _linha(self, extra=""):
        dec = time.time() - self.t0
        taxa = self._taxa()
        if self.total:
            pct = 100.0 * self.n / self.total
            corpo = (f"{self.n}/{self.total} {pct:5.1f}%  "
                     f"decorrido {hms(dec)}  restante ~{hms(self.restante())}")
        else:
            corpo = f"{self.n}  decorrido {hms(dec)}"
        if taxa >= 1.0:
            corpo += f"  {taxa:.1f}/s"
        elif taxa > 0:
            corpo += f"  {1.0/taxa:.1f} s cada"
        return f"   {self.rotulo} {corpo}" + (f"   {extra}" if extra else "")

    def _mostrar(self, extra="", forcar=False):
        agora = time.time()
        if not forcar and (agora - self.ultimo) < self.cada:
            return
        self.ultimo = agora
        linha = self._linha(extra)
        if self.log is not None:
            self.log(linha)
        elif self.tty:
            sys.stdout.write("\r" + linha.ljust(110)[:110])
            sys.stdout.flush()
        else:
            print(linha, flush=True)

    # ------------------------------------------------------------- publico
    def passo(self, n=1, extra=""):
        self.n += n
        self.marcas.append((time.time(), self.n))
        if len(self.marcas) > self.janela:
            del self.marcas[0]
        self._mostrar(extra)
        return self

    def fracao(self, f, extra=""):
        """Para trabalho que reporta fracao (0..1), como o GDAL."""
        if self.total:
            self.n = int(round(f * self.total))
        self.marcas.append((time.time(), self.n))
        if len(self.marcas) > self.janela:
            del self.marcas[0]
        self._mostrar(extra)
        return self

    def fim(self, extra=""):
        dec = time.time() - self.t0
        if self.log is None and self.tty:
            sys.stdout.write("\r" + " " * 110 + "\r")
            sys.stdout.flush()
        msg = f"   {self.rotulo} {self.n} em {hms(dec)}"
        (self.log or print)(msg + (f"   {extra}" if extra else ""))
        return dec


def callback_gdal(prog):
    """Adaptador para o callback de progresso do GDAL.

    Assinatura do GDAL: fn(fracao, mensagem, dado). Devolver False cancela; aqui
    nunca se cancela, so se reporta.
    """
    def fn(fracao, _msg="", _dado=None):
        prog.fracao(float(fracao))
        return 1
    return fn


class Etapas:
    """Cronometro do conjunto de passos, com estimativa do que falta.

    Precisa de um peso por passo, senao a estimativa fica sem sentido: no
    modelo do Vale o terreno pode levar horas e a escrita segundos, e tratar os
    dez passos como iguais daria um restante inutil.
    """

    def __init__(self, pesos, log=print):
        self.pesos = dict(pesos)
        self.total = sum(self.pesos.values()) or 1.0
        self.feito = 0.0
        self.t0 = time.time()
        self.log = log
        self.tempos = {}

    def inicia(self, nome):
        self._t = time.time()
        self._nome = nome
        peso = self.pesos.get(nome, 1.0)
        falta = self.total - self.feito
        dec = time.time() - self.t0
        est = (dec / self.feito * falta) if self.feito > 0 else None
        self.log(f"   [{100*self.feito/self.total:3.0f}% do conjunto]"
                 + (f"   restante estimado ~{hms(est)}" if est else ""))
        return peso

    def termina(self):
        dt = time.time() - self._t
        self.tempos[self._nome] = dt
        self.feito += self.pesos.get(self._nome, 1.0)
        return dt

    def resumo(self):
        L = [f"tempo total {hms(time.time() - self.t0)}"]
        for k, v in sorted(self.tempos.items(), key=lambda x: -x[1]):
            L.append(f"   {k:<12}{hms(v)}")
        return "\n".join(L)
