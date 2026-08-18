# -*- coding: utf-8 -*-
"""
Interface PySide6.

Divisao de responsabilidades: esta camada nao calcula nada. Ela carrega os
arquivos, chama qc/correction e desenha. Todo o julgamento sobre a geometria
mora nos modulos de analise, que rodam identicos no modo lote -- se a interface
e o lote discordassem, nao daria para confiar em nenhum dos dois.
"""
import os
import traceback

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Canvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as Toolbar
from matplotlib.figure import Figure
from PySide6 import QtCore, QtGui, QtWidgets

from . import correction, cross_sections, export, plotting, qc
from .dem import DEM
from .river_axis import EixoRio

FILTROS = ("todas", qc.OK, qc.ATENCAO, qc.INCERTO, qc.CRITICA)


class Tela(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HEC-RAS QC - controle de qualidade de secoes")
        self.resize(1500, 950)

        self.dem = None
        self.eixo = None
        self.secoes = []
        self.proposta = None
        self.candidatas = []
        self.desfazer = []          # pilha (posicao, secao_anterior)
        self.lim = qc.Limiares()
        self.log = None

        self._montar()

    # ------------------------------------------------------------- montagem
    def _montar(self):
        self._barra()

        self.fig_mapa = Figure(figsize=(6, 6))
        self.ax_mapa = self.fig_mapa.add_subplot(111)
        self.canvas_mapa = Canvas(self.fig_mapa)
        self.canvas_mapa.mpl_connect("button_press_event", self._clique_mapa)

        self.fig_perfil = Figure(figsize=(6, 3.2))
        self.ax_perfil = self.fig_perfil.add_subplot(111)
        self.canvas_perfil = Canvas(self.fig_perfil)

        self.fig_comp = Figure(figsize=(6, 3.0))
        self.ax_comp = self.fig_comp.add_subplot(111)
        self.canvas_comp = Canvas(self.fig_comp)

        esq = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(esq)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(Toolbar(self.canvas_mapa, self))
        lv.addWidget(self.canvas_mapa)

        dir_ = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        dir_.addWidget(self.canvas_perfil)
        dir_.addWidget(self.canvas_comp)
        dir_.addWidget(self._painel_info())
        dir_.setSizes([340, 260, 300])

        cima = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        cima.addWidget(esq)
        cima.addWidget(dir_)
        cima.setSizes([760, 740])

        tudo = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        tudo.addWidget(cima)
        tudo.addWidget(self._painel_tabela())
        tudo.setSizes([650, 300])
        self.setCentralWidget(tudo)
        self.statusBar().showMessage("abra o DEM, o eixo do rio e as secoes")

    def _barra(self):
        b = self.addToolBar("arquivos")
        b.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        for txt, fn in (("Abrir DEM", self.abrir_dem),
                        ("Abrir eixo do rio", self.abrir_eixo),
                        ("Abrir secoes", self.abrir_secoes)):
            a = QtGui.QAction(txt, self)
            a.triggered.connect(fn)
            b.addAction(a)
        b.addSeparator()
        a = QtGui.QAction("Analisar todas", self)
        a.triggered.connect(self.analisar_todas)
        b.addAction(a)
        a = QtGui.QAction("Desfazer", self)
        a.setShortcut("Ctrl+Z")
        a.triggered.connect(self.desfazer_ultima)
        b.addAction(a)
        b.addSeparator()
        for txt, fn in (("Exportar vetor", self.exportar_vetor),
                        ("Exportar perfis CSV", self.exportar_csv),
                        ("Exportar tabela QC", self.exportar_tabela)):
            a = QtGui.QAction(txt, self)
            a.triggered.connect(fn)
            b.addAction(a)

    def _painel_info(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)

        self.txt = QtWidgets.QTextEdit()
        self.txt.setReadOnly(True)
        self.txt.setFont(QtGui.QFont("Consolas", 9))
        lay.addWidget(self.txt, 3)

        cx = QtWidgets.QHBoxLayout()
        for txt, fn in (("Estender 30 m", lambda: self.estender(30.0)),
                        ("Estender 50 m", lambda: self.estender(50.0)),
                        ("Estender 100 m", lambda: self.estender(100.0)),
                        ("Perpendicular ao rio", self.gerar_perpendicular)):
            bt = QtWidgets.QPushButton(txt)
            bt.clicked.connect(fn)
            cx.addWidget(bt)
        lay.addLayout(cx)

        cy = QtWidgets.QHBoxLayout()
        self.bt_aceitar = QtWidgets.QPushButton("Aceitar proposta")
        self.bt_aceitar.clicked.connect(self.aceitar)
        self.bt_aceitar.setEnabled(False)
        bt_manter = QtWidgets.QPushButton("Manter original")
        bt_manter.clicked.connect(self.manter)
        cy.addWidget(self.bt_aceitar)
        cy.addWidget(bt_manter)
        lay.addLayout(cy)

        lay.addWidget(self._painel_limiares())
        return w

    def _painel_limiares(self):
        cx = QtWidgets.QGroupBox("limiares de QC")
        g = QtWidgets.QGridLayout(cx)
        self.campos = {}
        campos = [("pos_ok_min", "OK de", 0.0, 1.0, 0.01),
                  ("pos_ok_max", "ate", 0.0, 1.0, 0.01),
                  ("pos_atencao_min", "atencao de", 0.0, 1.0, 0.01),
                  ("pos_atencao_max", "ate", 0.0, 1.0, 0.01),
                  ("profundidade_min", "prof. min (m)", 0.0, 100.0, 0.1),
                  ("proeminencia_min", "proeminencia (m)", 0.0, 100.0, 0.1),
                  ("salto_talvegue", "salto talv. (m)", 0.0, 500.0, 0.5),
                  ("razao_largura", "razao largura", 1.1, 20.0, 0.1),
                  ("desvio_ortogonal", "desvio ang. (g)", 0.0, 90.0, 1.0),
                  ("espacamento", "amostragem (m)", 0.2, 50.0, 0.5)]
        for k, (nome, rot, mn, mx, passo) in enumerate(campos):
            sp = QtWidgets.QDoubleSpinBox()
            sp.setRange(mn, mx)
            sp.setSingleStep(passo)
            sp.setDecimals(2)
            sp.setValue(float(getattr(self.lim, nome)))
            self.campos[nome] = sp
            g.addWidget(QtWidgets.QLabel(rot), k // 2, 2 * (k % 2))
            g.addWidget(sp, k // 2, 2 * (k % 2) + 1)
        bt = QtWidgets.QPushButton("aplicar e reavaliar")
        bt.clicked.connect(self.aplicar_limiares)
        g.addWidget(bt, len(campos) // 2 + 1, 0, 1, 4)
        return cx

    def _painel_tabela(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("filtrar:"))
        self.cb_filtro = QtWidgets.QComboBox()
        self.cb_filtro.addItems(FILTROS)
        self.cb_filtro.currentTextChanged.connect(lambda _: self.encher_tabela())
        top.addWidget(self.cb_filtro)
        self.lb_contagem = QtWidgets.QLabel("")
        top.addWidget(self.lb_contagem, 1)
        lay.addLayout(top)

        self.tabela = QtWidgets.QTableWidget()
        cols = ["RS", "rio", "largura", "talvegue %", "prof. rel.",
                "orientacao", "QC", "status", "origem", "motivos"]
        self.tabela.setColumnCount(len(cols))
        self.tabela.setHorizontalHeaderLabels(cols)
        self.tabela.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tabela.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tabela.itemSelectionChanged.connect(self._selecao_tabela)
        self.tabela.horizontalHeader().setStretchLastSection(True)
        lay.addWidget(self.tabela)
        return w

    # --------------------------------------------------------------- abrir
    def _erro(self, e):
        traceback.print_exc()
        QtWidgets.QMessageBox.critical(self, "erro", str(e))

    def abrir_dem(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "DEM (GeoTIFF)", "", "Raster (*.tif *.tiff *.vrt);;todos (*)")
        if not p:
            return
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            self.dem = DEM(p)
        except Exception as e:                      # noqa: BLE001
            self._erro(e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        r = self.dem.resumo()
        self.statusBar().showMessage(
            f"DEM {r['dimensoes']} @ {r['resolucao']} m | {r['crs']} | "
            f"cotas {r['cotas']} | NoData {r['nodata']}")
        self.desenhar_mapa()

    def abrir_eixo(self):
        if self.dem is None:
            return self._aviso("abra o DEM primeiro: o CRS de trabalho vem dele")
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "eixo do rio", "", "vetor (*.shp *.geojson *.json *.gpkg)")
        if not p:
            return
        try:
            self.eixo = EixoRio.ler(p, self.dem.crs_metrico)
        except Exception as e:                      # noqa: BLE001
            return self._erro(e)
        self.statusBar().showMessage(
            f"eixo: {len(self.eixo.linhas)} linha(s), "
            f"{self.eixo.comprimento/1000:.1f} km")
        self.desenhar_mapa()

    def abrir_secoes(self):
        if self.dem is None:
            return self._aviso("abra o DEM primeiro")
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "secoes transversais", "",
            "vetor (*.shp *.geojson *.json *.gpkg)")
        if not p:
            return
        try:
            self.secoes = cross_sections.carregar(p, self.dem.crs_metrico,
                                                  self.eixo)
        except Exception as e:                      # noqa: BLE001
            return self._erro(e)
        self.log = os.path.splitext(p)[0] + "_qc_alteracoes.log"
        self.statusBar().showMessage(f"{len(self.secoes)} secoes carregadas")
        self.analisar_todas()

    def _aviso(self, msg):
        QtWidgets.QMessageBox.warning(self, "atencao", msg)

    # ------------------------------------------------------------- analise
    def _ler_limiares(self):
        for k, sp in self.campos.items():
            setattr(self.lim, k, float(sp.value()))
        return self.lim

    def aplicar_limiares(self):
        self._ler_limiares()
        self.analisar_todas()

    def analisar_todas(self):
        if not self.secoes or self.dem is None:
            return
        self._ler_limiares()
        prog = QtWidgets.QProgressDialog("extraindo perfis do terreno...",
                                         "cancelar", 0, len(self.secoes), self)
        prog.setWindowModality(QtCore.Qt.WindowModal)
        try:
            for i, s in enumerate(self.secoes):
                if prog.wasCanceled():
                    break
                s.extrair(self.dem, self.lim.espacamento, self.eixo,
                          self.lim.proeminencia_min)
                prog.setValue(i + 1)
        except Exception as e:                      # noqa: BLE001
            return self._erro(e)
        finally:
            prog.close()
        qc.avaliar_todas(self.secoes, self.lim)
        self.encher_tabela()
        self.desenhar_mapa()

    # -------------------------------------------------------------- tabela
    def encher_tabela(self):
        f = self.cb_filtro.currentText()
        vis = [s for s in self.secoes
               if f == "todas" or (s.qc and s.qc.status == f)]
        self._visiveis = vis
        self.tabela.setRowCount(len(vis))
        for r, s in enumerate(vis):
            def val(x, fmt="{:.1f}"):
                return fmt.format(x) if x is not None and np.isfinite(x) else "-"
            dados = [f"{s.rs:g}" if s.rs is not None else f"#{s.idx}",
                     s.rio,
                     val(s.largura, "{:.0f}"),
                     val(100 * s.posicao_relativa, "{:.1f}"),
                     val(s.profundidade_relativa, "{:.2f}"),
                     val(s.azimute, "{:.0f}"),
                     f"{s.qc.nota:.0f}" if s.qc else "-",
                     s.qc.status if s.qc else "-",
                     s.origem,
                     s.qc.resumo if s.qc else ""]
            for c, v in enumerate(dados):
                it = QtWidgets.QTableWidgetItem(str(v))
                if s.qc and c == 7:
                    it.setForeground(QtGui.QColor(
                        plotting.cor_status(s.qc.status)))
                    fonte = it.font()
                    fonte.setBold(True)
                    it.setFont(fonte)
                self.tabela.setItem(r, c, it)
        self.tabela.resizeColumnsToContents()
        c = qc.contagem(self.secoes)
        self.lb_contagem.setText(
            f"   OK {c[qc.OK]}   atencao {c[qc.ATENCAO]}   "
            f"incerto {c[qc.INCERTO]}   critica {c[qc.CRITICA]}"
            f"   (de {len(self.secoes)})")

    def _selecao_tabela(self):
        r = self.tabela.currentRow()
        if 0 <= r < len(getattr(self, "_visiveis", [])):
            self.selecionar(self._visiveis[r])

    def _clique_mapa(self, ev):
        if ev.xdata is None or not self.secoes:
            return
        from shapely.geometry import Point
        p = Point(ev.xdata, ev.ydata)
        s = min(self.secoes, key=lambda x: x.geom.distance(p))
        self.selecionar(s)

    # ------------------------------------------------------------ selecao
    def selecionar(self, s):
        self.atual = s
        self.proposta = None
        self.candidatas = []
        self.bt_aceitar.setEnabled(False)
        plotting.perfil(self.ax_perfil, s)
        self.canvas_perfil.draw_idle()
        self.ax_comp.clear()
        self.canvas_comp.draw_idle()
        self.desenhar_mapa(zoom=True)
        self.mostrar_info(s)

    def mostrar_info(self, s, proposta=None):
        t = s.talvegue or {}
        L = []
        L.append(f"{s.rotulo}")
        L.append("=" * 58)
        L.append(f"origem            : {s.origem}")
        L.append(f"largura           : {s.largura:.1f} m")
        L.append(f"pontos amostrados : {len(s.sta) if s.sta is not None else 0}"
                 f"  (espacamento {self.lim.espacamento:.1f} m)")
        if s.valida:
            nan = int((~np.isfinite(s.z)).sum())
            L.append(f"NoData no perfil  : {nan} ponto(s)")
            L.append(f"cota minima       : {np.nanmin(s.z):.2f} m")
            L.append(f"talvegue          : {s.z_talvegue:.2f} m em "
                     f"{s.sta[s.i_talvegue]:.1f} m")
            L.append(f"posicao relativa  : {100*s.posicao_relativa:.1f} % "
                     f"da largura")
            L.append(f"dist. a esquerda  : {s.dist_margem_esq:.1f} m")
            L.append(f"dist. a direita   : {s.dist_margem_dir:.1f} m")
            L.append(f"profundidade rel. : {s.profundidade_relativa:.2f} m")
            L.append(f"proeminencia      : {t.get('proeminencia', 0):.2f} m")
            ia, it = t.get("i_min_abs"), t.get("i_talvegue")
            if ia is not None and ia != it:
                L.append(f"minimo absoluto   : {s.z[ia]:.2f} m em "
                         f"{s.sta[ia]:.1f} m  (NAO e o talvegue escolhido)")
            if s.azimute is not None:
                L.append(f"orientacao        : {s.azimute:.1f} graus do eixo")
            if s.sta_eixo is not None:
                L.append(f"cruzamento c/eixo : {s.sta_eixo:.1f} m")
            else:
                L.append("cruzamento c/eixo : a secao NAO cruza o eixo")
        if s.qc:
            L.append("-" * 58)
            L.append(f"QC {s.qc.status}   nota {s.qc.nota:.0f}/100")
            for k in "ABCDE":
                d = s.qc.testes.get(k, {})
                if d:
                    marca = " " if d["status"] == qc.OK else ">"
                    L.append(f" {marca}{k}  {d['status']:<8} "
                             f"{d['motivo'] or 'ok'}")
            if s.qc.status in (qc.CRITICA, qc.INCERTO):
                L.append("")
                L.append("SECAO PROBLEMATICA -- use os botoes de correcao "
                         "abaixo.")
        if proposta is not None and proposta.qc:
            L.append("=" * 58)
            L.append(f"PROPOSTA: {proposta.origem}")
            L.append(f"  largura      {proposta.largura:.1f} m")
            L.append(f"  talvegue     {100*proposta.posicao_relativa:.1f} %")
            L.append(f"  profundidade {proposta.profundidade_relativa:.2f} m")
            L.append(f"  QC original {s.qc.nota:.0f}/100   ->   "
                     f"QC proposta {proposta.qc.nota:.0f}/100 "
                     f"[{proposta.qc.status}]")
        self.txt.setPlainText("\n".join(L))

    def desenhar_mapa(self, zoom=False):
        if self.dem is None:
            return
        lim = (self.ax_mapa.get_xlim(), self.ax_mapa.get_ylim())
        tinha = self.ax_mapa.has_data()
        plotting.mapa(self.ax_mapa, self.dem, self.eixo, self.secoes,
                      getattr(self, "atual", None))
        if zoom and getattr(self, "atual", None) is not None:
            plotting.zoom_secao(self.ax_mapa, self.atual)
        elif tinha:
            self.ax_mapa.set_xlim(*lim[0])
            self.ax_mapa.set_ylim(*lim[1])
        self.canvas_mapa.draw_idle()

    # ------------------------------------------------------------ correcao
    def _exige_selecao(self):
        if getattr(self, "atual", None) is None:
            self._aviso("selecione uma secao na tabela ou no mapa")
            return False
        return True

    def estender(self, d):
        if not self._exige_selecao():
            return
        self._ler_limiares()
        nova, cand = correction.estender(self.atual, self.dem, self.eixo,
                                         self.lim, passos=(d,))
        self._propor(nova, cand)

    def gerar_perpendicular(self):
        if not self._exige_selecao():
            return
        if self.eixo is None:
            return self._aviso("abra o eixo do rio para gerar a perpendicular")
        self._ler_limiares()
        nova = correction.perpendicular(self.atual, self.dem, self.eixo,
                                        self.lim)
        self._propor(nova, [nova] if nova else [])

    def _propor(self, nova, candidatas):
        if nova is None:
            return self._aviso("nao foi possivel gerar proposta")
        self.proposta = nova
        self.candidatas = candidatas
        plotting.comparacao(self.ax_comp, self.atual, nova)
        self.canvas_comp.draw_idle()
        self.mostrar_info(self.atual, nova)
        self.bt_aceitar.setEnabled(True)

    def aceitar(self):
        """Substitui a secao na camada de trabalho. O arquivo original nao e
        tocado -- a troca so chega ao disco quando o usuario exporta."""
        if self.proposta is None:
            return
        i = self.secoes.index(self.atual)
        self.desfazer.append((i, self.secoes[i]))
        if self.log:
            export.registrar(self.log, self.atual, self.proposta, True)
        self.secoes[i] = self.proposta
        qc.avaliar_todas(self.secoes, self.lim)
        self.atual = self.secoes[i]
        self.proposta = None
        self.bt_aceitar.setEnabled(False)
        self.encher_tabela()
        self.selecionar(self.atual)
        self.statusBar().showMessage(
            f"secao substituida; registrado em {os.path.basename(self.log)}"
            if self.log else "secao substituida")

    def manter(self):
        if self.proposta is not None and self.log:
            export.registrar(self.log, self.atual, self.proposta, False)
        self.proposta = None
        self.bt_aceitar.setEnabled(False)
        self.ax_comp.clear()
        self.canvas_comp.draw_idle()

    def desfazer_ultima(self):
        if not self.desfazer:
            return self.statusBar().showMessage("nada a desfazer")
        i, antiga = self.desfazer.pop()
        self.secoes[i] = antiga
        qc.avaliar_todas(self.secoes, self.lim)
        self.encher_tabela()
        self.selecionar(antiga)
        self.statusBar().showMessage("alteracao desfeita")

    # ----------------------------------------------------------- exportar
    def _destino(self, titulo, filtro, sufixo):
        p, _ = QtWidgets.QFileDialog.getSaveFileName(self, titulo, sufixo, filtro)
        return p

    def exportar_vetor(self):
        if not self.secoes:
            return
        p = self._destino("exportar secoes", "GeoJSON (*.geojson);;"
                          "Shapefile (*.shp)", "secoes_qc.geojson")
        if not p:
            return
        try:
            saidas = export.exportar_vetor(self.secoes, p,
                                           self.dem.crs_metrico)
        except Exception as e:                      # noqa: BLE001
            return self._erro(e)
        self.statusBar().showMessage("gravado: " + ", ".join(
            os.path.basename(s) for s in saidas))

    def exportar_csv(self):
        if not self.secoes:
            return
        p = self._destino("exportar perfis", "CSV (*.csv)", "perfis_hecras.csv")
        if p:
            self.statusBar().showMessage(
                "gravado: " + os.path.basename(
                    export.exportar_csv_perfis(self.secoes, p)))

    def exportar_tabela(self):
        if not self.secoes:
            return
        p = self._destino("exportar tabela QC", "CSV (*.csv)", "tabela_qc.csv")
        if p:
            self.statusBar().showMessage(
                "gravado: " + os.path.basename(
                    export.exportar_tabela(self.secoes, p)))


def rodar():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    t = Tela()
    t.show()
    return app.exec()
