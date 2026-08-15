"""
Módulo de Visualização Gráfica e Espaço-Temporal da Propagação de Cheias.
Gera gráficos estáticos de alta resolução (PNG) e interativos:
- Hidrogramas comparativos por trecho
- Análise de atenuação e atraso de pico (lag time)
- Diagrama espaço-temporal de propagação Q(x, t)
- Comparação de validação (Simulado vs Observado)
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from typing import Dict, Any, Optional

class FloodVisualizer:
    """
    Renderizador de gráficos hidrológicos e hidrodinâmicos para o modelo.
    """
    def __init__(self, output_dir: str = "output_plots"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Configuração de estilo limpo e profissional
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.edgecolor'] = '#333333'
        plt.rcParams['axes.linewidth'] = 0.8
        plt.rcParams['grid.color'] = '#e0e0e0'
        plt.rcParams['grid.linestyle'] = '--'

    def plot_reach_propagation(self, routing_results: Dict[str, Any], save_name: str = "1_propagacao_trechos.png") -> str:
        """
        Plota os hidrogramas de entrada, intermediários e saída em cada trecho.
        """
        t = routing_results['time_hours']
        outflows = routing_results['reach_outflows']
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        
        # Entrada
        q_in = routing_results['reach_inflows'][1]
        ax.plot(t, q_in, label='Entrada Montante (Trecho 1)', color='#1f77b4', linewidth=2.5, linestyle='--')
        
        # Trechos intermediários e saída
        colors = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        for rid, q_out in outflows.items():
            col = colors[(rid - 1) % len(colors)]
            is_final = (rid == max(outflows.keys()))
            lbl = f"Saída Trecho {rid}" + (" (Foz Final)" if is_final else "")
            lw = 3.0 if is_final else 1.8
            ax.plot(t, q_out, label=lbl, color=col, linewidth=lw)
            
        ax.set_title("Propagação da Onda de Cheia no Rio Itajaí-Mirim (Trecho a Trecho)", fontsize=13, fontweight='bold', pad=12)
        ax.set_xlabel("Tempo (horas)", fontsize=11)
        ax.set_ylabel("Vazão (m³/s)", fontsize=11)
        ax.grid(True, alpha=0.6)
        ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)
        
        plt.tight_layout()
        out_path = os.path.join(self.output_dir, save_name)
        plt.savefig(out_path)
        plt.close()
        return out_path

    def plot_attenuation_and_lag(self, routing_results: Dict[str, Any], save_name: str = "2_atenuacao_atraso_pico.png") -> str:
        """
        Destaca visualmente o atraso do pico (lag time), redução do pico e alargamento da onda.
        """
        t = routing_results['time_hours']
        q_in = routing_results['reach_inflows'][1]
        q_out = routing_results['final_outflow']
        m = routing_results['metrics']
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        
        ax.plot(t, q_in, label=f"Entrada (Pico: {m['peak_inflow_m3s']} m³/s em t={m['t_peak_inflow_h']}h)", color='#0284c7', linewidth=2.5)
        ax.plot(t, q_out, label=f"Saída Foz (Pico: {m['peak_outflow_m3s']} m³/s em t={m['t_peak_outflow_h']}h)", color='#dc2626', linewidth=2.5)
        
        # Linha indicando o atraso de pico (lag time)
        ax.annotate(
            f"Atraso (Lag): +{m['lag_time_hours']}h\nRedução: -{m['peak_reduction_pct']}%",
            xy=(m['t_peak_outflow_h'], m['peak_outflow_m3s']),
            xytext=(m['t_peak_inflow_h'] + (m['lag_time_hours'] / 2.0), (m['peak_inflow_m3s'] + m['peak_outflow_m3s']) / 2.0),
            arrowprops=dict(facecolor='#4b5563', shrink=0.08, width=1.5, headwidth=6),
            bbox=dict(boxstyle="round,pad=0.4", fc="#f3f4f6", ec="#9ca3af", lw=1.2),
            fontweight='bold', fontsize=10
        )
        
        ax.axvline(m['t_peak_inflow_h'], color='#0284c7', linestyle=':', alpha=0.7)
        ax.axvline(m['t_peak_outflow_h'], color='#dc2626', linestyle=':', alpha=0.7)
        
        ax.set_title("Análise Hidrológica de Atenuação, Atraso de Pico e Alargamento da Onda", fontsize=13, fontweight='bold', pad=12)
        ax.set_xlabel("Tempo (horas)", fontsize=11)
        ax.set_ylabel("Vazão (m³/s)", fontsize=11)
        ax.grid(True, alpha=0.6)
        ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)
        
        plt.tight_layout()
        out_path = os.path.join(self.output_dir, save_name)
        plt.savefig(out_path)
        plt.close()
        return out_path

    def plot_space_time_heatmap(self, routing_results: Dict[str, Any], reach_lengths_km: list, save_name: str = "3_diagrama_espaco_tempo.png") -> str:
        """
        Gera diagrama espaço-temporal da vazão Q(x, t):
        Eixo X: Tempo (h)
        Eixo Y: Posição ao longo do rio (km a jusante)
        """
        t = np.array(routing_results['time_hours'])
        outflows = routing_results['reach_outflows']
        
        # Posições acumuladas em km
        cum_km = [0.0]
        for l in reach_lengths_km:
            cum_km.append(cum_km[-1] + l)
            
        # Matriz de vazão [posicoes x tempo]
        matrix_q = [routing_results['reach_inflows'][1]]
        for rid in sorted(outflows.keys()):
            matrix_q.append(outflows[rid])
        matrix_q = np.array(matrix_q)
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        
        T, X = np.meshgrid(t, cum_km)
        contour = ax.contourf(T, X, matrix_q, levels=30, cmap='YlGnBu')
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label("Vazão Q(x, t) [m³/s]", fontsize=11)
        
        # Linhas de contorno
        ax.contour(T, X, matrix_q, levels=10, colors='white', alpha=0.3, linewidths=0.5)
        
        ax.set_title("Diagrama Espaço-Temporal da Propagação da Onda de Cheia (Rio Itajaí-Mirim)", fontsize=13, fontweight='bold', pad=12)
        ax.set_xlabel("Tempo (horas)", fontsize=11)
        ax.set_ylabel("Distância a partir de Montante (km)", fontsize=11)
        ax.grid(True, alpha=0.3, color='gray')
        
        plt.tight_layout()
        out_path = os.path.join(self.output_dir, save_name)
        plt.savefig(out_path)
        plt.close()
        return out_path

    def plot_validation_comparison(self, t_hours: list, q_sim: list, q_obs: list, metrics: Dict[str, Any], save_name: str = "4_validacao_simulado_vs_observado.png") -> str:
        """
        Gera gráfico de validação comparando hidrograma simulado x observado com tabela de métricas.
        """
        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        
        ax.plot(t_hours, q_obs, 'o-', color='#10b981', label='Vazão Observada (Brusque Centro)', linewidth=2.0, markersize=4)
        ax.plot(t_hours, q_sim, '-', color='#3b82f6', label='Vazão Simulada (Muskingum)', linewidth=2.5)
        
        # Caixa de texto com métricas
        textstr = '\n'.join((
            r'$\mathbf{Métricas\ de\ Validação:}$',
            f'RMSE: {metrics["rmse_m3s"]} m³/s',
            f'NSE: {metrics["nse"]:.3f}',
            f'Erro de Pico: {metrics["peak_error_m3s"]} m³/s ({metrics["peak_error_pct"]}%)',
            f'Erro de Horário: {metrics["t_peak_diff_h"]}h',
            f'Erro Volumétrico: {metrics["volume_error_pct"]}%'
        ))
        props = dict(boxstyle='round,pad=0.6', facecolor='#f8fafc', edgecolor='#cbd5e1', alpha=0.95)
        ax.text(0.68, 0.95, textstr, transform=ax.transAxes, fontsize=9.5, verticalalignment='top', bbox=props)
        
        ax.set_title("Validação Hidrológica: Hidrograma Simulado vs Observado em Brusque", fontsize=13, fontweight='bold', pad=12)
        ax.set_xlabel("Tempo (horas)", fontsize=11)
        ax.set_ylabel("Vazão (m³/s)", fontsize=11)
        ax.grid(True, alpha=0.6)
        ax.legend(loc='upper left', frameon=True, facecolor='white', framealpha=0.9)
        
        plt.tight_layout()
        out_path = os.path.join(self.output_dir, save_name)
        plt.savefig(out_path)
        plt.close()
        return out_path
