"""
================================================================================
SCRIPT DIDÁTICO: Visualização da Bacia do Itajaí e das 3 Barragens no HEC-RAS
================================================================================
Objetivo: Gerar figuras ilustrativas e explicativas para você entender 
como a água se comporta na bacia do Rio Itajaí-Açu e como as 3 barragens 
(Sul, Oeste e Norte) operam para conter cheias em Rio do Sul e Blumenau.

As figuras serão salvas na pasta 'figuras/'.
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def save_fig_safely(fig, filename):
    """Salva a figura com segurança no Windows mesmo se o arquivo estiver em uso pelo visualizador."""
    output_dir = Path("figuras").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = output_dir / filename
    
    try:
        if target_path.exists():
            try:
                target_path.unlink()
            except Exception:
                pass
        fig.savefig(str(target_path), dpi=300, bbox_inches='tight')
        print(f"-> Salvo com sucesso em: {target_path}")
    except OSError:
        alt_path = output_dir / f"novo_{filename}"
        fig.savefig(str(alt_path), dpi=300, bbox_inches='tight')
        print(f"-> Salvo em: {alt_path} (o arquivo original estava aberto em outro programa).")

# Configurar estilo dos gráficos
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#cccccc'
plt.rcParams['axes.linewidth'] = 1.2

# ==============================================================================
# FIGURA 1: Esquema Didático da Rede Hidrográfica e Localização das Barragens
# ==============================================================================
print("Gerando Figura 1: Diagrama da Bacia e Barragens...")

fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

# Desenhar Rios (linhas azuis)
ax.plot([1, 4], [8, 6], color='#1f77b4', linewidth=4, label='Rio Itajaí do Sul')
ax.plot([1, 4], [4, 6], color='#17becf', linewidth=4, label='Rio Itajaí do Oeste')
ax.plot([4, 7], [6, 6], color='#004c6d', linewidth=5, label='Rio Itajaí-Açu (Principal)')
ax.plot([5, 7], [9, 6], color='#2ca02c', linewidth=4, label='Rio Itajaí do Norte')
ax.plot([7.5, 8.5], [8.5, 7], color='#9467bd', linewidth=3, label='Rio Benedito')
ax.plot([7, 12], [6, 6], color='#004c6d', linewidth=6)
ax.plot([10, 11.5], [3, 5.8], color='#e377c2', linewidth=3.5, label='Rio Itajaí-Mirim')

# Desenhar Barragens (Ícones/Pontos)
ax.plot(2.2, 7.2, 's', color='#d62728', markersize=14, zorder=5) # Barragem Sul
ax.plot(2.2, 4.8, 's', color='#d62728', markersize=14, zorder=5) # Barragem Oeste
ax.plot(6.0, 7.5, 's', color='#d62728', markersize=14, zorder=5) # Barragem Norte

# Anotações das Barragens
ax.annotate('BARRAGEM SUL\n(Ituporanga)\n5 Comportas | Cap: 110M m³', (2.2, 7.2), xytext=(1.2, 8.3),
            arrowprops=dict(facecolor='#d62728', shrink=0.08, width=1.5, headwidth=6),
            fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#ffcccc", ec="#d62728"))

ax.annotate('BARRAGEM OESTE\n(Taió)\n7 Comportas | Cap: 110M m³', (2.2, 4.8), xytext=(1.2, 3.5),
            arrowprops=dict(facecolor='#d62728', shrink=0.08, width=1.5, headwidth=6),
            fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#ffcccc", ec="#d62728"))

ax.annotate('BARRAGEM NORTE\n(José Boiteux)\n2 Comportas | Cap: 357M m³', (6.0, 7.5), xytext=(5.0, 9.2),
            arrowprops=dict(facecolor='#d62728', shrink=0.08, width=1.5, headwidth=6),
            fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="#ffffcc", ec="#d62728"))

# Desenhar Cidades Chave (Pontos Amarelos)
cidades = [
    (4.0, 6.0, 'Rio do Sul\n(Junção Sul + Oeste)', (-30, -35)),
    (7.0, 6.0, 'Ibirama / Lontras\n(Junção Norte)', (0, -35)),
    (8.5, 6.0, 'Indaial', (0, -25)),
    (9.8, 6.0, 'Blumenau\n(Ponto Crítico de Nível)', (0, 15)),
    (12.0, 6.0, 'Foz em Itajaí / Navegantes\n(Nível do Mar / Maré)', (10, -10))
]

for cx, cy, nome, offset in cidades:
    ax.plot(cx, cy, 'o', color='#ff7f0e', markersize=10, markeredgecolor='black', zorder=6)
    ax.annotate(nome, (cx, cy), xytext=offset, textcoords='offset points',
                fontsize=9, fontweight='bold', color='#333333',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc", alpha=0.9))

ax.set_title('Rede Hidrográfica e Sistema de Barragens do Vale do Itajaí (HEC-RAS)', fontsize=14, fontweight='bold', pad=15)
ax.axis('off')
plt.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.9)
save_fig_safely(fig, "figura_1_rede_de_rios_e_barragens.png")
plt.close(fig)

# ==============================================================================
# FIGURA 2: Perfil Longitudinal de Altitude (Da Nascente ao Mar)
# ==============================================================================
print("Gerando Figura 2: Perfil Longitudinal do Rio...")

fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

distancia_km = np.linspace(250, 0, 250)

elevação_fundo = np.piecewise(distancia_km, 
    [distancia_km > 150, (distancia_km <= 150) & (distancia_km > 60), distancia_km <= 60],
    [lambda d: 340 + (d - 150) * 1.1, 
     lambda d: 12 + (d - 60) * (328/90), 
     lambda d: -15 + d * (27/60)])

ax.plot(distancia_km, elevação_fundo, color='#004c6d', linewidth=2.5, label='Fundo do Rio (Calha)')
ax.fill_between(distancia_km, elevação_fundo, -30, color='#e6f2ff', alpha=0.5)

ax.axvline(x=200, color='#d62728', linestyle='--', alpha=0.7)
ax.text(200, 380, ' Barragem Sul / Oeste\n (~350m a 390m)', color='#d62728', fontweight='bold')

ax.axvline(x=180, color='#2ca02c', linestyle='--', alpha=0.7)
ax.text(180, 280, ' Barragem Norte\n (~300m)', color='#2ca02c', fontweight='bold')

ax.plot(150, 340, 'o', color='#ff7f0e', markersize=8)
ax.text(150, 310, 'Rio do Sul (340m)', fontweight='bold', ha='center')

ax.plot(60, 12, 'o', color='#ff7f0e', markersize=8)
ax.text(60, 40, 'Blumenau (12m)', fontweight='bold', ha='center')

ax.plot(0, -15, 'o', color='#ff7f0e', markersize=8)
ax.text(0, -40, 'Foz Itajaí (0m)', fontweight='bold', ha='center')

ax.set_xlabel('Distância da Foz em Itajaí (km)', fontsize=11, fontweight='bold')
ax.set_ylabel('Altitude em Relação ao Nível do Mar (m)', fontsize=11, fontweight='bold')
ax.set_title('Perfil Longitudinal de Declividade da Bacia do Itajaí', fontsize=13, fontweight='bold')
ax.grid(True, linestyle=':', alpha=0.6)
ax.set_xlim(250, 0)
save_fig_safely(fig, "figura_2_perfil_longitudinal_cotas.png")
plt.close(fig)

# ==============================================================================
# FIGURA 3: Didática de Amortecimento de Cheia pelas Barragens
# ==============================================================================
print("Gerando Figura 3: Hidrograma de Amortecimento das Comportas...")

fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

horas = np.linspace(0, 48, 100)

vazao_entrada = 500 + 3500 * np.exp(-((horas - 18)/6)**2)
vazao_comportas_fechadas = 500 + 1200 * np.exp(-((horas - 28)/10)**2)

ax.plot(horas, vazao_entrada, color='#d62728', linestyle='--', linewidth=2.5, label='Vazão da Chuva Sem Barragem (Pico: 4000 m³/s)')
ax.plot(horas, vazao_comportas_fechadas, color='#1f77b4', linewidth=3, label='Vazão Retida com Comportas Fechadas (Pico: 1700 m³/s)')

ax.fill_between(horas, vazao_entrada, vazao_comportas_fechadas, where=(vazao_entrada > vazao_comportas_fechadas),
                color='#ffcccc', alpha=0.5, label='Volume de Água Acumulado no Reservatório (Amortecimento)')

ax.set_xlabel('Tempo Decorrido na Tempestade (Horas)', fontsize=11, fontweight='bold')
ax.set_ylabel('Vazão do Rio - Q (m³/s)', fontsize=11, fontweight='bold')
ax.set_title('Efeito Prático das Barragens: Redução do Pico de Cheia em Rio do Sul e Blumenau', fontsize=13, fontweight='bold')
ax.legend(loc='upper right', frameon=True, facecolor='white')
ax.grid(True, linestyle=':', alpha=0.6)
save_fig_safely(fig, "figura_3_hidrogramas_e_amortecimento.png")
plt.close(fig)

print("\n=========================================================================")
print("SUCESSO: Todas as figuras foram salvas com segurança na pasta 'figuras/'!")
print("=========================================================================")
