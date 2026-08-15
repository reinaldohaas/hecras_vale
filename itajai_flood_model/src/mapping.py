"""
Módulo de Mapeamento Geográfico dos Trechos Fluviais (RiverMapper).
Desacoplado do mecanismo matemático. Funciona 100% offline a partir de dados locais.
"""

import matplotlib.pyplot as plt
import pandas as pd
import json
import os
from typing import Dict, Any, Optional

class RiverMapper:
    """
    Gera mapas estáticos e esquemáticos dos trechos e estações do Rio Itajaí-Mirim.
    """
    def __init__(self, stations_csv_path: str, reaches_csv_path: str, geojson_path: Optional[str] = None):
        self.df_stations = pd.read_csv(stations_csv_path)
        self.df_reaches = pd.read_csv(reaches_csv_path)
        self.geojson_path = geojson_path
        
    def plot_river_map(self, output_path: str = "output_plots/5_mapa_trechos_itajai_mirim.png") -> str:
        """
        Plota o mapa georreferenciado dos trechos do Itajaí-Mirim e estações.
        """
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig, ax = plt.subplots(figsize=(9, 8), dpi=150)
        
        # Plotar estações e linha esquemática do rio
        lons = self.df_stations['longitude'].values
        lats = self.df_stations['latitude'].values
        names = self.df_stations['name'].values
        
        # Conectar nós em ordem
        ax.plot(lons, lats, color='#0284c7', linewidth=3.5, linestyle='-', label='Calha do Rio Itajaí-Mirim', zorder=2)
        
        # Destacar Canal Retificado (último trecho)
        if len(lons) >= 2:
            ax.plot(lons[-2:], lats[-2:], color='#dc2626', linewidth=4.0, linestyle='--', label='Canal Retificado Oficial (Alívio de Cheia)', zorder=3)
            
        # Plotar marcadores das estações
        ax.scatter(lons, lats, color='#f59e0b', edgecolor='#111827', s=90, zorder=4, label='Estações / Nós de Controle')
        
        # Rótulos das estações
        for lo, la, nm in zip(lons, lats, names):
            ax.annotate(
                nm,
                xy=(lo, la),
                xytext=(8, 4),
                textcoords='offset points',
                fontsize=8.5,
                fontweight='bold',
                color='#1f2937',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#cbd5e1', alpha=0.85)
            )
            
        ax.set_title("Esquema Geográfico dos Trechos do Rio Itajaí-Mirim (Vidal Ramos → Itajaí Foz)", fontsize=11.5, fontweight='bold', pad=12)
        ax.set_xlabel("Longitude (graus)", fontsize=10)
        ax.set_ylabel("Latitude (graus)", fontsize=10)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.9, fontsize=8.5)
        
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()
        return output_path
