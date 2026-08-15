"""
Módulo de Representação da Rede Fluvial e Trechos (River & Reach).
Gerencia a topologia, parâmetros físicos e carregamento a partir de arquivos CSV.
"""

import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Optional, Any
from .muskingum import MuskingumReach
from .muskingum_cunge import MuskingumCungeReach

@dataclass
class RiverReach:
    """
    Estrutura de dados para um trecho fluvial individual.
    """
    reach_id: int
    name: str
    upstream_node: str
    downstream_node: str
    length_km: float
    slope_m_km: float
    k_hours: float
    x_param: float
    roughness_manning: float
    initial_flow_m3s: float
    status: str # REAL_GEOM_PLACEHOLDER_ROUTING | REAL_OFFICIAL_CANAL | PLACEHOLDER

    def create_muskingum_solver(self, dt_hours: float = 1.0) -> MuskingumReach:
        return MuskingumReach(
            reach_id=self.reach_id,
            name=self.name,
            k_hours=self.k_hours,
            x_param=self.x_param,
            dt_hours=dt_hours
        )

    def create_cunge_solver(self, width_m: float = 35.0, reference_q: float = 200.0, dt_hours: float = 1.0) -> MuskingumCungeReach:
        return MuskingumCungeReach(
            reach_id=self.reach_id,
            name=self.name,
            length_km=self.length_km,
            slope_m_km=self.slope_m_km,
            width_m=width_m,
            manning_n=self.roughness_manning,
            reference_q_m3s=reference_q,
            dt_hours=dt_hours
        )


class RiverNetwork:
    """
    Representa o sistema fluvial composto por múltiplos trechos conectados em série.
    """
    def __init__(self, name: str = "Rio Itajaí-Mirim"):
        self.name = name
        self.reaches: List[RiverReach] = []
        self.stations: Dict[str, Dict[str, Any]] = {}
        
    def add_reach(self, reach: RiverReach):
        self.reaches.append(reach)
        self.reaches.sort(key=lambda r: r.reach_id)
        
    @classmethod
    def from_csv(cls, reaches_csv_path: str, stations_csv_path: Optional[str] = None, name: str = "Rio Itajaí-Mirim") -> 'RiverNetwork':
        """
        Carrega a rede a partir de arquivos CSV padronizados.
        """
        net = cls(name=name)
        df_reaches = pd.read_csv(reaches_csv_path)
        
        for _, row in df_reaches.iterrows():
            reach = RiverReach(
                reach_id=int(row['reach_id']),
                name=str(row['name']),
                upstream_node=str(row['upstream_node']),
                downstream_node=str(row['downstream_node']),
                length_km=float(row['length_km']),
                slope_m_km=float(row['slope_m_km']),
                k_hours=float(row['K_hours']),
                x_param=float(row['X_param']),
                roughness_manning=float(row['roughness_manning']),
                initial_flow_m3s=float(row['initial_flow_m3s']),
                status=str(row.get('status', 'PLACEHOLDER'))
            )
            net.add_reach(reach)
            
        if stations_csv_path:
            df_stations = pd.read_csv(stations_csv_path)
            for _, row in df_stations.iterrows():
                st_id = str(row['station_id'])
                net.stations[st_id] = {
                    'name': str(row['name']),
                    'location_type': str(row['location_type']),
                    'latitude': float(row['latitude']),
                    'longitude': float(row['longitude']),
                    'reach_id': int(row['reach_id']),
                    'river_km': float(row['river_km']),
                    'status': str(row.get('status', 'PLACEHOLDER'))
                }
                
        return net

    def get_total_length_km(self) -> float:
        return sum(r.length_km for r in self.reaches)

    def summary(self) -> pd.DataFrame:
        """
        Gera um resumo tabular dos trechos e parâmetros da rede.
        """
        data = []
        for r in self.reaches:
            data.append({
                'Trecho': r.reach_id,
                'Nome': r.name,
                'De': r.upstream_node,
                'Para': r.downstream_node,
                'Comprimento (km)': r.length_km,
                'Declividade (m/km)': r.slope_m_km,
                'K (h)': r.k_hours,
                'X': r.x_param,
                'Manning (n)': r.roughness_manning,
                'Status': r.status
            })
        return pd.DataFrame(data)
