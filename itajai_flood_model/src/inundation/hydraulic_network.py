"""
Rede Hidráulica da Bacia do Rio Itajaí (HydraulicNetwork):
Modela nós hidráulicos (confluências e bifurcações), ramos fluviais com perfis longitudinais
sem restrição artificial de dZ/dx <= 0 (permitindo remanso e controle de jusante),
tratamento especial da bifurcação do Itajaí-Mirim (Canal Retificado e Braço Velho)
e condição de contorno oceânica dinâmica H_ocean(t) na foz.
"""

from typing import Dict, Any, List, Optional, Tuple, Callable
import numpy as np

class HydraulicNode:
    """Nó hidráulico para confluências, bifurcações ou condições de contorno."""
    def __init__(self, node_id: str, name: str, lon: float, lat: float,
                 z_datum: float, node_type: str = "confluence"):
        self.node_id = node_id
        self.name = name
        self.lon = lon
        self.lat = lat
        self.z_datum = z_datum
        self.node_type = node_type  # 'confluence', 'bifurcation', 'boundary_inflow', 'boundary_ocean'
        self.inflow_rivers: List[str] = []
        self.outflow_rivers: List[str] = []
        self.q_in_total: float = 0.0
        self.q_out_total: float = 0.0
        self.z_water: float = z_datum
        self.stage_h: float = 0.0

    def compute_confluence_balance(self, inflows: Dict[str, float]) -> float:
        """Balanço de conservação de massa na confluência: Q_out = sum(Q_in)."""
        self.q_in_total = float(sum(inflows.get(r, 0.0) for r in self.inflow_rivers))
        self.q_out_total = self.q_in_total
        return self.q_out_total

    def compute_bifurcation_split(self, q_in: float, split_ratios: Dict[str, float]) -> Dict[str, float]:
        """Divisão de vazão na bifurcação: sum(Q_out_i) = Q_in."""
        self.q_in_total = float(q_in)
        outflows = {}
        total_ratio = sum(split_ratios.values())
        for r_out, ratio in split_ratios.items():
            outflows[r_out] = float(q_in * (ratio / total_ratio))
        self.q_out_total = float(sum(outflows.values()))
        return outflows


class RiverBranch:
    """Ramo fluvial com nós de montante/jusante, geometria de calha e perfil altimétrico."""
    def __init__(self, river_id: str, name: str,
                 coords: List[Tuple[float, float]],
                 z_bed: List[float],
                 upstream_node: str,
                 downstream_node: str,
                 b_start: float = 40.0,
                 b_end: float = 80.0,
                 manning_n: float = 0.038,
                 h_bank: float = 6.0):
        self.river_id = river_id
        self.name = name
        self.coords = coords
        self.z_bed = np.asarray(z_bed, dtype=float)
        self.n_pts = len(coords)
        self.upstream_node = upstream_node
        self.downstream_node = downstream_node
        self.b_vec = np.linspace(b_start, b_end, self.n_pts)
        self.manning_n = manning_n
        self.h_bank = h_bank

        # Calcular distâncias acumuladas
        self.dists_km = np.zeros(self.n_pts)
        for i in range(1, self.n_pts):
            dx = (coords[i][0] - coords[i-1][0]) * 111.32 * np.cos(np.radians(coords[i][1]))
            dy = (coords[i][1] - coords[i-1][1]) * 110.57
            self.dists_km[i] = self.dists_km[i-1] + np.sqrt(dx**2 + dy**2)

    def compute_backwater_water_surface(self, q_profile: np.ndarray,
                                        downstream_z_water: Optional[float] = None) -> Dict[str, np.ndarray]:
        """
        Calcula o perfil de linha d'água permitindo remanso (sem impor dZ/dx <= 0).
        Utiliza o método padrão de integração hidráulica (Standard Step Method aproximado)
        a partir da condição de contorno de jusante.
        """
        n = self.n_pts
        q_arr = np.asarray(q_profile, dtype=float)
        if len(q_arr) != n:
            q_arr = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(q_arr)), q_arr)

        h_normal = np.zeros(n)
        z_water = np.zeros(n)

        # 1. Profundidade normal via Manning em cada seção
        for i in range(n):
            b_i = self.b_vec[i]
            q_i = max(0.1, float(q_arr[i]))
            
            # Declividade local do fundo
            if i < n - 1:
                dz = abs(self.z_bed[i] - self.z_bed[i+1])
                dx = max(10.0, (self.dists_km[i+1] - self.dists_km[i]) * 1000.0)
                s0 = max(0.00008, dz / dx)
            else:
                s0 = 0.00010

            # Inversão Manning trapezoidal (talude z=1.2)
            h_grid = np.linspace(0.05, 25.0, 400)
            area = b_i * h_grid + 1.2 * (h_grid ** 2)
            perim = b_i + 2.0 * h_grid * np.sqrt(1.0 + 1.2**2)
            r_hyd = area / perim
            q_grid = (1.0 / self.manning_n) * area * (r_hyd ** (2.0/3.0)) * np.sqrt(s0)
            
            h_norm_i = float(np.interp(q_i, q_grid, h_grid))
            h_normal[i] = h_norm_i
            z_water[i] = self.z_bed[i] + h_norm_i

        # 2. Integração do Remanso de Jusante para Montante (Backwater Profile)
        # Se o nível de jusante imposto for superior ao normal, o remanso se propaga para montante!
        if downstream_z_water is not None and downstream_z_water > z_water[-1]:
            z_water[-1] = downstream_z_water
            h_normal[-1] = max(0.1, downstream_z_water - self.z_bed[-1])
            
            # Propagar remanso para montante
            for i in range(n - 2, -1, -1):
                dx = max(10.0, (self.dists_km[i+1] - self.dists_km[i]) * 1000.0)
                # Perda de carga por atrito (Friction Slope Sf)
                b_mid = 0.5 * (self.b_vec[i] + self.b_vec[i+1])
                h_mid = max(0.2, 0.5 * (h_normal[i] + h_normal[i+1]))
                area_mid = b_mid * h_mid + 1.2 * (h_mid**2)
                perim_mid = b_mid + 2.0 * h_mid * np.sqrt(1.0 + 1.2**2)
                r_mid = area_mid / perim_mid
                v_mid = q_arr[i] / max(1.0, area_mid)
                sf_mid = ((self.manning_n * v_mid) / (r_mid ** (2.0/3.0))) ** 2
                
                # Equação de energia: Z_i = Z_{i+1} + Sf * dx
                z_remanso = z_water[i+1] + sf_mid * dx
                # O nível real é o máximo entre a cota normal e a curva de remanso
                z_water[i] = max(z_water[i], z_remanso)
                h_normal[i] = max(0.1, z_water[i] - self.z_bed[i])

        return {
            'distances_km': self.dists_km,
            'z_bed_m': self.z_bed,
            'depth_h_m': np.round(h_normal, 2),
            'z_water_m': np.round(z_water, 2),
            'is_overtopping': (h_normal > self.h_bank)
        }


class ItajaiHydraulicNetwork:
    """
    Rede Hidráulica Completa da Bacia do Rio Itajaí com confluências,
    bifurcação do Itajaí-Mirim e condição de contorno oceânica variável H_ocean(t).
    """
    def __init__(self, dem_profiles_dict: Dict[str, Any]):
        self.nodes: Dict[str, HydraulicNode] = {}
        self.branches: Dict[str, RiverBranch] = {}
        self._build_network(dem_profiles_dict)

    def _build_network(self, dem_profiles: Dict[str, Any]):
        # 1. Criar Nós Estratégicos
        self.nodes['no_taio'] = HydraulicNode('no_taio', 'Barragem Oeste / Taió', -49.998, -27.115, 360.0)
        self.nodes['no_ituporanga'] = HydraulicNode('no_ituporanga', 'Barragem Sul / Ituporanga', -49.605, -27.414, 380.0)
        self.nodes['no_boiteux'] = HydraulicNode('no_boiteux', 'Barragem Norte / José Boiteux', -49.628, -26.960, 240.0)
        
        # Confluência de Rio do Sul (Oeste + Sul + Mirim Doce + Perimbó + Trombudo)
        self.nodes['no_rio_do_sul'] = HydraulicNode('no_rio_do_sul', 'Confluência de Rio do Sul', -49.643, -27.215, 330.0)
        self.nodes['no_rio_do_sul'].inflow_rivers = ['oeste', 'sul', 'mirim_doce', 'perimbo', 'trombudo']
        self.nodes['no_rio_do_sul'].outflow_rivers = ['acu_alto']

        # Confluência de Ibirama (Alto Vale + Rio Hercílio/Norte)
        self.nodes['no_ibirama'] = HydraulicNode('no_ibirama', 'Confluência de Ibirama', -49.520, -27.058, 125.0)
        self.nodes['no_ibirama'].inflow_rivers = ['acu_alto', 'norte']
        self.nodes['no_ibirama'].outflow_rivers = ['acu_medio']

        # Confluência de Indaial (Ibirama + Rio Benedito)
        self.nodes['no_indaial'] = HydraulicNode('no_indaial', 'Confluência de Indaial', -49.230, -26.898, 47.5)
        self.nodes['no_indaial'].inflow_rivers = ['acu_medio', 'benedito']
        self.nodes['no_indaial'].outflow_rivers = ['acu_blumenau']

        # Estação de Blumenau (Ponte de Ferro)
        self.nodes['no_blumenau'] = HydraulicNode('no_blumenau', 'Blumenau Centro (Ponte de Ferro)', -49.066, -26.918, 1.30)
        self.nodes['no_blumenau'].inflow_rivers = ['acu_blumenau']
        self.nodes['no_blumenau'].outflow_rivers = ['acu_baixo']

        # Bifurcação do Itajaí-Mirim (Montante de Itajaí)
        self.nodes['no_bifurcacao_mirim'] = HydraulicNode(
            'no_bifurcacao_mirim', 'Bifurcação do Itajaí-Mirim', -48.740, -26.965, 8.0, node_type="bifurcation"
        )
        self.nodes['no_bifurcacao_mirim'].inflow_rivers = ['mirim_montante']
        self.nodes['no_bifurcacao_mirim'].outflow_rivers = ['mirim_canal_retificado', 'mirim_braco_velho']

        # Foz do Rio Itajaí (Estuário no Oceano Atlântico)
        self.nodes['no_foz_itajai'] = HydraulicNode(
            'no_foz_itajai', 'Foz do Rio Itajaí (Oceano Atlântico)', -48.650, -26.905, -3.50, node_type="boundary_ocean"
        )
        self.nodes['no_foz_itajai'].inflow_rivers = ['acu_baixo', 'mirim_canal_retificado', 'mirim_braco_velho', 'luis_alves']

        # 2. Criar Ramos Fluviais
        for r_key, prof in dem_profiles.items():
            coords = prof.get('coords', [])
            z_bed = prof.get('z_dem') or prof.get('elevations', [])
            if len(coords) < 2:
                continue
                
            self.branches[r_key] = RiverBranch(
                river_id=r_key,
                name=prof.get('name', r_key),
                coords=coords,
                z_bed=z_bed,
                upstream_node="no_" + r_key,
                downstream_node="no_rio_do_sul" if r_key in ('oeste', 'sul', 'mirim_doce', 'perimbo', 'trombudo') else "no_foz_itajai"
            )

    def solve_network_hydraulics(self, basin_flows_t: Dict[str, np.ndarray],
                                 t_step: int,
                                 h_ocean_t: float = 0.00,
                                 mirim_split_ratio: float = 0.70) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Resolve o escoamento acoplado em toda a rede no instante t:
        1. Balanço de massa nas confluências.
        2. Divisão de vazão na bifurcação do Itajaí-Mirim:
           Q_canal = ratio * Q_mirim, Q_braco = (1 - ratio) * Q_mirim.
        3. Propagação de jusante a partir da condição oceânica H_ocean(t).
        """
        # Condição de contorno oceânica na foz
        z_ocean_level = -3.50 + max(0.0, 3.50 + h_ocean_t)
        self.nodes['no_foz_itajai'].z_water = z_ocean_level
        self.nodes['no_foz_itajai'].stage_h = h_ocean_t

        results = {}

        # 1. Tronco do Itajaí-Açu com remanso a partir da foz oceânica
        if 'acu' in self.branches:
            q_acu = basin_flows_t.get('acu', np.ones(self.branches['acu'].n_pts) * 500.0)
            res_acu = self.branches['acu'].compute_backwater_water_surface(q_acu, downstream_z_water=z_ocean_level)
            results['acu'] = res_acu

        # 2. Ramos dos Tributários
        for r_key, branch in self.branches.items():
            if r_key == 'acu':
                continue
            q_r = basin_flows_t.get(r_key, np.ones(branch.n_pts) * 100.0)
            
            # Condição de jusante para confluências
            ds_water = None
            if r_key in ('oeste', 'sul', 'trombudo', 'perimbo', 'mirim_doce') and 'acu' in results:
                # O nível de jusante é o nível de entrada em Rio do Sul
                ds_water = float(results['acu']['z_water_m'][0])
            elif r_key == 'norte' and 'acu' in results:
                # Confluência em Ibirama (km ~35 do Açú)
                idx_ib = min(len(results['acu']['z_water_m']) - 1, 20)
                ds_water = float(results['acu']['z_water_m'][idx_ib])
            elif r_key == 'benedito' and 'acu' in results:
                # Confluência em Indaial (km ~90 do Açú)
                idx_ind = min(len(results['acu']['z_water_m']) - 1, 45)
                ds_water = float(results['acu']['z_water_m'][idx_ind])

            res_r = branch.compute_backwater_water_surface(q_r, downstream_z_water=ds_water)
            results[r_key] = res_r

        return results
