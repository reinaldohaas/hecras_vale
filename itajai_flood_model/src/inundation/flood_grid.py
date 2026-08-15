"""
Módulo de Grade Topográfica e Filtro de Conectividade Hidráulica 2D (Flood-Fill):
Calcula a lâmina d'água h(x, y) = max(0, Z_water - Z_DEM) e aplica o algoritmo
de componentes conexos partindo das sementes da calha fluvial para eliminar depressões isoladas.
"""

from collections import deque
from typing import Dict, Any, List, Optional, Tuple, Set
import numpy as np

class HydraulicConnectivityFilter:
    """
    Executa o filtro de conectividade hidráulica 2D sobre uma matriz de lâmina d'água.
    Apenas células que possuem um caminho contínuo de células inundadas (h > 0)
    até pelo menos uma célula de calha fluvial (semente do rio) são mantidas.
    """

    @staticmethod
    def apply_connectivity_filter(depth_grid: np.ndarray,
                                  river_seed_mask: np.ndarray,
                                  min_depth_threshold: float = 0.05) -> np.ndarray:
        """
        depth_grid: matriz 2D (nrows, ncols) de lâminas d'água calculadas h = Z_water - Z_dem
        river_seed_mask: matriz 2D booleana (nrows, ncols) indicando onde passa a calha do rio
        min_depth_threshold: lâmina mínima para ser considerada inundada (ex: 5 cm)
        
        Retorna:
        filtered_depth_grid: matriz 2D onde apenas as áreas conectadas ao rio permanecem com h > 0.
        """
        nrows, ncols = depth_grid.shape
        is_wet = (depth_grid >= min_depth_threshold)
        
        visited = np.zeros((nrows, ncols), dtype=bool)
        connected_mask = np.zeros((nrows, ncols), dtype=bool)
        
        queue = deque()
        
        # Inserir todas as células de rio que estão molhadas como sementes iniciais
        seed_coords = np.argwhere(river_seed_mask & is_wet)
        for r, c in seed_coords:
            queue.append((int(r), int(c)))
            visited[r, c] = True
            connected_mask[r, c] = True

        # Se nenhuma semente de rio estiver molhada, mas houver células de rio, usar as células de rio como sementes
        if len(queue) == 0:
            seed_coords_all = np.argwhere(river_seed_mask)
            for r, c in seed_coords_all:
                if is_wet[r, c]:
                    queue.append((int(r), int(c)))
                    visited[r, c] = True
                    connected_mask[r, c] = True

        # Busca em Largura (BFS 2D - 8 vizinhos)
        neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
        
        while queue:
            curr_r, curr_c = queue.popleft()
            
            for dr, dc in neighbors:
                nr, nc = curr_r + dr, curr_c + dc
                if 0 <= nr < nrows and 0 <= nc < ncols:
                    if not visited[nr, nc] and is_wet[nr, nc]:
                        visited[nr, nc] = True
                        connected_mask[nr, nc] = True
                        queue.append((nr, nc))
                        
        filtered_depth = np.where(connected_mask, depth_grid, 0.0)
        return filtered_depth


class InundationGrid:
    """
    Representa uma grade espacial georreferenciada para cálculo de inundação.
    """
    def __init__(self, bounds: Tuple[float, float, float, float],
                 cell_size_m: float = 30.0,
                 nrows: int = 100, ncols: int = 100):
        """
        bounds: (min_lon, min_lat, max_lon, max_lat)
        cell_size_m: resolução da célula (ex: 30m para Copernicus DEM)
        """
        self.min_lon, self.min_lat, self.max_lon, self.max_lat = bounds
        self.cell_size_m = cell_size_m
        self.nrows = nrows
        self.ncols = ncols
        
        self.lons = np.linspace(self.min_lon, self.max_lon, ncols)
        self.lats = np.linspace(self.min_lat, self.max_lat, nrows)
        
        self.z_dem = np.zeros((nrows, ncols), dtype=float)
        self.river_mask = np.zeros((nrows, ncols), dtype=bool)

    def set_topography(self, z_dem_matrix: np.ndarray):
        """Define a matriz de elevações da grade."""
        if z_dem_matrix.shape != (self.nrows, self.ncols):
            raise ValueError(f"Dimensões incompatíveis: esperado {(self.nrows, self.ncols)}, recebido {z_dem_matrix.shape}")
        self.z_dem = z_dem_matrix.astype(float)

    def rasterize_river_line(self, coords: List[Tuple[float, float]]):
        """Marca as células por onde a calha do rio passa como sementes."""
        for lon, lat in coords:
            # Encontrar célula mais próxima
            c = int(np.clip(np.searchsorted(self.lons, lon), 0, self.ncols - 1))
            r = int(np.clip(np.searchsorted(self.lats, lat), 0, self.nrows - 1))
            self.river_mask[r, c] = True
            # Dilatar 1 célula ao redor para robustez
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < self.nrows and 0 <= nc < self.ncols:
                        self.river_mask[nr, nc] = True

    def compute_flood_depths(self, z_water_grid: np.ndarray,
                             apply_connectivity: bool = True) -> Dict[str, Any]:
        """
        Calcula as profundidades da lâmina d'água h(x,y) e métricas associadas.
        """
        raw_depths = np.maximum(0.0, z_water_grid - self.z_dem)
        
        if apply_connectivity:
            depths = HydraulicConnectivityFilter.apply_connectivity_filter(
                raw_depths, self.river_mask, min_depth_threshold=0.05
            )
        else:
            depths = raw_depths

        # Métricas
        wet_cells = np.sum(depths > 0.05)
        cell_area_km2 = (self.cell_size_m * self.cell_size_m) / 1e6
        total_area_km2 = float(wet_cells * cell_area_km2)
        
        total_volume_hm3 = float(np.sum(depths) * (self.cell_size_m * self.cell_size_m) / 1e6)
        max_depth_m = float(np.max(depths)) if wet_cells > 0 else 0.0
        mean_depth_m = float(np.mean(depths[depths > 0.05])) if wet_cells > 0 else 0.0

        return {
            'depths_grid': depths,
            'flooded_area_km2': np.round(total_area_km2, 2),
            'flooded_volume_hm3': np.round(total_volume_hm3, 2),
            'max_depth_m': np.round(max_depth_m, 2),
            'mean_depth_m': np.round(mean_depth_m, 2),
            'wet_cells_count': int(wet_cells)
        }
