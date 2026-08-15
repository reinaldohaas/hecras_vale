"""
Módulo de Superfície Raster 2D da Linha d'Água e Delimitação da Mancha (WaterSurfaceRasterEngine):
1. Interpolação da superfície da água Z_water(x,y,t) condicionada ao Corredor Fluvial (River Corridor Mask).
2. Cálculo da profundidade 2D: depth(x,y,t) = max(0, Z_water - Z_DEM).
3. Separação explícita entre Mapa Geométrico (DEM < H) e Mapa Hidraulicamente Conectado (Seed Fill BFS).
"""

from collections import deque
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

class WaterSurfaceRasterEngine:
    """
    Constrói a grade 2D de superfície livre d'água e cruza com a topografia do DEM Copernicus.
    """
    def __init__(self, bounds: Tuple[float, float, float, float] = (-50.2, -27.55, -48.55, -26.55),
                 grid_shape: Tuple[int, int] = (160, 220)):
        """
        bounds: (min_lon, min_lat, max_lon, max_lat) cobrindo toda a bacia do Itajaí
        grid_shape: (nrows, ncols)
        """
        self.min_lon, self.min_lat, self.max_lon, self.max_lat = bounds
        self.nrows, self.ncols = grid_shape
        
        self.lons = np.linspace(self.min_lon, self.max_lon, self.ncols)
        self.lats = np.linspace(self.min_lat, self.max_lat, self.nrows)
        self.lon_grid, self.lat_grid = np.meshgrid(self.lons, self.lats)
        
        # Dimensões da célula em metros (aproximadamente 30m para sub-regiões ou ~800m para bacia inteira)
        self.dx_m = (self.max_lon - self.min_lon) * 111320.0 * np.cos(np.radians(np.mean(self.lats))) / self.ncols
        self.dy_m = (self.max_lat - self.min_lat) * 110570.0 / self.nrows
        self.cell_area_km2 = (self.dx_m * self.dy_m) / 1e6

        self.z_dem = np.zeros(grid_shape, dtype=float)
        self.river_corridor_mask = np.zeros(grid_shape, dtype=bool)
        self.river_channel_seed_mask = np.zeros(grid_shape, dtype=bool)

    def load_dem_matrix(self, z_dem_matrix: np.ndarray):
        """Define o modelo digital de elevação."""
        if z_dem_matrix.shape != (self.nrows, self.ncols):
            # Interpolar para casar dimensões
            from scipy.ndimage import zoom
            zoom_r = self.nrows / z_dem_matrix.shape[0]
            zoom_c = self.ncols / z_dem_matrix.shape[1]
            self.z_dem = zoom(z_dem_matrix, (zoom_r, zoom_c), order=1)
        else:
            self.z_dem = np.asarray(z_dem_matrix, dtype=float)

    def build_synthetic_valley_dem(self):
        """Constrói um DEM realista da Bacia do Itajaí com declividades de montante a jusante."""
        # Gradiente regional Alto Vale (350m em Taió) -> Médio Vale (11m em Blumenau) -> Foz (0m em Itajaí)
        frac_lon = (self.lon_grid - self.min_lon) / (self.max_lon - self.min_lon) # 0 (Oeste) a 1 (Leste)
        frac_lat = (self.lat_grid - self.min_lat) / (self.max_lat - self.min_lat)
        
        # Calha principal e relevo de serras
        valley_axis_lat = -27.20 + 0.30 * frac_lon # Eixo de Rio do Sul a Itajaí
        dist_from_axis = np.abs(self.lat_grid - valley_axis_lat)
        
        # Elevação de fundo decrescente
        z_valley = 360.0 * (1.0 - frac_lon) ** 1.3 + 1.30
        # Encostas em "V" abrindo para planície larga em Blumenau/Gaspar
        valley_width_deg = 0.04 + 0.12 * frac_lon
        z_hills = 400.0 * np.minimum(1.0, (dist_from_axis / valley_width_deg) ** 1.5)
        
        self.z_dem = np.round(z_valley + z_hills, 1)

    def rasterize_river_corridor(self, river_branches: Dict[str, Any], max_corridor_km: float = 3.5):
        """
        Cria a máscara do corredor de inundação (River Corridor) ao redor das calhas dos rios
        para impedir interpolação de água fora dos vales hidraulicamente conectados.
        """
        self.river_corridor_mask.fill(False)
        self.river_channel_seed_mask.fill(False)
        
        max_corridor_deg = max_corridor_km / 111.0

        for r_key, branch in river_branches.items():
            coords = getattr(branch, 'coords', [])
            for lon, lat in coords:
                # Célula central do rio
                c = int(np.clip(np.searchsorted(self.lons, lon), 0, self.ncols - 1))
                r = int(np.clip(np.searchsorted(self.lats, lat), 0, self.nrows - 1))
                self.river_channel_seed_mask[r, c] = True
                
                # Raio do corredor do vale
                r_rad = int(np.ceil(max_corridor_deg / (self.dy_m / 110570.0)))
                c_rad = int(np.ceil(max_corridor_deg / (self.dx_m / 111320.0)))
                
                r_min, r_max = max(0, r - r_rad), min(self.nrows, r + r_rad + 1)
                c_min, c_max = max(0, c - c_rad), min(self.ncols, c + c_rad + 1)
                self.river_corridor_mask[r_min:r_max, c_min:c_max] = True

    def interpolate_2d_water_surface(self, river_water_levels: Dict[str, Dict[str, np.ndarray]]) -> np.ndarray:
        """
        Interpola as cotas absolutas Z_water(x,t) dos nós fluviais para a grade 2D
        condicionada estritamente ao corredor de vale.
        """
        z_water_grid = np.zeros((self.nrows, self.ncols), dtype=float)
        
        sample_pts = []
        sample_vals = []

        for r_key, res in river_water_levels.items():
            z_w = res.get('z_water_m', [])
            branch_coords = res.get('coords', [])
            if len(branch_coords) == 0 and r_key in ('acu', 'oeste', 'sul', 'norte', 'benedito', 'mirim', 'luis_alves'):
                # Usar coordenadas padrão da bacia
                pass

            for i in range(min(len(z_w), len(branch_coords))):
                lon, lat = branch_coords[i]
                sample_pts.append((lon, lat))
                sample_vals.append(float(z_w[i]))

        if len(sample_pts) < 3:
            # Fallback para gradiente padrão se amostras não fornecidas diretamente
            z_water_grid = 345.0 * (1.0 - (self.lon_grid - self.min_lon)/(self.max_lon - self.min_lon)) + 15.34
            return z_water_grid

        # Interpolação IDW (Inverse Distance Weighting) no corredor fluvial
        sample_pts_arr = np.array(sample_pts)
        sample_vals_arr = np.array(sample_vals)
        
        corridor_coords = np.argwhere(self.river_corridor_mask)
        for r_idx, c_idx in corridor_coords:
            pt_lon = self.lons[c_idx]
            pt_lat = self.lats[r_idx]
            
            dists = np.sqrt((sample_pts_arr[:, 0] - pt_lon)**2 + (sample_pts_arr[:, 1] - pt_lat)**2) + 1e-6
            weights = 1.0 / (dists ** 2.0)
            weights /= np.sum(weights)
            z_water_grid[r_idx, c_idx] = float(np.sum(weights * sample_vals_arr))

        return z_water_grid

    def compute_2d_inundation(self, z_water_2d: np.ndarray,
                              min_depth_m: float = 0.05) -> Dict[str, Any]:
        """
        Calcula a lâmina d'água 2D e produz:
        1. Mapa Geométrico Bruto: h_geom = max(0, Z_water - Z_DEM)
        2. Mapa Hidraulicamente Conectado: h_conn (via Flood-Fill BFS a partir do rio)
        """
        # 1. Mapa Geométrico (dentro do corredor)
        h_geom = np.where(
            self.river_corridor_mask,
            np.maximum(0.0, z_water_2d - self.z_dem),
            0.0
        )
        
        # 2. Conectividade Hidráulica 2D (Flood-Fill a partir das sementes da calha fluvial)
        is_wet = (h_geom >= min_depth_m)
        connected_mask = np.zeros((self.nrows, self.ncols), dtype=bool)
        visited = np.zeros((self.nrows, self.ncols), dtype=bool)
        
        queue = deque()
        seeds = np.argwhere(self.river_channel_seed_mask & is_wet)
        for r, c in seeds:
            queue.append((int(r), int(c)))
            visited[r, c] = True
            connected_mask[r, c] = True

        # Se calha não estiver molhada, usar todas as células de canal
        if len(queue) == 0:
            seeds = np.argwhere(self.river_channel_seed_mask)
            for r, c in seeds:
                if is_wet[r, c]:
                    queue.append((int(r), int(c)))
                    visited[r, c] = True
                    connected_mask[r, c] = True

        # BFS 8 vizinhos
        neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
        while queue:
            curr_r, curr_c = queue.popleft()
            for dr, dc in neighbors:
                nr, nc = curr_r + dr, curr_c + dc
                if 0 <= nr < self.nrows and 0 <= nc < self.ncols:
                    if not visited[nr, nc] and is_wet[nr, nc]:
                        visited[nr, nc] = True
                        connected_mask[nr, nc] = True
                        queue.append((nr, nc))

        h_connected = np.where(connected_mask, h_geom, 0.0)

        # 3. Métricas
        area_geom_km2 = float(np.sum(h_geom >= min_depth_m) * self.cell_area_km2)
        area_conn_km2 = float(np.sum(h_connected >= min_depth_m) * self.cell_area_km2)
        vol_conn_hm3 = float(np.sum(h_connected) * (self.dx_m * self.dy_m) / 1e6)
        max_depth_m = float(np.max(h_connected)) if area_conn_km2 > 0 else 0.0

        return {
            'geometric_depth_m': np.round(h_geom, 2),
            'connected_depth_m': np.round(h_connected, 2),
            'area_geometric_km2': round(area_geom_km2, 2),
            'area_connected_km2': round(area_conn_km2, 2),
            'volume_connected_hm3': round(vol_conn_hm3, 2),
            'max_depth_m': round(max_depth_m, 2),
            'isolated_pits_area_km2': round(max(0.0, area_geom_km2 - area_conn_km2), 2)
        }
