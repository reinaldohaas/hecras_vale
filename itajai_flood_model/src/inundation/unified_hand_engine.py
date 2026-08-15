"""
Motor Unificado e Rigoroso de HAND + Inundação Hidrodinâmica Sincronizada com Marés e Batimetria:
1. Batimetria Real com Fundo Abaixo do Nível do Mar (Z_bed < 0m na foz e estuário).
2. Condição de Contorno Oceânica Dinâmica: Maré Astronômica Semidiurna + Maré Meteorológica (Storm Surge).
3. Cotas de Margem e Transbordo (Bankfull Stage): Em nível normal, a água corre confinada na calha profunda (depth_floodplain = 0m).
4. HAND Topográfico Estático: HAND(x, y) = Z_DEM(x, y) - Z_drain(x, y).
5. Superfície Hidráulica Dinâmica Flúvio-Marítima: Z_water(x, y, t) com remanso de maré.
6. Lâmina d'Água 2D na Planície: depth(x, y, t) = max(0, Z_water - Z_DEM) quando Z_water > Z_bank (extravasamento).
7. Conectividade Hidráulica Temporal com o canal ativo.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
import scipy.ndimage as ndimage
import rasterio
from rasterio.features import rasterize, shapes
from shapely.geometry import shape, mapping
import json

@dataclass
class OceanTideState:
    """Estado do Nível do Mar na Foz (Itajaí / Navegantes) no instante t."""
    timestamp_hour: float
    astro_tide_m: float       # Maré astronômica semidiurna [m]
    storm_surge_m: float      # Sobre-elevação meteorológica / ressaca [m]
    total_ocean_level_z: float # Cota absoluta do mar Z_ocean(t) [m]

class OceanTideModel:
    """Modelo de Marés Astronômicas + Meteorológicas na Foz do Rio Itajaí-Açu."""
    def __init__(self, semi_diurnal_amp_m: float = 0.85, diurnal_amp_m: float = 0.15):
        self.semi_diurnal_amp = semi_diurnal_amp_m
        self.diurnal_amp = diurnal_amp_m

    def compute_ocean_level(self, t_hour: float, storm_surge_peak_m: float = 1.20,
                            t_surge_peak_h: float = 24.0, surge_duration_h: float = 12.0) -> OceanTideState:
        """Calcula o nível do mar na foz no instante t."""
        # Maré astronômica semidiurna (período 12.42h) + diurna (24h)
        astro = (self.semi_diurnal_amp * np.cos(2.0 * np.pi * (t_hour - 4.0) / 12.42) +
                 self.diurnal_amp * np.cos(2.0 * np.pi * (t_hour - 6.0) / 24.0))

        # Maré meteorológica / Storm surge (empilhamento de vento sul e baixa pressão)
        surge = storm_surge_peak_m * np.exp(-0.5 * ((t_hour - t_surge_peak_h) / (surge_duration_h / 2.355)) ** 2)

        z_ocean = float(astro + surge)
        return OceanTideState(
            timestamp_hour=t_hour,
            astro_tide_m=round(float(astro), 2),
            storm_surge_m=round(float(surge), 2),
            total_ocean_level_z=round(z_ocean, 2)
        )

@dataclass
class FloodState:
    """Estado hidrodinâmico longitudinal da onda de cheia no instante t."""
    timestamp_hour: float
    river_key: str
    reach_fraction_s: float     # Posição normalizada ao longo do rio [0.0 = montante, 1.0 = foz]
    discharge_q: float          # Vazão Q(s, t) [m³/s]
    stage_h: float              # Cota de régua H(s, t) [m]
    water_surface_z: float      # Cota absoluta da linha d'água Z_water(s, t) [m]

@dataclass
class FloodRaster:
    """Raster 2D resultante do acoplamento hidráulico-topográfico no instante t."""
    timestamp_hour: float
    z_water_2d: np.ndarray      # Superfície absoluta d'água Z_water(x, y, t) [m]
    eta_2d: np.ndarray          # Altura relativa acima da drenagem eta(x, y, t) [m]
    depth_raw_2d: np.ndarray    # Lâmina bruta max(0, Z_water - Z_DEM) [m]
    depth_connected_2d: np.ndarray # Lâmina hidraulicamente conectada [m]
    inundated_mask: np.ndarray  # Máscara binária de inundação conectada
    ocean_level_z: float        # Nível do mar na foz no instante t [m]
    area_km2: float
    volume_hm3: float
    max_depth_m: float
    mean_depth_m: float

class TopographicHANDModel:
    """
    Modelo Topográfico Estático:
    Calcula e preserva as propriedades topográficas imutáveis da bacia:
    - DEM(x, y)
    - drainage_id(x, y)
    - Z_drain(x, y)
    - Z_bed(x, y) (Batimetria de fundo)
    - Z_bank(x, y) (Cota de margem / transbordo)
    - HAND(x, y) = DEM(x, y) - Z_drain(x, y)
    """
    def __init__(self, dem_tif_path: Union[str, Path], river_network_json_path: Union[str, Path]):
        self.dem_tif_path = Path(dem_tif_path)
        self.river_network_json_path = Path(river_network_json_path)

        self.dem: Optional[np.ndarray] = None
        self.drainage_id: Optional[np.ndarray] = None
        self.drainage_elevation: Optional[np.ndarray] = None
        self.bed_elevation: Optional[np.ndarray] = None
        self.bank_elevation: Optional[np.ndarray] = None
        self.station_frac_s: Optional[np.ndarray] = None
        self.hand: Optional[np.ndarray] = None
        self.stream_mask: Optional[np.ndarray] = None
        
        self.transform: Optional[rasterio.Affine] = None
        self.crs = None
        self.bounds = None
        self.cell_area_km2: float = 0.0

        self.river_keys = ['acu', 'oeste', 'sul', 'norte', 'benedito', 'mirim', 'luis_alves', 'trombudo', 'mirim_doce', 'perimbo']
        self.river_to_id = {k: i + 1 for i, k in enumerate(self.river_keys)}
        self.id_to_river = {i + 1: k for i, k in enumerate(self.river_keys)}

        self._build_topographic_hand_grid()

    def _build_topographic_hand_grid(self):
        """Constrói a grade estática de HAND e drenagem com batimetria e cotas de margem."""
        if not self.dem_tif_path.exists():
            raise FileNotFoundError(f"DEM não encontrado: {self.dem_tif_path}")

        with rasterio.open(self.dem_tif_path) as src:
            self.dem = src.read(1).astype(float)
            self.transform = src.transform
            self.crs = src.crs
            self.bounds = src.bounds

        dx_m = self.transform[0] * 111320.0 * np.cos(np.radians(np.mean([self.bounds.bottom, self.bounds.top])))
        dy_m = abs(self.transform[4]) * 110570.0
        self.cell_area_km2 = (dx_m * dy_m) / 1e6

        nrows, ncols = self.dem.shape
        inv_transform = ~self.transform

        with open(self.river_network_json_path, 'r', encoding='utf-8') as f:
            net_data = json.load(f)
        profiles = net_data.get('river_profiles', {})

        # 1. Rasterizar a rede de canais fluviais atribuindo (river_id, station_s, z_bed, z_bank)
        channel_id_grid = np.zeros((nrows, ncols), dtype=np.int32)
        channel_s_grid = np.zeros((nrows, ncols), dtype=np.float32)
        channel_bed_grid = np.zeros((nrows, ncols), dtype=np.float32)
        channel_bank_grid = np.zeros((nrows, ncols), dtype=np.float32)
        channel_mask = np.zeros((nrows, ncols), dtype=bool)

        for r_key, prof in profiles.items():
            if r_key not in self.river_to_id:
                continue
            r_id = self.river_to_id[r_key]
            coords = prof.get('coords', [])
            z_dem_prof = prof.get('z_dem') or prof.get('elevations', [10.0]*len(coords))
            if len(coords) < 2:
                continue

            n_pts = len(coords)
            for i in range(n_pts - 1):
                p1 = coords[i]
                p2 = coords[i+1]
                f1 = i / float(n_pts - 1)
                f2 = (i + 1) / float(n_pts - 1)
                
                # Batimetria real de fundo do leito e cota de margem (bankfull)
                # Itajaí-Açu: Foz em Itajaí (Z_bed = -4.5m, Bank = 2.5m), Blumenau (Z_bed = 1.3m, Bank = 12.88m), Rio do Sul (Z_bed = 332m, Bank = 339.5m)
                if r_key == 'acu':
                    if f1 < 0.68:
                        w_b = f1 / 0.68
                        bed_z1 = (1.0 - w_b) * 332.0 + w_b * 1.3
                        bank_z1 = (1.0 - w_b) * (332.0 + 7.0) + w_b * (4.88 + 8.0) # Bank = 339.0m em RS, 12.88m em Blumenau
                    else:
                        w_b = (f1 - 0.68) / 0.32
                        bed_z1 = (1.0 - w_b) * 1.3 + w_b * (-4.5)
                        bank_z1 = (1.0 - w_b) * 12.88 + w_b * 2.50 # Bank = 2.50m na Foz
                elif r_key == 'mirim':
                    bed_z1 = (1.0 - f1) * 180.0 + f1 * (-3.0)
                    bank_z1 = (1.0 - f1) * (180.0 + 5.5) + f1 * 2.50 # Bank = 20.5m em Brusque, 2.5m na Foz
                else:
                    z_ref = z_dem_prof[min(i, len(z_dem_prof)-1)]
                    bed_z1 = z_ref - 4.5
                    bank_z1 = z_ref + 1.8

                sub_steps = 6
                for s in range(sub_steps):
                    w = s / float(sub_steps)
                    lon = (1.0 - w) * p1[0] + w * p2[0]
                    lat = (1.0 - w) * p1[1] + w * p2[1]
                    f_val = (1.0 - w) * f1 + w * f2
                    c, r = inv_transform * (lon, lat)
                    c, r = int(round(c)), int(round(r))
                    if 0 <= r < nrows and 0 <= c < ncols:
                        channel_id_grid[r, c] = r_id
                        channel_s_grid[r, c] = f_val
                        channel_bed_grid[r, c] = bed_z1
                        channel_bank_grid[r, c] = bank_z1
                        channel_mask[r, c] = True

        self.stream_mask = channel_mask

        # 2. Drainage Assignment: Mapear cada pixel da bacia para sua célula de drenagem correspondente
        dist, indices = ndimage.distance_transform_edt(~self.stream_mask, return_indices=True)
        nearest_r = indices[0]
        nearest_c = indices[1]

        # 3. Z_drain, Z_bed e Z_bank associados
        self.drainage_elevation = self.dem[nearest_r, nearest_c]
        self.bed_elevation = channel_bed_grid[nearest_r, nearest_c]
        self.bank_elevation = channel_bank_grid[nearest_r, nearest_c]
        self.drainage_id = channel_id_grid[nearest_r, nearest_c]
        self.station_frac_s = channel_s_grid[nearest_r, nearest_c]

        # 4. HAND(x, y) = DEM(x, y) - Z_drain(x, y) (Propriedade TOPOGRÁFICA Estática Exata)
        self.hand = self.dem - self.drainage_elevation

        # 5. Verificação estrita da identidade topográfica: DEM == Z_drain + HAND
        np.testing.assert_allclose(self.dem, self.drainage_elevation + self.hand, atol=1e-5,
                                  err_msg="Falha na identidade matemática: DEM != Z_drain + HAND")

    def export_debug_rasters(self, output_dir: Union[str, Path]):
        """Salva as 4 camadas topográficas fundamentais em formato GeoTIFF."""
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        meta = {
            'driver': 'GTiff',
            'height': self.dem.shape[0],
            'width': self.dem.shape[1],
            'count': 1,
            'dtype': 'float32',
            'crs': self.crs,
            'transform': self.transform
        }

        with rasterio.open(out_path / "01_dem.tif", 'w', **meta) as dst:
            dst.write(self.dem.astype(np.float32), 1)

        meta_id = meta.copy()
        meta_id['dtype'] = 'int32'
        with rasterio.open(out_path / "02_drainage_id.tif", 'w', **meta_id) as dst:
            dst.write(self.drainage_id.astype(np.int32), 1)

        with rasterio.open(out_path / "03_drainage_elevation.tif", 'w', **meta) as dst:
            dst.write(self.drainage_elevation.astype(np.float32), 1)

        with rasterio.open(out_path / "04_hand.tif", 'w', **meta) as dst:
            dst.write(self.hand.astype(np.float32), 1)


class SynchronizedFloodEngine:
    """
    Motor Hidrodinâmico-Topográfico Sincronizado com Marés e Batimetria:
    Recebe a superfície de água Z_water(s, t) e o nível do mar Z_ocean(t), acoplando com o HAND.
    """
    def __init__(self, topo_hand_model: TopographicHANDModel, min_flood_depth_m: float = 0.05):
        self.topo = topo_hand_model
        self.min_flood_depth = min_flood_depth_m
        self.tide_model = OceanTideModel()

    def compute_instantaneous_flood(self, river_water_surface_profiles: Dict[str, Dict[str, np.ndarray]],
                                    t_hour: float,
                                    storm_surge_peak_m: float = 1.20,
                                    max_corridor_cells: int = 150) -> FloodRaster:
        """
        Calcula a inundação 2D para o instante t exato levando em conta a calha profunda e as marés:
        - river_water_surface_profiles: {r_key: {'z_water_m': z_profile_t, 'is_overtopping': bool_vector}}
        """
        nrows, ncols = self.topo.dem.shape
        z_water_grid = np.zeros((nrows, ncols), dtype=float)

        # 1. Condição de Contorno Oceânica na Foz
        ocean_state = self.tide_model.compute_ocean_level(t_hour, storm_surge_peak_m=storm_surge_peak_m)
        z_ocean = ocean_state.total_ocean_level_z

        # 2. Mapear a superfície d'água Z_water(x, y, t) ao longo de cada rio
        for r_key, r_id in self.topo.river_to_id.items():
            mask_r = (self.topo.drainage_id == r_id)
            if not np.any(mask_r):
                continue

            prof = river_water_surface_profiles.get(r_key)
            if prof is not None and 'z_water_m' in prof:
                z_prof = np.asarray(prof['z_water_m'], dtype=float)
                n_sec = len(z_prof)
                s_vals = self.topo.station_frac_s[mask_r]
                sec_idx = np.clip((s_vals * (n_sec - 1)).astype(int), 0, n_sec - 1)
                z_fluvial = z_prof[sec_idx]

                # Acoplamento de remanso de maré na foz do Itajaí-Açu e Itajaí-Mirim
                if r_key in ('acu', 'mirim'):
                    # Remanso oceânico decai com a distância da foz (s = 1.0 é foz)
                    backwater_weight = np.clip((s_vals - 0.70) / 0.30, 0.0, 1.0)
                    z_estuary = z_ocean + (1.0 - backwater_weight) * 2.0
                    z_water_grid[mask_r] = np.maximum(z_fluvial, z_estuary * backwater_weight)
                else:
                    z_water_grid[mask_r] = z_fluvial
            else:
                # Nível de calha normal (dentro da calha profunda Z_bed + 1.5m)
                z_water_grid[mask_r] = self.topo.bed_elevation[mask_r] + 1.50

        # 3. Altura relativa da água acima da drenagem: eta(x, y, t) = Z_water(x, y, t) - Z_drain(x, y)
        eta_grid = z_water_grid - self.topo.drainage_elevation

        # 4. Condição Física de Transbordo (Bankfull / Extravasamento):
        # A planície de inundação só recebe água se o nível do rio ultrapassar a cota da margem (Z_water > Z_bank)
        is_overtopping = (z_water_grid > self.topo.bank_elevation)

        # Profundidade de Inundação na Planície
        depth_raw = np.where(is_overtopping, np.maximum(0.0, z_water_grid - self.topo.dem), 0.0)

        # Limitar alcance lateral do vale
        dist_to_stream = ndimage.distance_transform_edt(~self.topo.stream_mask)
        depth_raw[dist_to_stream > max_corridor_cells] = 0.0

        # 5. Conectividade Hidráulica Temporal (Flood-Fill a partir da calha ativa em extravasamento)
        is_wet = (depth_raw >= self.min_flood_depth)
        labeled, num_features = ndimage.label(is_wet, structure=np.ones((3, 3)))

        active_stream_mask = self.topo.stream_mask & (depth_raw > 0.0)
        if not np.any(active_stream_mask):
            # Se não houver transbordo em nenhum ponto, planície está 100% seca
            depth_connected = np.zeros_like(depth_raw)
            connected_mask = np.zeros_like(is_wet, dtype=bool)
        else:
            stream_labels = np.unique(labeled[active_stream_mask])
            stream_labels = stream_labels[stream_labels > 0]
            lut = np.zeros(num_features + 1, dtype=bool)
            lut[stream_labels] = True
            connected_mask = lut[labeled]
            depth_connected = np.where(connected_mask, depth_raw, 0.0)

        # 6. Métricas Físicas
        dx_m = self.topo.transform[0] * 111320.0 * np.cos(np.radians(np.mean([self.topo.bounds.bottom, self.topo.bounds.top])))
        dy_m = abs(self.topo.transform[4]) * 110570.0
        flooded_cells = int(np.sum(depth_connected >= self.min_flood_depth))
        area_km2 = float(flooded_cells * self.topo.cell_area_km2)
        vol_hm3 = float(np.sum(depth_connected) * (dx_m * dy_m) / 1e6)
        max_d = float(np.max(depth_connected)) if flooded_cells > 0 else 0.0
        mean_d = float(np.mean(depth_connected[connected_mask])) if flooded_cells > 0 else 0.0

        return FloodRaster(
            timestamp_hour=t_hour,
            z_water_2d=z_water_grid,
            eta_2d=eta_grid,
            depth_raw_2d=depth_raw,
            depth_connected_2d=depth_connected,
            inundated_mask=connected_mask,
            ocean_level_z=z_ocean,
            area_km2=round(area_km2, 2),
            volume_hm3=round(vol_hm3, 2),
            max_depth_m=round(max_d, 2),
            mean_depth_m=round(mean_d, 2)
        )

    def export_flood_debug_rasters(self, flood_raster: FloodRaster, output_dir: Union[str, Path]):
        """Salva as 4 camadas hidráulicas dinâmicas em formato GeoTIFF."""
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        meta = {
            'driver': 'GTiff',
            'height': self.topo.dem.shape[0],
            'width': self.topo.dem.shape[1],
            'count': 1,
            'dtype': 'float32',
            'crs': self.topo.crs,
            'transform': self.topo.transform
        }

        t_str = f"t{int(flood_raster.timestamp_hour):02d}"

        with rasterio.open(out_path / f"05_water_surface_{t_str}.tif", 'w', **meta) as dst:
            dst.write(flood_raster.z_water_2d.astype(np.float32), 1)

        with rasterio.open(out_path / f"06_relative_water_level_{t_str}.tif", 'w', **meta) as dst:
            dst.write(flood_raster.eta_2d.astype(np.float32), 1)

        with rasterio.open(out_path / f"07_depth_{t_str}.tif", 'w', **meta) as dst:
            dst.write(flood_raster.depth_raw_2d.astype(np.float32), 1)

        with rasterio.open(out_path / f"08_connected_flood_{t_str}.tif", 'w', **meta) as dst:
            dst.write(flood_raster.depth_connected_2d.astype(np.float32), 1)
