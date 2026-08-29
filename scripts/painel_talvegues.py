# -*- coding: utf-8 -*-
"""Painel interativo: talvegues dos rios + massas d'agua + mapa base.

    python scripts/painel_talvegues.py

Sai doc/painel/talvegues.html — autocontido, abre no navegador com
duplo clique. Leaflet via CDN; camadas base OSM / Google Satelite /
Esri Imagery; overlays: eixos dos rios do modelo (clicaveis), massas
d'agua e rios duplos da FBDS. Clicar num rio desenha o perfil do
talvegue embaixo; passar o mouse no perfil move um marcador no mapa.
"""
import glob
import json
import os

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
os.chdir(RAIZ)

G01_HDF = 'taha_ai.g01.hdf'
SAIDA = os.path.join('doc', 'painel', 'talvegues.html')


def talvegues():
    """Por rio/reach: [(RS, z_talvegue, lon, lat)] na ordem de jusante."""
    import h5py
    from pyproj import Transformer
    tr = Transformer.from_crs(31982, 4326, always_xy=True)
    f = h5py.File(G01_HDF, 'r')
    cs = f['Geometry/Cross Sections']
    at = cs['Attributes'][:]
    se_info = cs['Station Elevation Info'][:]
    se_val = cs['Station Elevation Values'][:]
    pl_info = cs['Polyline Info'][:]
    pl_pts = cs['Polyline Points'][:]
    rios = {}
    for k in range(len(at)):
        rio = at['River'][k].decode().strip()
        reach = at['Reach'][k].decode().strip()
        try:
            rs = float(at['RS'][k].decode())
        except ValueError:
            continue
        i0, n = se_info[k]
        if n < 2:
            continue
        z = float(se_val[i0:i0 + n, 1].min())
        j0, m = pl_info[k][0], pl_info[k][1]
        pts = pl_pts[j0:j0 + m]
        # ponto do talvegue ~ meio da cutline (bom o bastante p/ mapa)
        meio = pts[len(pts) // 2]
        lon, lat = tr.transform(float(meio[0]), float(meio[1]))
        rios.setdefault(rio, []).append(
            (rs, round(z, 2), round(lon, 6), round(lat, 6)))
    for nome in rios:
        rios[nome].sort(key=lambda t: t[0])
    return rios


def eixos_arquivados():
    """Rios arquivados (fora do modelo): Reach XY dos g01 de texto em
    modelos/*/ -- Luis Alves, Krauel."""
    import re
    from pyproj import Transformer
    tr = Transformer.from_crs(31982, 4326, always_xy=True)
    out = {}
    for g01 in glob.glob(os.path.join('modelos', '*', '*.g01')):
        try:
            txt = open(g01, encoding='latin-1').read()
        except OSError:
            continue
        nome = os.path.basename(os.path.dirname(g01))
        for m in re.finditer(
                r'River Reach=([^,]+),[^\n]*\n'
                r'Reach XY= *(\d+)\s*\n((?:[ \d.\-]+\n)+)', txt):
            n = int(m.group(2))
            vals = re.findall(r'[-\d.]+', m.group(3))
            vals = [float(v) for v in vals[:2 * n]]
            pts = list(zip(vals[0::2], vals[1::2]))
            if len(pts) < 2:
                continue
            lon, lat = tr.transform([p[0] for p in pts],
                                    [p[1] for p in pts])
            out.setdefault(nome, []).append(
                [[round(float(la), 6), round(float(lo), 6)]
                 for lo, la in zip(lon, lat)])
    return out


def centerlines():
    """Eixo suave de cada rio (River Centerlines do RAS), em 4326."""
    import h5py
    from pyproj import Transformer
    tr = Transformer.from_crs(31982, 4326, always_xy=True)
    f = h5py.File(G01_HDF, 'r')
    rc = f['Geometry/River Centerlines']
    at = rc['Attributes'][:]
    pl_info = rc['Polyline Info'][:]
    pl_pts = rc['Polyline Points'][:]
    eixos = {}
    for k in range(len(at)):
        rio = at['River Name'][k].decode().strip()
        j0, m = pl_info[k][0], pl_info[k][1]
        pts = pl_pts[j0:j0 + m]
        lon, lat = tr.transform(pts[:, 0], pts[:, 1])
        eixos.setdefault(rio, []).append(
            [[round(float(la), 6), round(float(lo), 6)]
             for lo, la in zip(lon, lat)])
    return eixos


def emendar_eixos(eixos):
    """Fecha os vaos visuais: cada ponta de reach e atada a juncao
    mais proxima (ou ao vertice mais proximo de OUTRO rio) ate 600 m.
    So visual -- a geometria do modelo nao muda."""
    import h5py
    import numpy as np
    from pyproj import Transformer
    tr = Transformer.from_crs(31982, 4326, always_xy=True)
    f = h5py.File(G01_HDF, 'r')
    jp = f['Geometry/Junctions/Points'][:]
    jlon, jlat = tr.transform(jp[:, 0].astype(float),
                              jp[:, 1].astype(float))
    juncoes = np.column_stack([jlat, jlon])
    todos = []
    for nome, partes in eixos.items():
        for p in partes:
            todos.extend([(nome, tuple(pt)) for pt in p])
    import math

    def d_m(a, b):
        return math.hypot((a[0] - b[0]) * 111000,
                          (a[1] - b[1]) * 99000)
    for nome, partes in eixos.items():
        for parte in partes:
            for idx in (0, -1):
                ponta = tuple(parte[idx])
                melhor, dm = None, 600.0
                for j in juncoes:
                    d = d_m(ponta, j)
                    if 1.0 < d < dm:
                        melhor, dm = [round(float(j[0]), 6),
                                      round(float(j[1]), 6)], d
                if melhor is None:
                    for outro, pt in todos:
                        if outro == nome:
                            continue
                        d = d_m(ponta, pt)
                        if 1.0 < d < dm:
                            melhor, dm = list(pt), d
                if melhor:
                    if idx == 0:
                        parte.insert(0, melhor)
                    else:
                        parte.append(melhor)
    return eixos


def bacia_prep():
    """Poligono oficial ANA (31982) preparado, com folga de 300 m."""
    from shapely.geometry import shape
    from shapely.prepared import prep
    bac = json.load(open('doc/qgis/bacia_itajai_ana.geojson',
                         encoding='utf-8'))
    pol = shape(bac['features'][0]['geometry']).buffer(300.0)
    return pol, prep(pol)


def bacia_4326():
    """O mesmo divisor, em graus (p/ dados que ja vem em 4326)."""
    from shapely.ops import transform as stransform
    from shapely.prepared import prep
    from pyproj import Transformer
    inv = Transformer.from_crs(31982, 4326, always_xy=True)
    pol, _ = bacia_prep()
    pol4326 = stransform(lambda x, y: inv.transform(x, y), pol)
    return pol4326, prep(pol4326)


def fbds_geojson(camada, tol=15.0, max_kb=2500):
    """Poligonos/linhas da FBDS DENTRO DA BACIA (divisor ANA),
    simplificados, em 4326."""
    import pyogrio
    from pyproj import Transformer
    from shapely.ops import transform as stransform
    tr = Transformer.from_crs(31982, 4326, always_xy=True)
    pol, prepped = bacia_prep()

    def p31982_4326(x, y):
        return tr.transform(x, y)
    feats = []
    for shp in sorted(glob.glob(f'doc/fbds/*/*_{camada}.shp')):
        try:
            g = pyogrio.read_dataframe(shp)
        except Exception:
            continue
        for geom in g.geometry:
            if geom is None:
                continue
            if not prepped.intersects(geom):
                continue                    # outra bacia: fora
            if not prepped.contains(geom):
                geom = geom.intersection(pol)
            s = geom.simplify(tol)
            if s.is_empty:
                continue
            s = stransform(p31982_4326, s)
            feats.append({'type': 'Feature', 'properties': {},
                          'geometry': s.__geo_interface__})
    gj = {'type': 'FeatureCollection', 'features': feats}
    txt = json.dumps(gj, separators=(',', ':'))
    print(f'  {camada}: {len(feats)} feicoes, {len(txt) // 1024} kB')
    if len(txt) > max_kb * 1024:
        print(f'  ({camada} acima de {max_kb} kB -- simplificando 3x)')
        return fbds_geojson(camada, tol * 3, max_kb * 10)
    return gj


HTML = """<!DOCTYPE html>
<html lang="pt-br"><head><meta charset="utf-8">
<title>Talvegues do Vale do Itajaí</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet"
 href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
 html,body{margin:0;height:100%;font-family:system-ui,sans-serif}
 #mapa{height:62%}
 #painel{height:38%;position:relative;background:#101418}
 #perfil{width:100%;height:100%;display:block}
 #titulo{position:absolute;top:6px;left:12px;color:#eee;font-size:14px;
         font-weight:600;pointer-events:none}
 #dica{position:absolute;top:6px;right:12px;color:#9ab;font-size:12px;
       pointer-events:none}
</style></head><body>
<div id="mapa"></div>
<div id="painel">
 <canvas id="perfil"></canvas>
 <div id="titulo">clique num rio (linha vermelha) para ver o talvegue</div>
 <div id="dica">passe o mouse no perfil para localizar no mapa</div>
</div>
<script>
const RIOS = @@RIOS@@;
const EIXOS = @@EIXOS@@;
const MASSAS = @@MASSAS@@;
const DUPLOS = @@DUPLOS@@;
const LAMINA = @@LAMINA@@;
const BARRAGENS = @@BARRAGENS@@;
const ARQUIVADOS = @@ARQUIVADOS@@;
const ANARIOS = @@ANARIOS@@;
const MUROS = @@MUROS@@;

const mapa = L.map('mapa');
const osm = L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',
  {maxZoom:19, attribution:'&copy; OpenStreetMap'});
const gsat = L.tileLayer('https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
  {maxZoom:20, attribution:'Google'});
const ghyb = L.tileLayer('https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
  {maxZoom:20, attribution:'Google'});
const esri = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  {maxZoom:19, attribution:'Esri'});
osm.addTo(mapa);

const camMassas = L.geoJSON(MASSAS, {style:{color:'#0a6cbf',weight:1,
  fillColor:'#3fa7ff',fillOpacity:0.45}});
const camDuplos = L.geoJSON(DUPLOS, {style:{color:'#0a6cbf',weight:1,
  fillColor:'#7fc4ff',fillOpacity:0.4}});
camMassas.addTo(mapa); camDuplos.addTo(mapa);

const camLamina = L.geoJSON(LAMINA, {style:{color:'#087f5b',weight:1,
  fillColor:'#12b886',fillOpacity:0.55}});

const camRios = L.layerGroup().addTo(mapa);
let marcador = null, atual = null;
const linhas = {};
for (const nome in RIOS) {
  const partes = EIXOS[nome] ||
    [RIOS[nome].map(p => [p[3], p[2]])];
  const l = L.polyline(partes, {color:'#d62828', weight:4,
                                opacity:0.95});
  l.bindTooltip(nome, {sticky:true});
  l.on('click', () => desenhar(nome));
  l.addTo(camRios);
  linhas[nome] = l;
}
setTimeout(() => {
  for (const n in linhas) linhas[n].bringToFront();
}, 300);
const camBarragens = L.geoJSON(BARRAGENS, {
  pointToLayer: (f, ll) => {
    const p = f.properties;
    const alt = p.BAR_NU_ALT_MAX_NIVEL_TERRENO || 0;
    if ((p.USO_PRINCIPAL||'').indexOf('inunda') >= 0)
      return L.circleMarker(ll, {radius:13, color:'#c92a2a',
        weight:4, fillColor:'#ffa8a8', fillOpacity:0.9});
    return L.circleMarker(ll, {radius: alt > 10 ? 9 : 5,
      color:'#5f3dc4', weight:2, fillColor:'#845ef7',
      fillOpacity:0.8});
  },
  onEachFeature: (f, l) => {
    const p = f.properties;
    l.bindPopup('<b>'+(p.BAR_NM_NOME||'?')+'</b><br>altura: '+
      (p.BAR_NU_ALT_MAX_NIVEL_TERRENO||'?')+' m<br>uso: '+
      (p.USO_PRINCIPAL||'?')+'<br>material: '+(p.TIPO_MATERIAL||'?'));
  }});
camBarragens.addTo(mapa);
const camAna = L.geoJSON(ANARIOS, {
  style: f => ({color:'#e8590c',
    weight: (f.properties.nome||'').indexOf('[')===0 ? 1.5 : 3,
    opacity:0.85}),
  onEachFeature: (f, l) => l.bindTooltip(f.properties.nome,
                                         {sticky:true})});
const CORMURO = {ponte:'#1864ab', represa:'#5f3dc4', '?':'#e67700'};
const camMuros = L.geoJSON(MUROS, {
  pointToLayer: (f, ll) => {
    const m = f.properties;
    return L.circleMarker(ll, {radius: Math.min(4 + m.altura_m, 14),
      color: CORMURO[m.classe] || '#e67700', weight: 2,
      fillOpacity: 0.75, fillColor: CORMURO[m.classe] || '#ffd43b'});
  },
  onEachFeature: (f, l) => {
    const m = f.properties;
    l.bindPopup('<b>muro ' + m.classe + '</b> ' + (m.nome || '') +
      '<br>degrau: ' + m.altura_m + ' m<br>área: ' + m.area_m2 + ' m²');
  }});
const camArquivados = L.layerGroup();
for (const nome in ARQUIVADOS) {
  const l = L.polyline(ARQUIVADOS[nome], {color:'#868e96', weight:3,
    dashArray:'8 6', opacity:0.9});
  l.bindTooltip(nome + ' (fora do modelo)', {sticky:true});
  l.addTo(camArquivados);
}
L.control.layers(
  {'OpenStreetMap':osm, 'Google Satélite':gsat,
   'Google Híbrido':ghyb, 'Esri Imagery':esri},
  {'Rios do modelo (centerline)':camRios,
   'Lâmina d\\'água SIG-SC 2010':camLamina,
   "Massas d'água (FBDS)":camMassas,
   'Rios duplos (FBDS)':camDuplos,
   'Barragens (SNISB)':camBarragens,
   'Rios arquivados (fora do modelo)':camArquivados,
   "Cursos d'água ANA (código Otto)":camAna,
   'Muros no MDT (degraus >1 m)':camMuros},
  {collapsed:false}).addTo(mapa);

const todos = Object.values(RIOS).flat();
mapa.fitBounds(todos.map(p => [p[3], p[2]]));

const cv = document.getElementById('perfil');
const tt = document.getElementById('titulo');
let dadosPerfil = null;
function desenhar(nome) {
  atual = nome;
  for (const n in linhas)
    linhas[n].setStyle({color: n===nome ? '#ffd60a' : '#d62828',
                        weight: n===nome ? 5 : 3});
  const v = RIOS[nome];
  const km = v.map(p => p[0]/1000), z = v.map(p => p[1]);
  dadosPerfil = {nome, v, km, z};
  tt.textContent = nome + ' — talvegue (' + v.length + ' seções)';
  pintar(null);
}
function pintar(ix) {
  const d = dadosPerfil; if (!d) return;
  const W = cv.width = cv.clientWidth * devicePixelRatio;
  const H = cv.height = cv.clientHeight * devicePixelRatio;
  const g = cv.getContext('2d');
  g.clearRect(0,0,W,H);
  const m = {l:70*devicePixelRatio, r:20*devicePixelRatio,
             t:34*devicePixelRatio, b:38*devicePixelRatio};
  const k0 = Math.min(...d.km), k1 = Math.max(...d.km);
  const z0 = Math.min(...d.z), z1 = Math.max(...d.z);
  const X = k => m.l + (k-k0)/(k1-k0||1)*(W-m.l-m.r);
  const Y = z => H-m.b - (z-z0)/(z1-z0||1)*(H-m.t-m.b);
  g.strokeStyle='#2c3b4a'; g.lineWidth=1; g.beginPath();
  for (let i=0;i<=5;i++){const y=m.t+i*(H-m.t-m.b)/5;
    g.moveTo(m.l,y); g.lineTo(W-m.r,y);}
  g.stroke();
  g.fillStyle='#9ab'; g.font = (12*devicePixelRatio)+'px sans-serif';
  for (let i=0;i<=5;i++){
    const z = z1 - i*(z1-z0)/5, y=m.t+i*(H-m.t-m.b)/5;
    g.fillText(z.toFixed(1)+' m', 8*devicePixelRatio, y+4);}
  for (let i=0;i<=8;i++){
    const k = k0 + i*(k1-k0)/8;
    g.fillText(k.toFixed(1)+' km', X(k)-16, H-10);}
  g.strokeStyle='#4dabf7'; g.lineWidth=2*devicePixelRatio;
  g.beginPath();
  d.km.forEach((k,i)=>{i?g.lineTo(X(k),Y(d.z[i])):g.moveTo(X(k),Y(d.z[i]));});
  g.stroke();
  if (ix!=null){
    g.fillStyle='#ffd60a';
    g.beginPath();
    g.arc(X(d.km[ix]), Y(d.z[ix]), 5*devicePixelRatio, 0, 7);
    g.fill();
    g.fillStyle='#ffd60a';
    g.fillText('RS '+d.v[ix][0]+'  z='+d.z[ix].toFixed(2)+' m',
               X(d.km[ix])+8, Y(d.z[ix])-8);
  }
}
cv.addEventListener('mousemove', e=>{
  const d = dadosPerfil; if (!d) return;
  const r = cv.getBoundingClientRect();
  const fx = (e.clientX-r.left)/r.width;
  const k0 = Math.min(...d.km), k1 = Math.max(...d.km);
  const alvo = k0 + fx*(k1-k0);
  let ix = 0, best = 1e18;
  d.km.forEach((k,i)=>{const q=Math.abs(k-alvo);
                       if(q<best){best=q;ix=i;}});
  pintar(ix);
  const p = d.v[ix];
  if (!marcador) marcador = L.circleMarker([p[3],p[2]],
    {radius:8, color:'#ffd60a', weight:3, fillOpacity:0.3}).addTo(mapa);
  else marcador.setLatLng([p[3],p[2]]);
});
window.addEventListener('resize', ()=>pintar(null));
</script></body></html>
"""


def main():
    print('extraindo talvegues...')
    rios = talvegues()
    for n, v in sorted(rios.items()):
        print(f'  {n}: {len(v)} secoes, z {v[0][1]}..{v[-1][1]} m')
    eixos = emendar_eixos(centerlines())
    print(f'centerlines: {sum(len(v) for v in eixos.values())} reaches')
    arquivados = eixos_arquivados()
    print(f'arquivados: {list(arquivados.keys())}')
    mur_arq = os.path.join('doc', 'painel', 'muros.geojson')
    muros = {'type': 'FeatureCollection', 'features': []}
    if os.path.exists(mur_arq):
        todos = json.load(open(mur_arq))['features']
        muros['features'] = [f for f in todos
                             if f['properties']['classe'] != '?'
                             or f['properties']['altura_m'] >= 2.0
                             ][:600]
        print(f'muros no painel: {len(muros["features"])}')
    ana_arq = os.path.join('doc', 'painel', 'ana_rios.geojson')
    anarios = {'type': 'FeatureCollection', 'features': []}
    if os.path.exists(ana_arq):
        anarios = json.load(open(ana_arq))
        print(f'cursos ANA: {len(anarios["features"])} tracos')
    print('recortando FBDS...')
    massas = fbds_geojson('MASSAS_DAGUA')
    duplos = fbds_geojson('RIOS_DUPLOS')
    lam_arq = os.path.join('doc', 'painel', 'lamina_sigsc.geojson')
    lamina = {'type': 'FeatureCollection', 'features': []}
    if os.path.exists(lam_arq):
        lamina = json.load(open(lam_arq))
        print(f'lamina SIG-SC: {len(lamina["features"])} poligonos')
    else:
        print('lamina SIG-SC ausente (rode lamina_do_sigsc.py)')
    # barragens do SNISB (cadastro oficial), recortadas pelo divisor
    bar_arq = os.path.join('doc', 'osm', 'snisb_barragens.geojson')
    barragens = {'type': 'FeatureCollection', 'features': []}
    if os.path.exists(bar_arq):
        from shapely.geometry import shape as sh2
        bruto = json.load(open(bar_arq, encoding='utf-8'))
        pol4326, prep4326 = bacia_4326()
        # SOMENTE as 3 barragens de contencao (ordem do professor)
        barragens['features'] = [
            f for f in bruto.get('features', [])
            if f.get('geometry')
            and 'inunda' in (f['properties'].get('USO_PRINCIPAL')
                             or '').lower()
            and prep4326.intersects(sh2(f['geometry']))]
        print(f'barragens SNISB na bacia: '
              f'{len(barragens["features"])}')
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    html = (HTML
            .replace('@@RIOS@@', json.dumps(rios,
                                            separators=(',', ':')))
            .replace('@@EIXOS@@', json.dumps(eixos,
                                             separators=(',', ':')))
            .replace('@@LAMINA@@', json.dumps(lamina,
                                              separators=(',', ':')))
            .replace('@@BARRAGENS@@', json.dumps(
                barragens, separators=(',', ':')))
            .replace('@@ARQUIVADOS@@', json.dumps(
                arquivados, separators=(',', ':')))
            .replace('@@ANARIOS@@', json.dumps(
                anarios, separators=(',', ':')))
            .replace('@@MUROS@@', json.dumps(
                muros, separators=(',', ':')))
            .replace('@@MASSAS@@', json.dumps(massas,
                                              separators=(',', ':')))
            .replace('@@DUPLOS@@', json.dumps(duplos,
                                              separators=(',', ':'))))
    with open(SAIDA, 'w', encoding='utf-8') as fh:
        fh.write(html)
    print(f'painel: {SAIDA} ({os.path.getsize(SAIDA) // 1024} kB)')


if __name__ == '__main__':
    main()
