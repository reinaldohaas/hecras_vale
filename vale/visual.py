# -*- coding: utf-8 -*-
"""
Visualizacao interativa em JavaScript, num arquivo HTML unico.

Sem servidor e sem biblioteca externa. O arquivo abre com duplo clique, roda
offline e pode ser mandado por e-mail -- o que importa quando o destinatario e
uma defesa civil e nao um laboratorio. Tudo (dados e codigo) vai embutido; o
desenho e Canvas puro.

O que a pagina mostra:
  MAPA        eixos dos rios e cutlines, coloridos pela profundidade no
              instante escolhido
  LINHA DO TEMPO  animacao da cheia, com play e arrasto
  PERFIL      perfil longitudinal do rio selecionado, com o leito, a lamina no
              instante e a envoltoria de maxima
  SECAO       a secao transversal clicada, com a lamina
  HIDROGRAMA  nivel e vazao no ponto selecionado

A escolha de cor e por PROFUNDIDADE, nao por cota. Mapa de cota mostra o
relevo (o rio desce 300 m) e a cheia some dentro dessa variacao -- foi
exatamente o que aconteceu com o primeiro mapa gerado, que parecia tudo
inundado a partir de 112 m ja no passo zero.
"""
import json
import os

import numpy as np

PAGINA = """<!doctype html>
<html lang="pt-BR"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITULO__</title>
<style>
:root{--bg:#0f1417;--pane:#161d21;--linha:#243036;--txt:#e6edf0;
      --fraco:#8fa3ad;--acento:#4fc3f7;--seco:#3a4a52}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--txt);
     font:13px/1.45 "Segoe UI",system-ui,sans-serif}
header{padding:10px 16px;border-bottom:1px solid var(--linha);
       display:flex;gap:16px;align-items:center;flex-wrap:wrap}
h1{font-size:15px;margin:0;font-weight:600;letter-spacing:.02em}
.sub{color:var(--fraco);font-size:12px}
main{display:grid;grid-template-columns:1.35fr 1fr;gap:1px;background:var(--linha);
     height:calc(100vh - 120px)}
section{background:var(--pane);position:relative;overflow:hidden}
.dir{display:grid;grid-template-rows:1fr 1fr 1fr;gap:1px;background:var(--linha)}
canvas{display:block;width:100%;height:100%}
.rot{position:absolute;top:8px;left:12px;font-size:11px;color:var(--fraco);
     text-transform:uppercase;letter-spacing:.08em;pointer-events:none}
footer{padding:10px 16px;border-top:1px solid var(--linha);
       display:flex;gap:14px;align-items:center}
button{background:var(--seco);color:var(--txt);border:0;border-radius:4px;
       padding:6px 14px;cursor:pointer;font-size:13px}
button:hover{background:#4a5d66}
input[type=range]{flex:1;accent-color:var(--acento)}
select{background:var(--seco);color:var(--txt);border:0;border-radius:4px;
       padding:5px 8px}
.leg{display:flex;gap:10px;align-items:center;font-size:11px;color:var(--fraco)}
.sw{width:14px;height:10px;border-radius:2px;display:inline-block}
#info{position:absolute;right:10px;top:8px;background:rgba(0,0,0,.55);
      padding:6px 9px;border-radius:4px;font-size:11px;color:var(--txt)}
</style></head><body>
<header>
  <h1>__TITULO__</h1>
  <span class="sub" id="quando"></span>
  <span class="leg" style="margin-left:auto">
    <span class="sw" style="background:#0d47a1"></span>0 m
    <span class="sw" style="background:#26a69a"></span>2
    <span class="sw" style="background:#fdd835"></span>5
    <span class="sw" style="background:#ef6c00"></span>8
    <span class="sw" style="background:#b71c1c"></span>12+
  </span>
</header>
<main>
  <section><span class="rot">mapa</span><div id="info"></div><canvas id="mapa"></canvas></section>
  <div class="dir">
    <section><span class="rot">perfil longitudinal</span><canvas id="perfil"></canvas></section>
    <section><span class="rot">secao transversal</span><canvas id="secao"></canvas></section>
    <section><span class="rot">hidrograma</span><canvas id="hidro"></canvas></section>
  </div>
</main>
<footer>
  <button id="play">play</button>
  <input type="range" id="t" min="0" value="0">
  <select id="rio"></select>
</footer>
<script>
const D = __DADOS__;
const cores = [[0,[13,71,161]],[2,[38,166,154]],[5,[253,216,53]],
               [8,[239,108,0]],[12,[183,28,28]]];
function cor(p){
  if(!(p>0)) return null;
  let a=cores[0], b=cores[cores.length-1];
  for(let i=0;i<cores.length-1;i++){
    if(p>=cores[i][0] && p<=cores[i+1][0]){a=cores[i];b=cores[i+1];break}}
  const f=Math.min(1,Math.max(0,(p-a[0])/((b[0]-a[0])||1)));
  const c=a[1].map((v,i)=>Math.round(v+(b[1][i]-v)*f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}
let t=0, tocando=false, rioSel=D.rios[0].nome, secSel=null;

// -------- enquadramento comum a todos os eixos
function limites(){
  let x0=1e18,y0=1e18,x1=-1e18,y1=-1e18;
  D.rios.forEach(r=>r.eixo.forEach(p=>{
    x0=Math.min(x0,p[0]);x1=Math.max(x1,p[0]);
    y0=Math.min(y0,p[1]);y1=Math.max(y1,p[1]);}));
  return {x0,y0,x1,y1};
}
const L = limites();
function ctx(id){const c=document.getElementById(id);
  const r=c.getBoundingClientRect();const dpr=devicePixelRatio||1;
  c.width=r.width*dpr;c.height=r.height*dpr;
  const g=c.getContext('2d');g.scale(dpr,dpr);
  g.clearRect(0,0,r.width,r.height);return [g,r.width,r.height];}

function desenhaMapa(){
  const [g,W,H]=ctx('mapa');
  const m=18, sx=(W-2*m)/(L.x1-L.x0), sy=(H-2*m)/(L.y1-L.y0);
  const s=Math.min(sx,sy);
  const px=p=>[m+(p[0]-L.x0)*s, H-m-(p[1]-L.y0)*s];
  D.rios.forEach(r=>{
    g.beginPath();
    r.eixo.forEach((p,i)=>{const q=px(p);i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])});
    g.strokeStyle=(r.nome===rioSel)?'#4fc3f7':'#37474f';
    g.lineWidth=(r.nome===rioSel)?2:1.2; g.stroke();
    r.secoes.forEach((sec,i)=>{
      const p = r.prof[t] ? r.prof[t][i] : 0;
      const c = cor(p);
      const a=px(sec.a), b=px(sec.b);
      g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]);
      g.strokeStyle = c || '#2b3940'; g.lineWidth = c?2.2:0.8; g.stroke();
      if(secSel && secSel.rio===r.nome && secSel.i===i){
        g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]);
        g.strokeStyle='#ffeb3b'; g.lineWidth=1; g.stroke();}
    });
  });
  desenhaMapa.px=px;
}

function desenhaPerfil(){
  const [g,W,H]=ctx('perfil');
  const r=D.rios.find(x=>x.nome===rioSel); if(!r) return;
  const m=34, n=r.rs.length;
  const x0=Math.min(...r.rs), x1=Math.max(...r.rs);
  const lo=Math.min(...r.leito), hi=Math.max(...r.maxws)+2;
  const px=(v,z)=>[m+(1-(v-x0)/((x1-x0)||1))*(W-m-8),
                   H-24-((z-lo)/((hi-lo)||1))*(H-40)];
  const linha=(vals,cor2,larg)=>{g.beginPath();
    for(let i=0;i<n;i++){const q=px(r.rs[i],vals[i]);
      i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])}
    g.strokeStyle=cor2;g.lineWidth=larg;g.stroke();};
  linha(r.maxws,'#37474f',1);
  const ws = r.ws[t] || r.leito;
  g.beginPath();
  for(let i=0;i<n;i++){const q=px(r.rs[i],ws[i]); i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])}
  for(let i=n-1;i>=0;i--){const q=px(r.rs[i],r.leito[i]); g.lineTo(q[0],q[1])}
  g.closePath(); g.fillStyle='rgba(79,195,247,.28)'; g.fill();
  linha(r.leito,'#8d6e63',1.6);
  linha(ws,'#4fc3f7',1.4);
  g.fillStyle='#8fa3ad';g.font='10px sans-serif';
  g.fillText('cota (m)',4,12); g.fillText('RS (km) — jusante à direita',m,H-6);
}

function desenhaSecao(){
  const [g,W,H]=ctx('secao');
  if(!secSel){g.fillStyle='#8fa3ad';g.font='12px sans-serif';
    g.fillText('clique numa seção no mapa',14,26);return;}
  const r=D.rios.find(x=>x.nome===secSel.rio); const s=r.secoes[secSel.i];
  const m=34, n=s.sta.length;
  const x1=s.sta[n-1], lo=Math.min(...s.z);
  const ws = (r.ws[t]||r.leito)[secSel.i];
  const hi=Math.max(Math.max(...s.z), ws+1);
  const px=(v,z)=>[m+(v/(x1||1))*(W-m-8), H-22-((z-lo)/((hi-lo)||1))*(H-38)];
  if(ws>lo){
    g.beginPath();
    for(let i=0;i<n;i++){const q=px(s.sta[i],Math.min(s.z[i],ws));
      i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])}
    for(let i=n-1;i>=0;i--){const q=px(s.sta[i],ws); g.lineTo(q[0],q[1])}
    g.closePath(); g.fillStyle='rgba(79,195,247,.30)'; g.fill();
  }
  g.beginPath();
  for(let i=0;i<n;i++){const q=px(s.sta[i],s.z[i]); i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])}
  g.strokeStyle='#cfd8dc'; g.lineWidth=1.3; g.stroke();
  g.fillStyle='#8fa3ad';g.font='10px sans-serif';
  g.fillText(`${secSel.rio}  RS ${s.rs.toFixed(0)}   lâmina ${(ws-lo).toFixed(2)} m`,m,14);
}

function desenhaHidro(){
  const [g,W,H]=ctx('hidro');
  const r=D.rios.find(x=>x.nome===rioSel); if(!r||!secSel) {
    g.fillStyle='#8fa3ad';g.font='12px sans-serif';
    g.fillText('selecione uma seção',14,26);return;}
  const serie=r.ws.map(v=>v[secSel.i]);
  const lo=Math.min(...serie), hi=Math.max(...serie);
  const m=34;
  const px=(i,v)=>[m+(i/((serie.length-1)||1))*(W-m-8),
                   H-20-((v-lo)/((hi-lo)||1))*(H-34)];
  g.beginPath();
  serie.forEach((v,i)=>{const q=px(i,v); i?g.lineTo(q[0],q[1]):g.moveTo(q[0],q[1])});
  g.strokeStyle='#4fc3f7';g.lineWidth=1.5;g.stroke();
  const q=px(t,serie[t]);
  g.beginPath();g.arc(q[0],q[1],3.5,0,7);g.fillStyle='#ffeb3b';g.fill();
  g.fillStyle='#8fa3ad';g.font='10px sans-serif';
  g.fillText(`nível ${serie[t].toFixed(2)} m   (máx ${hi.toFixed(2)})`,m,14);
}

function tudo(){
  desenhaMapa();desenhaPerfil();desenhaSecao();desenhaHidro();
  document.getElementById('quando').textContent =
    `${D.horas[t]}   —   passo ${t+1} de ${D.horas.length}`;
}
document.getElementById('t').max = D.horas.length-1;
document.getElementById('t').oninput = e=>{t=+e.target.value;tudo()};
document.getElementById('play').onclick = function(){
  tocando=!tocando; this.textContent = tocando?'pausa':'play';
  (function passo(){ if(!tocando) return;
    t=(t+1)%D.horas.length; document.getElementById('t').value=t; tudo();
    setTimeout(passo,140); })();
};
const sel=document.getElementById('rio');
D.rios.forEach(r=>{const o=document.createElement('option');
  o.value=o.textContent=r.nome; sel.appendChild(o)});
sel.onchange=e=>{rioSel=e.target.value;secSel=null;tudo()};
document.getElementById('mapa').onclick=ev=>{
  const c=ev.target.getBoundingClientRect();
  const mx=ev.clientX-c.left, my=ev.clientY-c.top;
  let melhor=null, dmin=1e9;
  D.rios.forEach(r=>r.secoes.forEach((s,i)=>{
    const a=desenhaMapa.px(s.a), b=desenhaMapa.px(s.b);
    const cx=(a[0]+b[0])/2, cy=(a[1]+b[1])/2;
    const d=(cx-mx)**2+(cy-my)**2;
    if(d<dmin){dmin=d;melhor={rio:r.nome,i:i}}}));
  if(melhor && dmin<3600){secSel=melhor;rioSel=melhor.rio;
    document.getElementById('rio').value=rioSel;tudo()}
};
addEventListener('resize',tudo); tudo();
</script></body></html>
"""


def _sub(x, n):
    """Reduz uma lista para no maximo n itens, preservando as pontas."""
    x = list(x)
    if len(x) <= n:
        return x
    idx = np.unique(np.linspace(0, len(x) - 1, n).astype(int))
    return [x[i] for i in idx]


def dados(trechos, wse=None, horas=None, max_secoes=140, max_pontos=60):
    """Empacota o que a pagina precisa.

    A reducao existe por um motivo pratico: um HTML com 1.400 secoes de 280
    pontos e 192 instantes passa de 100 MB e nenhum navegador abre com
    conforto. O que se perde e resolucao de desenho, nao de modelo.
    """
    rios = []
    for t in trechos:
        xs = _sub(t["xs"], max_secoes)
        leito = [float(np.min(d["z"])) for d in xs]
        secoes = []
        for d in xs:
            sta = np.asarray(d["sta"], float)
            z = np.asarray(d["z"], float)
            k = np.unique(np.linspace(0, len(sta) - 1, max_pontos).astype(int))
            secoes.append({"rs": float(d["rs"]),
                           "sta": [round(float(v), 1) for v in sta[k]],
                           "z": [round(float(v), 2) for v in z[k]],
                           "a": [round(d["cut"][0], 1), round(d["cut"][1], 1)],
                           "b": [round(d["cut"][2], 1), round(d["cut"][3], 1)]})
        n_t = len(horas) if horas else 1
        if wse is not None and t["rio"] in wse:
            w = np.asarray(wse[t["rio"]], float)   # (tempo, secao)
        else:
            w = np.tile(np.array(leito, float), (n_t, 1))
        ws = [[round(float(v), 2) for v in linha] for linha in w]
        prof = [[max(round(float(v - b), 2), 0.0) for v, b in zip(linha, leito)]
                for linha in w]
        rios.append({
            "nome": t["rio"],
            "eixo": [[round(x, 1), round(y, 1)]
                     for x, y in _sub(list(t["linha"].coords), 400)],
            "rs": [round(float(d["rs"]) / 1000.0, 3) for d in xs],
            "leito": [round(v, 2) for v in leito],
            "maxws": [round(float(np.max(w[:, i])), 2) for i in range(w.shape[1])],
            "ws": ws, "prof": prof, "secoes": secoes})
    return {"rios": rios,
            "horas": horas or ["condicao inicial"]}


def gerar(op, trechos, destino=None, wse=None, horas=None, titulo=None):
    d = dados(trechos, wse, horas)
    html = (PAGINA
            .replace("__TITULO__", titulo or f"{op.projeto} — Vale do Itajaí")
            .replace("__DADOS__", json.dumps(d, separators=(",", ":"))))
    destino = destino or op.caminho(f"{op.projeto}_visual.html")
    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    with open(destino, "w", encoding="utf-8") as f:
        f.write(html)
    return destino
