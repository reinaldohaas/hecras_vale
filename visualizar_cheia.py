# -*- coding: utf-8 -*-
"""
Sistema de visualizacao das cheias da bacia do Itajai.

Gera uma pagina unica, autocontida, que mostra a cheia de tres formas ligadas
pelo mesmo controle de tempo:

  1. PERFIL LONGITUDINAL animado -- leito, margens e lamina d'agua ao longo do
     rio, hora a hora, com as cidades marcadas.
  2. SECAO TRANSVERSAL na estaca escolhida, com a lamina do instante.
  3. HIDROGRAMA da estaca, com o instante atual destacado.

Le de DUAS fontes, com a mesma estrutura de saida:

  --fonte hecras  ->  <PROJETO>.p01.hdf     (solver do HEC-RAS)
  --fonte motor   ->  <PROJETO>_motor.npz   (motor quasi-permanente proprio)

Isso permite comparar os dois lado a lado no mesmo formato, que e como se
confere um contra o outro. Quando o HEC-RAS interrompe por instabilidade, a
pagina mostra ate onde ele foi e diz isso explicitamente, em vez de exibir um
resultado parcial como se fosse completo.

Uso:
  python visualizar_cheia.py Itajai_Rede_1983
  python visualizar_cheia.py Itajai_Rede_1983 --fonte motor
"""
import argparse
import datetime
import json
import os

import numpy as np

CIDADES = [
    ("Rio do Sul", 185.8), ("Ibirama", 148.0), ("Apiúna", 120.0),
    ("Indaial", 94.0), ("Blumenau", 68.0), ("Gaspar", 44.0),
    ("Ilhota", 28.0), ("Itajaí", 5.0),
]


# ------------------------------------------------------------------- LEITURA
def ler_hecras(projeto):
    import h5py
    base = ("Results/Unsteady/Output/Output Blocks/Base Output/"
            "Unsteady Time Series/Cross Sections")
    with h5py.File(f"{projeto}.p01.hdf", "r") as f:
        sol = f["Results/Unsteady"].attrs.get("Solution")
        sol = sol.decode() if isinstance(sol, bytes) else str(sol)
        g = f[base]
        ws, q = g["Water Surface"][:], g["Flow"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = np.array([x["River"].decode().strip() for x in at])
        rch = np.array([x["Reach"].decode().strip() for x in at])
        rs = np.array([float(x["RS"].decode()) for x in at])
        se = f["Geometry/Cross Sections/Station Elevation Values"][:]
        info = f["Geometry/Cross Sections/Station Elevation Info"][:]
    perfis = [se[a:a + n] for a, n in info]
    return dict(ws=ws, q=q, riv=riv, rch=rch, rs=rs, perfis=perfis,
                estado=sol, fonte="HEC-RAS 7.0.1")


def ler_motor(projeto):
    import motor_hidraulico as M
    d = np.load(f"{projeto}_motor.npz", allow_pickle=True)
    secs, _ = M.ler_geometria(projeto)
    ordem = {}
    for s in secs:
        ordem[(s["rio"], s["reach"], round(s["rs"], 2))] = s
    riv, rch, rs = d["river"], d["reach"], d["rs"]
    perfis = []
    for i in range(len(rs)):
        s = ordem.get((str(riv[i]), str(rch[i]), round(float(rs[i]), 2)))
        perfis.append(np.column_stack([s["sta"], s["z"]]) if s is not None
                      else np.zeros((2, 2)))
    return dict(ws=d["ws"], q=d["q"], riv=riv, rch=rch, rs=rs, perfis=perfis,
                estado="Concluido", fonte="Motor quasi-permanente")


# ------------------------------------------------------------------- MONTAGEM
def preparar(d, rio="Itajai_Acu", n_secoes=90, n_tempos=60):
    """Reduz a malha para caber na pagina sem perder a forma da cheia."""
    m = d["riv"] == rio
    idx = np.flatnonzero(m)
    idx = idx[np.argsort(-d["rs"][idx])]                 # montante -> jusante
    if len(idx) > n_secoes:
        idx = idx[np.linspace(0, len(idx) - 1, n_secoes).round().astype(int)]
    nt = d["ws"].shape[0]
    ts = np.arange(nt) if nt <= n_tempos else \
        np.linspace(0, nt - 1, n_tempos).round().astype(int)

    leito = [float(np.min(d["perfis"][j][:, 1])) for j in idx]
    topo = [float(np.max(d["perfis"][j][:, 1])) for j in idx]
    secs = []
    for j in idx:
        p = d["perfis"][j]
        k = np.linspace(0, len(p) - 1, min(len(p), 70)).round().astype(int)
        secs.append({"x": [round(float(v), 1) for v in p[k, 0]],
                     "z": [round(float(v), 2) for v in p[k, 1]]})
    return {
        "rio": rio,
        "rs": [round(float(d["rs"][j]) / 1000, 2) for j in idx],
        "leito": [round(v, 2) for v in leito],
        "topo": [round(v, 2) for v in topo],
        "horas": [int(t) for t in ts],
        "ws": [[round(float(d["ws"][t, j]), 2) for j in idx] for t in ts],
        "q": [[round(float(d["q"][t, j]), 1) for j in idx] for t in ts],
        "secoes": secs,
        "cidades": [{"nome": n, "rs": r} for n, r in CIDADES],
    }


def pagina(dados, projeto, meta):
    j = json.dumps(dados, ensure_ascii=False, separators=(",", ":"))
    aviso = ""
    if "Success" not in meta["estado"]:
        aviso = (f'<div class="aviso"><b>Simulação incompleta.</b> '
                 f'O solver parou com <code>{meta["estado"]}</code> após '
                 f'{meta["n_saidas"]} de {meta["n_previstas"]} saídas horárias. '
                 f'O que está abaixo é só o trecho calculado — não é a cheia '
                 f'inteira.</div>')
    return TEMPLATE.replace("__DADOS__", j) \
                   .replace("__PROJETO__", projeto) \
                   .replace("__FONTE__", meta["fonte"]) \
                   .replace("__INICIO__", meta["inicio"]) \
                   .replace("__AVISO__", aviso)


TEMPLATE = r"""<meta charset="utf-8">
<title>Cheia do Itajaí</title>
<style>
:root{
  --ceu:#eef3f7; --tinta:#16222c; --fraco:#5b7183; --linha:#c9d6e0;
  --agua:#2f7fb8; --agua-fraca:#9fc9e4; --terra:#8a7250; --leito:#3d3226;
  --alerta:#b3411f; --papel:#ffffff; --painel:#f7fafc;
}
:root:not([data-theme="light"]){}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --ceu:#0f1720; --tinta:#e7eef4; --fraco:#93a7b7; --linha:#25333f;
  --agua:#4aa3d8; --agua-fraca:#1d4e6d; --terra:#b09068; --leito:#c2ab8c;
  --alerta:#e2703f; --papel:#131c25; --painel:#18232e;
}}
:root[data-theme="dark"]{
  --ceu:#0f1720; --tinta:#e7eef4; --fraco:#93a7b7; --linha:#25333f;
  --agua:#4aa3d8; --agua-fraca:#1d4e6d; --terra:#b09068; --leito:#c2ab8c;
  --alerta:#e2703f; --papel:#131c25; --painel:#18232e;
}
*{box-sizing:border-box}
body{margin:0;background:var(--ceu);color:var(--tinta);
  font:15px/1.5 ui-sans-serif,system-ui,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:24px 20px 60px;
  display:flex;flex-direction:column;gap:20px}
h1{font-size:26px;margin:0;letter-spacing:-.01em}
.sub{color:var(--fraco);font-size:14px;margin:0}
.aviso{background:color-mix(in srgb,var(--alerta) 12%,var(--papel));
  border-left:4px solid var(--alerta);padding:12px 16px;border-radius:0 6px 6px 0}
.aviso code{font-size:13px}
.card{background:var(--papel);border:1px solid var(--linha);border-radius:10px;
  padding:16px 18px}
.card h2{font-size:13px;text-transform:uppercase;letter-spacing:.08em;
  color:var(--fraco);margin:0 0 12px;font-weight:600}
.ctrl{display:flex;align-items:center;gap:14px;flex-wrap:wrap}
.ctrl input[type=range]{flex:1;min-width:220px;accent-color:var(--agua)}
button{background:var(--agua);color:#fff;border:0;border-radius:6px;
  padding:7px 15px;font:inherit;font-weight:600;cursor:pointer}
button:hover{filter:brightness(1.1)}
.rel{font-variant-numeric:tabular-nums;color:var(--fraco);font-size:14px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:840px){.grid{grid-template-columns:1fr}}
svg{display:block;width:100%;height:auto}
.legenda{display:flex;gap:18px;flex-wrap:wrap;font-size:13px;
  color:var(--fraco);margin-top:10px}
.chave{display:inline-flex;align-items:center;gap:6px}
.sw{width:14px;height:3px;border-radius:2px;display:inline-block}
.dica{color:var(--fraco);font-size:13px;margin:8px 0 0}
</style>
<div class="wrap">
  <div>
    <h1>Cheia do Itajaí — __PROJETO__</h1>
    <p class="sub">__FONTE__ &middot; início em __INICIO__</p>
  </div>
  __AVISO__
  <div class="card">
    <div class="ctrl">
      <button id="play">▶ Animar</button>
      <input type="range" id="t" min="0" value="0">
      <span class="rel" id="rel"></span>
    </div>
  </div>
  <div class="card">
    <h2>Perfil longitudinal — clique para escolher a seção</h2>
    <svg id="perfil" viewBox="0 0 1100 380"></svg>
    <div class="legenda">
      <span class="chave"><i class="sw" style="background:var(--leito)"></i>leito escavado</span>
      <span class="chave"><i class="sw" style="background:var(--terra)"></i>topo da seção</span>
      <span class="chave"><i class="sw" style="background:var(--agua)"></i>lâmina d'água</span>
      <span class="chave"><i class="sw" style="background:var(--alerta)"></i>seção escolhida</span>
    </div>
  </div>
  <div class="grid">
    <div class="card">
      <h2 id="tit-sec">Seção transversal</h2>
      <svg id="secao" viewBox="0 0 540 300"></svg>
    </div>
    <div class="card">
      <h2 id="tit-hid">Hidrograma</h2>
      <svg id="hidro" viewBox="0 0 540 300"></svg>
    </div>
  </div>
  <p class="dica">Perfil de montante (esquerda) para a foz (direita). A escala
  vertical do perfil é logarítmica no desnível para caber a serra e a planície
  no mesmo gráfico.</p>
</div>
<script>
const D = __DADOS__;
const SVGNS = "http://www.w3.org/2000/svg";
let ti = 0, sel = Math.floor(D.rs.length/2), tocando = false, timer = null;

const el = id => document.getElementById(id);
const mk = (t,a) => { const e = document.createElementNS(SVGNS,t);
  for(const k in a) e.setAttribute(k,a[k]); return e; };
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

/* escala vertical comprimida: o vale vai de -14 m a 460 m */
const comp = z => Math.sign(z)*Math.log1p(Math.abs(z)/3);

function perfil(){
  const s = el("perfil"); s.replaceChildren();
  const W=1100,H=380,L=52,R=14,T=14,B=34;
  const n = D.rs.length;
  const zs = D.leito.concat(D.topo, D.ws[ti]);
  let lo=Math.min(...zs), hi=Math.max(...zs);
  const clo=comp(lo), chi=comp(hi);
  /* X por DISTANCIA, nao por indice da secao. O espacamento das secoes varia
     de 0,45 a 4,0 km -- ha refino na garganta do Salto Pilao -- entao plotar
     por indice dava a esse trecho 35% da largura para 21% do rio, e espremia
     todo o baixo vale. Ibirama aparecia em 0,35 quando deveria estar em 0,21. */
  const rsA = D.rs[0], rsB = D.rs[n-1], span = (rsA - rsB) || 1;
  const X = i => L + (W-L-R)*(rsA - D.rs[i])/span;
  const Y = z => T + (H-T-B)*(1-(comp(z)-clo)/((chi-clo)||1));

  /* grade */
  for(let k=0;k<=4;k++){
    const y=T+(H-T-B)*k/4;
    s.appendChild(mk("line",{x1:L,x2:W-R,y1:y,y2:y,stroke:css('--linha'),"stroke-width":1}));
  }
  const pol = (vals,fill,op) => {
    let d=`M ${X(0)} ${Y(vals[0])}`;
    for(let i=1;i<n;i++) d+=` L ${X(i)} ${Y(vals[i])}`;
    if(fill){ d+=` L ${X(n-1)} ${H-B} L ${X(0)} ${H-B} Z`; }
    return mk("path",{d:d,fill:fill||"none",stroke:fill?"none":"",opacity:op||1});
  };
  /* agua preenchida ate o leito */
  let d=`M ${X(0)} ${Y(D.ws[ti][0])}`;
  for(let i=1;i<n;i++) d+=` L ${X(i)} ${Y(D.ws[ti][i])}`;
  for(let i=n-1;i>=0;i--) d+=` L ${X(i)} ${Y(D.leito[i])}`;
  s.appendChild(mk("path",{d:d+" Z",fill:css('--agua'),opacity:.35}));

  const linha=(vals,cor,w)=>{ let p=`M ${X(0)} ${Y(vals[0])}`;
    for(let i=1;i<n;i++) p+=` L ${X(i)} ${Y(vals[i])}`;
    s.appendChild(mk("path",{d:p,fill:"none",stroke:cor,"stroke-width":w,
      "stroke-linejoin":"round"})); };
  linha(D.topo, css('--terra'), 1.2);
  linha(D.leito, css('--leito'), 1.8);
  linha(D.ws[ti], css('--agua'), 2.2);

  /* cidades */
  D.cidades.forEach(c=>{
    let i=0,melhor=1e9;
    D.rs.forEach((r,k)=>{ const dd=Math.abs(r-c.rs); if(dd<melhor){melhor=dd;i=k;} });
    if(melhor>8) return;
    s.appendChild(mk("line",{x1:X(i),x2:X(i),y1:T,y2:H-B,stroke:css('--linha'),
      "stroke-width":1,"stroke-dasharray":"3 4"}));
    const t=mk("text",{x:X(i),y:H-B+14,fill:css('--fraco'),"font-size":11,
      "text-anchor":"middle"}); t.textContent=c.nome; s.appendChild(t);
  });
  /* secao escolhida */
  s.appendChild(mk("line",{x1:X(sel),x2:X(sel),y1:T,y2:H-B,
    stroke:css('--alerta'),"stroke-width":2}));
  /* eixo Y */
  for(let k=0;k<=4;k++){
    const y=T+(H-T-B)*k/4;
    const zv=lo+(hi-lo)*(1-k/4);
    const t=mk("text",{x:L-8,y:y+4,fill:css('--fraco'),"font-size":11,
      "text-anchor":"end"}); t.textContent=zv.toFixed(0)+" m"; s.appendChild(t);
  }
  s.onclick = ev => {
    const b=s.getBoundingClientRect();
    const x=(ev.clientX-b.left)/b.width*W;
    /* inverte o eixo em DISTANCIA e procura a secao mais proxima */
    const rsAlvo = rsA - span*(x-L)/(W-L-R);
    let melhor=0, dmin=1e12;
    for(let i=0;i<n;i++){ const d=Math.abs(D.rs[i]-rsAlvo);
      if(d<dmin){dmin=d;melhor=i;} }
    sel=melhor; desenhar();
  };
}

function secao(){
  const s=el("secao"); s.replaceChildren();
  const W=540,H=300,L=48,R=12,T=12,B=30;
  const sc=D.secoes[sel], zw=D.ws[ti][sel];
  const xs=sc.x, zs=sc.z;
  const x0=Math.min(...xs), x1=Math.max(...xs);
  let lo=Math.min(...zs), hi=Math.max(Math.max(...zs), zw+1);
  const X=v=>L+(W-L-R)*(v-x0)/((x1-x0)||1);
  const Y=v=>T+(H-T-B)*(1-(v-lo)/((hi-lo)||1));
  /* agua */
  let dpol=null;
  for(let i=0;i<xs.length;i++){
    if(zs[i]<=zw){ if(!dpol) dpol=`M ${X(xs[i])} ${Y(zw)}`; }
  }
  if(dpol){
    let d=`M ${X(x0)} ${Y(zw)} L ${X(x1)} ${Y(zw)}`;
    for(let i=xs.length-1;i>=0;i--) d+=` L ${X(xs[i])} ${Y(Math.min(zs[i],zw))}`;
    s.appendChild(mk("path",{d:d+" Z",fill:css('--agua'),opacity:.4}));
  }
  /* terreno */
  let d=`M ${X(xs[0])} ${Y(zs[0])}`;
  for(let i=1;i<xs.length;i++) d+=` L ${X(xs[i])} ${Y(zs[i])}`;
  s.appendChild(mk("path",{d:d,fill:"none",stroke:css('--leito'),"stroke-width":2}));
  s.appendChild(mk("line",{x1:L,x2:W-R,y1:Y(zw),y2:Y(zw),stroke:css('--agua'),
    "stroke-width":1.5,"stroke-dasharray":"5 4"}));
  for(let k=0;k<=3;k++){
    const y=T+(H-T-B)*k/3, zv=lo+(hi-lo)*(1-k/3);
    const t=mk("text",{x:L-8,y:y+4,fill:css('--fraco'),"font-size":11,
      "text-anchor":"end"}); t.textContent=zv.toFixed(0); s.appendChild(t);
  }
  const prof=(zw-Math.min(...zs));
  el("tit-sec").textContent=`Seção transversal — RS ${D.rs[sel].toFixed(1)} km`
    +` · lâmina ${prof.toFixed(2)} m · cota ${zw.toFixed(2)} m`;
}

function hidro(){
  const s=el("hidro"); s.replaceChildren();
  const W=540,H=300,L=54,R=12,T=12,B=30;
  const q=D.ws.map((_,k)=>D.q[k][sel]);
  const n=q.length, hi=Math.max(...q,1), lo=Math.min(...q,0);
  const X=i=>L+(W-L-R)*i/((n-1)||1);
  const Y=v=>T+(H-T-B)*(1-(v-lo)/((hi-lo)||1));
  for(let k=0;k<=3;k++){
    const y=T+(H-T-B)*k/3;
    s.appendChild(mk("line",{x1:L,x2:W-R,y1:y,y2:y,stroke:css('--linha'),"stroke-width":1}));
    const t=mk("text",{x:L-8,y:y+4,fill:css('--fraco'),"font-size":11,
      "text-anchor":"end"}); t.textContent=(lo+(hi-lo)*(1-k/3)).toFixed(0); s.appendChild(t);
  }
  let d=`M ${X(0)} ${Y(q[0])}`;
  for(let i=1;i<n;i++) d+=` L ${X(i)} ${Y(q[i])}`;
  s.appendChild(mk("path",{d:d,fill:"none",stroke:css('--agua'),"stroke-width":2}));
  s.appendChild(mk("circle",{cx:X(ti),cy:Y(q[ti]),r:4,fill:css('--alerta')}));
  el("tit-hid").textContent=`Hidrograma — RS ${D.rs[sel].toFixed(1)} km`
    +` · ${q[ti].toFixed(0)} m³/s · pico ${hi.toFixed(0)} m³/s`;
}

function desenhar(){
  perfil(); secao(); hidro();
  el("rel").textContent=`hora ${D.horas[ti]} de ${D.horas[D.horas.length-1]}`;
}
el("t").max = D.horas.length-1;
el("t").oninput = e => { ti=+e.target.value; desenhar(); };
el("play").onclick = () => {
  tocando=!tocando;
  el("play").textContent = tocando ? "❚❚ Pausar" : "▶ Animar";
  if(tocando){ timer=setInterval(()=>{ ti=(ti+1)%D.horas.length;
    el("t").value=ti; desenhar(); },220); }
  else clearInterval(timer);
};
desenhar();
</script>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("projeto", nargs="?", default="Itajai_Rede_1983")
    ap.add_argument("--fonte", choices=["hecras", "motor"], default="hecras")
    ap.add_argument("--rio", default="Itajai_Acu")
    a = ap.parse_args()

    d = ler_motor(a.projeto) if a.fonte == "motor" else ler_hecras(a.projeto)
    dados = preparar(d, rio=a.rio)

    # data real do inicio, do proprio plano
    inicio = "?"
    p01 = f"{a.projeto}.p01"
    n_prev = d["ws"].shape[0]
    if os.path.exists(p01):
        for l in open(p01, encoding="utf-8", errors="ignore"):
            if l.startswith("Simulation Date="):
                p = l.split("=", 1)[1].strip().split(",")
                inicio = f"{p[0]} {p[1][:2]}:{p[1][2:]}"
                try:
                    d0 = datetime.datetime.strptime(p[0], "%d%b%Y")
                    d1 = datetime.datetime.strptime(p[2], "%d%b%Y")
                    n_prev = int((d1 - d0).total_seconds() // 3600
                                 + int(p[3][:2])) + 1
                except (ValueError, IndexError):
                    pass
    meta = {"estado": d["estado"], "fonte": d["fonte"], "inicio": inicio,
            "n_saidas": d["ws"].shape[0], "n_previstas": n_prev}

    # Grava DENTRO de app/, que e a pasta que rodar_modelo.sh sugere servir na
    # porta 8050. Antes a pagina saia na raiz do repositorio e o servidor nao
    # a enxergava -- a sugestao apontava para um lugar sem nada do que fora
    # gerado.
    os.makedirs("app", exist_ok=True)
    saida = os.path.join("app", f"{a.projeto}_cheia_{a.fonte}.html")
    with open(saida, "w", encoding="utf-8") as f:
        f.write(pagina(dados, a.projeto, meta))
    print(f"[OK] {saida}  ({os.path.getsize(saida)/1024:.0f} kB)")
    print(f"     {len(dados['rs'])} secoes x {len(dados['horas'])} instantes"
          f"  |  estado: {d['estado']}")


if __name__ == "__main__":
    main()
