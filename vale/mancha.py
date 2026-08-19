# -*- coding: utf-8 -*-
"""
Mancha de inundacao animada, em planta, a partir do resultado do solver.

DE ONDE SAI A MANCHA. Nao ha mapa de inundacao pronto no resultado 1D: o que o
solver devolve e a COTA DA LAMINA em cada secao, a cada hora. A mancha e
construida aqui, e a construcao e o conteudo:

    para cada secao e cada instante, procura-se ate onde a cota da lamina
    alcanca o perfil transversal, a partir do talvegue e para os dois lados.
    Isso da duas ESTACAS (esquerda e direita) que, projetadas na linha de
    corte da secao, viram dois pontos no terreno. Costurando esses pontos
    entre secoes vizinhas sai o poligono molhado.

PARAR NA PRIMEIRA SUBIDA, e nao pegar tudo que esta abaixo da cota. Um vale
vizinho mais baixo que a lamina apareceria molhado sem ter ligacao com o rio --
a agua nao pula a encosta. A varredura sai do talvegue e para no primeiro ponto
acima da lamina, que e o que a hidraulica 1D de fato representa.

O QUE ESTA MANCHA NAO E. Ela e a extensao lateral que o modelo 1D admite, e a
1D so sabe da secao onde ha secao: entre duas cutlines a borda e interpolada em
linha reta. Nao substitui mapeamento 2D, e sobre Copernicus nem chega perto --
o terreno tem o dossel dentro. Serve para ver POR ONDE a cheia passa e QUANDO,
que e a pergunta que a pagina antiga nao respondia.
"""
import json
import os

import numpy as np


def _por_secao(hdf):
    """Lamina por secao e por instante, do .p01.hdf. Devolve (chaves, wse, t)."""
    import h5py

    f = h5py.File(hdf, "r")
    base = ("Results/Unsteady/Output/Output Blocks/Base Output/"
            "Unsteady Time Series")
    g = f[base + "/Cross Sections"]
    wse = np.asarray(g["Water Surface"][:], float)
    t = [s.decode("latin-1", "replace") for s in
         f[base + "/Time Date Stamp"][:]]
    att = g["Cross Section Attributes"][:]
    chaves = []
    for a in att:
        v = [x.decode("latin-1", "replace").strip() if isinstance(x, bytes)
             else str(x).strip() for x in a]
        rio = v[0] if v else ""
        rs = next((x for x in reversed(v) if x.replace(".", "").isdigit()), "")
        chaves.append((rio, rs))
    f.close()
    return chaves, wse, t


def _molhado(sta, z, i_thal, cota):
    """Estacas alcancadas pela lamina, saindo do talvegue para os dois lados."""
    n = len(z)
    if cota <= z[i_thal]:
        return None
    e = i_thal
    while e > 0 and z[e - 1] <= cota:
        e -= 1
    dd = i_thal
    while dd < n - 1 and z[dd + 1] <= cota:
        dd += 1
    if e == dd:
        return None
    # interpola a borda entre o ultimo ponto seco e o primeiro molhado
    def borda(i, j):
        if i == j:
            return float(sta[i])
        dz = z[i] - z[j]
        f = 0.0 if abs(dz) < 1e-9 else (cota - z[j]) / dz
        f = float(np.clip(f, 0.0, 1.0))
        return float(sta[j] + f * (sta[i] - sta[j]))
    return (borda(e - 1, e) if e > 0 else float(sta[0]),
            borda(dd + 1, dd) if dd < n - 1 else float(sta[-1]))


def extrair(estado, hdf, max_secoes=600, log=print):
    """Bordas molhadas por secao e por hora, prontas para a pagina."""
    chaves, wse, t = _por_secao(hdf)
    idx = {(r, s): i for i, (r, s) in enumerate(chaves)}
    secoes, faltando = [], 0
    for rio, v in (estado.get("xs_pronto") or {}).items():
        w = sorted(v, key=lambda d: -d["rs"])
        passo = max(1, len(w) // max(1, max_secoes // max(len(estado["xs_pronto"]), 1)))
        for d in w[::passo]:
            i = idx.get((rio, f"{d['rs']:.2f}")) or idx.get((rio, f"{d['rs']:g}"))
            if i is None:
                faltando += 1
                continue
            secoes.append((d, i))
    if faltando:
        log(f"      {faltando} secoes sem correspondencia no resultado")
    if not secoes:
        return None

    n_t = wse.shape[0]
    saida = []
    for d, i in secoes:
        sta = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        it = int(d.get("i_thal", int(np.argmin(z))))
        x0, y0, x1, y1 = [float(c) for c in d["cut"]]
        larg = float(sta[-1] - sta[0]) or 1.0
        esq, dir_ = [], []
        for k in range(n_t):
            m = _molhado(sta, z, it, float(wse[k, i]))
            if m is None:
                esq.append(-1)
                dir_.append(-1)
            else:
                esq.append(int(round(1000 * (m[0] - sta[0]) / larg)))
                dir_.append(int(round(1000 * (m[1] - sta[0]) / larg)))
        saida.append({"rio": d["rio"], "rs": round(float(d["rs"]), 1),
                      "p0": [round(x0, 1), round(y0, 1)],
                      "p1": [round(x1, 1), round(y1, 1)],
                      "z": round(float(z[it]), 2),
                      "e": esq, "d": dir_})
    log(f"      mancha: {len(saida)} secoes x {n_t} instantes")
    return {"secoes": saida, "tempos": t,
            "wse_max": float(np.nanmax(wse)) if wse.size else 0.0}


# ---------------------------------------------------------------- PAGINA
_JS = r"""
const D = DADOS;
const S = D.secoes, T = D.tempos, NT = T.length;
const porRio = {};
S.forEach((s,i)=>{ (porRio[s.rio] ||= []).push(i); });
let bb=[1e18,1e18,-1e18,-1e18];
S.forEach(s=>{ [s.p0,s.p1].forEach(p=>{
  bb[0]=Math.min(bb[0],p[0]); bb[1]=Math.min(bb[1],p[1]);
  bb[2]=Math.max(bb[2],p[0]); bb[3]=Math.max(bb[3],p[1]); }); });
const cv=document.getElementById('mapa'), cx=cv.getContext('2d');
let esc=1, ox=0, oy=0;
function ajustar(){
  const r=cv.getBoundingClientRect(), dpr=window.devicePixelRatio||1;
  cv.width=r.width*dpr; cv.height=r.height*dpr; cx.setTransform(dpr,0,0,dpr,0,0);
  const w=bb[2]-bb[0], h=bb[3]-bb[1], m=18;
  esc=Math.min((r.width-2*m)/w,(r.height-2*m)/h);
  ox=m+((r.width-2*m)-w*esc)/2; oy=m+((r.height-2*m)-h*esc)/2;
}
const PX=p=>ox+(p[0]-bb[0])*esc, PY=p=>cv.getBoundingClientRect().height-(oy+(p[1]-bb[1])*esc);
function ponto(s,f){ // f em milesimos da largura da cutline
  const t=f/1000; return [s.p0[0]+t*(s.p1[0]-s.p0[0]), s.p0[1]+t*(s.p1[1]-s.p0[1])];
}
function desenhar(k){
  const r=cv.getBoundingClientRect();
  cx.clearRect(0,0,r.width,r.height);
  cx.fillStyle=getComputedStyle(document.body).getPropertyValue('--mapa');
  cx.fillRect(0,0,r.width,r.height);
  for(const rio in porRio){
    const idx=porRio[rio];
    // eixo do rio: meio de cada cutline
    cx.beginPath();
    idx.forEach((i,j)=>{ const s=S[i], p=ponto(s,500);
      j?cx.lineTo(PX(p),PY(p)):cx.moveTo(PX(p),PY(p)); });
    cx.strokeStyle=getComputedStyle(document.body).getPropertyValue('--rio');
    cx.lineWidth=1; cx.stroke();
    // mancha: costura as bordas esquerdas na ida e as direitas na volta
    let run=[];
    const fechar=()=>{ if(run.length<2){run=[];return;}
      cx.beginPath();
      run.forEach((i,j)=>{ const s=S[i], p=ponto(s,s.e[k]);
        j?cx.lineTo(PX(p),PY(p)):cx.moveTo(PX(p),PY(p)); });
      for(let j=run.length-1;j>=0;j--){ const s=S[run[j]], p=ponto(s,s.d[k]);
        cx.lineTo(PX(p),PY(p)); }
      cx.closePath();
      cx.fillStyle=getComputedStyle(document.body).getPropertyValue('--agua');
      cx.fill();
      cx.strokeStyle=getComputedStyle(document.body).getPropertyValue('--borda');
      cx.lineWidth=0.6; cx.stroke(); run=[]; };
    idx.forEach(i=>{ (S[i].e[k]>=0) ? run.push(i) : fechar(); });
    fechar();
  }
}
let k=0, tocando=false, timer=null;
const sl=document.getElementById('t'), rot=document.getElementById('quando'),
      bt=document.getElementById('play');
sl.max=NT-1;
function ir(v){ k=Math.max(0,Math.min(NT-1,v|0)); sl.value=k;
  rot.textContent=T[k]+'   —   hora '+k+' de '+(NT-1); desenhar(k); }
sl.addEventListener('input',e=>ir(+e.target.value));
function parar(){ tocando=false; clearInterval(timer); bt.textContent='▶ tocar'; }
bt.addEventListener('click',()=>{ if(tocando){parar();return;}
  tocando=true; bt.textContent='‖ pausar';
  timer=setInterval(()=>{ ir(k>=NT-1?0:k+1); }, 90); });
addEventListener('resize',()=>{ajustar();desenhar(k);});
addEventListener('keydown',e=>{ if(e.key==='ArrowRight')ir(k+1);
  if(e.key==='ArrowLeft')ir(k-1); if(e.key===' '){e.preventDefault();bt.click();} });
ajustar(); ir(0);
"""


def pagina(dados, destino, titulo="", resumo=""):
    """Grava a pagina com a mancha animada."""
    css = """
:root{--bg:#f7f7f5;--tinta:#1a1a18;--fraco:#6b6b66;--linha:#dcdcd6;
      --mapa:#eeeee9;--rio:#3a6ea5;--agua:rgba(58,110,165,.42);--borda:#2c5f8f}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#16161a;--tinta:#ececeb;--fraco:#9a9a94;--linha:#2e2e34;
  --mapa:#1e1e24;--rio:#7fb2e5;--agua:rgba(127,178,229,.38);--borda:#9ec8ef}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--tinta);
  font:15px/1.5 ui-sans-serif,system-ui,-apple-system,Segoe UI,sans-serif}
.env{max-width:1120px;margin:0 auto;padding:28px 22px 40px}
h1{font-size:1.45rem;margin:0 0 4px;letter-spacing:-.01em}
.sub{color:var(--fraco);font-size:.9rem;margin:0 0 20px}
#mapa{width:100%;height:min(64vh,620px);border:1px solid var(--linha);
  border-radius:8px;display:block}
.ctrl{display:flex;gap:14px;align-items:center;margin-top:14px;flex-wrap:wrap}
button{font:inherit;padding:7px 16px;border:1px solid var(--linha);
  border-radius:6px;background:var(--mapa);color:var(--tinta);cursor:pointer}
button:hover{border-color:var(--rio)}
input[type=range]{flex:1;min-width:220px;accent-color:var(--rio)}
#quando{font-variant-numeric:tabular-nums;color:var(--fraco);font-size:.88rem;
  min-width:250px}
.nota{margin-top:22px;padding-top:16px;border-top:1px solid var(--linha);
  color:var(--fraco);font-size:.85rem;max-width:70ch}
"""
    html = (
        "<!doctype html>\n<html lang=\"pt-BR\"><head><meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n"
        f"<title>{titulo or 'Mancha de inundacao'}</title>\n"
        f"<style>{css}</style></head><body>\n<div class=\"env\">\n"
        f"<h1>{titulo or 'Mancha de inundacao'}</h1>\n"
        f"<p class=\"sub\">{resumo}</p>\n"
        "<canvas id=\"mapa\"></canvas>\n"
        "<div class=\"ctrl\"><button id=\"play\">▶ tocar</button>"
        "<input type=\"range\" id=\"t\" min=\"0\" value=\"0\">"
        "<span id=\"quando\"></span></div>\n"
        "<p class=\"nota\">A mancha e construida da cota da lamina calculada em "
        "cada secao: a varredura sai do talvegue para os dois lados e para no "
        "primeiro ponto acima da agua, entao vale vizinho mais baixo nao "
        "aparece molhado sem ligacao com o rio. Entre duas secoes a borda e "
        "reta &mdash; e o que a hidraulica 1D representa. Setas movem uma hora; "
        "espaco toca e pausa.</p>\n</div>\n<script>\nconst DADOS = "
        + json.dumps(dados, separators=(",", ":")) + ";\n"
        + _JS.replace("DADOS;", "DADOS;", 1) + "\n</script>\n</body></html>\n")
    os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
    with open(destino, "w", encoding="utf-8") as f:
        f.write(html)
    return destino


def gerar(estado, hdf, destino, titulo="", resumo="", log=print):
    d = extrair(estado, hdf, log=log)
    if d is None:
        log("      sem dados para a mancha")
        return None
    return pagina(d, destino, titulo, resumo)
