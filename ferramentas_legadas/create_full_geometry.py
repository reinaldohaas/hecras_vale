"""
Gerador do modelo HEC-RAS 7.0.1 da Bacia do Rio Itajai-Acu (Vale do Itajai).

Rede 1D conectada por UMA juncao (confluencia):

    Itajai_Sul   (Trecho_Sul)   \\
    Itajai_Oeste (Trecho_Oeste)  >--[Confluencia]--> Itajai_Acu (Trecho_Principal) --> Mar
    Itajai_Norte (Trecho_Norte) /

As 3 barragens de contencao (Sul, Oeste, Norte) entram como Inline Structures
(vertedouro + comportas) quando INCLUDE_DAMS = True.

FORMATOS validados contra os projetos-exemplo oficiais do HEC-RAS
(neeraip/hecras-example-models):
  - series temporais em COLUNAS FIXAS de 8 caracteres, 10 por linha
    (causa raiz por que o modelo antigo nunca computava);
  - Boundary Location com 6 campos e padding fixo (Rio=16, Reach=16, RS=8);
  - Normal Depth = "Friction Slope=<decl>" (um unico valor);
  - condicao inicial unsteady = "Use Restart= 0" + "Initial Flow Loc=" por trecho;
  - juncao = "Up River,Reach="/"Dn River,Reach="/"Junc L&A=" (liga por NOME).
"""
import numpy as np

PROJECT = "Itajai_Bacia_Completa"
INCLUDE_DAMS = True         # rede conectada + 3 barragens (inline structures)
DAM_GATES    = True         # True = vertedouro + comportas; False = so vertedouro
NHOURS = 49                 # 48 h + hora 0
JUNCTION_BED = 50.0
DS_SLOPE = 0.0005           # declividade p/ Normal Depth de jusante
SEC_STEP = 2000.0
EDIT_TIME = "Node Last Edited Time= Aug/03/2026 00:00:00"

def p16(s):  return f"{s:<16}"
def f8(v):   return f"{v:8.3f}"

# --- afluentes: river, reach, L, desnivel(m), n_comportas, titulo, coord XY ---
TRIBS = [
    ("Itajai_Sul",   "Trecho_Sul",   60000.0,  90.0, 5, "Barragem Sul (Ituporanga)",   650000, 6970000),
    ("Itajai_Oeste", "Trecho_Oeste", 60000.0,  80.0, 7, "Barragem Oeste (Taio)",       600000, 7010000),
    ("Itajai_Norte", "Trecho_Norte", 60000.0, 100.0, 2, "Barragem Norte (Jose Boiteux)",640000, 7050000),
]
MAIN = ("Itajai_Acu", "Trecho_Principal", 100000.0)   # river, reach, L
JX, JY = 660000, 7000000                               # ponto da confluencia

TRIB_FLOW = {  # (base, pico) m3/s
    "Itajai_Sul":   (200.0, 1200.0),
    "Itajai_Oeste": (250.0, 1500.0),
    "Itajai_Norte": (300.0, 2000.0),
}
DAM_FRAC = 0.5   # barragem no meio de cada afluente
GATE_OPEN = 4.0  # abertura fixa das comportas (m)
GATE_NAME = "Comportas   "  # 12 chars (deve casar geometria x fluxo)


def dam_rs_of(L):
    """RS da barragem: meio de grade menos meio passo (fica entre 2 secoes)."""
    return round((L * DAM_FRAC) / SEC_STEP) * SEC_STEP - SEC_STEP / 2


def series(vals):
    """Serie temporal em colunas fixas de 8 caracteres, 10 por linha."""
    return "\n".join("".join(f"{v:8.2f}" for v in vals[i:i+10])
                     for i in range(0, len(vals), 10))


def hydrograph(base, peak, n=NHOURS, tp=18, te=40):
    v = []
    for h in range(n):
        if h <= tp:   q = base + (peak-base)*(h/tp)
        elif h <= te: q = peak - (peak-base)*((h-tp)/(te-tp))
        else:         q = base
        v.append(q)
    return v


def xsec(rs, zb, length):
    L = []
    L.append(f"Type RM Length L Ch R = 1 ,{rs:>10.2f} ,{length:>10.2f},{length:>10.2f},{length:>10.2f}")
    L.append(EDIT_TIME)
    pts = [(-400.0, zb+30), (-100.0, zb), (0.0, zb), (100.0, zb), (400.0, zb+30)]
    L.append("#Sta/Elev= 5 ")
    L.append("".join(f8(x)+f8(y) for x, y in pts))
    L.append("Bank Sta=-100,100")
    L.append("#Mann= 3 , -1 , 0 ")
    L.append(f8(-400)+f8(0.06)+f8(0)+f8(-100)+f8(0.035)+f8(0)+f8(100)+f8(0.06)+f8(0))
    return L


def dam(rs, zb, ngates, title):
    """Inline Structure (vertedouro + comportas), formato oficial HEC-RAS (SI)."""
    crest = round(zb + 15.0, 2)     # crista do vertedouro
    notch = round(zb + 12.0, 2)     # rebaixo central (secao vertente)
    ginv  = round(zb + 2.0, 2)      # soleira das comportas (acima do leito)
    gw, gh = 6.0, 8.0               # largura e altura de cada comporta
    spill_ht = round(crest - zb, 2)   # altura do vertedouro acima do leito
    # perfil espelha o exemplo oficial (Inline_3Gates): crista com rebaixo central
    weir_se = [(-400, crest), (-60, crest), (-50, notch),
               (50, notch), (60, crest), (400, crest)]
    L = []
    L.append(f"Type RM Length L Ch R = 5 ,{rs:<10.2f},,,")
    L.append("BEGIN DESCRIPTION:")
    L.append(title)
    L.append("END DESCRIPTION:")
    L.append(EDIT_TIME)
    L.append("IW Pilot Flow=0")
    L.append(f"#Inline Weir SE= {len(weir_se)} ")
    # perfil em COLUNAS FIXAS de 8 caracteres, max 10 valores por linha
    se_vals = [v for pair in weir_se for v in pair]
    for i in range(0, len(se_vals), 10):
        L.append("".join(f8(v) for v in se_vals[i:i+10]))
    L.append("IW Dist,WD,Coef,Skew,MaxSub,Min_El,Is_Ogee,SpillHt,DesHd")
    L.append(f"10,20,1.7,0,0.95,,-1 ,{spill_ht},3,2,2,")
    if DAM_GATES:
        L.append("IW Gate Name Wd,H,Inv,GCoef,Exp_T,Exp_O,Exp_H,Type,WCoef,Is_Ogee,SpillHt,DesHd,#Openings")
        L.append(f"{GATE_NAME},{gw},{gh},{ginv},0.68,0,1,0.5, 1 ,1.7,-1 ,{gh},3, {ngates} ,10,0.8, 0 ")
        xs = np.linspace(-((ngates-1)*7)/2.0, ((ngates-1)*7)/2.0, ngates)
        L.append("".join(f"{x:8.2f}" for x in xs))
    return L


def build_reach(river, reach, L, z_up, z_dn, x0, y0, x1, y1,
                dam_rs=None, ngates=0, title=""):
    out = [f"River Reach={p16(river)},{p16(reach)}"]
    out.append("Reach XY= 2 ")
    out.append(f"{x0:16.4f}{y0:16.4f}{x1:16.4f}{y1:16.4f}")
    def zb_at(station):
        frac = (station - SEC_STEP) / (L - SEC_STEP) if L != SEC_STEP else 0.0
        return z_dn + frac * (z_up - z_dn)
    rs = L
    while rs >= SEC_STEP - 0.1:
        ln = SEC_STEP if rs > SEC_STEP else 0.0
        out += xsec(rs, zb_at(rs), ln)
        # insere a barragem no RS unico entre esta secao e a proxima (jusante)
        if (INCLUDE_DAMS and dam_rs is not None
                and (rs - SEC_STEP) < dam_rs < rs):
            out += dam(dam_rs, zb_at(dam_rs), ngates, title)
        rs -= SEC_STEP
    return out


def write_geometry():
    g = ["Geom Title=Bacia do Itajai-Acu (rede 1D, 3 afluentes + confluencia)",
         "Program Version=7.01"]
    for river, reach, L, drop, ng, title, x0, y0 in TRIBS:
        dam_rs = dam_rs_of(L) if INCLUDE_DAMS else None
        g += build_reach(river, reach, L, JUNCTION_BED+drop, JUNCTION_BED,
                         x0, y0, JX, JY, dam_rs, ng, title)
    mr, mrc, mL = MAIN
    g += build_reach(mr, mrc, mL, JUNCTION_BED, 0.0, JX, JY, 730000, 7020000)

    # juncao (liga por NOME de rio/reach)
    g.append(f"Junct Name={p16('Confluencia')}")
    g.append("Junct Desc=Confluencia dos 3 afluentes no Itajai-Acu, 0 , 0 , 0 ,0")
    g.append(f"Junct X Y & Text X Y={JX},{JY},{JX+3000},{JY+3000}")
    for river, reach, *_ in TRIBS:
        g.append(f"Up River,Reach={p16(river)},{p16(reach)}")
    g.append(f"Dn River,Reach={p16(mr)},{p16(mrc)}")
    for _ in TRIBS:
        g.append("Junc L&A=1000,")

    with open(f"{PROJECT}.g01", "w") as f:
        f.write("\n".join(g) + "\n")
    print(f"[OK] {PROJECT}.g01  (dams={INCLUDE_DAMS})")


def write_unsteady():
    def bl(river, reach, rs):
        return (f"Boundary Location={p16(river)},{p16(reach)},{rs:<8}"
                f",        ,                ,                ")
    u = ["Flow Title=Cenario_Cheia_Bacia_Completa", "Program Version=7.01",
         "Use Restart= 0 "]
    # condicoes iniciais (uma por trecho, no RS de montante)
    for river, reach, L, *_ in TRIBS:
        base, _ = TRIB_FLOW[river]
        u.append(f"Initial Flow Loc={p16(river)},{p16(reach)},{L:<8.0f},{base:.0f}")
    mr, mrc, mL = MAIN
    u.append(f"Initial Flow Loc={p16(mr)},{p16(mrc)},{mL:<8.0f},750")
    # contornos de montante (hidrogramas de cheia)
    for river, reach, L, *_ in TRIBS:
        base, peak = TRIB_FLOW[river]
        u.append(bl(river, reach, f"{L:.2f}"))
        u.append("Interval=1HOUR")
        u.append(f"Flow Hydrograph= {NHOURS} ")
        u.append(series(hydrograph(base, peak)))
    # controle das comportas de cada barragem (abertura fixa)
    if INCLUDE_DAMS and DAM_GATES:
        for river, reach, L, *_ in TRIBS:
            u.append(bl(river, reach, f"{dam_rs_of(L):.2f}"))
            u.append(f"Gate Name={GATE_NAME}")
            u.append("Gate DSS Path=")
            u.append("Gate Use DSS=False")
            u.append("Gate Time Interval=1HOUR")
            u.append("Gate Use Fixed Start Time=False")
            u.append("Gate Fixed Start Date/Time=,")
            u.append(f"Gate Openings= {NHOURS} ")
            u.append(series([GATE_OPEN] * NHOURS))
    # contorno de jusante (profundidade normal)
    u.append(bl(mr, mrc, f"{SEC_STEP:.2f}"))
    u.append(f"Friction Slope={DS_SLOPE}")
    with open(f"{PROJECT}.u01", "w") as f:
        f.write("\n".join(u) + "\n")
    print(f"[OK] {PROJECT}.u01")


def write_plan():
    p = ["Plan Title=Simulacao_Bacia_Completa", "Program Version=7.01",
         "Short Identifier=BaciaFull",
         "Simulation Date=01SEP2008,0000,03SEP2008,0000",
         "Geom File=g01", "Flow File=u01", "Subcritical Flow",
         "Computation Interval=1MIN", "Output Interval=1HOUR",
         "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
         "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1", "Run RASMapper=0"]
    with open(f"{PROJECT}.p01", "w") as f:
        f.write("\n".join(p) + "\n")
    print(f"[OK] {PROJECT}.p01")


def write_prj():
    pr = [f"Proj Title={PROJECT}", "Current Plan=p01", "Default Exp/Contr=0.3,0.1",
          "SI Units", "Geom File=g01", "Unsteady File=u01", "Plan File=p01",
          "Y Axis Title=Elevation", "X Axis Title(PR)=Distance",
          "X Axis Title(CS)=Station"]
    with open(f"{PROJECT}.prj", "w") as f:
        f.write("\n".join(pr) + "\n")
    print(f"[OK] {PROJECT}.prj")


if __name__ == "__main__":
    write_prj(); write_geometry(); write_unsteady(); write_plan()
    print("Modelo da bacia completa gerado.")
