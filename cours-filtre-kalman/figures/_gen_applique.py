#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures des chapitres appliqués « vitesse depuis un inclinomètre ».
Toutes les trajectoires sont issues de vrais filtres de Kalman simulés
(pas de dessin à main levée). Unités : angle en degrés, vitesse en deg/s.

Palette cohérente avec style.css :
  a priori / modèle A   bleu   #3987e5
  mesure / brut         ambre  #f0a44f
  estimé / modèle B     vert   #35c98b
  vérité terrain        gris   #c3ccd8
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
C_A     = "#3987e5"   # premier filtre / a priori
C_MEAS  = "#f0a44f"   # mesure / dérivée brute
C_B     = "#35c98b"   # deuxième filtre / estimé retenu
C_TRUTH = "#c3ccd8"
C_GRID  = "#242b34"
C_AXIS  = "#3a4653"
C_DIM   = "#9aa6b6"
C_ACC   = "#34d1bf"
C_WARN  = "#f4c15b"
FONT = "-apple-system, Segoe UI, sans-serif"
MONO = "JetBrains Mono, monospace"

# ---------- SVG helpers ----------
def esc(s): return s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
def svg_open(w,h): return (f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" '
                           f'font-family="{FONT}" role="img">')
def line(x1,y1,x2,y2,st=C_GRID,w=1,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{st}" stroke-width="{w}"{d}/>'
def poly(pts,st,w=2.2,fill="none",op=1.0,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    p=" ".join(f"{x:.2f},{y:.2f}" for x,y in pts)
    return (f'<polyline points="{p}" fill="{fill}" stroke="{st}" stroke-width="{w}" '
            f'opacity="{op:.2f}" stroke-linejoin="round" stroke-linecap="round"{d}/>')
def polygon(pts,fill,op=1.0):
    p=" ".join(f"{x:.2f},{y:.2f}" for x,y in pts)
    return f'<polygon points="{p}" fill="{fill}" opacity="{op:.2f}"/>'
def txt(x,y,s,fill=C_DIM,size=13,anc="middle",wt="400",mono=False,style=""):
    fam=f' font-family="{MONO}"' if mono else ""
    st=f' font-style="{style}"' if style else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" fill="{fill}" font-size="{size}" '
            f'text-anchor="{anc}" font-weight="{wt}"{fam}{st}>{esc(s)}</text>')
def circ(x,y,r,fill,op=1.0,st="none",sw=0):
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" opacity="{op:.2f}" stroke="{st}" stroke-width="{sw}"/>'
def rrect(x,y,w,h,fill="none",st=C_AXIS,sw=1,rx=8,op=1.0):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" '
            f'fill="{fill}" stroke="{st}" stroke-width="{sw}" opacity="{op:.2f}"/>')
def arrow(x1,y1,x2,y2,st=C_DIM,w=1.8,head=6):
    ang=math.atan2(y2-y1,x2-x1); hx,hy=x2-head*math.cos(ang),y2-head*math.sin(ang)
    a1=ang+math.radians(150); a2=ang-math.radians(150)
    p=(f'{x2:.1f},{y2:.1f} {x2+head*math.cos(a1):.1f},{y2+head*math.sin(a1):.1f} '
       f'{x2+head*math.cos(a2):.1f},{y2+head*math.sin(a2):.1f}')
    return line(x1,y1,hx,hy,st,w)+f'<polygon points="{p}" fill="{st}"/>'
def save(name,body):
    with open(os.path.join(HERE,name),"w") as f: f.write(body+"</svg>\n")
    print("écrit",name)

# ---------- mini algèbre linéaire (listes) ----------
def eye(n): return [[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]
def mm(A,B):
    n,k,m=len(A),len(B),len(B[0])
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def mv(A,x): return [sum(A[i][j]*x[j] for j in range(len(x))) for i in range(len(A))]
def T(A): return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]
def add(A,B): return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

# LCG reproductible -> bruit ~[-1,1]
class RNG:
    def __init__(self,seed=12345): self.s=seed & ((1<<64)-1)
    def u(self):
        self.s=(self.s*6364136223846793005+1)&((1<<64)-1)
        return (self.s>>33)/(1<<31)-1.0
    def n(self):  # gaussien approx (somme de 3 uniformes)
        return (self.u()+self.u()+self.u())/1.732

# ---------- pas de filtre (mesure scalaire) ----------
def kf_predict(x,P,F,Q,B=None,u=None):
    x2=mv(F,x)
    if B is not None:
        x2=[x2[i]+sum(B[i][j]*u[j] for j in range(len(u))) for i in range(len(x2))]
    P2=add(mm(mm(F,P),T(F)),Q)
    return x2,P2
def kf_update(x,P,Hrow,R,z):
    n=len(x)
    Ph=[sum(P[i][j]*Hrow[j] for j in range(n)) for i in range(n)]     # P Hᵀ (colonne)
    S=sum(Hrow[i]*Ph[i] for i in range(n))+R
    K=[Ph[i]/S for i in range(n)]
    yv=z-sum(Hrow[i]*x[i] for i in range(n))
    x=[x[i]+K[i]*yv for i in range(n)]
    # Joseph : (I-KH)P(I-KH)ᵀ + K R Kᵀ
    IKH=[[ (1.0 if i==j else 0.0)-K[i]*Hrow[j] for j in range(n)] for i in range(n)]
    P=add(mm(mm(IKH,P),T(IKH)), [[K[i]*R*K[j] for j in range(n)] for i in range(n)])
    P=[[(P[i][j]+P[j][i])*0.5 for j in range(n)] for i in range(n)]
    return x,P,yv,S

# ---------- cadre de tracé temporel ----------
def frame(s,L,Rr,T_,B,tmax,ymin,ymax,xlabel,ylabel,tticks,yticks):
    def X(t): return L+t/tmax*(Rr-L)
    def Y(v): return B-(v-ymin)/(ymax-ymin)*(B-T_)
    for tt in tticks:
        s.append(line(X(tt),T_,X(tt),B,C_GRID,1)); s.append(txt(X(tt),B+18,f"{tt:g}",C_DIM,11,mono=True))
    for yy in yticks:
        s.append(line(L,Y(yy),Rr,Y(yy),C_GRID,1)); s.append(txt(L-8,Y(yy)+4,f"{yy:g}",C_DIM,11,anc="end",mono=True))
    s.append(line(L,Y(0) if ymin<0<ymax else B,Rr,Y(0) if ymin<0<ymax else B,C_AXIS,1.2))
    s.append(line(L,T_,L,B,C_AXIS,1.2))
    s.append(txt((L+Rr)/2,B+34,xlabel,C_DIM,12.5))
    s.append(txt(L-34,T_-10,ylabel,C_DIM,12,anc="start"))
    return X,Y

def legend(s,x,y,items):  # items: (color, dash, label)
    for (c,dash,lab) in items:
        s.append(line(x,y,x+22,y,c,2.6,dash=dash)); s.append(txt(x+28,y+4,lab,c,12,anc="start"))
        x+=34+len(lab)*7.0
    return x

# ===========================================================
# Figure A — dérivée brute vs Kalman
# ===========================================================
def fig_diff():
    W,H=720,360; L,Rr,T_,B=60,690,60,300
    dt=0.01; N=400; tmax=N*dt
    rng=RNG(7)
    # vérité : vitesse sinusoïdale
    th=0.0; om=lambda t: 30.0*math.sin(2*math.pi*0.35*t)
    truth_th=[]; truth_om=[]
    for k in range(N):
        t=k*dt; truth_th.append(th); truth_om.append(om(t)); th+=om(t)*dt
    sig=0.4
    meas=[truth_th[k]+sig*rng.n() for k in range(N)]
    # dérivée brute
    diff=[0.0]+[ (meas[k]-meas[k-1])/dt for k in range(1,N)]
    # KF ordre 1
    F=[[1,dt],[0,1]]; Hrow=[1,0]; R=sig*sig
    sa=90.0  # deg/s^2 white accel
    Q=[[sa*sa*dt**4/4, sa*sa*dt**3/2],[sa*sa*dt**3/2, sa*sa*dt**2]]
    x=[meas[0],0.0]; P=[[4,0],[0,400]]; kf_om=[]
    for k in range(N):
        x,P=kf_predict(x,P,F,Q); x,P,_,_=kf_update(x,P,Hrow,R,meas[k]); kf_om.append(x[1])
    s=[svg_open(W,H)]
    X,Y=frame(s,L,Rr,T_,B,tmax,-90,90,"temps (s)","deg/s",[0,1,2,3,4],[-60,0,60])
    s.append(poly([(X(k*dt),Y(max(-90,min(90,diff[k])))) for k in range(N)],C_MEAS,1.0,op=0.5))
    s.append(poly([(X(k*dt),Y(truth_om[k])) for k in range(N)],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(k*dt),Y(kf_om[k])) for k in range(N)],C_B,2.8))
    legend(s,L+6,T_-26,[(C_MEAS,None,"dérivée brute Δθ/Δt"),(C_TRUTH,"6 4","vitesse vraie"),(C_B,None,"Kalman")])
    save("fig-incl-diff.svg","".join(s))

# ===========================================================
# Figure B — ordre 1 (NCV) vs ordre 2 (NCA) sur une rampe de vitesse
# ===========================================================
def fig_order():
    W,H=720,360; L,Rr,T_,B=60,690,60,300
    dt=0.01; N=400; tmax=N*dt
    rng=RNG(11)
    # vérité : vitesse en rampe puis plateau
    def omf(t):
        if t<1.0: return 0.0
        if t<2.0: return 55.0*(t-1.0)
        return 55.0
    th=0.0; truth_th=[]; truth_om=[]
    for k in range(N):
        t=k*dt; truth_th.append(th); truth_om.append(omf(t)); th+=omf(t)*dt
    sig=0.4; meas=[truth_th[k]+sig*rng.n() for k in range(N)]; R=sig*sig
    # ordre 1
    F1=[[1,dt],[0,1]]; H1=[1,0]; sa=60.0
    Q1=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
    x1=[meas[0],0.0]; P1=[[4,0],[0,400]]; om1=[]
    # ordre 2
    F2=[[1,dt,dt*dt/2],[0,1,dt],[0,0,1]]; H2=[1,0,0]; sj=500.0
    Q2=[[sj*sj*dt**5/20,sj*sj*dt**4/8,sj*sj*dt**3/6],
        [sj*sj*dt**4/8, sj*sj*dt**3/3,sj*sj*dt**2/2],
        [sj*sj*dt**3/6, sj*sj*dt**2/2, sj*sj*dt]]
    x2=[meas[0],0.0,0.0]; P2=[[4,0,0],[0,400,0],[0,0,4000]]; om2=[]
    for k in range(N):
        x1,P1=kf_predict(x1,P1,F1,Q1); x1,P1,_,_=kf_update(x1,P1,H1,R,meas[k]); om1.append(x1[1])
        x2,P2=kf_predict(x2,P2,F2,Q2); x2,P2,_,_=kf_update(x2,P2,H2,R,meas[k]); om2.append(x2[1])
    s=[svg_open(W,H)]
    X,Y=frame(s,L,Rr,T_,B,tmax,-10,75,"temps (s)","deg/s",[0,1,2,3,4],[0,25,50])
    s.append(poly([(X(k*dt),Y(truth_om[k])) for k in range(N)],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(k*dt),Y(om1[k])) for k in range(N)],C_A,2.6))
    s.append(poly([(X(k*dt),Y(om2[k])) for k in range(N)],C_B,2.6))
    legend(s,L+6,T_-26,[(C_TRUTH,"6 4","vitesse vraie"),(C_A,None,"ordre 1 (θ,ω)"),(C_B,None,"ordre 2 (θ,ω,α)")])
    # zone rampe
    s.append(txt(X(1.5),B-6,"rampe",C_WARN,11,style="italic"))
    save("fig-incl-order.svg","".join(s))

# ===========================================================
# Figure C — retard capteur : modèle naïf vs modèle avec retard
# ===========================================================
def fig_lag():
    W,H=720,360; L,Rr,T_,B=60,690,60,300
    dt=0.01; N=500; tmax=N*dt
    rng=RNG(23)
    f0=0.35
    truth_om=[45.0*math.sin(2*math.pi*f0*k*dt) for k in range(N)]
    th=0.0; truth_th=[]
    for k in range(N): truth_th.append(th); th+=truth_om[k]*dt
    tau=0.20
    # sortie capteur : retard 1er ordre
    thm=truth_th[0]; sens=[]
    for k in range(N):
        thm=thm+dt*(truth_th[k]-thm)/tau; sens.append(thm)
    sig=0.35; meas=[sens[k]+sig*rng.n() for k in range(N)]; R=sig*sig
    # KF naïf : modèle [θ,ω], suppose z = θ
    F1=[[1,dt],[0,1]]; H1=[1,0]; sa=200.0
    Q1=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
    x1=[meas[0],0.0]; P1=[[4,0],[0,400]]; om1=[]
    # KF avec retard : état [θ,ω,θ_m], z = θ_m
    sl=350.0
    F2=[[1,dt,0],[0,1,0],[dt/tau,0,1-dt/tau]]; H2=[0,0,1]
    Q2=[[sl*sl*dt**4/4,sl*sl*dt**3/2,0],[sl*sl*dt**3/2,sl*sl*dt**2,0],[0,0,1e-4]]
    x2=[meas[0],0.0,meas[0]]; P2=[[4,0,0],[0,400,0],[0,0,4]]; om2=[]
    for k in range(N):
        x1,P1=kf_predict(x1,P1,F1,Q1); x1,P1,_,_=kf_update(x1,P1,H1,R,meas[k]); om1.append(x1[1])
        x2,P2=kf_predict(x2,P2,F2,Q2); x2,P2,_,_=kf_update(x2,P2,H2,R,meas[k]); om2.append(x2[1])
    s=[svg_open(W,H)]
    X,Y=frame(s,L,Rr,T_,B,tmax,-62,62,"temps (s)","deg/s",[0,1,2,3,4,5],[-40,0,40])
    s.append(poly([(X(k*dt),Y(truth_om[k])) for k in range(N)],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(k*dt),Y(om1[k])) for k in range(N)],C_A,2.4))
    s.append(poly([(X(k*dt),Y(om2[k])) for k in range(N)],C_B,2.6))
    legend(s,L+6,T_-26,[(C_TRUTH,"6 4","vitesse vraie"),(C_A,None,"sans modèle de retard"),(C_B,None,"retard modélisé")])
    save("fig-incl-lag.svg","".join(s))

# ===========================================================
# Figure D — commande en anticipation (feedforward)
# ===========================================================
def fig_command():
    W,H=720,360; L,Rr,T_,B=60,690,60,300
    dt=0.01; N=400; tmax=N*dt
    rng=RNG(31)
    # commande = accélération angulaire connue (créneaux)
    def uf(t):
        if t<0.8: return 0.0
        if t<1.4: return 120.0
        if t<2.2: return 0.0
        if t<2.8: return -120.0
        return 0.0
    th=0.0; om=0.0; truth_th=[]; truth_om=[]
    for k in range(N):
        t=k*dt; truth_th.append(th); truth_om.append(om)
        a=uf(t)+8.0*rng.n()  # perturbation légère
        th+=om*dt+0.5*a*dt*dt; om+=a*dt
    sig=0.4; meas=[truth_th[k]+sig*rng.n() for k in range(N)]; R=sig*sig
    F=[[1,dt],[0,1]]; Hrow=[1,0]
    Bm=[[dt*dt/2],[dt]]
    # sans commande : gros Q pour tolérer l'accélération
    sa=180.0; Qn=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
    xn=[meas[0],0.0]; Pn=[[4,0],[0,400]]; omn=[]
    # avec commande : petit Q (le mouvement est expliqué par u)
    sc=15.0; Qc=[[sc*sc*dt**4/4,sc*sc*dt**3/2],[sc*sc*dt**3/2,sc*sc*dt**2]]
    xc=[meas[0],0.0]; Pc=[[4,0],[0,400]]; omc=[]
    for k in range(N):
        t=k*dt
        xn,Pn=kf_predict(xn,Pn,F,Qn); xn,Pn,_,_=kf_update(xn,Pn,Hrow,R,meas[k]); omn.append(xn[1])
        xc,Pc=kf_predict(xc,Pc,F,Qc,Bm,[uf(t)]); xc,Pc,_,_=kf_update(xc,Pc,Hrow,R,meas[k]); omc.append(xc[1])
    s=[svg_open(W,H)]
    X,Y=frame(s,L,Rr,T_,B,tmax,-80,80,"temps (s)","deg/s",[0,1,2,3,4],[-50,0,50])
    s.append(poly([(X(k*dt),Y(truth_om[k])) for k in range(N)],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(k*dt),Y(omn[k])) for k in range(N)],C_A,2.4))
    s.append(poly([(X(k*dt),Y(omc[k])) for k in range(N)],C_B,2.6))
    legend(s,L+6,T_-26,[(C_TRUTH,"6 4","vitesse vraie"),(C_A,None,"sans commande"),(C_B,None,"commande B·u")])
    save("fig-incl-command.svg","".join(s))

# ===========================================================
# Figure E — cohérence NIS pour trois réglages de Q
# ===========================================================
def fig_nis():
    W,H=720,340; L,Rr,T_,B=60,690,45,275
    dt=0.01; N=500; tmax=N*dt
    rng=RNG(41)
    def omf(t): return 25.0*math.sin(2*math.pi*0.3*t)+10.0*math.sin(2*math.pi*0.11*t)
    th=0.0; truth_th=[]
    for k in range(N):
        t=k*dt; truth_th.append(th); th+=omf(t)*dt
    sig=0.4; meas=[truth_th[k]+sig*rng.n() for k in range(N)]; R=sig*sig
    F=[[1,dt],[0,1]]; Hrow=[1,0]
    def run(sa):
        Q=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
        x=[meas[0],0.0]; P=[[4,0],[0,400]]; nis=[]
        for k in range(N):
            x,P=kf_predict(x,P,F,Q); x,P,yv,S=kf_update(x,P,Hrow,R,meas[k]); nis.append(yv*yv/S)
        # moyenne glissante W
        Wd=45; avg=[]
        for k in range(N):
            a=max(0,k-Wd+1); avg.append(sum(nis[a:k+1])/(k-a+1))
        return avg
    small=run(20.0); good=run(33.0); big=run(60.0)
    s=[svg_open(W,H)]
    X,Y=frame(s,L,Rr,T_,B,tmax,0,5.0,"temps (s)","NIS (moy. glissante)",[0,1,2,3,4,5],[0,1,2,3,4])
    # bande d'acceptation chi2(W)/W ~ [0.63, 1.45] pour W=45, 95%
    lo,hi=0.63,1.45
    s.append(polygon([(X(0),Y(lo)),(X(tmax),Y(lo)),(X(tmax),Y(hi)),(X(0),Y(hi))],C_ACC,0.12))
    s.append(line(X(0),Y(1.0),X(tmax),Y(1.0),C_ACC,1.2,dash="4 4"))
    s.append(txt(Rr-6,Y(1.0)-6,"cible = 1",C_ACC,11,anc="end"))
    s.append(txt(L+10,Y(hi)-5,"bande de cohérence 95 %",C_ACC,10.5,anc="start",style="italic"))
    s.append(poly([(X(k*dt),Y(min(4.9,small[k]))) for k in range(N)],C_MEAS,2.4))
    s.append(poly([(X(k*dt),Y(good[k])) for k in range(N)],C_B,2.6))
    s.append(poly([(X(k*dt),Y(big[k])) for k in range(N)],C_A,2.4))
    legend(s,L+6,T_-22,[(C_MEAS,None,"Q trop petit"),(C_B,None,"Q correct"),(C_A,None,"Q trop grand")])
    save("fig-incl-nis.svg","".join(s))

# ===========================================================
# Figure F — structure sur PLC : observateur à gain constant
# ===========================================================
def fig_plc():
    W,H=720,360
    s=[svg_open(W,H)]
    # OB cyclique (cadre englobant)
    s.append(rrect(30,50,660,210,"none",C_AXIS,1.4,12))
    s.append(txt(46,40,"OB cyclique (période dt fixe)  —  appelé 1×/cycle",C_DIM,12.5,anc="start",wt="600"))
    # FB
    fx,fy,fw,fh=250,80,300,138
    s.append(rrect(fx,fy,fw,fh,"#10201d",C_B,1.6,10))
    s.append(txt(fx+fw/2,fy+24,'FB "IncliVitesse"',C_B,14,wt="700"))
    s.append(txt(fx+fw/2,fy+50,"x⁻ = F·x + B·u",  "#cfe3ff",13,mono=True))
    s.append(txt(fx+fw/2,fy+72,"y  = z − H·x⁻",   "#d6f5e6",13,mono=True))
    s.append(txt(fx+fw/2,fy+94,"x  = x⁻ + K∞·y",  "#d6f5e6",13,mono=True))
    s.append(txt(fx+fw/2,fy+120,"gains K∞ = constantes (α, β)","#9fe9cf",11.5,style="italic"))
    # entrées
    s.append(arrow(70,130,fx,130,C_MEAS,2,8)); s.append(txt(72,118,"z  inclinomètre",C_MEAS,12,anc="start",wt="600"))
    s.append(arrow(70,185,fx,185,C_ACC,2,8));  s.append(txt(72,205,"u  commande",C_ACC,12,anc="start",wt="600"))
    # sortie
    s.append(arrow(fx+fw,150,660,150,C_B,2,8)); s.append(txt(598,138,"ω estimé",C_B,12.5,wt="700"))
    s.append(txt(598,166,"→ régulation",C_DIM,11,style="italic"))
    # bouclage état (par le bas, à l'intérieur du cadre OB)
    yb=fy+fh+20
    s.append(line(fx+fw-45,fy+fh,fx+fw-45,yb,C_DIM,1.6))
    s.append(line(fx+fw-45,yb,fx+45,yb,C_DIM,1.6))
    s.append(arrow(fx+45,yb,fx+45,fy+fh,C_DIM,1.6,7))
    s.append(txt(fx+fw/2,yb+16,"état x mémorisé (statique du FB) → cycle suivant",C_DIM,10.5,style="italic"))
    # encart offline (sous le cadre OB, pleine largeur)
    s.append(rrect(30,290,660,40,"#151a22",C_WARN,1.2,8))
    s.append(txt(50,308,"Hors ligne (une seule fois) :",C_WARN,12,anc="start",wt="700"))
    s.append(txt(50,324,"résoudre l'équation de Riccati (DARE) → gain constant K∞. Aucune matrice à inverser en ligne.",C_DIM,11,anc="start"))
    save("fig-incl-plc.svg","".join(s))

if __name__=="__main__":
    fig_diff(); fig_order(); fig_lag(); fig_command(); fig_nis(); fig_plc()
    print("Figures appliquées générées.")
