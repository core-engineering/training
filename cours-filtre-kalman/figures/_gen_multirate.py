#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures du chapitre « fusion multi-fréquences ».
Réutilise les briques de _gen_applique.py (KF simulé, helpers SVG).
Toutes les courbes sont issues de vrais filtres de Kalman simulés.
"""
import math
from _gen_applique import (svg_open, line, poly, polygon, txt, circ, rrect, arrow,
                           save, frame, legend, RNG, kf_predict, kf_update,
                           C_A, C_MEAS, C_B, C_TRUTH, C_GRID, C_AXIS, C_DIM, C_ACC, C_WARN)

C_SLOW = "#c07cf0"   # capteur lent / accompagnement (violet)

# ------- scénario commun (vérité + mesures) -------
def scenario():
    dt=0.01; N=600
    rng=RNG(101)
    f0=0.30
    truth_th=[20.0*math.sin(2*math.pi*f0*k*dt) for k in range(N)]
    truth_om=[20.0*2*math.pi*f0*math.cos(2*math.pi*f0*k*dt) for k in range(N)]
    sig_f=0.6; sig_s=0.06
    drop=set(range(250,330))                    # capteur rapide absent
    fast=[(None if k in drop else truth_th[k]+sig_f*rng.n()) for k in range(N)]
    slow=[(truth_th[k]+sig_s*rng.n() if k%25==0 else None) for k in range(N)]
    return dt,N,truth_th,truth_om,fast,slow,sig_f,sig_s,drop

def run_filter(use_fast, use_slow):
    dt,N,truth_th,truth_om,fast,slow,sig_f,sig_s,drop = scenario()
    F=[[1,dt],[0,1]]; H=[1,0]
    sa=120.0
    Q=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
    x=[0.0,0.0]; P=[[25.0,0],[0,900.0]]
    om=[]; th=[]; sig=[]
    for k in range(N):
        x,P=kf_predict(x,P,F,Q)
        if use_fast and fast[k] is not None:
            x,P,_,_=kf_update(x,P,H,sig_f*sig_f,fast[k])
        if use_slow and slow[k] is not None:
            x,P,_,_=kf_update(x,P,H,sig_s*sig_s,slow[k])
        om.append(x[1]); th.append(x[0]); sig.append(math.sqrt(max(P[0][0],0)))
    return om,th,sig

# ===========================================================
# Figure A — le rythme : prédire toujours, corriger si dispo
# ===========================================================
def fig_schedule():
    W,H=720,300; L,Rr=175,688
    s=[svg_open(W,H)]
    tmax=6.0
    def X(t): return L+t/tmax*(Rr-L)
    rows=[(70,"prédiction","chaque cycle",C_ACC),
          (140,"capteur rapide","100 Hz · bruité",C_A),
          (210,"capteur lent","4 Hz · précis",C_SLOW)]
    for (yy,name,desc,c) in rows:
        s.append(line(L,yy,Rr,yy,C_AXIS,1.2))
        s.append(txt(L-14,yy-4,name,c,12.5,anc="end",wt="700"))
        s.append(txt(L-14,yy+13,desc,C_DIM,10.5,anc="end"))
    # prédiction : ticks denses
    k=0
    while k<=120:
        xx=X(k*0.05)
        s.append(line(xx,70-6,xx,70+6,C_ACC,1.4)); k+=1
    # capteur rapide : un tick par (petit) pas, trou de dropout t in [2.5,3.3]
    t=0.0
    while t<=6.0:
        if not (2.5<=t<3.3):
            s.append(circ(X(t),140,2.4,C_A))
        t+=0.05
    s.append(rrect(X(2.5),128,X(3.3)-X(2.5),24,"#2a1f14",C_WARN,1.2,5,op=0.5))
    s.append(txt((X(2.5)+X(3.3))/2,120,"absent (panne / trame perdue)",C_WARN,10.5))
    # capteur lent : tous les 0.25 s
    t=0.0
    while t<=6.0:
        s.append(circ(X(t),210,3.6,C_SLOW)); t+=0.25
    # axe temps
    for gt in range(0,7):
        s.append(line(X(gt),250,X(gt),256,C_AXIS,1)); s.append(txt(X(gt),272,f"{gt}",C_DIM,11,mono=True))
    s.append(txt((L+Rr)/2,292,"temps (s)",C_DIM,12))
    # annotations
    s.append(txt(X(0.0),44,"prédire TOUJOURS · corriger SI une mesure est présente",C_ACC,12.5,anc="start",wt="600"))
    save("fig-multi-schedule.svg","".join(s))

# ===========================================================
# Figure B — la covariance respire au rythme des mesures (sawtooth)
# ===========================================================
def fig_band():
    dt,N,truth_th,truth_om,fast,slow,sig_f,sig_s,drop = scenario()
    _,_,sig = run_filter(True,True)
    W,H=720,340; L,Rr,T_,B=62,690,45,275
    s=[svg_open(W,H)]
    tmax=N*dt
    X,Y=frame(s,L,Rr,T_,B,tmax,0,1.2,"temps (s)","incertitude √P₀₀  (°)",[0,1,2,3,4,5,6],[0,0.5,1.0])
    # zone dropout
    s.append(polygon([(X(250*dt),T_),(X(330*dt),T_),(X(330*dt),B),(X(250*dt),B)],C_WARN,0.10))
    s.append(txt((X(250*dt)+X(330*dt))/2,T_-6,"capteur rapide absent",C_WARN,11))
    # marques des mises à jour lentes
    for k in range(N):
        if k%25==0: s.append(line(X(k*dt),B,X(k*dt),B-5,C_SLOW,1.4))
    s.append(poly([(X(k*dt),Y(min(1.2,sig[k]))) for k in range(N)],C_B,2.6))
    s.append(txt(X(0.15),Y(0.9),"√P₀₀ = écart-type de l'estimation d'angle",C_B,11.5,anc="start"))
    # légende marque lente
    s.append(line(L+8,T_-24,L+30,T_-24,C_SLOW,2)); s.append(txt(L+36,T_-20,"mise à jour lente (précise) → chute de P",C_SLOW,11,anc="start"))
    save("fig-multi-band.svg","".join(s))

# ===========================================================
# Figure C — fusion vs source unique (vitesse estimée)
# ===========================================================
def fig_fusion():
    dt,N,truth_th,truth_om,fast,slow,sig_f,sig_s,drop = scenario()
    om_fu,_,_ = run_filter(True,True)
    om_fa,_,_ = run_filter(True,False)
    om_sl,_,_ = run_filter(False,True)
    W,H=720,350; L,Rr,T_,B=62,690,58,290
    s=[svg_open(W,H)]
    tmax=N*dt
    X,Y=frame(s,L,Rr,T_,B,tmax,-60,60,"temps (s)","vitesse ω  (°/s)",[0,1,2,3,4,5,6],[-40,0,40])
    s.append(polygon([(X(250*dt),T_),(X(330*dt),T_),(X(330*dt),B),(X(250*dt),B)],C_WARN,0.08))
    s.append(txt((X(250*dt)+X(330*dt))/2,B-8,"dropout rapide",C_WARN,10.5))
    s.append(poly([(X(k*dt),Y(truth_om[k])) for k in range(N)],C_TRUTH,2.2,dash="6 4"))
    s.append(poly([(X(k*dt),Y(max(-60,min(60,om_sl[k])))) for k in range(N)],C_SLOW,1.8,op=0.9))
    s.append(poly([(X(k*dt),Y(max(-60,min(60,om_fa[k])))) for k in range(N)],C_A,1.8,op=0.9))
    s.append(poly([(X(k*dt),Y(om_fu[k])) for k in range(N)],C_B,2.8))
    legend(s,L+4,T_-26,[(C_TRUTH,"6 4","vraie"),(C_SLOW,None,"lent seul"),(C_A,None,"rapide seul"),(C_B,None,"fusion")])
    save("fig-multi-fusion.svg","".join(s))

# ===========================================================
# Figure D — structure sur PLC : prédire + corrections conditionnelles
# ===========================================================
def fig_plc():
    W,H=720,360
    s=[svg_open(W,H)]
    s.append(rrect(28,48,664,232,"none",C_AXIS,1.4,12))
    s.append(txt(44,40,"OB cyclique rapide — 1×/cycle",C_DIM,12.5,anc="start",wt="600"))
    # bloc prédiction
    px,py,pw,ph=60,80,180,70
    s.append(rrect(px,py,pw,ph,"#122033",C_ACC,1.5,10))
    s.append(txt(px+pw/2,py+27,"PRÉDIRE",C_ACC,13,wt="700"))
    s.append(txt(px+pw/2,py+48,"x⁻=F·x+B·u ; P⁻=…","#cfe3ff",11,mono=True))
    s.append(txt(px+pw/2,py+62,"(toujours)",C_DIM,10.5,style="italic"))
    # deux corrections conditionnelles
    def corr(x,y,name,rate,col):
        s.append(rrect(x,y,190,58,"#10201d",col,1.4,10))
        s.append(txt(x+95,y+20,name,col,12.5,wt="700"))
        s.append(txt(x+95,y+38,"si dispo → update scalaire",C_DIM,10.5,style="italic"))
        s.append(txt(x+95,y+52,rate,col,10.5))
    corr(300,72,"CORRIGER · rapide","100 Hz, R_f",C_A)
    corr(300,150,"CORRIGER · lent","4 Hz, R_s",C_SLOW)
    # flèches
    s.append(arrow(px+pw,115,300,101,C_DIM,1.8,7))
    s.append(arrow(px+pw,120,300,179,C_DIM,1.8,7))
    s.append(arrow(495,101,540,101,C_DIM,1.8,7))
    s.append(arrow(495,179,540,140,C_DIM,1.8,7))
    s.append(rrect(540,86,120,70,"#10201d",C_B,1.5,10))
    s.append(txt(600,116,"x, P",C_B,14,wt="700",mono=True))
    s.append(txt(600,138,"ω estimé",C_B,11.5))
    # entrées capteurs
    s.append(txt(300,64,"drapeaux de disponibilité par capteur",C_DIM,10.5,anc="start",style="italic"))
    # note périodique
    s.append(rrect(28,296,664,44,"#151a22",C_WARN,1.2,8))
    s.append(txt(46,314,"Rapport de fréquences entier et fixe :",C_WARN,12,anc="start",wt="700"))
    s.append(txt(46,330,"le gain devient périodique → on peut précalculer un jeu de gains constants (un par phase du cycle).",C_DIM,11,anc="start"))
    save("fig-multi-plc.svg","".join(s))

if __name__=="__main__":
    fig_schedule(); fig_band(); fig_fusion(); fig_plc()
    # petits diagnostics chiffrés
    om_fu,_,sg_fu=run_filter(True,True)
    om_fa,_,_=run_filter(True,False)
    dt,N,tth,tom,fa,sl,sf,ss,drop=scenario()
    def rmse(a): return (sum((a[k]-tom[k])**2 for k in range(120,N))/(N-120))**0.5
    print("RMSE vitesse  fusion=%.2f  rapide_seul=%.2f  (°/s)"%(rmse(om_fu),rmse(om_fa)))
    print("√P nominal≈%.3f  pic dropout≈%.3f (°)"%(sg_fu[240],max(sg_fu[250:330])))
    print("Figures multi-fréquences générées.")
