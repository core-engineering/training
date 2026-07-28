#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures du chapitre « choisir le repère d'état » : caméra eye-in-hand
(lente, asynchrone) + cinématique effecteur connue (rapide, chaque cycle),
cible lente. KF linéaire sur la pose cible en repère BASE.
Tout est simulé (vrai KF + vraies transformations de repère).
Modèle 2D : position (X,Y) + angle traité à part ; rotation R(θ_eff).
"""
import math
from _gen_applique import (svg_open, line, poly, polygon, txt, circ, rrect, arrow,
                           save, frame, legend, RNG, kf_predict, kf_update,
                           C_A, C_MEAS, C_B, C_TRUTH, C_GRID, C_AXIS, C_DIM, C_ACC, C_WARN)
C_BASE = C_B       # repère base (vert)
C_REL  = C_A       # repère effecteur / relatif (bleu)
C_CAM  = "#f0a44f" # caméra (ambre)

def R(th):  return ((math.cos(th),-math.sin(th)),(math.sin(th),math.cos(th)))
def Rt(th): return ((math.cos(th), math.sin(th)),(-math.sin(th),math.cos(th)))
def mul2(M,v): return (M[0][0]*v[0]+M[0][1]*v[1], M[1][0]*v[0]+M[1][1]*v[1])
def sub(a,b): return (a[0]-b[0],a[1]-b[1])

# ---- vérité (fonctions du temps) ----
def p_tgt(t):  return (1.60+0.05*math.sin(2*math.pi*0.08*t), 0.70+0.04*math.sin(2*math.pi*0.06*t+1))
def v_tgt(t):  return (0.05*2*math.pi*0.08*math.cos(2*math.pi*0.08*t), 0.04*2*math.pi*0.06*math.cos(2*math.pi*0.06*t+1))
def p_eff(t):  return (0.55*math.sin(2*math.pi*0.45*t), 0.45*math.sin(2*math.pi*0.33*t+0.5))
def v_eff(t):  return (0.55*2*math.pi*0.45*math.cos(2*math.pi*0.45*t), 0.45*2*math.pi*0.33*math.cos(2*math.pi*0.33*t+0.5))
def th_eff(t): return 0.7*math.sin(2*math.pi*0.40*t)
def om_eff(t): return 0.7*2*math.pi*0.40*math.cos(2*math.pi*0.40*t)

def p_rel_true(t):
    return mul2(Rt(th_eff(t)), sub(p_tgt(t),p_eff(t)))
def v_rel_true(t):
    pr=p_rel_true(t); w=om_eff(t)
    base=mul2(Rt(th_eff(t)), sub(v_tgt(t),v_eff(t)))
    return (base[0]-w*(-pr[1]), base[1]-w*(pr[0]))   # − ω×p_rel

def simulate():
    dt=0.01; N=600; rng=RNG(77)
    # caméra : arrivées asynchrones (~7 Hz) avec gigue, + latence de 4 cycles
    lat=4; cam_k=[]; k=8
    while k<N:
        cam_k.append(k)
        k+=13+int(round(2*rng.u()))     # ~0.13 s ± gigue
    cam_set=set(cam_k)
    sc=0.010    # bruit caméra (m) sur la position relative
    # KF base : état [X,Y,VX,VY]
    F=[[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]]
    sa=0.05     # accel cible très faible (cible lente) -> Q petit
    q=[[sa*sa*dt**4/4,0,sa*sa*dt**3/2,0],[0,sa*sa*dt**4/4,0,sa*sa*dt**3/2],
       [sa*sa*dt**3/2,0,sa*sa*dt**2,0],[0,sa*sa*dt**3/2,0,sa*sa*dt**2]]
    x=[1.0,0.3,0.0,0.0]; P=[[0.5,0,0,0],[0,0.5,0,0],[0,0,0.2,0],[0,0,0,0.2]]
    rec=dict(t=[],prel_true=[],prel_est=[],prel_cam=[],vrelx_true=[],vrelx_ok=[],vrelx_naive=[],
             basex_true=[],basex_est=[],sigb=[],cam_t=[],cam_prel=[])
    for k in range(N):
        t=k*dt
        x,P=kf_predict(x,P,F,q)
        # correction caméra (à l'arrivée), avec la pose effecteur à l'instant de capture
        if k in cam_set:
            tc=(k-lat)*dt
            zrel=(p_rel_true(tc)[0]+sc*rng.n(), p_rel_true(tc)[1]+sc*rng.n())
            zbase=(p_eff(tc)[0]+mul2(R(th_eff(tc)),zrel)[0], p_eff(tc)[1]+mul2(R(th_eff(tc)),zrel)[1])
            # deux updates scalaires (X puis Y), H=[1,0,0,0]/[0,1,0,0]
            x,P,_,_=kf_update(x,P,[1,0,0,0], sc*sc, zbase[0])
            x,P,_,_=kf_update(x,P,[0,1,0,0], sc*sc, zbase[1])
            rec['cam_t'].append(t); rec['cam_prel'].append(zrel[0])
        # SORTIE chaque cycle : recomposition en repère effecteur (courant)
        pe=p_eff(t); ve=v_eff(t); th=th_eff(t); w=om_eff(t)
        prel_est=mul2(Rt(th), (x[0]-pe[0], x[1]-pe[1]))
        vbase=mul2(Rt(th), (x[2]-ve[0], x[3]-ve[1]))
        vrel_ok    = (vbase[0]-w*(-prel_est[1]), vbase[1]-w*(prel_est[0]))  # correct
        vrel_naive = (vbase[0], vbase[1])                                   # sans −ω×p_rel
        rec['t'].append(t)
        rec['prel_true'].append(p_rel_true(t)[0]); rec['prel_est'].append(prel_est[0])
        rec['vrelx_true'].append(v_rel_true(t)[0]); rec['vrelx_ok'].append(vrel_ok[0]); rec['vrelx_naive'].append(vrel_naive[0])
        rec['basex_true'].append(p_tgt(t)[0]); rec['basex_est'].append(x[0]); rec['sigb'].append(math.sqrt(max(P[0][0],0)))
    rec['dt']=dt; rec['N']=N
    return rec

# ===========================================================
# Figure A — même cible, deux dynamiques
# ===========================================================
def fig_two(rec):
    W,H=720,350; L,Rr,T_,B=62,690,55,290
    s=[svg_open(W,H)]; N=rec['N']; dt=rec['dt']; tmax=N*dt
    X,Y=frame(s,L,Rr,T_,B,tmax,-1.4,2.2,"temps (s)","position X (m)",[0,1,2,3,4,5,6],[-1,0,1,2])
    s.append(poly([(X(k*dt),Y(rec['basex_true'][k])) for k in range(N)],C_BASE,2.8))
    s.append(poly([(X(k*dt),Y(rec['prel_true'][k])) for k in range(N)],C_REL,2.0))
    s.append(txt(X(3.2),Y(1.85),"cible en repère BASE — lente, Q minuscule",C_BASE,11.5,anc="start",wt="600"))
    s.append(txt(X(3.2),Y(-1.1),"cible en repère EFFECTEUR (ce que voit la caméra) — rapide",C_REL,11.5,anc="end",wt="600"))
    save("fig-frame-two.svg","".join(s))

# ===========================================================
# Figure B — structure (schéma)
# ===========================================================
def fig_structure():
    W,H=720,360; s=[svg_open(W,H)]
    s.append(rrect(26,54,668,214,"none",C_AXIS,1.4,12))
    s.append(txt(42,46,"OB cyclique rapide — 1×/cycle",C_DIM,12.5,anc="start",wt="600"))
    # KF base (centre)
    kx,ky,kw,kh=300,92,150,120
    s.append(rrect(kx,ky,kw,kh,"#10201d",C_BASE,1.6,10))
    s.append(txt(kx+kw/2,ky+24,"KF linéaire",C_BASE,13,wt="700"))
    s.append(txt(kx+kw/2,ky+42,"état = pose cible",C_BASE,11,wt="600"))
    s.append(txt(kx+kw/2,ky+58,"en repère BASE",C_BASE,11,wt="600"))
    s.append(txt(kx+kw/2,ky+80,"prédire chaque cycle",C_DIM,10,style="italic"))
    s.append(txt(kx+kw/2,ky+96,"corriger si image",C_DIM,10,style="italic"))
    # entrée caméra (haut) via déprojection
    s.append(rrect(300,20,190,30,"#1a1410",C_CAM,1.3,7))
    s.append(txt(395,39,"déprojection z_base",C_CAM,11,wt="700"))
    s.append(arrow(395,50,kx+kw/2,ky,C_CAM,1.8,7))
    s.append(arrow(150,26,300,30,C_CAM,1.8,7)); s.append(txt(60,22,"caméra (async)",C_CAM,11,anc="start",wt="600"))
    # entrée effecteur (gauche bas) : CD + Jacobien
    s.append(rrect(40,150,150,64,"#0f1622",C_REL,1.3,8))
    s.append(txt(115,172,"CD + Jacobien",C_REL,11.5,wt="700"))
    s.append(txt(115,190,"p_eff, v_eff, θ, ω",C_REL,10.5,mono=True))
    s.append(txt(115,206,"(connu, chaque cycle)",C_DIM,9.5,style="italic"))
    s.append(arrow(190,168,300,26,C_REL,1.5,7)); s.append(txt(210,150,"pose @ t_capture",C_REL,9.5,anc="start"))
    s.append(arrow(190,196,540,196,C_REL,1.5,7))
    # sortie : recomposition
    s.append(rrect(540,150,140,74,"#10201d",C_BASE,1.5,9))
    s.append(txt(610,170,"recomposition",C_BASE,11,wt="700"))
    s.append(txt(610,188,"p_rel, v_rel",C_BASE,11,mono=True))
    s.append(txt(610,206,"− ω×p_rel !",C_WARN,10.5,wt="700"))
    s.append(arrow(kx+kw,ky+kh/2,540,175,C_BASE,1.8,7))
    s.append(arrow(680,187,700,187,C_BASE,1.8,7)); s.append(txt(612,238,"→ régulation (lisse, chaque cycle)",C_DIM,10.5))
    # buffer
    s.append(rrect(300,300,300,40,"#151a22",C_WARN,1.2,8))
    s.append(txt(320,318,"Buffer horodaté des poses effecteur",C_WARN,11,anc="start",wt="700"))
    s.append(txt(320,334,"→ associer chaque image à p_eff(t_capture) malgré la latence.",C_DIM,10.5,anc="start"))
    save("fig-frame-structure.svg","".join(s))

# ===========================================================
# Figure C — sortie lisse chaque cycle malgré une caméra rare et gigueuse
# ===========================================================
def fig_output(rec):
    W,H=720,340; L,Rr,T_,B=62,690,52,275
    s=[svg_open(W,H)]; N=rec['N']; dt=rec['dt']; tmax=N*dt
    # on zoome sur 1.5..4.5 s pour bien voir
    t0,t1=1.5,4.5
    def X(t): return L+(t-t0)/(t1-t0)*(Rr-L)
    ys=[rec['prel_true'][k] for k in range(N) if t0<=k*dt<=t1]
    ymin,ymax=min(ys)-0.15,max(ys)+0.15
    def Y(v): return B-(v-ymin)/(ymax-ymin)*(B-T_)
    for gt in [2,3,4]:
        s.append(line(X(gt),T_,X(gt),B,C_GRID,1)); s.append(txt(X(gt),B+16,str(gt),C_DIM,11,mono=True))
    s.append(line(L,B,Rr,B,C_AXIS,1.2)); s.append(line(L,T_,L,B,C_AXIS,1.2))
    s.append(txt((L+Rr)/2,B+34,"temps (s)",C_DIM,12.5))
    s.append(txt(L-30,T_-10,"p_rel X (m) — ce que reçoit la régulation",C_DIM,11.5,anc="start")
             )
    # vérité relative
    idx=[k for k in range(N) if t0<=k*dt<=t1]
    s.append(poly([(X(k*dt),Y(rec['prel_true'][k])) for k in idx],C_TRUTH,2.4,dash="6 4"))
    # sortie filtrée (chaque cycle)
    s.append(poly([(X(k*dt),Y(rec['prel_est'][k])) for k in idx],C_BASE,2.8))
    # mesures caméra brutes (rares, gigueuses)
    for (tc,pr) in zip(rec['cam_t'],rec['cam_prel']):
        if t0<=tc<=t1: s.append(circ(X(tc),Y(pr),3.0,C_CAM))
    legend(s,L+6,T_-24,[(C_TRUTH,"6 4","vérité relative"),(C_CAM,None,"images caméra (rares)"),(C_BASE,None,"sortie filtrée (chaque cycle)")])
    save("fig-frame-output.svg","".join(s))

# ===========================================================
# Figure D — le terme −ω×p_rel dans la vitesse relative
# ===========================================================
def fig_velocity(rec):
    W,H=720,340; L,Rr,T_,B=62,690,55,275
    s=[svg_open(W,H)]; N=rec['N']; dt=rec['dt']; tmax=N*dt
    vs=rec['vrelx_true']+rec['vrelx_ok']+rec['vrelx_naive']
    lim=max(abs(min(vs)),abs(max(vs)))*1.05
    X,Y=frame(s,L,Rr,T_,B,tmax,-lim,lim,"temps (s)","vitesse relative Vx (m/s)",[0,1,2,3,4,5,6],[-2,-1,0,1,2])
    s.append(poly([(X(k*dt),Y(rec['vrelx_naive'][k])) for k in range(N)],C_WARN,2.0,op=0.9))
    s.append(poly([(X(k*dt),Y(rec['vrelx_true'][k])) for k in range(N)],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(k*dt),Y(rec['vrelx_ok'][k])) for k in range(N)],C_BASE,2.6))
    legend(s,L+6,T_-26,[(C_TRUTH,"6 4","vraie"),(C_BASE,None,"correcte (avec −ω×p_rel)"),(C_WARN,None,"naïve (sans le terme)")])
    save("fig-frame-velocity.svg","".join(s))

if __name__=="__main__":
    rec=simulate()
    fig_two(rec); fig_structure(); fig_output(rec); fig_velocity(rec)
    N=rec['N']
    rp=(sum((rec['prel_est'][k]-rec['prel_true'][k])**2 for k in range(100,N))/(N-100))**0.5
    rv_ok=(sum((rec['vrelx_ok'][k]-rec['vrelx_true'][k])**2 for k in range(100,N))/(N-100))**0.5
    rv_na=(sum((rec['vrelx_naive'][k]-rec['vrelx_true'][k])**2 for k in range(100,N))/(N-100))**0.5
    print("RMSE p_rel sortie = %.4f m"%rp)
    print("RMSE v_rel : correcte=%.3f  naïve=%.3f  (m/s)"%(rv_ok,rv_na))
    print("Figures repère générées.")
