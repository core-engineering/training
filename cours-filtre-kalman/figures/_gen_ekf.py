#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures du chapitre EKF : suivi de la position relative d'une cible par
deux inclinomètres (rapides) + une caméra (lente, mesure non-linéaire).
Tout est simulé par un vrai EKF (pas de dessin à main levée).
État x = [px, py, vx, vy, θ, φ].
"""
import math
from _gen_applique import (svg_open, line, poly, polygon, txt, circ, rrect, arrow,
                           save, frame, legend, RNG, mm, mv, T, add,
                           C_A, C_MEAS, C_B, C_TRUTH, C_GRID, C_AXIS, C_DIM, C_ACC, C_WARN)
from _gen_figures import ellipse_pts
C_CAM = "#f0a44f"   # caméra (ambre)
C_INC = "#3987e5"   # inclinomètres (bleu)

# ---------- EKF (mesure scalaire, innovation fournie) ----------
def predict(x,P,F,Q):
    x2=mv(F,x); P2=add(mm(mm(F,P),T(F)),Q); return x2,P2
def update_scalar(x,P,Hrow,R,y):
    n=len(x)
    Ph=[sum(P[i][j]*Hrow[j] for j in range(n)) for i in range(n)]
    S=sum(Hrow[i]*Ph[i] for i in range(n))+R
    K=[Ph[i]/S for i in range(n)]
    x=[x[i]+K[i]*y for i in range(n)]
    IKH=[[(1.0 if i==j else 0.0)-K[i]*Hrow[j] for j in range(n)] for i in range(n)]
    P=add(mm(mm(IKH,P),T(IKH)),[[K[i]*R*K[j] for j in range(n)] for i in range(n)])
    P=[[(P[i][j]+P[j][i])*0.5 for j in range(n)] for i in range(n)]
    return x,P,S

Hd=3.0  # profondeur / portée nominale connue (m)

def simulate():
    dt=0.02; N=400; rng=RNG(2024)
    # vérité
    px,py=2.2,-1.6; vx,vy=-0.12,0.30
    def th_t(t): return 0.16*math.sin(2*math.pi*0.13*t)
    def ph_t(t): return 0.11*math.sin(2*math.pi*0.09*t+1.0)
    truth=[]; incl=[]; cam=[]
    for k in range(N):
        t=k*dt
        truth.append((px,py,th_t(t),ph_t(t)))
        # inclinomètres (chaque cycle)
        si=math.radians(0.6)
        incl.append((th_t(t)+si*rng.n(), ph_t(t)+si*rng.n()))
        # caméra (tous les 25 cycles)
        if k%25==0:
            sc=math.radians(0.4)
            zx=math.atan2(px,Hd)-th_t(t)+sc*rng.n()
            zy=math.atan2(py,Hd)-ph_t(t)+sc*rng.n()
            cam.append((k,zx,zy))
        px+=vx*dt; py+=vy*dt
    # EKF
    F=[[1,0,dt,0,0,0],[0,1,0,dt,0,0],[0,0,1,0,0,0],
       [0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]
    sa=0.25; sr=math.radians(6.0)   # bruit accel cible ; bruit taux attitude
    Qp=[[sa*sa*dt**4/4,sa*sa*dt**3/2],[sa*sa*dt**3/2,sa*sa*dt**2]]
    Q=[[0.0]*6 for _ in range(6)]
    Q[0][0]=Qp[0][0]; Q[0][2]=Qp[0][1]; Q[2][0]=Qp[0][1]; Q[2][2]=Qp[1][1]
    Q[1][1]=Qp[0][0]; Q[1][3]=Qp[0][1]; Q[3][1]=Qp[0][1]; Q[3][3]=Qp[1][1]
    Q[4][4]=sr*sr*dt*dt; Q[5][5]=sr*sr*dt*dt
    # init : position mal connue (grand P), attitude ok, vitesse inconnue
    x=[1.2,-0.8,0.0,0.0,incl[0][0],incl[0][1]]
    P=[[0.0]*6 for _ in range(6)]
    for d,val in ((0,4.0),(1,4.0),(2,1.0),(3,1.0),(4,0.02),(5,0.02)): P[d][d]=val
    Ri=math.radians(0.6)**2; Rc=math.radians(0.4)**2
    est=[]; ell=[]; perr=[]; aerr=[]; nis=[]; sigp=[]; siga=[]
    ci=0
    for k in range(N):
        x,P=predict(x,P,F,Q)
        # inclinomètres : deux updates scalaires linéaires
        zi=incl[k]
        Hth=[0,0,0,0,1,0]; Hph=[0,0,0,0,0,1]
        x,P,_=update_scalar(x,P,Hth,Ri, zi[0]-x[4])
        x,P,_=update_scalar(x,P,Hph,Ri, zi[1]-x[5])
        # caméra : EKF (deux composantes non-linéaires)
        if ci<len(cam) and cam[ci][0]==k:
            _,zx,zy=cam[ci]; ci+=1
            # composante x : h = atan2(px,Hd) - θ
            hx=math.atan2(x[0],Hd)-x[4]
            Hx=[Hd/(Hd*Hd+x[0]*x[0]),0,0,0,-1,0]
            x,P,Sx=update_scalar(x,P,Hx,Rc, zx-hx)
            yx=zx-hx
            hy=math.atan2(x[1],Hd)-x[5]
            Hy=[0,Hd/(Hd*Hd+x[1]*x[1]),0,0,0,-1]
            x,P,Sy=update_scalar(x,P,Hy,Rc, zy-hy)
            nis.append((k, yx*yx/Sx))
        tp=truth[k]
        est.append((x[0],x[1]))
        ell.append((x[0],x[1],P[0][0],P[0][1],P[1][1]))
        perr.append(math.hypot(x[0]-tp[0],x[1]-tp[1]))
        aerr.append(math.degrees(abs(x[4]-tp[2])))
        sigp.append(math.sqrt(max(P[0][0],0)))          # √P position px (m)
        siga.append(math.degrees(math.sqrt(max(P[4][4],0))))  # √P attitude θ (°)
    return dict(dt=dt,N=N,truth=truth,cam=cam,est=est,ell=ell,perr=perr,aerr=aerr,
                nis=nis,sigp=sigp,siga=siga)

# ===========================================================
# Figure A — la géométrie (schéma vue de côté)
# ===========================================================
def fig_geometry():
    W,H=720,380; s=[svg_open(W,H)]
    camx,camy=330,66
    # plateforme + caméra
    s.append(rrect(camx-72,camy-26,144,30,"#151a22",C_AXIS,1.3,6))
    s.append(txt(camx,camy-8,"plateforme + caméra",C_DIM,11.5,wt="600"))
    # inclinomètres
    s.append(circ(camx+40,camy-11,4,C_INC)); s.append(circ(camx+54,camy-11,4,C_INC))
    s.append(txt(camx+86,camy-7,"2 inclinos → θ, φ",C_INC,10.5,anc="start"))
    cy0=camy+6
    # verticale gravité
    s.append(line(camx,cy0,camx,350,C_DIM,1.2,dash="4 5"))
    s.append(txt(camx+6,344,"verticale (gravité)",C_DIM,10.5,anc="start"))
    # axe optique incliné de θ
    ang=math.radians(22)
    ox,oy=camx+250*math.sin(ang),cy0+250*math.cos(ang)
    s.append(line(camx,cy0,ox,oy,C_ACC,1.6,dash="2 3"))
    s.append(txt(ox+6,oy,"axe optique (tilt θ)",C_ACC,10.5,anc="start"))
    # arc theta (petit)
    s.append(f'<path d="M {camx:.0f} {cy0+34:.0f} A 34 34 0 0 0 {camx+34*math.sin(ang):.1f} {cy0+34*math.cos(ang):.1f}" fill="none" stroke="{C_ACC}" stroke-width="1.3"/>')
    s.append(txt(camx+11,cy0+30,"θ",C_ACC,12,mono=True,wt="700"))
    # cible
    tx,ty=camx+205,320
    s.append(circ(tx,ty,7,C_B)); s.append(txt(tx+14,ty-2,"cible",C_B,12,anc="start",wt="700"))
    s.append(txt(tx+14,ty+14,"(px, py) relative",C_B,10.5,anc="start"))
    # ligne de visée
    s.append(line(camx,cy0,tx,ty,C_TRUTH,1.6))
    bang=math.atan2(tx-camx,ty-cy0)
    s.append(f'<path d="M {camx:.0f} {cy0+56:.0f} A 56 56 0 0 0 {camx+56*math.sin(bang):.1f} {cy0+56*math.cos(bang):.1f}" fill="none" stroke="{C_TRUTH}" stroke-width="1.2"/>')
    s.append(txt(camx+66,cy0+92,"β = atan(px / H)",C_TRUTH,11,anc="start"))
    # profondeur H
    s.append(line(camx-100,cy0,camx-100,ty,C_DIM,1,dash="2 3"))
    s.append(txt(camx-106,(cy0+ty)/2,"H (portée connue)",C_DIM,10.5,anc="end"))
    # équation mesure
    s.append(rrect(28,300,270,64,"#10201d",C_CAM,1.2,8))
    s.append(txt(44,322,"Mesure caméra (non-linéaire) :",C_CAM,11.5,anc="start",wt="700"))
    s.append(txt(44,342,"zx = atan(px/H) − θ",C_CAM,11,anc="start",mono=True))
    s.append(txt(44,358,"zy = atan(py/H) − φ",C_CAM,11,anc="start",mono=True))
    save("fig-ekf-geometry.svg","".join(s))

# ===========================================================
# Figure B — linéarisation : atan et sa tangente
# ===========================================================
def fig_lin():
    W,H=720,330; L,Rr,T_,B=64,690,40,270
    s=[svg_open(W,H)]
    p0,p1=-4.0,4.0
    def X(p): return L+(p-p0)/(p1-p0)*(Rr-L)
    def Y(v): return B-(v+1.3)/(2.6)*(B-T_)   # v in [-1.3,1.3] rad
    # axes
    s.append(line(L,Y(0),Rr,Y(0),C_AXIS,1.2)); s.append(line(X(0),T_,X(0),B,C_AXIS,1.2))
    for gp in (-4,-2,2,4): s.append(txt(X(gp),Y(0)+16,str(gp),C_DIM,11,mono=True))
    s.append(txt(Rr-4,Y(0)-8,"px / H",C_DIM,11,anc="end"))
    s.append(txt(X(0)+8,T_+8,"atan(px/H)  (rad)",C_DIM,11,anc="start"))
    # courbe atan
    s.append(poly([(X(p0+i*(p1-p0)/300),Y(math.atan(p0+i*(p1-p0)/300))) for i in range(301)],C_B,2.8))
    # point de linéarisation p̂
    ph=1.0
    s.append(circ(X(ph),Y(math.atan(ph)),5,"#fff"))
    s.append(txt(X(ph)+8,Y(math.atan(ph))+18,"p̂ (estimation courante)","#fff",11,anc="start"))
    # tangente : atan(ph) + (1/(1+ph²))(p-ph)
    sl=1.0/(1+ph*ph)
    s.append(poly([(X(p),Y(math.atan(ph)+sl*(p-ph))) for p in (p0,p1)],C_INC,2.2,dash="6 4"))
    s.append(txt(X(-3.4),Y(math.atan(ph)+sl*(-3.4-ph))-10,"tangente = linéarisation EKF (jacobienne)",C_INC,11.5,anc="start"))
    # zone d'erreur loin du point
    pf=3.4
    s.append(line(X(pf),Y(math.atan(pf)),X(pf),Y(math.atan(ph)+sl*(pf-ph)),C_WARN,2))
    s.append(txt(X(pf)-6,Y((math.atan(pf)+math.atan(ph)+sl*(pf-ph))/2),"erreur de",C_WARN,10.5,anc="end"))
    s.append(txt(X(pf)-6,Y((math.atan(pf)+math.atan(ph)+sl*(pf-ph))/2)+14,"linéarisation",C_WARN,10.5,anc="end"))
    save("fig-ekf-linearization.svg","".join(s))

# ===========================================================
# Figure C — trajectoire suivie (vue de dessus) + ellipses
# ===========================================================
def fig_track(sim):
    W,H=720,430; s=[svg_open(W,H)]
    truth=sim['truth']; est=sim['est']; ell=sim['ell']; cam=sim['cam']; N=sim['N']
    xs=[t[0] for t in truth]+[e[0] for e in est]
    ys=[t[1] for t in truth]+[e[1] for e in est]
    xmin,xmax=min(xs)-0.3,max(xs)+0.3; ymin,ymax=min(ys)-0.4,max(ys)+0.4
    L,Rr,T_,B=60,690,40,390
    sxp=(Rr-L)/(xmax-xmin); syp=(B-T_)/(ymax-ymin)
    def X(p): return L+(p-xmin)*sxp
    def Y(p): return B-(p-ymin)*syp
    # grille
    import math as _m
    for gx in range(int(_m.ceil(xmin)),int(xmax)+1):
        s.append(line(X(gx),T_,X(gx),B,C_GRID,1)); s.append(txt(X(gx),B+16,str(gx),C_DIM,10,mono=True))
    for gy in range(int(_m.ceil(ymin)),int(ymax)+1):
        s.append(line(L,Y(gy),Rr,Y(gy),C_GRID,1)); s.append(txt(L-8,Y(gy)+4,str(gy),C_DIM,10,anc="end",mono=True))
    s.append(txt((L+Rr)/2,B+34,"px  (m, repère stabilisé)",C_DIM,12))
    s.append(txt(L-40,T_-12,"py (m)",C_DIM,12,anc="start"))
    # caméra origine
    s.append(circ(X(0) if xmin<0<xmax else L, Y(0) if ymin<0<ymax else B, 0.1,"none"))
    # ellipses 1σ à intervalles + aux fixes caméra
    camk=set(c[0] for c in cam)
    for k in range(N):
        if k%50==0 or k in camk:
            cx,cy,a,b,c=ell[k]
            col=C_CAM if k in camk else C_ACC
            wv=2.0 if k in camk else 1.1
            s.append(poly(ellipse_pts(X(cx),Y(cy),a,b,c,1.0,sxp,-syp), col, wv, op=0.85 if k in camk else 0.4))
    # trajectoires
    s.append(poly([(X(t[0]),Y(t[1])) for t in truth],C_TRUTH,2.4,dash="6 4"))
    s.append(poly([(X(e[0]),Y(e[1])) for e in est],C_B,2.6))
    # marqueurs fixes caméra
    for (k,_,_) in cam:
        s.append(circ(X(est[k][0]),Y(est[k][1]),3.2,C_CAM))
    # départ
    s.append(circ(X(est[0][0]),Y(est[0][1]),4,"#fff"))
    s.append(txt(X(est[0][0])-10,Y(est[0][1])+16,"départ (P grand)","#fff",10.5,anc="end"))
    # légende
    lx,ly=L+10,T_+12
    s.append(line(lx,ly,lx+22,ly,C_TRUTH,2.4,dash="6 4")); s.append(txt(lx+28,ly+4,"cible vraie",C_TRUTH,11,anc="start"))
    s.append(line(lx+130,ly,lx+152,ly,C_B,2.6)); s.append(txt(lx+158,ly+4,"EKF estimé",C_B,11,anc="start"))
    s.append(circ(lx+270,ly,3.2,C_CAM)); s.append(txt(lx+280,ly+4,"fix caméra (ellipse resserrée)",C_CAM,11,anc="start"))
    save("fig-ekf-track.svg","".join(s))

# ===========================================================
# Figure D — l'incertitude respire au rythme de la caméra
#   √P position : dent de scie (baisse seulement aux fix caméra)
#   √P attitude : plate et basse (inclinomètres à chaque cycle)
# ===========================================================
def fig_uncert(sim):
    W,H=720,340; L,Rr,T_,B=62,690,50,275
    s=[svg_open(W,H)]; dt=sim['dt']; N=sim['N']; sigp=sim['sigp']; siga=sim['siga']; cam=sim['cam']
    tmax=N*dt
    # axe : position √P en m [0..0.6] ; attitude en ° affichée à l'échelle 1°↔0.6m
    X,Y=frame(s,L,Rr,T_,B,tmax,0,0.6,"temps (s)","√P position  (m)",[0,2,4,6,8],[0,0.2,0.4,0.6])
    for (k,_,_) in cam:
        s.append(line(X(k*dt),B,X(k*dt),B-6,C_CAM,1.4))
    # position : dent de scie
    s.append(poly([(X(k*dt),Y(min(0.6,sigp[k]))) for k in range(N)],C_B,2.6))
    s.append(txt(X(3.6),Y(0.55),"position (px,py) : ne baisse qu'aux fix caméra →",C_B,11,anc="start"))
    # attitude : √P (°) à l'échelle 1° -> 0.6 m
    s.append(poly([(X(k*dt),Y(min(0.6,siga[k]*0.6))) for k in range(N)],C_INC,2.0,op=0.95))
    s.append(txt(X(0.2),Y(siga[200]*0.6)+16,"attitude θ : √P plate et basse (inclinos chaque cycle)",C_INC,11,anc="start"))
    # marque cadence caméra
    s.append(line(L+8,T_-22,L+30,T_-22,C_CAM,1.6)); s.append(txt(L+36,T_-18,"fix caméra (2 Hz)",C_CAM,11,anc="start"))
    save("fig-ekf-uncert.svg","".join(s))

if __name__=="__main__":
    sim=simulate()
    fig_geometry(); fig_lin(); fig_track(sim); fig_uncert(sim)
    N=sim['N']
    rp=(sum(e*e for e in sim['perr'][100:])/(N-100))**0.5
    nism=sum(v for _,v in sim['nis'])/max(1,len(sim['nis']))
    # erreur juste avant / après un fix caméra
    ks=[c[0] for c in sim['cam'] if c[0]>150][:1]
    print("RMSE position (régime) = %.3f m"%rp)
    print("NIS caméra moy = %.2f (cible ~1)"%nism)
    print("perr avant/après 1er fix>150 : %.3f -> %.3f"%(sim['perr'][ks[0]-1],sim['perr'][ks[0]+1]) if ks else "n/a")
    print("Figures EKF générées.")
