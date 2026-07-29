#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures du chapitre « estimer la variance d'un capteur (R) ».
- histogramme de mesures au repos + gaussienne ajustée (biais, σ)
- variance d'Allan (séparer bruit blanc / instabilité de biais / dérive)
- piège de la fenêtre : variance naïve qui gonfle vs estimateur lag-1 robuste
Tout est simulé (échantillons tirés d'un modèle bruit blanc + dérive).
"""
import math
from _gen_applique import (svg_open, line, poly, polygon, txt, circ, rrect, arrow,
                           save, frame, legend, C_A, C_MEAS, C_B, C_TRUTH,
                           C_GRID, C_AXIS, C_DIM, C_ACC, C_WARN)

# ---- gaussien EXACT (moyenne 0, variance 1) — indispensable pour retrouver σ² ----
class G:
    def __init__(self, seed=1): self.s = seed & ((1<<64)-1)
    def u01(self):
        self.s = (self.s*6364136223846793005 + 1) & ((1<<64)-1)
        return (self.s >> 11) / (1<<53)          # [0,1)
    def n(self):
        u1 = max(self.u01(), 1e-12); u2 = self.u01()
        return math.sqrt(-2*math.log(u1))*math.cos(2*math.pi*u2)

# ---- signal capteur au repos : biais + dérive (marche aléatoire) + bruit blanc ----
def signal(N, sw, srw, seed=5):
    g=G(seed); bias0=0.15; rw=0.0; out=[]
    for k in range(N):
        rw += srw*g.n()
        out.append(bias0 + rw + sw*g.n())
    return out

# ===========================================================
# Figure A — histogramme de mesures au repos + gaussienne
# ===========================================================
def fig_hist():
    W,H=720,360; L,Rr,T_,B=60,690,45,300
    s=[svg_open(W,H)]
    g=G(3); N=4000; sw=0.40; bias=0.15
    data=[bias+sw*g.n() for _ in range(N)]
    mu=sum(data)/N
    var=sum((d-mu)**2 for d in data)/(N-1); sig=math.sqrt(var)
    lo,hi=mu-4*sig,mu+4*sig; nb=40; bw=(hi-lo)/nb
    counts=[0]*nb
    for d in data:
        i=int((d-lo)/bw)
        if 0<=i<nb: counts[i]+=1
    cmax=max(counts)
    def X(v): return L+(v-lo)/(hi-lo)*(Rr-L)
    def Y(c): return B-c/(cmax*1.12)*(B-T_)
    # axe x
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    for gx in (-1,-0.5,0,0.5,1,1.5):
        if lo<=gx<=hi: s.append(txt(X(gx),B+18,f"{gx:g}",C_DIM,11,mono=True))
    s.append(txt((L+Rr)/2,H-8,"mesure (° — capteur immobile)",C_DIM,12.5))
    # barres
    for i in range(nb):
        x0=X(lo+i*bw); x1=X(lo+(i+1)*bw)
        s.append(rrect(x0+0.5,Y(counts[i]),max(1,x1-x0-1),B-Y(counts[i]),"#1e3a44",C_ACC,0.8,2,op=0.9))
    # gaussienne ajustée
    def g(v): return math.exp(-0.5*((v-mu)/sig)**2)/(sig*math.sqrt(2*math.pi))
    gmax=g(mu)
    pts=[(X(lo+j*(hi-lo)/240), Y(g(lo+j*(hi-lo)/240)/gmax*cmax)) for j in range(241)]
    s.append(poly(pts,C_B,2.8))
    # moyenne (biais) + ±σ
    s.append(line(X(mu),Y(cmax*1.05),X(mu),B,C_MEAS,1.6,dash="4 4"))
    s.append(txt(X(mu),T_-2,"moyenne = biais (≈ %.2f°) — À RETIRER, ce n'est pas du bruit"%mu,C_MEAS,11,wt="600"))
    for k in (1,2):
        for sgn in (1,-1):
            s.append(line(X(mu+sgn*k*sig),Y(g(mu+sgn*k*sig)/gmax*cmax),X(mu+sgn*k*sig),B,C_ACC,1,dash="2 3"))
    s.append(txt(X(mu+sig),Y(cmax*0.55),"±σ",C_ACC,12,anc="start",wt="700"))
    s.append(rrect(L+6,T_+6,232,40,"#10201d",C_B,1.2,7))
    s.append(txt(L+18,T_+24,"σ = √s²  = %.3f°"%sig,C_B,12,anc="start",mono=True))
    s.append(txt(L+18,T_+40,"R = s² = %.4f deg²"%var,C_B,12,anc="start",mono=True,wt="700"))
    save("fig-var-hist.svg","".join(s))

# ===========================================================
# Figure B — variance d'Allan (log-log)
# ===========================================================
def allan(y, tau0):
    N=len(y); out=[]
    m=1
    while m< N//8:
        K=N//m
        avg=[sum(y[i*m:(i+1)*m])/m for i in range(K)]
        sd=sum((avg[i+1]-avg[i])**2 for i in range(K-1))/(2*(K-1))
        out.append((m*tau0, math.sqrt(sd)))
        m=int(m*1.6)+1
    return out

def fig_allan():
    W,H=720,360; L,Rr,T_,B=70,690,45,285
    s=[svg_open(W,H)]
    tau0=0.01; N=40000
    y=signal(N, sw=0.02, srw=0.0018, seed=9)
    pts=allan(y,tau0)
    xs=[math.log10(t) for (t,_) in pts]; ys=[math.log10(a) for (_,a) in pts]
    xmin,xmax=min(xs)-0.1,max(xs)+0.1; ymin,ymax=min(ys)-0.15,max(ys)+0.15
    def X(lx): return L+(lx-xmin)/(xmax-xmin)*(Rr-L)
    def Y(ly): return B-(ly-ymin)/(ymax-ymin)*(B-T_)
    # grille décades
    gx=math.floor(xmin)
    while gx<=xmax:
        s.append(line(X(gx),T_,X(gx),B,C_GRID,1)); s.append(txt(X(gx),B+18,f"10{_sup(gx)}",C_DIM,11)); gx+=1
    gy=math.floor(ymin)
    while gy<=ymax:
        s.append(line(L,Y(gy),Rr,Y(gy),C_GRID,1)); s.append(txt(L-8,Y(gy)+4,f"10{_sup(gy)}",C_DIM,11,anc="end")); gy+=1
    s.append(line(L,B,Rr,B,C_AXIS,1.3)); s.append(line(L,T_,L,B,C_AXIS,1.3))
    s.append(txt((L+Rr)/2,H-8,"temps d'intégration  τ  (s)",C_DIM,12.5))
    s.append(txt(L-38,T_-12,"σ_Allan (°)",C_DIM,12,anc="start")
             )
    s.append(poly([(X(lx),Y(ly)) for lx,ly in zip(xs,ys)],C_B,2.8))
    for lx,ly in zip(xs,ys): s.append(circ(X(lx),Y(ly),2.2,C_B,op=0.8))
    # pentes indicatives
    s.append(txt(X(xmin+0.5),Y(ys[2]+0.35),"pente −½ : bruit BLANC",C_ACC,11.5,anc="start",wt="600"))
    s.append(txt(X(xmin+0.5),Y(ys[2]+0.18),"→ lire R ici (τ = τ₀)",C_ACC,11,anc="start"))
    imin=ys.index(min(ys))
    s.append(circ(X(xs[imin]),Y(ys[imin]),4,C_MEAS))
    s.append(txt(X(xs[imin]),Y(ys[imin])+22,"plancher : instabilité de biais",C_MEAS,11,wt="600"))
    s.append(txt(X(xmax-0.1),Y(ys[-1]-0.05),"pente +½ : dérive",C_WARN,11.5,anc="end",wt="600"))
    # point R
    s.append(circ(X(xs[0]),Y(ys[0]),4,C_B))
    save("fig-var-allan.svg","".join(s))

def _sup(n):
    m={'-':'⁻','0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'}
    return ''.join(m[c] for c in str(n))

# ===========================================================
# Figure C — piège de la fenêtre : variance naïve vs estimateur lag-1
# ===========================================================
def fig_window():
    W,H=720,340; L,Rr,T_,B=64,690,50,275
    s=[svg_open(W,H)]
    tau0=0.01; N=6000; sw=0.30; srw=0.011
    y=signal(N, sw=sw, srw=srw, seed=4)
    Rtrue=sw*sw
    # à mesure que la fenêtre grandit
    Ws=list(range(50,N,50)); naive=[]; lag1=[]
    for w in Ws:
        seg=y[:w]; mu=sum(seg)/w
        naive.append(sum((d-mu)**2 for d in seg)/(w-1))
        lag1.append(sum((seg[i+1]-seg[i])**2 for i in range(w-1))/(2*(w-1)))
    tmax=N*tau0
    def X(w): return L+(w*tau0)/tmax*(Rr-L)
    ymax=max(max(naive),Rtrue)*1.15
    def Y(v): return B-v/ymax*(B-T_)
    for gt in [0,10,20,30,40,50,60]:
        s.append(line(X(gt/tau0),T_,X(gt/tau0),B,C_GRID,1)); s.append(txt(X(gt/tau0),B+16,str(gt),C_DIM,10,mono=True))
    s.append(line(L,B,Rr,B,C_AXIS,1.3)); s.append(line(L,T_,L,B,C_AXIS,1.3))
    s.append(txt((L+Rr)/2,B+34,"durée de la fenêtre de mesure (s)",C_DIM,12.5))
    s.append(txt(L-40,T_-12,"variance estimée (deg²)",C_DIM,12,anc="start"))
    # R vrai (blanc)
    s.append(line(L,Y(Rtrue),Rr,Y(Rtrue),C_TRUTH,1.6,dash="6 4"))
    s.append(txt(Rr-6,Y(Rtrue)-6,"R vrai (bruit blanc)",C_TRUTH,11,anc="end"))
    s.append(poly([(X(w),Y(naive[i])) for i,w in enumerate(Ws)],C_WARN,2.6))
    s.append(poly([(X(w),Y(lag1[i])) for i,w in enumerate(Ws)],C_B,2.6))
    legend(s,L+6,T_-26,[(C_WARN,None,"variance naïve (gonfle avec la dérive)"),(C_B,None,"estimateur lag-1 (robuste)")])
    save("fig-var-window.svg","".join(s))

if __name__=="__main__":
    fig_hist(); fig_allan(); fig_window()
    # diagnostics
    y=signal(40000, sw=0.02, srw=0.0018, seed=9); a=allan(y,0.01)
    print("Allan : σ_A(τ0=%.3f)=%.4f  (σ_w visé=0.02)"%(a[0][0],a[0][1]))
    y2=signal(6000, sw=0.30, srw=0.006, seed=4)
    def naive(w):
        seg=y2[:w]; mu=sum(seg)/w; return sum((d-mu)**2 for d in seg)/(w-1)
    def lag1(w):
        seg=y2[:w]; return sum((seg[i+1]-seg[i])**2 for i in range(w-1))/(2*(w-1))
    print("Fenêtre 60 s : naïve=%.4f  lag1=%.4f  (R vrai=%.4f)"%(naive(5999),lag1(5999),0.09))
    print("Figures variance générées.")
