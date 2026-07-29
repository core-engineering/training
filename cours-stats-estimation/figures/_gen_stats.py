#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures du cours « Statistiques pour l'estimation — du hasard aux capteurs ».
Autonome. Tout est calculé exactement : échantillons tirés d'un générateur
gaussien exact (Box-Muller sur LCG 64 bits), courbes analytiques exactes.

Palette (cohérente avec style.css du cours) :
  accent / loi (ensemble)   violet #b98bf5
  réalisation observée      vert   #35c98b
  mesure / donnée           ambre  #f0a44f
  autres réalisations       gris   #6b7686
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
C_LOI   = "#b98bf5"   # la loi / le niveau ensemble
C_REA   = "#35c98b"   # LA réalisation observée
C_MEAS  = "#f0a44f"   # données / mesures
C_GRAY  = "#6b7686"   # autres réalisations (spaghetti)
C_TRUTH = "#c3ccd8"
C_GRID  = "#242b34"
C_AXIS  = "#3a4653"
C_DIM   = "#9aa6b6"
C_ACC   = "#34d1bf"   # sarcelle (calculs)
C_WARN  = "#f4c15b"
C_BLUE  = "#3987e5"
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
def circ(x,y,r,fill,op=1.0):
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" opacity="{op:.2f}"/>'
def rrect(x,y,w,h,fill="none",st=C_AXIS,sw=1,rx=6,op=1.0):
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

def gauss(x,mu,sig): return math.exp(-0.5*((x-mu)/sig)**2)/(sig*math.sqrt(2*math.pi))

# ---------- gaussien exact ----------
class G:
    def __init__(self, seed=1): self.s = seed & ((1<<64)-1)
    def u01(self):
        self.s = (self.s*6364136223846793005 + 1) & ((1<<64)-1)
        return (self.s >> 11) / (1<<53)
    def n(self):
        u1=max(self.u01(),1e-12); u2=self.u01()
        return math.sqrt(-2*math.log(u1))*math.cos(2*math.pi*u2)

# ===========================================================
# Fig 1 — les deux visages : réalisations (temps) vs loi (ensemble)
# ===========================================================
def fig_deux_niveaux():
    W,H=720,380; s=[svg_open(W,H)]
    # panneau gauche : 15 réalisations temporelles d'un capteur immobile
    L,Rr,T_,B=55,470,60,320
    Nt=120; mu=5.0; sig=0.5
    def X(k): return L+k/(Nt-1)*(Rr-L)
    def Y(v): return B-(v-3.2)/(6.8-3.2)*(B-T_)
    s.append(line(L,B,Rr,B,C_AXIS,1.2)); s.append(line(L,T_,L,B,C_AXIS,1.2))
    s.append(txt((L+Rr)/2,B+22,"temps (échantillons)",C_DIM,12))
    s.append(txt(L-8,Y(mu)+4,"μ",C_LOI,13,anc="end",mono=True,wt="700"))
    s.append(line(L,Y(mu),Rr,Y(mu),C_LOI,1.2,dash="5 4"))
    trajs=[]
    for r in range(15):
        g=G(100+r); tr=[mu+sig*g.n() for _ in range(Nt)]; trajs.append(tr)
        if r>0: s.append(poly([(X(k),Y(tr[k])) for k in range(Nt)],C_GRAY,1.0,op=0.35))
    s.append(poly([(X(k),Y(trajs[0][k])) for k in range(Nt)],C_REA,2.0))
    s.append(txt(L+8,T_-28,"NIVEAU RÉALISATION — des signaux dans le temps",C_REA,12.5,anc="start",wt="700"))
    s.append(txt(L+8,T_-12,"la verte : celle que VOTRE capteur produit ; les grises : celles qu'il aurait pu produire",C_DIM,10.5,anc="start",style="italic"))
    # coupe verticale à k*
    ks=86
    s.append(line(X(ks),T_,X(ks),B,C_MEAS,1.6,dash="3 3"))
    s.append(txt(X(ks),B+38,"coupe à un instant donné →",C_MEAS,11,wt="600"))
    # panneau droit : histogramme vertical des valeurs à k* + gaussienne
    L2,R2=500,690
    vals=[]
    for r in range(400):
        g=G(100+r)
        v=mu
        for _ in range(ks+1): v=mu+sig*g.n()
        vals.append(v)
    nb=22; lo,hi=3.2,6.8; bw=(hi-lo)/nb
    counts=[0]*nb
    for v in vals:
        i=int((v-lo)/bw)
        if 0<=i<nb: counts[i]+=1
    cmax=max(counts)
    def Y2(v): return B-(v-lo)/(hi-lo)*(B-T_)
    def X2(c): return L2+c/cmax*(R2-L2-40)
    s.append(line(L2,T_,L2,B,C_AXIS,1.2))
    for i in range(nb):
        y0=Y2(lo+(i+1)*bw); y1=Y2(lo+i*bw)
        s.append(rrect(L2+1,y0,max(1,X2(counts[i])-L2),max(1,y1-y0-1),"#251a35",C_LOI,0.8,2,op=0.9))
    # gaussienne verticale
    gmax=gauss(mu,mu,sig)
    s.append(poly([(X2(gauss(lo+j*(hi-lo)/160,mu,sig)/gmax*cmax),Y2(lo+j*(hi-lo)/160)) for j in range(161)],C_LOI,2.6))
    s.append(txt((L2+R2)/2,T_-28,"NIVEAU ENSEMBLE — la loi",C_LOI,12.5,wt="700"))
    s.append(txt((L2+R2)/2,T_-12,"N(μ, σ²) : où tombent les valeurs possibles",C_DIM,10.5,style="italic"))
    save("fig-deux-niveaux.svg","".join(s))

# ===========================================================
# Fig 2 — la densité gaussienne et ses aires (68/95/99.7)
# ===========================================================
def fig_pdf():
    W,H=720,340; L,Rr,T_,B=55,690,50,280
    s=[svg_open(W,H)]
    mu,sig=5.0,0.5
    lo,hi=mu-4*sig,mu+4*sig
    def X(v): return L+(v-lo)/(hi-lo)*(Rr-L)
    gmax=gauss(mu,mu,sig)
    def Y(d): return B-d/gmax*(B-T_)*0.92
    # aires
    for k,col,op in ((3,C_LOI,0.08),(2,C_LOI,0.14),(1,C_LOI,0.24)):
        pts=[(X(mu-k*sig+j*(2*k*sig)/120),Y(gauss(mu-k*sig+j*(2*k*sig)/120,mu,sig))) for j in range(121)]
        s.append(polygon(pts+[(X(mu+k*sig),B),(X(mu-k*sig),B)],col,op))
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    for v in (mu-3*sig,mu-2*sig,mu-sig,mu,mu+sig,mu+2*sig,mu+3*sig):
        s.append(line(X(v),B,X(v),B+5,C_AXIS,1))
    for lab,v in (("μ−3σ",mu-3*sig),("μ−σ",mu-sig),("μ",mu),("μ+σ",mu+sig),("μ+3σ",mu+3*sig)):
        s.append(txt(X(v),B+20,lab,C_DIM,11.5,mono=True))
    s.append(poly([(X(lo+j*(hi-lo)/240),Y(gauss(lo+j*(hi-lo)/240,mu,sig))) for j in range(241)],C_LOI,2.8))
    s.append(line(X(mu),Y(gmax),X(mu),B,C_LOI,1.2,dash="4 4"))
    s.append(txt(X(mu),Y(gmax*0.52)-2,"68 %",'#fff',14,wt="700"))
    s.append(txt(X(mu),Y(gmax*0.52)+16,"dans ±1σ",C_DIM,10.5))
    s.append(txt(X(mu-1.55*sig),Y(gmax*0.16),"95 % dans ±2σ",C_LOI,11,anc="end"))
    s.append(txt(X(mu+2.45*sig),Y(gmax*0.035),"99,7 % dans ±3σ",C_LOI,10.5,anc="start"))
    s.append(txt((L+Rr)/2,H-6,"mesure de l'inclinomètre (°) — capteur immobile",C_DIM,12))
    save("fig-pdf.svg","".join(s))

# ===========================================================
# Fig 3 — jointe / conditionnelle : nuage (erreur angle, erreur vitesse)
# ===========================================================
def fig_jointe():
    W,H=720,380; s=[svg_open(W,H)]
    cx,cy=330,195; sx,sy=105,-95
    # covariance corrélée (erreurs angle/vitesse d'un filtre)
    a,b,c=1.0,0.65,0.8
    l1=(a+c)/2+math.sqrt(((a-c)/2)**2+b*b); l2=(a+c)/2-math.sqrt(((a-c)/2)**2+b*b)
    th=0.5*math.atan2(2*b,a-c); ct,stt=math.cos(th),math.sin(th)
    s.append(line(70,cy,600,cy,C_AXIS,1.2)); s.append(line(cx,35,cx,355,C_AXIS,1.2))
    s.append(txt(604,cy+4,"erreur d'angle δθ",C_DIM,12,anc="start"))
    s.append(txt(cx+8,44,"erreur de vitesse δω",C_DIM,12,anc="start"))
    g=G(42)
    for _ in range(300):
        e1,e2=g.n(),g.n()
        wx=ct*math.sqrt(l1)*e1-stt*math.sqrt(l2)*e2
        wy=stt*math.sqrt(l1)*e1+ct*math.sqrt(l2)*e2
        s.append(circ(cx+wx*sx*0.55,cy+wy*sy*0.55,2.0,C_LOI,op=0.5))
    # tranche conditionnelle en δθ = z
    zx=0.9
    xs=cx+zx*sx*0.55
    s.append(line(xs,40,xs,355,C_MEAS,2,dash="6 4"))
    s.append(txt(xs+6,56,"on observe δθ = z",C_MEAS,12,anc="start",wt="600"))
    # loi conditionnelle : mu2 = (b/a) z ; var2 = c - b²/a
    mu2=b/a*zx; var2=c-b*b/a; s2=math.sqrt(var2)
    scale=110
    pts=[]
    for j in range(121):
        w=-2.6+5.2*j/120
        pts.append((xs+gauss(w,mu2,s2)*scale,cy+w*sy*0.55))
    s.append(polygon([(xs,cy+(-2.6)*sy*0.55)]+pts+[(xs,cy+2.6*sy*0.55)],C_REA,0.15))
    s.append(poly(pts,C_REA,2.6))
    ym=cy+mu2*sy*0.55
    s.append(circ(xs,ym,4,C_REA))
    s.append(txt(xs+130,ym-10,"loi de δω sachant δθ = z :",C_REA,12,anc="start",wt="700"))
    s.append(txt(xs+130,ym+8,"recentrée, resserrée",C_REA,11,anc="start",style="italic"))
    s.append(txt(90,340,"corrélation > 0 : les deux erreurs penchent du même côté →",C_DIM,11,anc="start"))
    s.append(txt(90,356,"observer l'une renseigne sur l'autre. C'est le moteur du filtre de Kalman.",C_DIM,11,anc="start"))
    save("fig-jointe.svg","".join(s))

# ===========================================================
# Fig 4 — théorème central limite : somme d'uniformes
# ===========================================================
def fig_tcl():
    W,H=720,330; s=[svg_open(W,H)]
    panels=[(1,"n = 1 (uniforme)"),(2,"n = 2"),(12,"n = 12 → gaussienne")]
    pw=200; gap=30; x0=45; T_,B=70,270
    for pi,(n,titre) in enumerate(panels):
        Lp=x0+pi*(pw+gap); Rp=Lp+pw
        g=G(7+pi)
        Ns=6000
        # somme de n uniformes centrées réduites : (u-0.5)*sqrt(12/n) sommées
        vals=[]
        for _ in range(Ns):
            ssum=0.0
            for _ in range(n): ssum+=(g.u01()-0.5)
            vals.append(ssum*math.sqrt(12.0/n))   # variance 1 quelle que soit n
        nb=30; lo,hi=-3.4,3.4; bw=(hi-lo)/nb
        counts=[0]*nb
        for v in vals:
            i=int((v-lo)/bw)
            if 0<=i<nb: counts[i]+=1
        cmax=max(counts)
        def X(v,Lp=Lp,Rp=Rp): return Lp+(v-lo)/(hi-lo)*(Rp-Lp)
        def Y(c): return B-c/cmax*(B-T_)*0.92
        s.append(line(Lp,B,Rp,B,C_AXIS,1.2))
        for i in range(nb):
            s.append(rrect(X(lo+i*bw)+0.5,Y(counts[i]),max(1,(Rp-Lp)/nb-1),B-Y(counts[i]),"#251a35",C_MEAS if n==1 else C_ACC if n==2 else C_LOI,0.7,1,op=0.85))
        # gaussienne N(0,1) superposée
        gm=gauss(0,0,1)
        s.append(poly([(X(lo+j*(hi-lo)/160),Y(gauss(lo+j*(hi-lo)/160,0,1)/gm*cmax)) for j in range(161)],C_TRUTH,1.8,dash="4 3"))
        s.append(txt((Lp+Rp)/2,T_-24,titre,C_LOI if n==12 else C_DIM,12.5,wt="700" if n==12 else "600"))
    s.append(txt(W/2,T_-46,"somme de n bruits indépendants (normalisée) — l'histogramme tend vers la cloche",C_DIM,12.5))
    s.append(txt(W/2,B+30,"tirets gris : gaussienne N(0,1) — le point de convergence universel",C_TRUTH,11.5))
    save("fig-tcl.svg","".join(s))

# ===========================================================
# Fig 5 — processus : ensemble vs moyenne temporelle (ergodicité)
# ===========================================================
def fig_processus():
    W,H=720,400; s=[svg_open(W,H)]
    L,Rr,T_,B=55,690,60,300
    Nt=300; mu=2.0; sig=0.45
    # AR(1) stationnaire : x_{k+1} = a x_k + e, a=0.9, var stationnaire sig²
    aco=0.90; se=sig*math.sqrt(1-aco*aco)
    def X(k): return L+k/(Nt-1)*(Rr-L)
    def Y(v): return B-(v-0.2)/(3.8-0.2)*(B-T_)
    s.append(line(L,B,Rr,B,C_AXIS,1.2)); s.append(line(L,T_,L,B,C_AXIS,1.2))
    s.append(txt((L+Rr)/2,B+24,"temps (échantillons)",C_DIM,12))
    trajs=[]
    for r in range(12):
        g=G(500+r); x=mu+sig*g.n(); tr=[]
        for _ in range(Nt):
            tr.append(x); x=mu+aco*(x-mu)+se*g.n()
        trajs.append(tr)
        if r>0: s.append(poly([(X(k),Y(tr[k])) for k in range(Nt)],C_GRAY,1.0,op=0.32))
    s.append(poly([(X(k),Y(trajs[0][k])) for k in range(Nt)],C_REA,1.9))
    # bande ensemble μ ± σ (stationnaire : constante)
    s.append(line(L,Y(mu),Rr,Y(mu),C_LOI,1.6,dash="6 4"))
    s.append(polygon([(L,Y(mu-sig)),(Rr,Y(mu-sig)),(Rr,Y(mu+sig)),(L,Y(mu+sig))],C_LOI,0.10))
    s.append(txt(Rr-4,T_-10,"bande μ ± σ (ENSEMBLE — sur toutes les trajectoires)",C_LOI,11.5,anc="end",wt="600"))
    # moyenne temporelle cumulée de LA trajectoire verte
    csum=0.0; tavg=[]
    for k in range(Nt):
        csum+=trajs[0][k]; tavg.append(csum/(k+1))
    s.append(poly([(X(k),Y(tavg[k])) for k in range(Nt)],C_MEAS,2.4))
    s.append(txt(L+6,B+40,"orange : moyenne TEMPORELLE cumulée de la trajectoire verte → converge vers μ (ERGODICITÉ) :",C_MEAS,11.5,anc="start",wt="600"))
    s.append(txt(L+6,B+56,"une seule (longue) trajectoire suffit pour estimer les statistiques d'ensemble.",C_MEAS,11,anc="start",style="italic"))
    save("fig-processus.svg","".join(s))

# ===========================================================
# Fig 6 — autocorrélation : bruit blanc vs bruit corrélé
# ===========================================================
def fig_autocorr():
    W,H=720,330; s=[svg_open(W,H)]
    Nt=4000; lagmax=25
    def autocorr(y,lag):
        mu=sum(y)/len(y); v=sum((d-mu)**2 for d in y)/len(y)
        return sum((y[i]-mu)*(y[i+lag]-mu) for i in range(len(y)-lag))/((len(y)-lag)*v)
    g=G(21); white=[g.n() for _ in range(Nt)]
    aco=0.85; se=math.sqrt(1-aco*aco); x=0.0; col=[]
    for _ in range(Nt): x=aco*x+se*g.n(); col.append(x)
    panels=[(white,"BRUIT BLANC — aucune mémoire",C_ACC,45),(col,"bruit CORRÉLÉ (filtré) — mémoire",C_WARN,395)]
    T_,B=70,270
    for (y,titre,colr,Lp) in panels:
        Rp=Lp+280
        def X(l,Lp=Lp,Rp=Rp): return Lp+l/lagmax*(Rp-Lp)
        def Y(v): return B-(v+0.25)/1.25*(B-T_)
        s.append(line(Lp,Y(0),Rp,Y(0),C_AXIS,1.2))
        s.append(line(Lp,T_,Lp,B,C_AXIS,1.2))
        for l in range(lagmax+1):
            r=autocorr(y,l)
            s.append(line(X(l),Y(0),X(l),Y(r),colr,4))
            s.append(circ(X(l),Y(r),2.6,colr))
        s.append(txt((Lp+Rp)/2,T_-24,titre,colr,12.5,wt="700"))
        s.append(txt((Lp+Rp)/2,B+24,"décalage temporel (lags)",C_DIM,11.5))
        s.append(txt(Lp-6,Y(1.0)+4,"1",C_DIM,10.5,anc="end",mono=True))
        s.append(txt(Lp-6,Y(0.0)+4,"0",C_DIM,10.5,anc="end",mono=True))
    s.append(txt(W/2,T_-46,"autocorrélation ρ(lag) : « la valeur d'un instant prédit-elle la suivante ? »",C_DIM,12.5))
    save("fig-autocorr.svg","".join(s))

# ===========================================================
# Fig 7 — marche aléatoire : l'enveloppe en √N
# ===========================================================
def fig_marche():
    W,H=720,360; s=[svg_open(W,H)]
    L,Rr,T_,B=55,690,45,300
    Nt=400; sig=0.30
    lim=sig*math.sqrt(Nt)*2.9
    def X(k): return L+k/Nt*(Rr-L)
    def Y(v): return B-(v+lim)/(2*lim)*(B-T_)
    s.append(line(L,Y(0),Rr,Y(0),C_AXIS,1.2)); s.append(line(L,T_,L,B,C_AXIS,1.2))
    s.append(txt((L+Rr)/2,B+24,"nombre de pas N",C_DIM,12))
    for r in range(28):
        g=G(900+r); x=0.0; tr=[0.0]
        for _ in range(Nt): x+=sig*g.n(); tr.append(x)
        s.append(poly([(X(k),Y(tr[k])) for k in range(0,Nt+1,2)],C_GRAY,0.9,op=0.35))
    # une réalisation en vert
    g=G(899); x=0.0; tr=[0.0]
    for _ in range(Nt): x+=sig*g.n(); tr.append(x)
    s.append(poly([(X(k),Y(tr[k])) for k in range(Nt+1)],C_REA,1.9))
    # enveloppe ±σ√N et ±2σ√N
    for m,dsh,opq in ((1,None,1.0),(2,"5 4",0.75)):
        s.append(poly([(X(k),Y( m*sig*math.sqrt(k))) for k in range(Nt+1)],C_LOI,2.2 if m==1 else 1.5,op=opq,dash=dsh))
        s.append(poly([(X(k),Y(-m*sig*math.sqrt(k))) for k in range(Nt+1)],C_LOI,2.2 if m==1 else 1.5,op=opq,dash=dsh))
    s.append(txt(X(330),Y(sig*math.sqrt(330))-10,"±σ·√N",C_LOI,13,wt="700",mono=True))
    s.append(txt(X(300),Y(2*sig*math.sqrt(300))-10,"±2σ·√N",C_LOI,11.5,mono=True))
    s.append(txt(L+10,T_+16,"chaque pas ajoute un petit bruit ⇒ les variances S'ADDITIONNENT : Var = N·σ²",C_DIM,12,anc="start"))
    s.append(txt(L+10,T_+33,"l'étalement croît en √N — c'est la dérive (et l'origine du Q du filtre de Kalman)",C_DIM,11.5,anc="start",style="italic"))
    save("fig-marche.svg","".join(s))

# ===========================================================
# Fig 8 — PSD : blanc plat vs filtré, aire = variance
# ===========================================================
def fig_psd():
    W,H=720,330; s=[svg_open(W,H)]
    L,Rr,T_,B=65,690,55,270
    fmax=100.0; fc=20.0; d2=1.0   # densité de puissance du blanc (unité²/Hz)
    def X(f): return L+f/fmax*(Rr-L)
    def Y(p): return B-p/1.15*(B-T_)
    s.append(line(L,B,Rr,B,C_AXIS,1.3)); s.append(line(L,T_,L,B,C_AXIS,1.3))
    for gf in (0,20,40,60,80,100):
        s.append(line(X(gf),B,X(gf),B+5,C_AXIS,1)); s.append(txt(X(gf),B+20,str(gf),C_DIM,11,mono=True))
    s.append(txt((L+Rr)/2,B+40,"fréquence (Hz)",C_DIM,12.5))
    s.append(txt(L-40,T_-10,"densité de puissance d²(f)  (unité²/Hz)",C_DIM,11.5,anc="start"))
    # blanc : plat (aire jusqu'à la bande capteur)
    s.append(poly([(X(0),Y(d2)),(X(fmax),Y(d2))],C_ACC,2.6))
    s.append(txt(X(78),Y(d2)-10,"bruit BLANC : plat — d² constant",C_ACC,12,wt="600"))
    # filtré 1er ordre fc : d²/(1+(f/fc)²), aire ombrée = variance
    pts=[(X(f),Y(d2/(1+(f/fc)**2))) for f in [fmax*j/300 for j in range(301)]]
    s.append(polygon(pts+[(X(fmax),B),(X(0),B)],C_LOI,0.20))
    s.append(poly(pts,C_LOI,2.6))
    s.append(txt(X(30),Y(0.30),"filtré (passe-bas f_c)",C_LOI,12,anc="start",wt="600"))
    s.append(txt(X(16),Y(0.13),"aire = VARIANCE σ²",'#fff',12.5,anc="start",wt="700"))
    s.append(line(X(fc),B,X(fc),Y(0.5),C_LOI,1.2,dash="3 3"))
    s.append(txt(X(fc),Y(0.55)-6,"f_c",C_LOI,11.5,mono=True))
    s.append(txt(X(48),Y(0.62),"σ² = ∫ d²(f) df ≈ d² × (π/2)·f_c",C_TRUTH,12.5,anc="start"))
    s.append(txt(X(48),Y(0.50),"→ σ = d·√BW  (la formule datasheet !)",C_TRUTH,12,anc="start"))
    save("fig-psd.svg","".join(s))

# ===========================================================
# Fig 9 — la précision d'un estimateur : loi de la moyenne selon N
# ===========================================================
def fig_estimateur():
    W,H=720,330; s=[svg_open(W,H)]
    L,Rr,T_,B=55,690,55,270
    mu=5.0; sig=0.5
    lo,hi=mu-1.8*sig,mu+1.8*sig
    def X(v): return L+(v-lo)/(hi-lo)*(Rr-L)
    gmax=gauss(mu,mu,sig/math.sqrt(100))
    def Y(d): return B-d/gmax*(B-T_)*0.92
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    for lab,v in (("μ−σ",mu-sig),("μ",mu),("μ+σ",mu+sig)):
        s.append(line(X(v),B,X(v),B+5,C_AXIS,1)); s.append(txt(X(v),B+20,lab,C_DIM,11.5,mono=True))
    cases=[(1,C_GRAY,"N = 1  →  écart-type σ"),(10,C_ACC,"N = 10  →  σ/√10"),(100,C_LOI,"N = 100  →  σ/10")]
    ly=T_+18
    for n,colr,lab in cases:
        sn=sig/math.sqrt(n)
        s.append(poly([(X(lo+j*(hi-lo)/300),Y(gauss(lo+j*(hi-lo)/300,mu,sn))) for j in range(301)],colr,2.6))
        s.append(line(L+14,ly,L+38,ly,colr,2.8)); s.append(txt(L+46,ly+4,lab,colr,12,anc="start",wt="600",mono=True))
        ly+=20
    s.append(line(X(mu),Y(gmax),X(mu),B,C_TRUTH,1,dash="4 4"))
    s.append(txt((L+Rr)/2,T_-24,"la MOYENNE de N mesures est elle-même une variable aléatoire — d'écart-type σ/√N",C_DIM,12.5))
    s.append(txt((L+Rr)/2,H-4,"moyenne empirique μ̂ (°)",C_DIM,12))
    save("fig-estimateur.svg","".join(s))

# ===========================================================
# Fig 10 — fusion de deux gaussiennes (préfiguration Kalman)
# ===========================================================
def fig_fusion():
    W,H=720,320; L,Rr,T_,B=55,690,40,265
    s=[svg_open(W,H)]
    x0,x1=-1.0,9.0
    def X(x): return L+(x-x0)/(x1-x0)*(Rr-L)
    ma,sa=3.0,1.4; mb,sb=5.5,0.9
    vc=1/(1/sa**2+1/sb**2); mc=vc*(ma/sa**2+mb/sb**2); sc=math.sqrt(vc)
    ymax=gauss(mc,mc,sc)*1.06
    def Y(v): return B-v/ymax*(B-T_)
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    for gx in range(0,10,2): s.append(txt(X(gx),B+18,str(gx),C_DIM,11,mono=True))
    def curve(m,sg): return [(X(x0+j*(x1-x0)/240),Y(gauss(x0+j*(x1-x0)/240,m,sg))) for j in range(241)]
    s.append(polygon(curve(mc,sc)+[(X(x1),B),(X(x0),B)],C_REA,0.15))
    s.append(poly(curve(ma,sa),C_BLUE,2.4)); s.append(poly(curve(mb,sb),C_MEAS,2.4))
    s.append(poly(curve(mc,sc),C_REA,2.8))
    s.append(txt(X(ma)-8,Y(gauss(ma,ma,sa))-10,"capteur A (incertain)",C_BLUE,12,anc="middle",wt="600"))
    s.append(txt(X(mb)+60,Y(gauss(mb,mb,sb))-10,"capteur B (plus sûr)",C_MEAS,12,wt="600"))
    s.append(txt(X(mc),T_+6,"fusion : penche vers le plus sûr, plus étroite que les deux",C_REA,12.5,wt="700"))
    s.append(txt((L+Rr)/2,H-4,"grandeur mesurée",C_DIM,12))
    save("fig-fusion-stats.svg","".join(s))

if __name__=="__main__":
    fig_deux_niveaux(); fig_pdf(); fig_jointe(); fig_tcl()
    fig_processus(); fig_autocorr(); fig_marche(); fig_psd()
    fig_estimateur(); fig_fusion()
    print("Toutes les figures stats générées.")
