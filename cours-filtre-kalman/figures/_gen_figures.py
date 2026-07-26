#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Génère les figures SVG du cours « Filtre de Kalman — fondements théoriques ».
Toutes les courbes (gaussiennes, ellipses de covariance, trajectoire filtrée)
sont calculées exactement — jamais dessinées à main levée.

Palette (cohérente avec style.css) :
  a priori / prédiction   bleu   #3987e5
  mesure / vraisemblance  ambre  #f0a44f
  a posteriori / estimé   vert   #35c98b
  vérité terrain          gris   #c3ccd8
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- palette ---------------------------------------------------------------
C_PRIOR = "#3987e5"
C_MEAS  = "#f0a44f"
C_POST  = "#35c98b"
C_TRUTH = "#c3ccd8"
C_GRID  = "#242b34"
C_AXIS  = "#3a4653"
C_DIM   = "#9aa6b6"
C_ACC   = "#34d1bf"
FONT    = "-apple-system, Segoe UI, sans-serif"
MONO    = "JetBrains Mono, monospace"

# ---- petits utilitaires SVG ------------------------------------------------
def esc(s): return (s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;"))

def svg_open(w, h):
    return (f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" '
            f'font-family="{FONT}" role="img">')

def line(x1,y1,x2,y2,stroke=C_GRID,w=1,dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{w}"{d}/>')

def polyline(pts, stroke, w=2.2, fill="none", opacity=1.0, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    p = " ".join(f"{x:.2f},{y:.2f}" for x,y in pts)
    return (f'<polyline points="{p}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{w}" opacity="{opacity:.2f}" '
            f'stroke-linejoin="round" stroke-linecap="round"{d}/>')

def polygon(pts, fill, opacity=1.0, stroke="none", w=0):
    p = " ".join(f"{x:.2f},{y:.2f}" for x,y in pts)
    return (f'<polygon points="{p}" fill="{fill}" opacity="{opacity:.2f}" '
            f'stroke="{stroke}" stroke-width="{w}"/>')

def text(x,y,s,fill=C_DIM,size=13,anchor="middle",weight="400",mono=False,style=""):
    fam = f' font-family="{MONO}"' if mono else ""
    st  = f' font-style="{style}"' if style else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" fill="{fill}" font-size="{size}" '
            f'text-anchor="{anchor}" font-weight="{weight}"{fam}{st}>{esc(s)}</text>')

def circle(x,y,r,fill,stroke="none",w=0,opacity=1.0):
    return (f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{w}" opacity="{opacity:.2f}"/>')

def rect(x,y,w,h,fill="none",stroke=C_AXIS,sw=1,rx=6,opacity=1.0):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" '
            f'opacity="{opacity:.2f}"/>')

def arrow(x1,y1,x2,y2,stroke=C_DIM,w=1.8,head=6):
    ang = math.atan2(y2-y1, x2-x1)
    hx, hy = x2 - head*math.cos(ang), y2 - head*math.sin(ang)
    a1 = ang + math.radians(150); a2 = ang - math.radians(150)
    p = (f'{x2:.1f},{y2:.1f} '
         f'{x2+head*math.cos(a1):.1f},{y2+head*math.sin(a1):.1f} '
         f'{x2+head*math.cos(a2):.1f},{y2+head*math.sin(a2):.1f}')
    return (line(x1,y1,hx,hy,stroke,w) +
            f'<polygon points="{p}" fill="{stroke}"/>')

def save(name, body):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(body + "</svg>\n")
    print("écrit", name)

def gauss(x, mu, sig):
    return math.exp(-0.5*((x-mu)/sig)**2) / (sig*math.sqrt(2*math.pi))

# ---------------------------------------------------------------------------
# Figure 1 — Fusion de deux gaussiennes (a priori × mesure → a posteriori)
# ---------------------------------------------------------------------------
def fig_fusion():
    W,H = 720, 340
    L,R,T,B = 55, 690, 30, 285
    s = [svg_open(W,H)]
    # domaine
    x0, x1 = -6.0, 10.0
    def X(x): return L + (x-x0)/(x1-x0)*(R-L)
    # deux gaussiennes : prédiction large, mesure plus sûre
    mp, sp = 1.0, 2.0     # a priori
    mm, sm = 5.0, 1.3     # mesure
    # posterior = produit (fusion optimale scalaire)
    vp, vm = sp*sp, sm*sm
    spost2 = 1.0/(1.0/vp + 1.0/vm)
    mpost  = spost2*(mp/vp + mm/vm)
    spost  = math.sqrt(spost2)
    # échelle verticale commune
    ymax = max(gauss(mpost,mpost,spost), gauss(mm,mm,sm), gauss(mp,mp,sp))*1.08
    def Y(v): return B - v/ymax*(B-T)
    # axe
    s.append(line(L,B,R,B,C_AXIS,1.4))
    for gx in range(-6,11,2):
        s.append(line(X(gx),B,X(gx),B+5,C_AXIS,1))
        s.append(text(X(gx),B+20,str(gx),C_DIM,12,mono=True))
    s.append(text((L+R)/2, H-6, "état  x", C_DIM, 13))
    # courbes
    def curve(mu,sig):
        return [(X(x0+i*(x1-x0)/240), Y(gauss(x0+i*(x1-x0)/240,mu,sig))) for i in range(241)]
    # remplissage posterior
    fillpts = curve(mpost,spost) + [(X(x1),B),(X(x0),B)]
    s.append(polygon(fillpts, C_POST, 0.16))
    s.append(polyline(curve(mp,sp), C_PRIOR, 2.4))
    s.append(polyline(curve(mm,sm), C_MEAS, 2.4))
    s.append(polyline(curve(mpost,spost), C_POST, 2.8))
    # moyennes
    for mu,c in ((mp,C_PRIOR),(mm,C_MEAS),(mpost,C_POST)):
        s.append(line(X(mu),Y(gauss(mu,mu,{mp:sp,mm:sm,mpost:spost}[mu])),X(mu),B,c,1,dash="3 3"))
    # étiquettes
    s.append(text(X(mp)-4, Y(gauss(mp,mp,sp))-8, "a priori  N(x⁻, P⁻)", C_PRIOR, 13, anchor="middle", weight="600"))
    s.append(text(X(mm)+70, Y(gauss(mm,mm,sm))-8, "mesure  N(z, R)", C_MEAS, 13, anchor="middle", weight="600"))
    s.append(text(X(mpost), T-2, "a posteriori — plus étroite que les deux", C_POST, 13, anchor="middle", weight="700"))
    save("fig-fusion.svg", "".join(s))

# ---------------------------------------------------------------------------
# Figure 2 — La boucle bayésienne récursive (schéma prédiction / correction)
# ---------------------------------------------------------------------------
def fig_bayes_loop():
    W,H = 720, 300
    s = [svg_open(W,H)]
    # deux boîtes
    bw,bh = 210, 92
    yb = 70
    xp, xc = 90, 420
    s.append(rect(xp,yb,bw,bh, "#122033", C_PRIOR, 1.6, 10))
    s.append(rect(xc,yb,bw,bh, "#12271f", C_POST, 1.6, 10))
    s.append(text(xp+bw/2, yb+26, "PRÉDICTION", C_PRIOR, 15, weight="700"))
    s.append(text(xp+bw/2, yb+50, "x⁻ = F·x + B·u", "#cfe3ff", 14, mono=True))
    s.append(text(xp+bw/2, yb+72, "P⁻ = F·P·Fᵀ + Q", "#cfe3ff", 14, mono=True))
    s.append(text(xc+bw/2, yb+26, "CORRECTION", C_POST, 15, weight="700"))
    s.append(text(xc+bw/2, yb+50, "K = P⁻Hᵀ (H P⁻Hᵀ+R)⁻¹", "#d6f5e6", 13.5, mono=True))
    s.append(text(xc+bw/2, yb+72, "x = x⁻ + K (z − H x⁻)", "#d6f5e6", 13.5, mono=True))
    # flèche prédiction -> correction (a priori)
    s.append(arrow(xp+bw, yb+bh/2, xc, yb+bh/2, C_DIM, 2, 8))
    s.append(text((xp+bw+xc)/2, yb+bh/2-10, "a priori", C_DIM, 12, style="italic"))
    s.append(text((xp+bw+xc)/2, yb+bh/2+22, "x⁻, P⁻", C_PRIOR, 12, mono=True))
    # boucle correction -> prédiction (a posteriori) par le bas
    ybot = yb+bh+55
    s.append(line(xc+bw/2, yb+bh, xc+bw/2, ybot, C_DIM, 2))
    s.append(line(xc+bw/2, ybot, xp+bw/2, ybot, C_DIM, 2))
    s.append(arrow(xp+bw/2, ybot, xp+bw/2, yb+bh, C_DIM, 2, 8))
    s.append(text((xp+bw/2+xc+bw/2)/2, ybot+20, "a posteriori  x, P  → devient l'état du cycle suivant", C_POST, 12.5, style="italic"))
    # entrées : mesure z arrive dans correction
    s.append(arrow(xc+bw/2, 22, xc+bw/2, yb, C_MEAS, 2, 8))
    s.append(text(xc+bw/2+70, 18, "mesure z (capteur)", C_MEAS, 12.5, weight="600"))
    # entrée commande u dans prédiction
    s.append(arrow(xp+bw/2, 22, xp+bw/2, yb, C_ACC, 2, 8))
    s.append(text(xp+bw/2-58, 18, "commande u", C_ACC, 12.5, weight="600"))
    save("fig-bayes-loop.svg", "".join(s))

# ---------------------------------------------------------------------------
#  outils ellipse de covariance (2×2 symétrique définie positive)
# ---------------------------------------------------------------------------
def eig2(a,b,c):
    """valeurs/vecteurs propres de [[a,b],[b,c]]."""
    tr = a+c; det = a*c - b*b
    disc = math.sqrt(max((tr/2)**2 - det, 0.0))
    l1, l2 = tr/2+disc, tr/2-disc
    th = 0.5*math.atan2(2*b, a-c)
    return l1, l2, th

def ellipse_pts(cx, cy, a, b, c, k, sx, sy, n=120):
    """points pixel de l'ellipse k-sigma de covariance [[a,b],[b,c]] centrée (cx,cy).
       sx, sy : facteurs d'échelle pixel/unité (sy négatif car y vers le bas)."""
    l1,l2,th = eig2(a,b,c)
    ax, ay = k*math.sqrt(l1), k*math.sqrt(l2)
    ct, stt = math.cos(th), math.sin(th)
    pts=[]
    for i in range(n+1):
        t = 2*math.pi*i/n
        ex, ey = ax*math.cos(t), ay*math.sin(t)
        wx = ct*ex - stt*ey
        wy = stt*ex + ct*ey
        pts.append((cx + wx*sx, cy + wy*sy))
    return pts

# ---------------------------------------------------------------------------
# Figure 3 — Gaussienne bivariée : moyenne, ellipses 1σ / 2σ, échantillons
# ---------------------------------------------------------------------------
def fig_gaussian_2d():
    W,H = 720, 380
    s = [svg_open(W,H)]
    cx, cy = 360, 195
    sx, sy = 26, -26            # pixels par unité
    # covariance corrélée
    a,b,c = 4.0, 2.3, 2.2
    # grille + axes
    for gx in range(-6,7,2):
        s.append(line(cx+gx*sx, 40, cx+gx*sx, 330, C_GRID, 1))
    for gy in range(-4,5,2):
        s.append(line(90, cy+gy*sy, 630, cy+gy*sy, C_GRID, 1))
    s.append(line(90, cy, 630, cy, C_AXIS, 1.3))
    s.append(line(cx, 40, cx, 330, C_AXIS, 1.3))
    s.append(text(636, cy+4, "x₁", C_DIM, 13, anchor="start", mono=True))
    s.append(text(cx+8, 48, "x₂", C_DIM, 13, anchor="start", mono=True))
    # échantillons pseudo-aléatoires (Box-Muller, LCG déterministe)
    seed = 987654321
    def rnd():
        nonlocal seed
        seed = (seed*6364136223846793005 + 1) & ((1<<64)-1)
        return ((seed>>11)/(1<<53))
    l1,l2,th = eig2(a,b,c)
    ct,stt = math.cos(th), math.sin(th)
    for _ in range(140):
        u1,u2 = max(rnd(),1e-9), rnd()
        g1 = math.sqrt(-2*math.log(u1))*math.cos(2*math.pi*u2)
        g2 = math.sqrt(-2*math.log(u1))*math.sin(2*math.pi*u2)
        ex, ey = math.sqrt(l1)*g1, math.sqrt(l2)*g2
        wx = ct*ex - stt*ey; wy = stt*ex + ct*ey
        s.append(circle(cx+wx*sx, cy+wy*sy, 2.1, C_POST, opacity=0.55))
    # ellipses
    s.append(polyline(ellipse_pts(cx,cy,a,b,c,2,sx,sy), C_ACC, 1.6, dash="5 4", opacity=0.8))
    s.append(polyline(ellipse_pts(cx,cy,a,b,c,1,sx,sy), C_ACC, 2.6))
    s.append(circle(cx,cy,4,"#fff"))
    s.append(text(cx-10, cy-10, "μ", "#fff", 15, anchor="end", mono=True, weight="700"))
    # axes propres
    for sign in (1,-1):
        s.append(line(cx, cy, cx+sign*math.sqrt(l1)*ct*sx, cy+sign*math.sqrt(l1)*stt*sy, C_MEAS, 1.4, dash="2 3"))
    s.append(text(cx+120, 58, "1σ", C_ACC, 13, weight="700"))
    s.append(text(cx+165, 44, "2σ", C_ACC, 12, weight="600", style="italic"))
    save("fig-gaussian-2d.svg", "".join(s))

# ---------------------------------------------------------------------------
# Figure 4 — Conditionnement : observer x₁ = z resserre la loi de x₂
#            (c'est exactement la mise à jour de Kalman)
# ---------------------------------------------------------------------------
def fig_conditioning():
    W,H = 720, 380
    s=[svg_open(W,H)]
    cx, cy = 300, 200
    sx, sy = 30, -30
    a,b,c = 4.0, 2.6, 3.0
    # axes
    s.append(line(70, cy, 545, cy, C_AXIS, 1.3))
    s.append(line(cx, 40, cx, 350, C_AXIS, 1.3))
    s.append(text(548, cy+4, "x₁ (mesuré)", C_DIM, 12.5, anchor="start"))
    s.append(text(cx+8, 48, "x₂ (caché)", C_DIM, 12.5, anchor="start"))
    # ellipse jointe
    s.append(polyline(ellipse_pts(cx,cy,a,b,c,2,sx,sy), C_PRIOR, 1.4, dash="5 4", opacity=0.7))
    s.append(polyline(ellipse_pts(cx,cy,a,b,c,1,sx,sy), C_PRIOR, 2.4))
    s.append(circle(cx,cy,3.5,C_PRIOR))
    s.append(text(cx-70, cy-58, "loi jointe  p(x₁, x₂)", C_PRIOR, 13, weight="600"))
    # mesure observée x1 = z
    zx = 2.4
    xz = cx + zx*sx
    s.append(line(xz, 45, xz, 350, C_MEAS, 2, dash="6 4"))
    s.append(text(xz+6, 60, "x₁ = z observé", C_MEAS, 13, anchor="start", weight="600"))
    # loi conditionnelle p(x2 | x1=z) : moyenne = c*... , var = c - b^2/a
    mu2 = b/a * zx                 # E[x2|x1=z] (moyennes nulles)
    var2 = c - b*b/a               # variance conditionnelle (réduite !)
    sig2 = math.sqrt(var2)
    # tracer cette gaussienne 1D le long de la verticale x1=z, couchée vers la droite
    scaleg = 120
    gpts=[]
    for i in range(121):
        x2 = -4 + 8*i/120
        val = gauss(x2, mu2, sig2)
        gpts.append((xz + val*scaleg, cy - x2*sy*0/1 + ( - x2)*0))  # placeholder
    # correct mapping: vertical axis is x2 -> y = cy - x2*sy? sy negative so cy + x2*sy
    gpts=[]
    for i in range(121):
        x2 = -4 + 8*i/120
        val = gauss(x2, mu2, sig2)
        y = cy + x2*sy
        gpts.append((xz + val*scaleg, y))
    s.append(polyline(gpts, C_POST, 2.6))
    # remplissage
    fill = gpts + [(xz, cy + (-4)*sy)]
    s.append(polygon([(xz, cy+(-4)*sy)]+gpts+[(xz, cy+4*sy)], C_POST, 0.14))
    # moyenne conditionnelle
    ym = cy + mu2*sy
    s.append(circle(xz, ym, 4, C_POST))
    s.append(line(xz, ym, xz+ (gauss(mu2,mu2,sig2))*scaleg, ym, C_POST, 1, dash="2 2"))
    s.append(text(xz+150, ym-8, "p(x₂ | x₁=z)", C_POST, 13, anchor="start", weight="700"))
    s.append(text(xz+150, ym+12, "moyenne recalée,", C_POST, 11.5, anchor="start", style="italic"))
    s.append(text(xz+150, ym+28, "variance réduite", C_POST, 11.5, anchor="start", style="italic"))
    save("fig-conditioning.svg", "".join(s))

# ---------------------------------------------------------------------------
# Figure 5 — Cycle de la covariance : la prédiction gonfle, la correction resserre
# ---------------------------------------------------------------------------
def fig_predict_update():
    W,H = 720, 330
    s=[svg_open(W,H)]
    sx, sy = 20, -20
    cys = 175
    centers = [150, 360, 570]
    titles = [("P au cycle k−1", C_POST), ("P⁻ prédit  (F P Fᵀ + Q)", C_PRIOR), ("P corrigé  (I−KH) P⁻", C_POST)]
    # trois états de covariance
    covs = [(2.2, 0.6, 1.6),      # posterior précédent
            (5.2, 2.8, 3.2),      # prédiction : plus grand + corrélé (cisaillé par F)
            (1.5, 0.4, 0.9)]      # corrigé : nettement plus petit
    for (cx, (a,b,c), (tt,tc)) in zip(centers, covs, titles):
        # petite grille repère
        for gx in (-2,0,2):
            s.append(line(cx+gx*sx, 60, cx+gx*sx, 290, C_GRID, 1))
        for gy in (-2,0,2):
            s.append(line(cx-60, cys+gy*sy, cx+60, cys+gy*sy, C_GRID, 1))
        s.append(polyline(ellipse_pts(cx,cys,a,b,c,1,sx,sy), tc, 2.6))
        s.append(polyline(ellipse_pts(cx,cys,a,b,c,2,sx,sy), tc, 1.4, dash="4 4", opacity=0.6))
        s.append(circle(cx,cys,3,"#fff"))
        s.append(text(cx, 45, tt, tc, 13, weight="700"))
    # flèches
    s.append(arrow(centers[0]+95, cys, centers[1]-95, cys, C_PRIOR, 2, 8))
    s.append(text((centers[0]+centers[1])/2, cys-95, "prédire", C_PRIOR, 12.5, weight="600"))
    s.append(text((centers[0]+centers[1])/2, cys-78, "↑ incertitude", C_PRIOR, 11.5, style="italic"))
    s.append(arrow(centers[1]+95, cys, centers[2]-95, cys, C_POST, 2, 8))
    s.append(text((centers[1]+centers[2])/2, cys-95, "corriger", C_POST, 12.5, weight="600"))
    s.append(text((centers[1]+centers[2])/2, cys-78, "↓ incertitude", C_POST, 11.5, style="italic"))
    save("fig-predict-update.svg", "".join(s))

# ---------------------------------------------------------------------------
# Figure 6 — Un filtre en action : suivi 1D position/vitesse
#            vérité, mesures bruitées, estimé filtré + bande ±1σ
# ---------------------------------------------------------------------------
def fig_kalman_run():
    W,H = 720, 360
    L,Rr,T,B = 58, 690, 30, 300
    s=[svg_open(W,H)]
    N = 42
    dt = 1.0
    # vérité : départ 3, vitesse 1.1 puis léger virage
    truth=[]
    p, v = 3.0, 1.1
    for k in range(N):
        truth.append(p)
        if k==20: v = 0.35    # la vraie vitesse change : on teste la poursuite
        p += v*dt
    # mesures bruitées (LCG déterministe)
    seed = 2024
    def noise(scale):
        nonlocal seed
        seed = (seed*6364136223846793005 + 1) & ((1<<64)-1)
        u = (seed>>33)/ (1<<31) - 1.0     # ~[-1,1]
        return u*scale
    Rvar = 3.0
    meas = [truth[k] + noise(2.6) for k in range(N)]
    # filtre de Kalman 2D (p,v), mesure = position
    # F, H
    F = [[1,dt],[0,1]]
    Q = [[0.02,0],[0,0.02]]
    x = [meas[0], 0.0]
    P = [[10.0,0],[0,10.0]]
    est=[]; band=[]
    def mul(A,B):
        return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    def matvec(A,vv):
        return [A[i][0]*vv[0]+A[i][1]*vv[1] for i in range(2)]
    def T2(A): return [[A[0][0],A[1][0]],[A[0][1],A[1][1]]]
    for k in range(N):
        # prédiction
        x = matvec(F,x)
        FP = mul(F,P); P = mul(FP, T2(F))
        P = [[P[0][0]+Q[0][0],P[0][1]],[P[1][0],P[1][1]+Q[1][1]]]
        # correction (H = [1 0])
        Spp = P[0][0] + Rvar
        Kk = [P[0][0]/Spp, P[1][0]/Spp]
        y = meas[k] - x[0]
        x = [x[0]+Kk[0]*y, x[1]+Kk[1]*y]
        # P = (I-KH)P
        P = [[(1-Kk[0])*P[0][0], (1-Kk[0])*P[0][1]],
             [P[1][0]-Kk[1]*P[0][0], P[1][1]-Kk[1]*P[0][1]]]
        est.append(x[0]); band.append(math.sqrt(max(P[0][0],0)))
    # échelle
    allv = truth+meas+est
    ymin, ymax = min(allv)-2, max(allv)+2
    def X(k): return L + k/(N-1)*(Rr-L)
    def Y(val): return B - (val-ymin)/(ymax-ymin)*(B-T)
    # grille
    for gk in range(0,N,7):
        s.append(line(X(gk),T,X(gk),B,C_GRID,1))
        s.append(text(X(gk),B+18,str(gk),C_DIM,11,mono=True))
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    s.append(text((L+Rr)/2, H-6, "cycle  k", C_DIM, 13))
    s.append(text(L-6, T+6, "position", C_DIM, 12, anchor="end"))
    # bande ±1σ
    up = [(X(k),Y(est[k]+band[k])) for k in range(N)]
    dn = [(X(k),Y(est[k]-band[k])) for k in range(N)]
    s.append(polygon(up+dn[::-1], C_POST, 0.16))
    # vérité
    s.append(polyline([(X(k),Y(truth[k])) for k in range(N)], C_TRUTH, 2.2, dash="6 4"))
    # mesures
    for k in range(N):
        s.append(circle(X(k),Y(meas[k]),2.6,C_MEAS,opacity=0.85))
    # estimé
    s.append(polyline([(X(k),Y(est[k])) for k in range(N)], C_POST, 2.8))
    # marqueur du changement de vitesse
    s.append(line(X(20),T,X(20),B,C_MEAS,1,dash="2 4"))
    s.append(text(X(20), T-6, "la vraie vitesse change", C_MEAS, 11.5))
    # légende
    lx, ly = L+8, T+16
    s.append(line(lx,ly,lx+22,ly,C_TRUTH,2.2,dash="6 4")); s.append(text(lx+28,ly+4,"vérité",C_TRUTH,12,anchor="start"))
    s.append(circle(lx+120,ly,2.6,C_MEAS)); s.append(text(lx+130,ly+4,"mesures",C_MEAS,12,anchor="start"))
    s.append(line(lx+215,ly,lx+237,ly,C_POST,2.8)); s.append(text(lx+243,ly+4,"estimé ±1σ",C_POST,12,anchor="start"))
    save("fig-kalman-run.svg", "".join(s))

# ---------------------------------------------------------------------------
# Figure 7 — Le gain de Kalman arbitre : K en fonction du bruit de mesure
# ---------------------------------------------------------------------------
def fig_gain():
    W,H = 720, 330
    L,Rr,T,B = 60, 685, 35, 265
    s=[svg_open(W,H)]
    # cas scalaire régime : K = P⁻ / (P⁻ + R), P⁻ fixé
    Pm = 1.0
    r0, r1 = 0.0, 6.0
    def X(r): return L + (r-r0)/(r1-r0)*(Rr-L)
    def Y(k): return B - k*(B-T)     # k in [0,1]
    # grille horizontale
    for kk in (0,0.25,0.5,0.75,1.0):
        s.append(line(L,Y(kk),Rr,Y(kk),C_GRID,1))
        s.append(text(L-8,Y(kk)+4,f"{kk:.2f}",C_DIM,11,anchor="end",mono=True))
    s.append(line(L,B,Rr,B,C_AXIS,1.3))
    s.append(line(L,T,L,B,C_AXIS,1.3))
    for gr in range(0,7):
        s.append(text(X(gr),B+18,str(gr),C_DIM,11,mono=True))
    s.append(text((L+Rr)/2,H-8,"bruit de mesure  R   (P⁻ fixé = 1)",C_DIM,13))
    s.append(text(L-30,T-12,"gain K",C_DIM,12,anchor="start"))
    # courbe
    pts=[(X(r0+i*(r1-r0)/300), Y(Pm/(Pm+ (r0+i*(r1-r0)/300)))) for i in range(301)]
    s.append(polyline(pts, C_ACC, 3))
    # régimes
    s.append(circle(X(0.06),Y(Pm/(Pm+0.06)),4,C_MEAS))
    s.append(text(X(0.06)+8,Y(Pm/(Pm+0.06))-8,"R→0 : K→1", C_MEAS,12.5,anchor="start",weight="600"))
    s.append(text(X(0.06)+8,Y(Pm/(Pm+0.06))+10,"on croit la mesure", C_MEAS,11,anchor="start",style="italic"))
    s.append(circle(X(5.4),Y(Pm/(Pm+5.4)),4,C_PRIOR))
    s.append(text(X(5.4),Y(Pm/(Pm+5.4))-14,"R→∞ : K→0", C_PRIOR,12.5,anchor="middle",weight="600"))
    s.append(text(X(5.4),Y(Pm/(Pm+5.4))-30,"on croit le modèle", C_PRIOR,11,anchor="middle",style="italic"))
    # point d'équilibre R=P
    s.append(line(X(1),B,X(1),Y(0.5),C_ACC,1,dash="3 3"))
    s.append(text(X(1),B+34,"R = P⁻ ⇒ K = ½",C_ACC,11.5))
    save("fig-gain.svg", "".join(s))

if __name__ == "__main__":
    fig_fusion()
    fig_bayes_loop()
    fig_gaussian_2d()
    fig_conditioning()
    fig_predict_update()
    fig_kalman_run()
    fig_gain()
    print("Toutes les figures ont été générées.")
