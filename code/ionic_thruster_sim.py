"""
Reduced-order electrohydrodynamic (EHD) simulation of a laboratory-scale ionic thruster
with variable emitter-to-collector spacing d.  Version 2 (field-line space-charge model).

Model chain (SI units):
 1. Charge Simulation Method (CSM): Laplace field of a thin emitter wire facing a 7-tube collector row.
 2. Positive-DC Peek/Whitehead onset field; Laplace field lines are traced from the emitter to the collector (Deutsch approximation).
 3. Along every field-line tube the unipolar space-charge equations are solved with the Kaptzov condition
    (field at the emitter surface = Peek onset field), giving the corona current I(V).
 4. Ion-drag thrust T = eta * I * d / mu  and power P = V * I.
 5. Monte Carlo propagation of parameter uncertainty (parallel over spacings).
An equivalent-coaxial mapping (version 1 of the model) is kept for comparison in verification_checks.py.
"""
import numpy as np, csv, sys, json, os
from pathlib import Path
from multiprocessing import Pool
from scipy.optimize import brentq
from scipy.stats import spearmanr

EPS0 = 8.854187817e-12
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
PARAMS = dict(
    re=1.0e-4,      # emitter wire radius (m)
    Rt=2.0e-3,      # collector tube radius (m)
    pitch=7.0e-3,   # tube pitch (m)
    ntube=7,
    Lw=0.048,       # emitter wire length (m)
    mu=1.71e-4,     # positive-ion mobility (m^2 V^-1 s^-1), measured mean in DC corona in air
    delta=0.985,    # relative air density
    ms=0.965,       # Peek surface factor (nominal), solid wire between smooth (1.0) and weathered (0.93)
    eta=0.8,        # momentum-transfer efficiency (nominal)
    Vb_fit=(0.2514, 0.756, 0.350),  # V_b[kV] = A * p[kPa]^b * d[mm]^c, fitted to 30 measured breakdown voltages of a wire-to-grid thruster (Usuelli et al. dataset)
    p_kPa=100.0,    # ambient pressure used to evaluate the fit (kPa)
    Vcap=30e3,      # supply limit (V), as in the cited test rigs
)
SD = dict(re=0.05e-4, mu=0.12e-4, delta=0.01, d_err=0.5e-3, bf=0.10)   # bf: scatter of sparkover voltage (fit residual ~3%, plus extrapolation and geometry)
ETA_RANGE = (0.6, 1.0); MS_RANGE = (0.93, 1.0)

# ------------------------------------------------------------------ CSM electrostatics
def csm_setup(d, re, Rt, pitch, nt, Me=8, Mc=16, Ne=16, Nc=32):
    cx = d + Rt; ys = (np.arange(nt) - (nt - 1) / 2) * pitch
    Q, Cp, Vt = [], [], []
    for i in range(Me):
        a = 2*np.pi*i/Me; Q.append((0.5*re*np.cos(a), 0.5*re*np.sin(a)))
    for i in range(Ne):
        a = 2*np.pi*(i+0.5)/Ne; Cp.append((re*np.cos(a), re*np.sin(a))); Vt.append(1.0)
    for y in ys:
        for i in range(Mc):
            a = 2*np.pi*i/Mc; Q.append((cx+0.6*Rt*np.cos(a), y+0.6*Rt*np.sin(a)))
        for i in range(Nc):
            a = 2*np.pi*(i+0.5)/Nc; Cp.append((cx+Rt*np.cos(a), y+Rt*np.sin(a))); Vt.append(0.0)
    return np.array(Q), np.array(Cp), np.array(Vt)

def csm_full(d, re, Rt, pitch, nt, **kw):
    Q, Cp, Vt = csm_setup(d, re, Rt, pitch, nt, **kw)
    r = np.linalg.norm(Cp[:, None, :] - Q[None, :, :], axis=2)
    A = np.hstack([-np.log(r), np.ones((len(Cp), 1))]); A = np.vstack([A, np.r_[np.ones(len(Q)), 0.0]])
    x, *_ = np.linalg.lstsq(A, np.r_[Vt, 0.0], rcond=None)
    res = np.max(np.abs(A[:-1] @ x - Vt))
    return x[:-1], x[-1], Q, res

def csm_solve(d, re, Rt, pitch, nt, **kw):
    lam, c, Q, res = csm_full(d, re, Rt, pitch, nt, **kw)
    th = np.linspace(0, 2*np.pi, 360, endpoint=False); P = np.c_[re*np.cos(th), re*np.sin(th)]
    dv = P[:, None, :] - Q[None, :, :]; r2 = (dv**2).sum(2)
    E = (lam[None, :, None]*dv/r2[:, :, None]).sum(1)
    return lam, Q, res, np.linalg.norm(E, axis=1).max()

def field_map(d, V, nx=420, ny=320, xlim=(-0.01, 0.075), ylim=(-0.035, 0.035)):
    p = PARAMS; lam, Q, res, beta = csm_solve(d, p['re'], p['Rt'], p['pitch'], p['ntube'])
    xs = np.linspace(*xlim, nx); ys = np.linspace(*ylim, ny); X, Y = np.meshgrid(xs, ys)
    dx = X[..., None]-Q[:, 0]; dy = Y[..., None]-Q[:, 1]; r2 = dx**2+dy**2+1e-14
    return xs, ys, -(lam*0.5*np.log(r2)).sum(-1)*V, (lam*dx/r2).sum(-1)*V, (lam*dy/r2).sum(-1)*V

def peek_Ec(re, delta, ms):
    # Peek/Whitehead onset field for POSITIVE DC polarity: E = 33.7 m delta (1 + 0.24/sqrt(delta r)) kV/cm, r in cm
    return ms*3.37e6*delta*(1+0.24/np.sqrt(delta*re*100.0))      # V/m

# ------------------------------------------------------------------ Laplace field lines
def trace_lines(d, re, Rt, pitch, nt, nlines=120):
    lam, c, Q, _ = csm_full(d, re, Rt, pitch, nt)
    cx = d+Rt; ys = (np.arange(nt)-(nt-1)/2)*pitch; cent = np.c_[np.full(nt, cx), ys]
    def field(P):
        dv = P-Q; r2 = (dv**2).sum(1); return (lam[:, None]*dv/r2[:, None]).sum(0)
    def psi(P):
        dv = P-Q; return c-(lam*0.5*np.log((dv**2).sum(1))).sum()
    lines = []
    for j in range(nlines):
        th = 2*np.pi*(j+0.5)/nlines; P = np.array([re*np.cos(th)*1.0005, re*np.sin(th)*1.0005])
        S_, G, PS = [0.0], [np.linalg.norm(field(P))], [psi(P)]; s = 0.0; hit = False
        for it in range(6000):
            r = np.linalg.norm(P); ds = min(2e-3, 0.03*max(r, re))
            def f(Pp): e = field(Pp); return e/np.linalg.norm(e)
            k1 = f(P); k2 = f(P+0.5*ds*k1); k3 = f(P+0.5*ds*k2); k4 = f(P+ds*k3)
            Pn = P+ds*(k1+2*k2+2*k3+k4)/6; s += ds
            if np.linalg.norm(cent-Pn, axis=1).min() < Rt:
                hit = True; S_.append(s); G.append(np.linalg.norm(field(Pn))); PS.append(0.0); break
            P = Pn; S_.append(s); G.append(np.linalg.norm(field(P))); PS.append(psi(P))
            if np.linalg.norm(P) > 0.4: break
        lines.append((np.array(S_), np.array(G), np.array(PS), hit))
    return lines, 2*np.pi/nlines

class FieldLineModel:
    """Deutsch-Kaptzov unipolar space-charge model along Laplace field-line tubes."""
    def __init__(self, d, re=1.0e-4, nlines=120, M=160):
        p = PARAMS; lines, dth = trace_lines(d, re, p['Rt'], p['pitch'], p['ntube'], nlines)
        L = len(lines); self.dth = dth; self.d = d
        self.s = np.zeros((L, M)); self.g = np.zeros((L, M)); self.dpsi = np.zeros(L); self.Phi = np.zeros(L)
        tau = np.linspace(0, 1, M)
        for l, (s, g, ps, hit) in enumerate(lines):
            idx = tau*(len(s)-1); k = np.arange(len(s))
            self.s[l] = np.interp(idx, k, s); self.g[l] = np.interp(idx, k, g)
            self.dpsi[l] = ps[0]-ps[-1]; self.Phi[l] = g[0]*re*dth
        self.A = self.Phi[:, None]/self.g; self.ds = np.diff(self.s, axis=1)
        self.C = np.concatenate([np.zeros((L, 1)), np.cumsum(0.5*(self.A[:, 1:]+self.A[:, :-1])*self.ds, axis=1)], axis=1)
        self.invA = 1.0/self.A
    def _vline(self, I, u0):
        u = np.sqrt(u0[:, None]**2+2*I[:, None]*self.C/EPS0); f = u*self.invA
        return (0.5*(f[:, 1:]+f[:, :-1])*self.ds).sum(1)
    def onset(self, Ec, re):
        u0 = np.full(len(self.Phi), Ec*re*self.dth); return np.min(self._vline(np.zeros(len(u0)), u0)/self.dpsi)
    def current(self, V, Ec, re, mu, Lw):
        L = len(self.Phi); u0 = np.full(L, Ec*re*self.dth); dV = V*self.dpsi
        act = self._vline(np.zeros(L), u0) < dV
        if not act.any(): return 0.0
        lo = np.zeros(L); hi = np.full(L, 1e-10)
        for _ in range(60):
            need = act & (self._vline(hi, u0) < dV)
            if not need.any(): break
            hi[need] *= 4
        for _ in range(45):
            mid = 0.5*(lo+hi); big = self._vline(mid, u0) > dV
            hi = np.where(big, mid, hi); lo = np.where(big, lo, mid)
        return mu*np.where(act, 0.5*(lo+hi), 0.0).sum()*Lw      # current scales linearly with mobility

def _build(d): return FieldLineModel(float(d), PARAMS['re'], 120)
BANK = {}
def build_bank(dmin=0.0085, dmax=0.0425, step=0.0005):
    grid = np.round(np.arange(dmin, dmax+1e-9, step), 5)
    with Pool() as pool: models = pool.map(_build, list(grid))
    BANK.update({g: m for g, m in zip(grid, models)}); return grid
def bank_get(d):
    g = np.array(list(BANK.keys())); return BANK[g[np.argmin(abs(g-d))]]

def Vb_nominal(d, p=PARAMS):
    A, b, c = p['Vb_fit']; return 1e3*A*p['p_kPa']**b*(d*1e3)**c
def Vop_A(d, p=PARAMS): return min(0.8*Vb_nominal(d, p), p['Vcap'])

def evaluate(model, V, re, mu, delta, ms, eta, Lw, d=None):
    Ec = peek_Ec(re, delta, ms); I = model.current(V, Ec, re, mu, Lw)
    return model.onset(Ec, re), I, eta*I*(d if d is not None else model.d)/mu

# ------------------------------------------------------------------ Monte Carlo (one spacing)
FIXEDV = (10e3, 20e3, 30e3); P_FIXED = 1.0
def mc_spacing(args):
    d, N, seed = args; rng = np.random.default_rng(seed); p = PARAMS
    acc = {k: [] for k in ['V0','IA','PA','TA','spA','VC','TC','spC','IC']+[f'{a}{k}' for k in range(3) for a in ('TB','IB','sp')]}; X = []
    VA = Vop_A(d)
    for _ in range(N):
        re = rng.normal(p['re'], SD['re']); mu = rng.normal(p['mu'], SD['mu']); delta = rng.normal(p['delta'], SD['delta'])
        ms = rng.uniform(*MS_RANGE); eta = rng.uniform(*ETA_RANGE); dact = d+rng.normal(0, SD['d_err']); bf = rng.normal(1, SD['bf'])
        mod = bank_get(dact); Vb = Vb_nominal(dact)*bf
        V0, I, T = evaluate(mod, VA, re, mu, delta, ms, eta, p['Lw'], dact)
        acc['V0'].append(V0); acc['IA'].append(I); acc['PA'].append(VA*I); acc['TA'].append(T); acc['spA'].append(VA >= Vb)
        for k, Vf in enumerate(FIXEDV):
            _, Ik, Tk = evaluate(mod, Vf, re, mu, delta, ms, eta, p['Lw'], dact)
            acc[f'TB{k}'].append(Tk); acc[f'IB{k}'].append(Ik); acc[f'sp{k}'].append(Vf >= Vb)
        Ec = peek_Ec(re, delta, ms)
        g = lambda V: V*mod.current(V, Ec, re, mu, p['Lw'])-P_FIXED
        lo, hi = V0*1.0001, p['Vcap']*1.5
        if g(hi) > 0:
            Vc = brentq(g, lo, hi, xtol=5.0, maxiter=40); Ic = mod.current(Vc, Ec, re, mu, p['Lw'])
            acc['VC'].append(Vc); acc['IC'].append(Ic); acc['TC'].append(eta*Ic*dact/mu); acc['spC'].append(Vc >= Vb)
        else:
            for k in ('VC','IC','TC'): acc[k].append(np.nan)
            acc['spC'].append(True)
        X.append([re, dact-d, mu, delta, ms, eta, bf, T])
    A = {k: np.array(v, float) for k, v in acc.items()}; okA = ~A['spA'].astype(bool)
    def ms_(x, m=None):
        x = x if m is None else x[m]; x = x[~np.isnan(x)]
        return (np.mean(x), np.std(x, ddof=1), *np.percentile(x, [2.5, 97.5])) if len(x) > 1 else (np.nan,)*4
    row = dict(d_cm=d*100, V0_kV=ms_(A['V0'])[0]/1e3, V0_sd=ms_(A['V0'])[1]/1e3, VA_kV=VA/1e3,
               IA_mA=ms_(A['IA'], okA)[0]*1e3, IA_sd=ms_(A['IA'], okA)[1]*1e3, PA_W=ms_(A['PA'], okA)[0], PA_sd=ms_(A['PA'], okA)[1],
               TA_mN=ms_(A['TA'], okA)[0]*1e3, TA_sd=ms_(A['TA'], okA)[1]*1e3, TA_lo=ms_(A['TA'], okA)[2]*1e3, TA_hi=ms_(A['TA'], okA)[3]*1e3,
               TA_gf=ms_(A['TA'], okA)[0]/9.80665*1e3, TPA_mNW=np.nanmean(A['TA'][okA]/np.maximum(A['PA'][okA], 1e-12))*1e3, spark_A=A['spA'].mean(),
               TC_mN=ms_(A['TC'])[0]*1e3, TC_sd=ms_(A['TC'])[1]*1e3, VC_kV=ms_(A['VC'])[0]/1e3, spark_C=np.nanmean(A['spC'].astype(float)))
    for k, Vf in enumerate(FIXEDV):
        okk = ~A[f'sp{k}'].astype(bool)
        row[f'TB{k}_mN'] = ms_(A[f'TB{k}'], okk)[0]*1e3 if okk.sum() > 5 else np.nan
        row[f'TB{k}_sd'] = ms_(A[f'TB{k}'], okk)[1]*1e3 if okk.sum() > 5 else np.nan
        row[f'IB{k}_mA'] = ms_(A[f'IB{k}'], okk)[0]*1e3 if okk.sum() > 5 else np.nan
        row[f'spark_B{k}'] = A[f'sp{k}'].mean()
    sens = None
    if abs(d-0.025) < 1e-9:
        Xa = np.array(X); names = ['wire radius','spacing error','ion mobility','air density','surface roughness','momentum efficiency','sparkover scatter']
        sens = {n: float(spearmanr(Xa[okA, i], Xa[okA, -1])[0]) for i, n in enumerate(names)}
    return row, sens

def run(N=300, seed=20261015):
    spacings = np.round(np.arange(1.0, 4.0001, 0.25), 2)*1e-2
    build_bank(); print('bank built', len(BANK), flush=True)
    with Pool() as pool: out = pool.map(mc_spacing, [(float(d), N, seed+i) for i, d in enumerate(spacings)])
    rows = [o[0] for o in out]; sens = [o[1] for o in out if o[1] is not None][0]
    return rows, sens

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rows, sens = run(N=N)
    RESULTS.mkdir(exist_ok=True)
    with open(RESULTS / 'results_main.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(RESULTS / 'sensitivity.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['parameter', 'spearman_rho']); [w.writerow([k, v]) for k, v in sens.items()]
    print('done', sens)
