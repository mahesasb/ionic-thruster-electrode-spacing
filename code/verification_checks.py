"""Verification checks of the model (paper Section 3.1, Figure 3a-b).

Writes results/mapping_check.csv, results/iv_check.csv and results/validation_summary.json.
Run: python code/verification_checks.py
"""
import numpy as np, json
from pathlib import Path
from scipy.optimize import brentq
import ionic_thruster_sim as S
EPS0 = S.EPS0; p = S.PARAMS; mu = p['mu']
RESULTS = S.RESULTS; RESULTS.mkdir(exist_ok=True)

def coax_I(V, Ec, r0, R, mu, Lw, n=2500):
    V0 = Ec*r0*np.log(R/r0)
    if V <= V0: return 0.0
    s = np.linspace(0, np.log(R/r0), n); r = r0*np.exp(s)
    f = lambda K: np.trapezoid(np.sqrt((Ec*r0/r)**2+K*(1-(r0/r)**2))*r, s)-V
    hi = (V/(R-r0))**2*4+1.0
    while f(hi) < 0: hi *= 4
    K = brentq(f, 0.0, hi, xtol=1e-12*hi, rtol=1e-10)
    return 2*np.pi*EPS0*mu*K*Lw

if __name__ == '__main__':
    out = {}
    print("A) field-line model vs equivalent-coaxial mapping (nominal parameters, mu=%.2e)" % mu)
    rows = []
    for dcm in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        d = dcm/100; m = S.FieldLineModel(d, p['re'], 120)
        lam, Q, res, beta = S.csm_solve(d, p['re'], p['Rt'], p['pitch'], p['ntube']); Ec = S.peek_Ec(p['re'], p['delta'], p['ms']); V0 = Ec/beta
        R = p['re']*np.exp(V0/(Ec*p['re']))
        for fac in (1.3, 1.8, 2.5):
            V = fac*V0
            if V > 45e3: continue
            Ic = coax_I(V, Ec, p['re'], R, mu, p['Lw']); Id = m.current(V, Ec, p['re'], mu, p['Lw'])
            rows.append((dcm, V/1e3, Ic*1e6, Id*1e6, Id/Ic)); print("%.1f  %5.1f  coax %8.2f  field-line %8.2f  ratio %.2f" % rows[-1])
    r = np.array([x[4] for x in rows]); print("ratio field-line/coax: mean %.2f, min %.2f, max %.2f (n=%d)" % (r.mean(), r.min(), r.max(), len(r)))
    np.savetxt(RESULTS / 'mapping_check.csv', np.array(rows), delimiter=',', header='d_cm,V_kV,I_coax_uA,I_fieldline_uA,ratio', comments='')
    out['ratio'] = dict(mean=float(r.mean()), min=float(r.min()), max=float(r.max()), n=len(r))
    print("\nB) field-line density convergence (d=2.5 cm, V=20 kV)")
    conv = []
    for nl in (60, 120, 240):
        m = S.FieldLineModel(0.025, p['re'], nl); Ec = S.peek_Ec(p['re'], p['delta'], p['ms'])
        I = m.current(20e3, Ec, p['re'], mu, p['Lw'])*1e6; V0 = m.onset(Ec, p['re'])/1e3; conv.append((nl, I, V0)); print(nl, "lines: I = %.3f uA ; onset %.3f kV" % (I, V0))
    out['conv'] = conv
    print("\nC) scaling checks (field-line model, nominal)")
    m2 = S.FieldLineModel(0.02, p['re'], 120); Ec = S.peek_Ec(p['re'], p['delta'], p['ms']); V0 = m2.onset(Ec, p['re'])
    Vs = np.linspace(V0*1.05, 30e3, 25); Is = np.array([m2.current(V, Ec, p['re'], mu, p['Lw']) for V in Vs]); X = Vs*(Vs-V0)
    k = (X*Is).sum()/(X*X).sum(); R2 = 1-((Is-k*X)**2).sum()/((Is-Is.mean())**2).sum()
    print("I = k V (V - V0) at d=2 cm: R2 = %.4f (V0 = %.2f kV)" % (R2, V0/1e3)); out['R2'] = float(R2)
    np.savetxt(RESULTS / 'iv_check.csv', np.c_[Vs/1e3, Is*1e6, k*X*1e6], delimiter=',', header='V_kV,I_fieldline_uA,I_fit_uA', comments='')
    slopes = {}
    for Vf in (20e3, 30e3):
        ds = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])/100; Is2 = []; Ts = []
        for dv in ds:
            mm = S.FieldLineModel(dv, p['re'], 120); I = mm.current(Vf, Ec, p['re'], mu, p['Lw']); Is2.append(I); Ts.append(I*dv/mu)
        Is2 = np.array(Is2); Ts = np.array(Ts)
        sI = np.polyfit(np.log(ds), np.log(Is2), 1)[0]; sT = np.polyfit(np.log(ds), np.log(Ts), 1)[0]; slopes[int(Vf/1e3)] = (float(sI), float(sT))
        print("V=%d kV: dlnI/dlnd = %.2f ; dlnT/dlnd = %.2f" % (Vf/1e3, sI, sT))
    out['slopes'] = slopes
    print("\nD) reference condition 20 kV, d = 2 cm")
    I = m2.current(20e3, Ec, p['re'], mu, p['Lw']); out['ref'] = dict(I_uA=I*1e6, P_W=20e3*I)
    for eta in (0.6, 0.8, 1.0): print("eta=%.1f: I=%.1f uA P=%.2f W T=%.2f mN T/P=%.2f N/kW" % (eta, I*1e6, 20e3*I, eta*I*0.02/mu*1e3, eta*0.02/(mu*20e3)*1e3))
    for mu_ in (1.54e-4, 1.71e-4, 1.93e-4): print("ideal ion-drag bound T/P at mu=%.2e: %.2f N/kW" % (mu_, 0.02/(mu_*20e3)*1e3))
    print("\nE) Peek onset for 0.1 mm wire at d=2 cm: %.2f kV (nominal ms)" % (V0/1e3)); out['V0_2cm'] = V0/1e3
    for ms in (0.93, 1.0): print("   ms=%.2f -> %.2f kV" % (ms, m2.onset(S.peek_Ec(p['re'], p['delta'], ms), p['re'])/1e3))
    json.dump(out, open(RESULTS / 'validation_summary.json', 'w'))
