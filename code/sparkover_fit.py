"""Power-law fit of the breakdown voltage and cross-validation of its extrapolation.

Data: Fig. 4 of Usuelli et al., J. Phys. D 59 (2026) 015203 (dataset doi:10.5281/zenodo.17085847).
Columns of data/input/breakdown_data.csv: d_mm, pd (Pa m), Vbr (kV).
Model: Vb = a * p^b * d^c (Vb in kV, p in kPa, d in mm), fitted by least squares on ln(Vb).
"""
from pathlib import Path
import numpy as np

DATA = Path(__file__).resolve().parents[1] / 'data' / 'input'

d_mm, pd_, vbr = np.loadtxt(DATA / 'breakdown_data.csv', delimiter=',', skiprows=1, unpack=True)
p = pd_ / (d_mm / 1000.0) / 1000.0          # kPa
X = np.column_stack([np.ones_like(p), np.log(p), np.log(d_mm)])
y = np.log(vbr)

def fit(Xm, ym):
    return np.linalg.lstsq(Xm, ym, rcond=None)[0]

b = fit(X, y)
res = vbr / np.exp(X @ b) - 1
n, k = X.shape
s2 = np.sum((y - X @ b) ** 2) / (n - k)
se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
print(f'n = {n}; Vb = {np.exp(b[0]):.4f} p^{b[1]:.3f} d^{b[2]:.3f}')
print(f'rms residual {100*np.sqrt(np.mean(res**2)):.2f} %, max {100*np.max(np.abs(res)):.2f} %')
print(f'standard errors: ln a {se[0]:.3f}, p exponent {se[1]:.3f}, d exponent {se[2]:.3f}')

print('leave-one-gap-out:')
for g in np.unique(d_mm):
    tr = d_mm != g
    bb = fit(X[tr], y[tr])
    e = vbr[~tr] / np.exp(X[~tr] @ bb) - 1
    print(f'  d = {g:.0f} mm: rms {100*np.sqrt(np.mean(e**2)):.1f} %, max {100*np.max(np.abs(e)):.1f} %')

tr = p <= 70
bb = fit(X[tr], y[tr])
e = vbr[~tr] / np.exp(X[~tr] @ bb) - 1
print(f'train p <= 70 kPa, predict {np.sum(~tr)} points above: mean {100*e.mean():+.1f} %, max {100*np.max(np.abs(e)):.1f} %')

X2 = np.column_stack([np.ones_like(p), np.log(pd_)])
b2 = fit(X2, y)
r2 = vbr / np.exp(X2 @ b2) - 1
print(f'pd-only power law: rms {100*np.sqrt(np.mean(r2**2)):.1f} %, max {100*np.max(np.abs(r2)):.1f} %')
