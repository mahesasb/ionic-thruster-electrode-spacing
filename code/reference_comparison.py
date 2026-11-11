"""Comparison of the model with the reference measurement of a wire-to-grid thruster (paper Table 6, Figure 3c).

Model: nominal parameters, 20 kV, d = 2 cm. Measured values are those reported by Usuelli et al., J. Phys. D 59
(2026) 015203 for 20 kV, 2 cm, 101 kPa: thrust density 5.4 N/m2, thrust-to-power 6.39 N/kW, emitter spacing 17.5 mm.
Run: python code/reference_comparison.py
"""
import ionic_thruster_sim as S

p = S.PARAMS
V, d = 20e3, 0.02
mod = S.FieldLineModel(d, p['re'], 120)
Ec = S.peek_Ec(p['re'], p['delta'], p['ms'])
I = mod.current(V, Ec, p['re'], p['mu'], p['Lw'])
A = p['Lw'] * p['ntube'] * p['pitch']                       # frontal area of the modelled collector row (m2)

meas = dict(T_A=5.4, T_P=6.39, spacing=17.5e-3)             # N/m2, N/kW, m
meas['P_A'] = meas['T_A'] / (meas['T_P'] * 1e-3)            # W/m2
meas['j'] = meas['P_A'] / V * 1e3                           # mA/m2
meas['I_L'] = meas['j'] * meas['spacing']                   # mA/m

print(f'I = {I*1e6:.1f} uA, P = {V*I:.2f} W, frontal area = {A*1e3:.3f} 1e-3 m2')
rows = [
    ('Power per frontal area (W/m2)', V * I / A, meas['P_A']),
    ('Current density (mA/m2)', I / A * 1e3, meas['j']),
    ('Current per emitter length (mA/m)', I / p['Lw'] * 1e3, meas['I_L']),
]
for name, m_, x in rows:
    print(f'{name:<38} model {m_:8.2f}  measured {x:8.2f}  ratio {m_/x:.2f}')
for eta in (0.6, 1.0):
    TA = eta * I * d / p['mu'] / A
    TP = eta * d / (p['mu'] * V) * 1e3
    print(f'eta = {eta}: thrust density {TA:.2f} N/m2 (ratio {TA/meas["T_A"]:.2f}), thrust-to-power {TP:.2f} N/kW (ratio {TP/meas["T_P"]:.2f})')
for mu in (1.54e-4, 1.71e-4, 1.93e-4):
    print(f'ideal ion-drag bound T/P at mu = {mu:.2e}: {d / (mu * V) * 1e3:.2f} N/kW')
