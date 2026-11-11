"""Convergence of the charge simulation method (paper Table 2): d = 2.5 cm, nominal parameters.

Prints, for four numbers of line charges and contour points, the maximum boundary-condition residual
(relative to the applied voltage), the maximum surface field per volt beta, and the onset voltage V0 = Ec / beta.
Run: python code/csm_convergence.py
"""
import ionic_thruster_sim as S

p = S.PARAMS
d = 0.025
Ec = S.peek_Ec(p['re'], p['delta'], p['ms'])
settings = [(4, 8, 8, 16), (6, 12, 12, 24), (8, 16, 16, 32), (12, 24, 24, 48)]   # Me, Mc, Ne, Nc
print('line charges (wire/tube)  contour points (wire/tube)  max residual  beta (1/m)  V0 (kV)')
for Me, Mc, Ne, Nc in settings:
    lam, Q, res, beta = S.csm_solve(d, p['re'], p['Rt'], p['pitch'], p['ntube'], Me=Me, Mc=Mc, Ne=Ne, Nc=Nc)
    print(f'{Me:>3} / {Mc:<3}                  {Ne:>3} / {Nc:<3}                   {res:.1e}      {beta:7.0f}   {Ec / beta / 1e3:6.2f}')
