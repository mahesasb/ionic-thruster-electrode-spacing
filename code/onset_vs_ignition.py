"""Compare the model onset voltage with measured ignition voltages and test the thrust sensitivity to the onset.

Data: Fig. 3 of Usuelli et al., J. Phys. D 59 (2026) 015203 (dataset doi:10.5281/zenodo.17085847);
columns of data/input/ignition_data.csv: gap_mm, p_kPa, Vign_kV.
Run: python code/onset_vs_ignition.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import ionic_thruster_sim as m

df = pd.read_csv(Path(__file__).resolve().parents[1] / 'data' / 'input' / 'ignition_data.csv')
P = m.PARAMS

def model_onset(re, ms=P['ms']):
    mods = {g: m.FieldLineModel(g / 1000.0, re, 120) for g in sorted(df.gap_mm.unique())}
    return np.array([mods[g].onset(m.peek_Ec(re, pk / 101.325, ms), re) / 1e3
                     for g, pk in zip(df.gap_mm, df.p_kPa)])

for label, re in (('measured wire radius 0.05 mm', 0.05e-3), ('simulated wire radius 0.10 mm', 0.10e-3)):
    v = model_onset(re)
    r100 = df.p_kPa == 100
    ratio = df.Vign_kV / v
    print(label)
    print('  measured/model at 100 kPa by gap:', [(int(g), round(x, 2)) for g, x in zip(df.gap_mm[r100], ratio[r100])])
    n_meas = np.polyfit(np.log(df.gap_mm[r100]), np.log(df.Vign_kV[r100]), 1)[0]
    n_mod = np.polyfit(np.log(df.gap_mm[r100]), np.log(v[r100]), 1)[0]
    print(f'  gap exponent at 100 kPa: measured {n_meas:.2f}, model {n_mod:.2f}')

print('case A thrust change when the onset field is scaled (nominal parameters):')
for d in (0.010, 0.025, 0.040):
    mod = m.FieldLineModel(d, P['re'], 120)
    V = min(0.8 * 1e3 * P['Vb_fit'][0] * P['p_kPa'] ** P['Vb_fit'][1] * (d * 1e3) ** P['Vb_fit'][2], P['Vcap'])
    T0 = m.evaluate(mod, V, P['re'], P['mu'], P['delta'], P['ms'], P['eta'], P['Lw'], d)[2]
    for k in (1.1, 1.2, 1.5):
        T = m.evaluate(mod, V, P['re'], P['mu'], P['delta'], P['ms'] * k, P['eta'], P['Lw'], d)[2]
        print(f'  d = {d*100:.1f} cm, onset x{k}: {100*(T/T0-1):+.1f} %')
