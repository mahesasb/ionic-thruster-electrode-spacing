# Ionic thruster electrode spacing: simulation code, input data, and results

Supplementary material for the paper *Monte Carlo-based uncertainty evaluation of electrode spacing effects on corona current, sparkover margin, and thrust of a laboratory-scale ionic thruster* (M. S. Bhuwana, R. S. Z. Aini, E. L. Talakua).

The code is a reduced-order electrohydrodynamic model of a wire emitter facing a seven-tube collector row: charge simulation method (CSM), Peek onset, field-line space-charge current (Deutsch approximation with the Kaptzov condition), ion-drag thrust, and Monte Carlo uncertainty propagation. The results are model predictions and have not been validated experimentally for the modelled device.

## Layout

```
code/            model and analysis scripts
data/input/      published breakdown and ignition voltages used by the scripts (see data/input/README.md)
results/         result tables written by the scripts
```

## Which file produces what

| Paper element | Script | Output |
| --- | --- | --- |
| Table 1 (parameters and ranges) | `code/ionic_thruster_sim.py` (`PARAMS`, `SD`, `ETA_RANGE`, `MS_RANGE`) | no output file |
| Table 2 (CSM convergence) | `code/csm_convergence.py` | printed |
| Section 2.7, sparkover fit and its extrapolation | `code/sparkover_fit.py` | printed |
| Section 3.1 checks, Figure 3a-b | `code/verification_checks.py` | `results/iv_check.csv`, `results/mapping_check.csv`, `results/validation_summary.json` |
| Tables 4 and 5, Figures 5 and 6a, Sections 3.2-3.5 | `code/ionic_thruster_sim.py` | `results/results_main.csv` |
| Figure 6b, Section 3.6 | `code/ionic_thruster_sim.py` | `results/sensitivity.csv` |
| Table 6, Figure 3c, Section 3.7 | `code/reference_comparison.py` | printed |
| Table 3 (onset rows), Limitations | `code/onset_vs_ignition.py` | printed |

Figures 1 and 2 are schematic drawings. Figure 4 uses `field_map()` in `ionic_thruster_sim.py`; the plotting scripts for Figures 3-6 are not included.

Column names of `results/results_main.csv`: `d_cm` spacing; case A (80 % of the sparkover voltage): `V0_kV` onset, `VA_kV` operating voltage, `IA_mA` current, `PA_W` power, `TA_mN` thrust (`TA_sd`, `TA_lo`, `TA_hi` spread), `TA_gf` thrust in gram-force, `TPA_mNW` thrust-to-power, `spark_A` sparkover probability; case B (fixed 10, 20, 30 kV, index 0, 1, 2): `TB*_mN`, `IB*_mA`, `spark_B*`; case C (fixed 1 W): `TC_mN`, `VC_kV`, `spark_C`.

## Usage

```bash
pip install -r requirements.txt
python code/csm_convergence.py
python code/sparkover_fit.py
python code/reference_comparison.py
python code/onset_vs_ignition.py
python code/verification_checks.py      # about 1 minute on one CPU
python code/ionic_thruster_sim.py 300   # Monte Carlo, 300 draws per spacing; slow on few CPUs
```

The scripts find their input and output folders from their own location, so they can be run from any directory. `ionic_thruster_sim.py` overwrites `results/results_main.csv` and `results/sensitivity.csv`; the Monte Carlo seeds are fixed. A run with fewer draws (for example `python code/ionic_thruster_sim.py 3`) is only a smoke test and does not reproduce the tables.

## Reproduced values (checked)

- `csm_convergence.py`: residual 1.5e-4 to 3.1e-8, onset 6.74 to 7.19 kV (Table 2).
- `sparkover_fit.py`: Vb = 0.2514 p^0.756 d^0.350, rms residual 3.23 %, maximum 9.35 %; leave-one-gap-out error at 20 mm 4.6 % rms.
- `reference_comparison.py`: 20 kV, 2 cm: 82.5 uA, 1.65 W; model/measured ratios 0.83, 0.83, 2.3, 0.46-0.76, 0.55-0.92 (Table 6).
- `onset_vs_ignition.py`: measured/model onset 1.06-1.26 at 100 kPa for the 0.1 mm radius; gap exponent 0.39 measured, 0.23 modelled.
- `verification_checks.py`: field-line to coaxial current ratio 1.25-3.15 (n = 21); R2 = 1.0000.

## License and citation

MIT license (see `LICENSE`). Cite the paper and, if you use the code, `CITATION.cff`.
