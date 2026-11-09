# Ionic Thruster Electrode Spacing: simulation code and results

Supplementary material for the paper *Effect of Electrode Spacing on Corona Current, Sparkover Margin, and Thrust of a Laboratory-Scale Ionic Thruster: A High-Voltage Simulation Study* (M. S. Bhuwana, R. S. Z. Aini, E. L. Talakua).

The code is a reduced-order electrohydrodynamic model: charge simulation method, Peek onset, field-line space-charge current (Deutsch approximation with the Kaptzov condition), ion-drag thrust, and Monte Carlo uncertainty. The results are model predictions and have not been validated experimentally for the modelled device.

## Files

| File | Content |
| --- | --- |
| `ionic_thruster_sim.py` | Main model, operating cases A, B, C, Monte Carlo, and sensitivity analysis |
| `validation2.py` | Verification checks: current–voltage fit, field-line convergence, coaxial comparison, scaling exponents |
| `results_main.csv` | Results per spacing (Tables 4 and 5 of the paper) |
| `sensitivity.csv` | Spearman rank correlations at d = 2.5 cm (Figure 6b) |
| `mapping_check.csv` | Field-line to coaxial current ratio (Figure 3b) |
| `iv_check.csv` | Current–voltage curve at d = 2 cm and its quadratic fit (Figure 3a) |
| `validation_summary.json` | Summary numbers of the verification checks |

## Usage

```bash
pip install -r requirements.txt
python ionic_thruster_sim.py   # writes results_main.csv and sensitivity.csv
python validation2.py          # writes mapping_check.csv, iv_check.csv, validation_summary.json
```

Run time is several minutes; the Monte Carlo step uses multiprocessing.

## Data

No experimental data were generated. The sparkover relation and the reference measurement (20 kV, 2 cm, 101 kPa) come from the dataset deposited with Usuelli et al., J. Phys. D 59, 015203 (2026), doi: 10.5281/zenodo.17085847.

## Citation and license

Please cite the paper above. Add a license file of your choice (for example MIT for code and CC BY 4.0 for data) before making the repository public.
