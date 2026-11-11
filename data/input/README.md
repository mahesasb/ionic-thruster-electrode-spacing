# Input data

The two tables are values read from Fig. 3 and Fig. 4 of the open-access article and its deposited dataset:

D. Usuelli, A. Favalli, S. Trovato, R. Terenzi, M. Belan, "Atmospheric EAD grid thruster performance: a theoretical-experimental approach on geometry and pressure effects", J. Phys. D: Appl. Phys. 59 (2026) 015203, doi: 10.1088/1361-6463/ae251d. Dataset: doi: 10.5281/zenodo.17085847.

| File | Content | Columns |
| --- | --- | --- |
| `breakdown_data.csv` | Breakdown voltage of a wire-to-grid thruster, 30 points, gaps 20-50 mm, 20-90 kPa (Fig. 4) | `d_mm`, `pd` (Pa m), `Vbr` (kV) |
| `ignition_data.csv` | Ignition (corona onset) voltage, 45 points, gaps 20-60 mm, 20-100 kPa (Fig. 3) | `gap_mm`, `p_kPa`, `Vign_kV` |

Please cite the article above when you reuse these values.
