[README.md](https://github.com/user-attachments/files/25310736/README.md)
# Orthogonal Differential Adherence

Reproducibility materials for:

> Rodriguez AM. Differential Adherence: Orthogonal Decomposition of Medication Nonadherence Reveals Clinically Significant Heterogeneity Destroyed by Sum-Score Instruments. *Submitted*, 2026.

## Data

Download `NeuroGerAd_Data_OSF.xlsx` from the NeuroGerAd study repository:

- **OSF:** [https://osf.io/kuaph](https://osf.io/kuaph)
- **Citation:** Prell T, Schönenberg A, Mendorf S, Mühlhammer HM, Grosskreutz J, Teschner U. Data on medication adherence in adults with neurological disorders: The NeuroGerAd study. *Sci Data*. 2022;9(1):734. doi:10.1038/s41597-022-01847-9

Place the file in the same directory as the analysis script.

## Usage

```bash
pip install -r requirements.txt
python differential_adherence_analysis.py
python sensitivity_analysis.py          # Table 2 sensitivity analysis
```

Or open `differential_adherence_analysis.ipynb` in Jupyter and run all cells.

## Outputs

| File | Description |
|------|-------------|
| `results_summary.txt` | Full statistical results: descriptives, ANOVA, pairwise contrasts |
| `fig1_3d_vs_sumscore.png` | Orthogonal decomposition vs sum-score reduction |
| `fig2_group_heterogeneity.png` | Group-level heterogeneity in differential adherence |
| `fig3_information_loss.png` | Orthogonality and between-group signal comparison |
| sensitivity_results.txt | Table 2: omnibus ANOVA and key pairwise contrast under 4 missingness thresholds |

## License

Analysis code is released under the MIT License. The NeuroGerAd dataset is governed by its own terms (CC BY 4.0; noncommercial scientific use). See the [OSF repository](https://osf.io/kuaph) for data license details.
