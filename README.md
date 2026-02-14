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

# Figures that will generate:

<img width="4171" height="1617" alt="fig3_information_loss" src="https://github.com/user-attachments/assets/b1e8871f-2432-47f6-b38b-a4d8c5186919" />
<img width="5368" height="1612" alt="fig2_group_heterogeneity" src="https://github.com/user-attachments/assets/70775be4-8fb5-4cda-bba3-59a9f4cf72c2" />
<img width="4507" height="2068" alt="fig1_3d_vs_sumscore" src="https://github.com/user-attachments/assets/aa0c4b52-c5e8-46f3-bc28-910331f4759e" />


## License

Analysis code is released under the MIT License. The NeuroGerAd dataset is governed by its own terms (CC BY 4.0; noncommercial scientific use). See the [OSF repository](https://osf.io/kuaph) for data license details.
