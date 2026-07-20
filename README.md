# scoring-rules-2024

This repository contains code, data, and results from a simulation study on the properness of scoring rules in survival analysis, including supplementary experiments and visualizations.

See [HTML report](https://survival-org.github.io/scoring-rules-2024/).
Instructions for running each simulation are provided in the corresponding sections of the report.

## Contents of the HTML Report

- **Weibull distribution variability** [[link](https://survival-org.github.io/scoring-rules-2024/#weib)] (main manuscript; Figure 1)
- **Empirical violations of properness** [[link](https://survival-org.github.io/scoring-rules-2024/#tab1)] (main manuscript; Table 1)
- **SBS improperness: Bias and Tail Effects** [[link](https://survival-org.github.io/scoring-rules-2024/#sbs-bias)] (main manuscript; Figures 2,3)
- **SBS improperness: Administrative Censoring** [[link](https://survival-org.github.io/scoring-rules-2024/#sbs-admin)] (main manuscript; investigation related to the experiments section)
- **Empirical violations of properness using G(t)** [[link](https://survival-org.github.io/scoring-rules-2024/#tabc1)] (main manuscript; Table C1)
- **Sensitivity to Model Misspecification** [[link](https://survival-org.github.io/scoring-rules-2024/#sen)] (main manuscript; Figures 4,5,6,7,C1; Tables 3,C2)
- **Real-World Benchmark Evaluation** [[link](https://survival-org.github.io/scoring-rules-2024/#bench)]
  An exploratory comparison of survival metrics including 3 scoring rules (RCLL, RCLL\*, ISBS), the C-index and D-calibration, on a few real-world survival datasets.

## Reproducibility

All scripts used for the analyses were run with the R package versions listed in the [session info](https://survival-org.github.io/scoring-rules-2024/#r-session-info).
The most important packages, depending on the analysis script, are:

- `mlr3proba` $\to$ `0.8.9`
- `mlr3` $\to$ `1.7.1`
- `mlr3extralearners` $\to$ `1.5.1`
- `mlr3pipelines` $\to$ `0.10.0`
- `survdistr` $\to$ `0.0.3`

## Citation

> Sonabend, R., Zobolas, J., De Bin, R., Kopper, P., Burk, L., & Bender, A. (2024). Examining properness in the external validation of survival models with squared and logarithmic losses. https://arxiv.org/abs/2212.05260v3

