# scoring-rules-2024

This repository contains the code, data, and results from a simulation study on the properness of scoring rules in survival analysis, alongside related experiments and visualizations.

Check the [supplementary HTML report](https://survival-org.github.io/scoring-rules-2024/), which summarizes the results and includes reproducibility instructions.
It also links code outputs to manuscript figures and tables for clarity.

See also [here the R script](https://github.com/survival-org/scoring-rules-2024/blob/main/ISBS_counterexample.R) that contains a numerical evaluation that complements the analytical proof of ISBS improperness.

## HTML Report Contents

The report includes:

- **High variability for generated survival and censoring distributions** [[link](https://survival-org.github.io/scoring-rules-2024/#high-variability-for-generated-survival-and-censoring-distributions)]  
  (corresponds to *Figure 1* in the main manuscript)
  
- **Appendix D experiment results:**  
  - *Table D1*: Empirical violations of properness [[link](https://survival-org.github.io/scoring-rules-2024/#table-d1)]
  - *Table D2*: Empirical violations of properness using G(t) [[link](https://survival-org.github.io/scoring-rules-2024/#table-d2)]
  
- **Degenerate Model Exploits ISBS Scoring Rule** [[link](https://survival-org.github.io/scoring-rules-2024/#degenerate-model-exploits-isbs-scoring-rule)]  
  A minimal model that outperforms established methods under ISBS, illustrating its vulnerability to gaming.

- **Real-World Benchmark Evaluation** [[link](https://survival-org.github.io/scoring-rules-2024/#real-world-benchmark-evaluation)]  
  An exploratory comparison of scoring rules (RCLL, RCLL*, ISBS, C-index, D-calibration) on validation tasks using real-world survival datasets.

## Citation

> Sonabend, R., Zobolas, J., Kopper, P., Burk LMU Munich, L., & Bender, A. (2022). Examining properness in the external validation of survival models with squared and logarithmic losses. https://arxiv.org/abs/2212.05260v3

