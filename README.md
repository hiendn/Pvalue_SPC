# P-value SPC manuscript code and supplementary tables

This repository contains the R code, numerical settings, figures, and
supplementary tables for **Statistical process control via p-values**.
The reported comparison uses five charts and Phase II sample sizes
`n = 3, 5, 8, 20`.

## Supplementary tables

- `Table_S1a_IC_calibration_and_validation.csv`: five-chart IC calibration
  and independent validation.
- `Table_S1b_full_OOC_performance.csv`: all 120 five-chart OOC
  configurations, including ECD, PTS, uncertainty, and survivor counts.
- `Table_S2_VSI_performance.csv`: VSI calibration and independent
  elapsed-time evaluation.
- `Table_S3_implementation_example.csv`: the nine-inspection numerical trace.
- `Table_S4_elementary_bound_results.csv`: reusable-baseline and stationary
  AR(1) theorem-bound studies.
- `Table_S5_localisation_results.csv`: localisation performance and
  stopped-run diagnostics.

The tables are in [`supplementary_tables/`](supplementary_tables/). The
committed source summaries used to construct them are in [`data/`](data/).

## Code

The core simulations use base R. From the repository root:

```sh
Rscript R/run_core_simulations.R full
Rscript R/postprocess_vsi_example.R full
Rscript R/make_supplementary_tables.R
```

Use `pilot` instead of `full` in the first two commands for a short
structural check. The elementary-bound and localisation generators are in
`R/auxiliary/`; their run instructions are given in the file headers.
Simulation reruns are written under `R/output/{mode}`; the cited
supplementary tables are rebuilt from the committed numerical snapshots in
`data/`.

Corresponding uniform-$p$-value figures can be generated with:

```sh
Rscript scripts/make_figures.R
```

They are written to `generated_figures/`, leaving the supplied manuscript
graphics unchanged.

The reported calculations were run with R 4.4.0. One of the 15,000 `n=8`
raw exact-KS IC validation runs was unsignalled at inspection 5000 and was
recorded as 5001; the IC table therefore reports recorded, capped run
lengths. None of the 120 retained OOC configurations was censored. The
localisation simulations use a Bonferroni global alarm followed by Holm
localisation; three weak Cauchy configurations contain capped run lengths,
as identified by their censoring-rate column.
