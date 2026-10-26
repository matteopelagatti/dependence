# Replication-materials audit — 26 October 2026

**Conclusion: `inst/paper` is not yet a complete replication package for the current paper and supplement.** Most central simulation code is present. Several supplementary studies and all three retained applications have source code elsewhere. Some numerical checks still lack an identified generator. Presence of a script is not evidence that it reproduces the published numbers.

This is a source and saved-output audit, not a full rerun. No analysis code, manuscript text, or research data was changed or copied during this audit. The accompanying `replication_source_inventory.csv` records inspected source candidates, their locations relative to the source roots below, and SHA-256 hashes. It does not designate every candidate as the final version.

Validation: all 42 inventoried R source files parsed successfully under R 4.5.2, and all three Python source files passed syntax parsing. This checks syntax only; no analysis scripts were executed. The inventory has 50 entries, including the two manuscript sources and three timing CSVs outside `inst/paper`.

## Source roots

- `repository`: the `dependence` Git repository.
- `All2`: the author's `autodep/R/All2` directory.
- `May2026`: the author's `autodep/R/Code/May2026` directory.
- `old_script`: the author's `autodep/R/Code/old_script` directory.
- `manuscript`: the synced Overleaf project, specifically `paper_jasa.tex` and `supplement_jasa.tex`.

The audit covered the requested repository and All2 directories, then followed relevant newer source files in May2026 and older figure generators in old_script. Cardiomyopathy is excluded because that application has been removed from the manuscript.

## Coverage of the current manuscript

Paths in the second column are relative to `repository/inst/paper`, unless a different source root is stated.

| Material | Located source | Status / remaining work |
|---|---|---|
| Continuous power study; main power figure and tables; supplementary full profiles | `02_power_continuous.R`, `02b_topup_statistics.R`, `results/power_continuous_figures.py` | Present, with aggregate CSVs and simulation caches. Establish a single canonical workflow: the resumable and older plotting scripts are alternatives, not interchangeable steps. Verify regeneration against the final aggregates. |
| Polynomial and spline order calibration | `01_null_calibration.R`, `01b_null_calibration_spline.R` | Present with calibration CSVs. |
| Nominal contingency study | `03_contingency.R` | Present with aggregate CSV and figure. |
| Order-scan calibration | `04_scan_calibration.R` | Present. Its output directory depends on the working directory/editor context; make this deterministic in the runbook. |
| Statistic timings and resampling timings | Scaling section of `02_power_continuous.R`; `06_timing_resampling.R`, `06b_timing_scan.R`, `06c_timing_hhg_fast.R` | Code present. Three timing CSVs are under `repository/results`, outside `inst/paper`. Record hardware, software and threading; timing values are machine-dependent. |
| S2.2 local-efficiency table (`tab:lp_eff`) | No matching final generator identified | Add the calculation of the noncentrality thresholds and the simulated Cramer–von Mises column, including truncation, seed and Monte Carlo settings. The manuscript specifies 100,000 replications and indices up to 30. |
| S3.2 numerical leverage claims and S3.3 orthogonality table (`tab:orth`) | No matching final experiment driver identified | Package basis implementations are not the experiment driver. Add scripts for the reported grids, normalization, error criterion and summaries. The older `quick_orthonormal_basis.R` is not a verified generator for these results. |
| S4.5 six null-density panels | `old_script/sim_for_pq_size.R`, density block at lines 76–102 | Missing from target. Candidate uses a historical `indeptest.R`, a hardcoded working directory and old return fields. Extract the relevant block and identify its exact dependency before adapting it. |
| S4.6 order sensitivity | `May2026/sim_sensitivity.R` | Missing from target. This is the ten-scheme source matching `sensitivity_Bn_plot.pdf`; the four-type variant is different. Isolate caches because variants reuse filenames. |
| S4.7 natural-label contingency comparison | `May2026/sim_contingency.R` | Missing from target. Candidate for the earlier, unadjusted experiment described in the text; it assumes `indeptest` is available without loading `dependence`. Confirm its outputs and API version. |
| S4.8 discrete-margin Genest simulation design | `May2026/dependencetest/exec/simulation_contingency_genest.R` | Missing from target. Generates synthetic data; this study does not require the restricted brain-injury dataset. |
| Main noisy-circle figure and S5 deformation figures | `old_script/example_of_deformation_plots.R` and `old_script/indeptest_with_cat.R` | Missing from target. Identified generator for all twelve `ex1`, `ex2`, `ex3` panels. Preserve the historical basis orders and RNG sequence when adapting; current defaults need not match. |
| S6 feature-screening simulations and tables | `May2026/sim_varsel_rho.R`, `summarize_varsel_rho.R`; also `sim_varsel.R`, `summarize_varsel.R` | Missing from target. The rho driver covers 0.2, 0.5 and 0.8. Check summary input paths: the single-rho driver and its summary currently use different directory conventions. Specify ties explicitly for discrete responses. |
| Main brain-injury application and S7.1 | `May2026/Brain application.R`; earlier `All2/BrainInjury_mod.R` | Missing from target. Include analysis/preprocessing code and access instructions, **not the restricted data**. Resolve the tie-handling mismatch below before claiming replication. |
| S7.2 yeast application | `May2026/Hughes application.R` | Missing from target. Loads `HHG::Yeast_hughes`; provide package version and loading instructions. No private data file is required. Verify tie convention and final output. |
| S7.3 housing application | `May2026/House prices topk screening.R`, `House prices topk screening logRMSE.R`, `House prices topk ranks.R`, `House prices topk replot.R` | Missing from target. The saved best-k results and fold-based rank summaries match the two manuscript tables at their printed precision. This verifies saved summaries, not a clean rerun. Data access/redistribution and cache provenance remain to be documented. |

`05_scan_K_calibration.R`, the pilot scripts, `basis_and_bic_experiment.py` and `grid_scan_experiment.py` include exploratory work. Label such work separately from the scripts needed to reproduce published outputs. Do not include deleted cardiomyopathy analyses or the unrelated applications in the large historical `All2/examples.R` file merely to increase coverage.

## Findings requiring resolution

1. **Brain-injury ties disagree with the manuscript.** `May2026/Brain application.R`, line 116, explicitly requests `ties = "first"`. The main paper, line 501, and the supplementary table note, line 1524, explicitly say mid-ranks. These are different statistics for tied observations. Locate a later mid-rank run, or rerun with `ties = "average"` and update all affected numbers and conclusions. Do not silently change the code while retaining the old results. No restricted data was copied for this audit.
2. **Other application scripts rely on defaults.** The current package's numeric `indeptest` method defaults to random tie-breaking. The yeast script omits `ties`; other numerical/discrete analyses also need an explicit audit of this setting. The current manuscript says its tied-data construction retains mid-ranks. This matters for bounded, ordinal and integer-valued variables, even when their R storage type is numeric.
3. **Historical code cannot simply be dropped into the current package.** All2 contains legacy `order`, `basis = "orthopoly"`, `ties.method` and `$pvalue` usage, local paths and assumptions about an existing workspace. Use it as provenance, not as an executable final package without adaptation and validation.
4. **Housing has alternative analysis branches and reusable caches.** The q=1 and “finalize Bn” scripts include an exploratory comparison of Bn variants. Identify the exact branch behind the reported q=7 results. Pin data, folds, parameters, package versions and seeds in cache metadata so a changed analysis cannot reuse stale screening results.
5. **Environment and execution instructions are incomplete.** A replication README, an explicit script-to-table/figure map, dependency versions, run order, expected outputs and runtime estimates are needed. The six package versions supplied for the manuscript do not cover all dependencies (e.g. minerva, copula, mvtnorm, data.table, e1071, xgboost and Python plotting libraries). Preserve the historical environment separately from any modern rerun environment.

## Data disposition

| Data | Recommended treatment |
|---|---|
| Brain injury supplied through Genest/coauthors | Exclude raw data and unreviewed workspace dumps. Supply preprocessing/analysis code, variable schema, lawful access instructions and, if needed, clearly labelled synthetic demonstration data. A synthetic demonstration is not replication of the empirical results. |
| Yeast Hughes dataset | Load via `HHG` and document the exact package/data version. Do not copy unrelated yeast/proteomics files from All2. |
| Kaggle Ames housing training data | Scripts use `data/house_train.csv`. Record the exact original source, file identity and preprocessing. Public download availability alone does not establish permission to redistribute: use official download instructions until the applicable terms are confirmed. The competition test file is not needed for the reported cross-validation. |
| Simulated datasets | Include generators and seeds. Retain compact final summaries and enough intermediate results to rebuild figures without rerunning expensive simulations; full multi-gigabyte generated CSVs are unnecessary if regeneration is specified. |

The original local housing training file can be fingerprinted without placing its contents in the repository. No research data was added during this audit.

## Recommended completion sequence

1. Resolve the brain-injury discrepancy and identify final versions of the other tied-data analyses.
2. Bring the verified missing scripts into an organized `inst/paper` tree; adapt paths and APIs explicitly, recording provenance and any changes. Add the still-unlocated numerical experiment drivers.
3. Add permitted data or precise retrieval instructions; keep restricted files outside the distributable tree. Add compact final reference outputs and a manifest linking them to manuscript labels.
4. Write a runbook and pin the environment. Separate quick reconstruction from saved summaries, a small pipeline check, and full simulation runs. Run from a fresh session and compare results against the manuscript, using Monte Carlo tolerances where appropriate.
5. Archive a fixed replication release alongside the manuscript. `.Rbuildignore` excludes both `inst/paper` and root `results`: a CRAN source tarball alone will **not** distribute these materials. Keeping replication materials out of the CRAN package is compatible with releasing a separate archive.

The [JASA Applications and Case Studies reproducibility guidance](https://jasa-acs.github.io/repro-guide/pages/author-guidelines.html) is a useful checklist for the documented workflow and restricted-data treatment. A completed source inventory by itself does not establish successful computational reproduction.
