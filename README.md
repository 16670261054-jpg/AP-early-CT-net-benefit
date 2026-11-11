# AP-early-CT-net-benefit

Analysis code for: **Routine admission data identify a third of patients with acute pancreatitis in whom early CT adds little net benefit** (Mi Y, Han J; submitted).

Six routine admission variables (age, calcium, lymphocyte percentage, red blood cell count, creatinine, alanine aminotransferase) form a continuous risk and four admission phenotypes. In a Chinese surgical-ward cohort with CT on admission in every patient (n = 1,343), the code estimates the incremental value of the modified CT severity index (MCTSI) over this admission stratification and the net benefit of an MCTSI-based rule within an admission-defined low-risk stratum. Associations with the late course are also estimated in MIMIC-IV v3.1 (n = 4,439) and eICU-CRD v2.0 (n = 941).

## Layout

```
scripts/   Python analysis pipeline
  01_build_local.py           ward cohort; continuous risk (cross-fitted in the derivation period)
  02_build_mimic.py           MIMIC-IV cohort, ICD-10 subsets, coded outcomes
  03_build_eicu.py            eICU cohort, diagnosis-string outcomes
  m01_mimic_lipase_dod.py     MIMIC-IV lipase and date-of-death extract
  04_gradients.py             gradients, Wilson CIs, trend tests, adjusted ORs, likelihood-ratio tests, continuous risk
  05_ct_incremental.py        10 x 5-fold out-of-fold AUC, DeLong, NRI, IDI, calibration, decision curves
  05c_lowrisk_final.py        admission-defined low-risk stratum and stratum-specific decision curves
  05d_tableS3.py              Supplementary Tables S3-S4
  06_secondary.py             infection sites, competing-risk incidence, death timing, LOS, organ support
  07_extras.py                E-values, likelihood ratios, subgroup analyses, classifier agreement
  08_sens_lipase_death90.py   lipase-restricted and long-term mortality sensitivity analyses
  09_figdata_profile.py       phenotype profiles and continuous-risk figure data
rfig/      R (ggplot2) figure scripts: Figures 1-6 and Supplementary Figures S1-S7
models/    fitted six-variable GBDT, imputation medians, phenotype classifier (no patient-level data)
```

Run the Python scripts in numeric order (`m01` before `08`), then the R scripts from inside `rfig/`.

## Data access

- MIMIC-IV and eICU-CRD are available to credentialed users at [PhysioNet](https://physionet.org). Patient-level extracts are not redistributed, in line with the data use agreements.
- The ward cohort dataset is available from the corresponding author on reasonable request.
- Input locations are set with environment variables (`AP_PROBES`, `AP_SAP_INPUTS`, `PHENO_INPUTS`, `MIMIC_HOSP`, `EICU_DIR`) and default to `inputs/` in the repository root. Scripts write patient-level working files to `data/` and aggregate statistics to `results/`, from which all figures and tables are produced.

## Dependencies

Python 3.11 (pandas, numpy, scikit-learn, statsmodels, scipy, joblib, xgboost); R 4.6 (ggplot2, data.table, patchwork, ggdist, scales, stringr).
