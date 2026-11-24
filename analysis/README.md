# Study analysis scripts

These scripts produce the study's results: logistic regression panels of one to five salivary biomarkers, each fitted with age, ranked by balanced accuracy at the Youden threshold under five-fold cross-validation repeated five times; nested cross-validation of the whole procedure; and the calculator's single-biomarker equation. Their settings are those of the study's analysis: seed 777, high DMFT defined as DMFT > 13.9, panels of one to five biomarkers, and two covariate sets, biomarker panels with age (`Biomarkers_Age`) and age alone (`Age`). With the study data, script 01 reports 111 participants (56 LD, 55 HD).

The participant data are not in this repository. Each file name keeps the script's version number.

| Script | What it does | Results |
|---|---|---|
| `01_Preprocessing_DMFT.V8.R` | Reads and merges the three source Excel files; k-nearest-neighbour imputation; log2(raw + 5 × 10⁻⁷); standardization to the LD group; LD and HD groups; the covariate sets (`add_set()` lines) | Cohort |
| `02_LogReg_DMFT.V8.R` | Panel search: every single biomarker and all 111,156 pairs of the 472 candidates, then beam search keeping the 50 best panels of each size, to five biomarkers; age alone on the same folds | Supplementary Table S1 |
| `09_Nested_CV_DMFT.V8.R` | Nested cross-validation: in each of five outer training sets, the whole search (inner five-fold cross-validation repeated five times) and the threshold; the best panel of each size, age alone, lasso, elastic net (mixing parameter 0.5) and random forest (500 trees) | Tables 2 and 3 |
| `04b_Export_DMFT_Equation_Single_Biomarker.V9.R` | Refits the best single biomarker with age to all participants, checks the refit against script 02's workbook, and writes the raw-input equation, `dmft_age_model.csv` and the `PUBLISHED_MODEL` block for `app/app.R` | Section 3.4, the calculator |
| `05b_Calibration_10rep_DMFT.V11.R` | Rebuilds the repeated cross-validation predictions of the best panel of each size on script 02's folds, checks them against the workbooks, and computes calibration intercept and slope and the Brier score | Supplementary File S1 |
| `05c_Calibration_Repeated_Nested_DMFT.V16.R` | ROC curves of each panel size and calibration of the best panel under repeated (05b) and nested (09) cross-validation | Figure 1 (nested), Supplementary Figure S1 (repeated) |
| `07_Table1_DMFT_Age_BMI.R` | Descriptive statistics of the cohort by group | Table 1 (a subset of its rows) |
| `06_Shiny_DMFT.R` | The desktop calculator (see the main README) | |

## Run order

Use a fresh R session and run the scripts in this order, in the same session, using their full paths: script 01 changes the working directory.

1. `01_Preprocessing_DMFT.V8.R`. Set `DATA_DIR_CHOICES` (section 1) to the folder that holds the three source Excel files and check their file names and sheets. The script then works in that folder's `Shared with Qi` subfolder, which must exist; the other scripts read and write there.
2. `02_LogReg_DMFT.V8.R`. Writes the results workbook of each covariate set (`Results_LogReg_HD_vs_LD_<set>_<registry_key>.xlsx`), the summary files and the fold indices.
3. `09_Nested_CV_DMFT.V8.R`. Reads script 02's `DMFT_Covariate_Sets_Summary_<registry_key>.xlsx` and writes `DMFT_Nested_CV_<registry_key>.xlsx`. This is the longest step.
4. `04b_Export_DMFT_Equation_Single_Biomarker.V9.R`. Writes to `DMFT_single_biomarker_equation_Biomarkers_Age`.
5. `05b_Calibration_10rep_DMFT.V11.R`. Writes to `DMFT_10rep_outputs_Biomarkers_Age`.
6. `05c_Calibration_Repeated_Nested_DMFT.V16.R`. Reads 05b's folder and 09's workbook and writes to `DMFT_calibration_repeated_nested_Biomarkers_Age`; `Figure_ROC_Calibration_nested_max5.png` is Figure 1 and `Figure_ROC_Calibration_repeated_max5.png` is Supplementary Figure S1 (ROC curves beside the calibration of the best panel).
7. `07_Table1_DMFT_Age_BMI.R`, at any point after script 01.

The registry key, which ends every output file name, is `DMFT_cutoff13.9_conf_AgeDOS_allBiomarkers_noRatios_noVif_exhaustive12_beam3` with these settings. Scripts 02 and 09 use 10 cores (`numCores`, `NESTED_CORES`).

## Packages

- 01: readxl, dplyr, VIM
- 02: caret, pROC, foreach, doSNOW, dplyr, writexl; fastglm is used when installed
- 09: readxl, writexl, glmnet and ranger (rstudioapi only for choosing a file in a browser)
- 04b, 05c: readxl; 05b: readxl and caret (fastglm when installed)
- 07: dplyr, gtsummary, gt, flextable, officer

## Checks

The scripts check their inputs against each other: 04b refits the model and compares it with script 02's odds ratios and fitted log-odds; 05b reproduces script 02's out-of-fold predictions and script 09's nested estimates before using them; 05c reproduces 05b's and 09's numbers. `tests/test_selected_model.R` runs 04b on a synthetic workbook and checks that the calculator imports its model file with identical predictions; see the main README.

Cross-validation after panel and threshold selection describes development performance and can be optimistic; the nested estimates describe the selection procedure. Independent validation is still required.
