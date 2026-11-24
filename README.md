# Salivary biomarker DMFT calculator

[Open the calculator](https://qiyan14.github.io/dmft-calculator/)

**Model status:** the calculator bundles the single-biomarker study model, *Bacteroides heparinolyticus* (BH) with age in years, with the full-precision coefficients of its final fit to the 111 study participants (`dmft_age_model.csv`, written by `analysis/04b_Export_DMFT_Equation_Single_Biomarker.V9.R`). Under repeated cross-validation, balanced accuracy rose with each added biomarker, but when panel selection was repeated within each training set (nested cross-validation), the best single biomarker generalized better (balanced accuracy 82.0% versus 70.4% to 78.4% with two to five biomarkers) and had a higher balanced accuracy than age alone (78.4%) and than lasso, elastic-net and random forest models that used more biomarkers (66.6% to 77.5%). The classification threshold is 0.3834346393. In the development cohort (111 adults: 56 with DMFT ≤ 13.9 and 55 with DMFT > 13.9), five-fold cross-validation repeated 5 times gave balanced accuracy 83.29%, sensitivity 88.73%, specificity 77.86% and AUC 0.840. Panel and threshold selection used these predictions, so this performance is optimistic; independent validation is required.

Every study panel is a logistic regression of one to five biomarkers together with age in years; BMI and sex are not covariates. The search scored every single biomarker and all 111,156 pairs of the 472 candidates, then built larger panels by beam search, adding each remaining biomarker to each of the 50 highest-ranked panels of the preceding size; it does not evaluate every possible panel of three or more biomarkers. Panels were capped at five biomarkers so that, with age, a model has at most six predictors.

The scripts in `analysis/` run this analysis with the study's settings (seed 777, five-fold cross-validation repeated five times, one to five biomarkers); [analysis/README.md](analysis/README.md) describes what each one does and writes.

## Run the study analysis

The participant data are not in this repository. Use a fresh R session and run these scripts from `analysis/` in order, in the same session and with their full paths (script 01 changes the working directory):

1. `01_Preprocessing_DMFT.V8.R`: preprocessing and the two covariate sets, biomarker panels with age and age alone. Set `DATA_DIR_CHOICES` to the folder that holds the three source Excel files and check their file names and sheets; the outputs of every script go to its `Shared with Qi` subfolder, which must exist. With the study data it reports 111 participants (56 LD, 55 HD). The all-biomarker DMFT analysis does not need `script1_panels.R`.
2. `02_LogReg_DMFT.V8.R`: panel search and repeated cross-validation.
3. `09_Nested_CV_DMFT.V8.R`: nested cross-validation of each panel size, age alone, lasso, elastic net and random forest.
4. `04b_Export_DMFT_Equation_Single_Biomarker.V9.R`: the calculator's equation, `dmft_age_model.csv` and a `PUBLISHED_MODEL` block for `app/app.R`.
5. `05b_Calibration_10rep_DMFT.V11.R`: calibration under repeated cross-validation.
6. `05c_Calibration_Repeated_Nested_DMFT.V15.R`: ROC curves and calibration under repeated and nested cross-validation.
7. `07_Table1_DMFT_Age_BMI.R`: descriptive statistics of the cohort, at any point after script 01.

## Shiny calculator

Run `shiny::runApp("analysis/06_Shiny_DMFT.R")` from the repository folder, or open the file in RStudio and click Run App. Shiny is required; readxl and writexl are optional for Excel uploads/downloads. The online calculator and the desktop app contain the same bundled selected model. If a `dmft_age_model.csv` is placed in the desktop app's working directory, it is used only when it matches the bundled coefficients and threshold; use **Saved model** for any other model.

Choose **Study DMFT model (trained on the 111 study participants)**, enter the raw *B. heparinolyticus* (BH) value and set age with the slider (years), then click **Calculate risk**. **Fill example values** enters illustrative values that do not come from a participant (0 for BH; age 40; estimated risk 55.8%, above the 38.3% threshold). The app applies `log2(raw + 5e-7)` exactly once; do not standardize the inputs, and enter 0 for a biomarker that was not detected. After a calculation, the gauge, equation and calculation table update as the inputs change.

Training on uploaded data remains available, with age and BMI selected by default when present and exhaustive sizes 1–2 followed by beam search. The browser refuses pair searches above 250,000 pairs rather than silently narrowing them; use the desktop analysis script for larger studies.

## Verification and limits

`tests/test_selected_model.R` compares the app's calculation with the published equation (BH and age) for 10,000 random inputs (including zero values and ages 18–100), checks that the biomarker slopes equal the logarithms of the workbook's per-doubling odds ratios, saves and reloads the model, and checks the example values. It also checks that all 111,156 pairs are evaluated and that beam extensions start at size three, and it runs `04b_Export_DMFT_Equation_Single_Biomarker.V9.R` on a synthetic results workbook: the model file it writes must give the same predictions in the calculator as the fitted model, the desktop app must refuse it in place of the bundled model, and stale or incomplete model files must be rejected (this part needs readxl and writexl). These are implementation checks on synthetic data; they do not estimate study performance.

Cross-validation after panel and threshold selection is developmental performance and may be optimistic; independent validation is still required.

## Repository paths

The numbered analysis scripts are in `analysis/`; the website source is `app/app.R`. `docs/app.json` contains the matching browser app.

Run the browser-source checks from the repository root with:

```r
# In a terminal:
# Rscript tests/test_selected_model.R app/app.R analysis/04b_Export_DMFT_Equation_Single_Biomarker.V9.R
```

The calculator runs in your browser through Shinylive. Entered data are processed locally. No participant data are included in the repository. See [PUBLISHING.md](PUBLISHING.md) for website updates.
