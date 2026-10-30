# Salivary biomarker DMFT calculator

[Open the calculator](https://qiyan14.github.io/dmft-calculator/)

The latest figure, equation and calibration scripts read the new age/BMI summary workbook. See [analysis/README.md](analysis/README.md) for the revised run order and the added held-out permutation-importance figure.

**Model status:** the calculator bundles the full-precision coefficients of the selected six-biomarker age/BMI model from the final fit (`DMFT_equation.txt`): *Lachnoanaerobaculum orale* (LO), *Olsenella profusa* (OP), unclassified Bacteroidetes G5 (BG5), *Bacteroides heparinolyticus* (BH), Weeksellaceae G1 sp. HMT 900 (W900) and Peptostreptococcaceae G3 bacterium HMT 950 (P950), with age in years and BMI in kg/m². The classification threshold is 0.4303368621. In the development cohort (111 adults: 56 with DMFT ≤ 13 and 55 with DMFT > 13), five-fold cross-validation repeated 100 times gave balanced accuracy 88.07%, sensitivity 90.33%, specificity 85.82% and AUC 0.899. Panel and threshold selection used these predictions, so performance is exploratory and may be optimistic.

This version adjusts every study biomarker panel for **Age_DOS and BMI**, both continuous in their original units. Sex is no longer a model covariate or a complete-case eligibility requirement. Table 1 can still describe sex and other participant characteristics.

The study search evaluates every single biomarker and every pair, then uses beam search from three through six biomarkers. With 472 candidates there are 111,156 pairs. Beam search extends retained panels at the preceding size; it does not evaluate all possible panels of size three or larger. The existing sweep, finalist scoring, cross-validation settings, ranking, and repeat-percentile outputs are preserved. Runtime increases because every pair is evaluated.

## Run the study analysis

Use a fresh R session. Unzip these files into one folder. In `01_Preprocessing_DMFT_Age_BMI.R`, set `DATA_DIR_CHOICES` to the directory containing your three source Excel files and check the input filenames/sheets. The script retains your existing `Shared with Qi` output directory; ensure that directory exists. The default all-biomarker DMFT analysis does not need `script1_panels.R`; that file is required only if you change to selected/publication panels.

Run these scripts in the same session using their full paths (preprocessing changes the working directory):

1. `01_Preprocessing_DMFT_Age_BMI.R`: preprocessing, LD reference means/SDs, age/BMI registry.
2. `02_LogReg_DMFT_Age_BMI.R`: adjusted panel search; saves the workbook and exact resampling indices.
3. `03_Figures_DMFT_Age_BMI.R`: selects the best panel at each size, reproduces out-of-fold predictions and checks them against the workbook; writes ROC, radar, heatmap, performance plots and the selected-model CSV.
4. `04_Export_DMFT_Equation_Age_BMI.R`: selects the best overall panel from the exported Excel workbook by balanced accuracy, then AUC, with smaller panels breaking exact ties. Refits using matching preprocessed participant data and writes the full-precision risk equation and `dmft_age_bmi_model.csv`.
5. `05_Calibration_DMFT_Age_BMI.R`: calibration plots, Brier score and participant-bootstrap intervals from the new out-of-fold predictions.
6. `07_Table1_DMFT_Age_BMI.R`: descriptive table for the same age/BMI complete-case cohort. Smoking and diabetes rows remain excluded.

The figure and equation scripts derive the workbook name automatically. Set their `RESULTS_FILE` only when using another location. New outputs contain the model tag `Biomarkers_Age_BMI` and registry suffix `_exhaustive12_beam3`; old age/sex workbooks are rejected.

The heatmap retains its original colors, black row outlines, and mismatch highlighting confined to the non-biomarker annotations. The ROC omits the covariates-only line and 90%-specificity annotation.

## Shiny calculator

Run `shiny::runApp("06_Shiny_DMFT_Age_BMI.R")`, or open the file in RStudio and click Run App. Shiny is required; readxl and writexl are optional for Excel uploads/downloads. The online calculator and the desktop app contain the same bundled selected model. If a `dmft_age_bmi_model.csv` is placed in the desktop app's working directory, it is used only when it matches the bundled coefficients and threshold; use **Load a saved model** for any other model.

Choose **Selected DMFT model**, enter the six raw biomarker values, set age with the slider (years) and enter BMI (kg/m²), then click **Calculate risk**. **Fill example values** enters illustrative values that do not come from a participant (1, 1, 1, 0, 0, 0 for LO, OP, BG5, BH, W900, P950; age 40; BMI 25; estimated risk 73.8%). The app applies `log2(raw + 5e-7)` exactly once; do not standardize the inputs, and enter 0 for a biomarker that was not detected. After a calculation, the gauge, equation and calculation table update as the inputs change. There is no age comparison chart, because BMI is continuous rather than categorical.

Training on uploaded data remains available, with age and BMI selected by default when present and exhaustive sizes 1–2 followed by beam search. The browser refuses pair searches above 250,000 pairs rather than silently narrowing them; use the desktop analysis script for larger studies.

## Verification and limits

`tests/test_selected_model.R` compares the app's calculation with the published age/BMI equation for 10,000 random inputs (including zero values, ages 18–100 and BMI 12–70), checks that the biomarker slopes equal the logarithms of the workbook's per-doubling odds ratios, saves and reloads the model, and checks the example values. It also uses synthetic fixtures to check all 111,156 pairs, beam extensions from size three, workbook panel selection, raw-equation equivalence, CSV import/export and rejection of stale models. These are implementation checks; they do not estimate study performance.

Cross-validation after panel and threshold selection is developmental performance and may be optimistic; independent validation is still required.

## Repository paths

The numbered study scripts are in `analysis/`; the website source is `app/app.R`. `docs/app.json` contains the matching browser app. The previous figure and calibration filenames are retained as compatible copies of the new code.

Run the browser-source checks from the repository root with:

```r
# In a terminal:
# Rscript tests/test_selected_model.R app/app.R analysis/04_Export_DMFT_Equation_Age_BMI.R
```

The calculator runs in your browser through Shinylive. Entered data are processed locally. No participant data are included in the repository. See [PUBLISHING.md](PUBLISHING.md) for website updates.
