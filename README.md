# Salivary biomarker DMFT calculator

[Open the calculator](https://qiyan14.github.io/dmft-calculator/)

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

Run `shiny::runApp("06_Shiny_DMFT_Age_BMI.R")`, or open the file in RStudio and click Run App. Shiny is required; readxl and writexl are optional for Excel uploads/downloads. In Selected DMFT model, upload the newly exported `dmft_age_bmi_model.csv`. For automatic desktop loading, put that CSV in the app's working directory before starting it.

Enter raw biomarker values, age in years and BMI in kg/m². The app applies `log2(raw + 5e-7)` exactly once; do not standardize the inputs. The age slider and overlapping comparison curves show the entered BMI versus BMI five units higher, holding biomarkers fixed. These are model scenarios, not causal effects.

The selected panel, intercept, coefficients and classification threshold are read from the new export; no panel size is hardcoded. Previous age/sex coefficients and performance are not valid for the revised model. Until you run the new analysis and load its CSV, selected-model risk calculation is unavailable. Training on uploaded data remains available, with age and BMI selected by default when present and exhaustive sizes 1–2 followed by beam search. The browser refuses pair searches above 250,000 pairs rather than silently narrowing them; use the desktop analysis script for larger studies.

## Verification and limits

Run `Rscript test_age_bmi.R` from this folder. The base-R tests use synthetic fixtures to check all 111,156 pairs, beam extensions from size three, workbook panel selection, cohort exclusions, raw-equation equivalence, CSV import/export, age/BMI curves, and rejection of stale models. They do not estimate study performance.

All scripts were syntax-checked. The export and prediction tests passed on synthetic data. The complete study has not been rerun because the original participant files are not present. New coefficients, figures, cohort counts and performance must come from that run. Cross-validation after panel/threshold selection is developmental performance and may be optimistic; these changes do not provide independent validation.

## Repository paths

The numbered study scripts are in `analysis/`; the website source is `app/app.R`. `docs/app.json` contains the matching browser app. The previous figure and calibration filenames are retained as compatible copies of the new code.

Run the browser-source checks from the repository root with:

```r
# In a terminal:
# Rscript tests/test_selected_model.R app/app.R analysis/04_Export_DMFT_Equation_Age_BMI.R
```

The calculator runs in your browser through Shinylive. Entered data are processed locally. No participant data are included in the repository. See [PUBLISHING.md](PUBLISHING.md) for website updates.
