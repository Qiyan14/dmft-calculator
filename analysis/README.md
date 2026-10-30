# DMFT figures, importance, equation and calibration — new results

This bundle is updated for the supplied age/BMI summary workbook. It includes the workbook and the original preprocessing, search and Shiny scripts for continuity.

The new selected panel has six biomarkers: Lachnoanaerobaculum orale; Olsenella profusa; unclassified Bacteroidetes G5; Bacteroides heparinolyticus; Weeksellaceae G1 sp. HMT 900; and Peptostreptococcaceae G3 bacterium HMT 950. Age_DOS and BMI remain continuous adjustment variables. The summary reports 111 participants (56 LD, 55 HD), balanced accuracy 0.8807435065, AUC 0.8987005195, sensitivity 0.9032727273, specificity 0.8582142857 and threshold 0.4303368621. The BA endpoints 0.8469926948–0.9058238636 are repeat percentiles, not independent-validation confidence intervals.

## What changed

- `00_DMFT_Result_Reader.R` reads the new Best_By_Size and Best_Overall sheets, checks their agreement and preserves the actual Excel row provenance. Keep this helper beside the updated scripts.
- `03_Figures_DMFT_Age_BMI.R` reads all six winning panels directly from the summary and regenerates the ROC, radar, heatmap and performance figures. It checks that repeated-CV refits reproduce the saved BA, AUC and threshold before exporting them. The selected-panel coefficient CSV is regenerated from the matching fit. Original heatmap colors, non-biomarker mismatch highlighting, row outlines and ROC settings are retained.
- `08_Variable_Importance_DMFT_Age_BMI.R` adds held-out permutation importance for the six selected biomarkers plus age and BMI. It runs automatically after script 03. Points show mean decreases in balanced accuracy in percentage points at the fixed saved threshold; whiskers show the 2.5th–97.5th percentiles across 20 shuffles, not confidence intervals. It also exports AUC decreases and individual shuffle results. Negative importance is retained. Correlated variables can substitute for one another; this measures model reliance, not causality.
- `04_Export_DMFT_Equation_Age_BMI.R` takes the selected panel and threshold from the summary, refits its coefficients with the matching participant data, and checks raw-input predictions. For summary inputs, it requires the verified refits from script 03. The selected record is Best_By_Size, Excel row 7. The Shiny output remains `dmft_age_bmi_model.csv`.
- `05_Calibration_DMFT_Age_BMI.R` reads the new summary and requires the matching participant-level out-of-fold predictions. It stops on mismatched cohort, threshold or performance, then calculates calibration intercept, slope, Brier score, participant-bootstrap intervals and repeat percentiles. The calibration plot and results text use the summary's winning panel.

## Run order

1. Keep your matching age/BMI preprocessing session open, with `data.dfPH`, `registry_key` and the LD reference statistics. If you restarted R, rerun `01_Preprocessing_DMFT_Age_BMI.R` with exactly the settings used to create this summary. You do not need to repeat the full panel search if the saved resampling file exists.
2. Put the supplied summary workbook and its matching `DMFT_Age_BMI_Resampling_<registry_key>.rds` in your R working directory. The resampling file was saved by script 02. If the summary is elsewhere, set `RESULTS_FILE` in scripts 03/04 and `CAL_RESULTS_FILE` in script 05 to its full path; keep the resampling RDS beside the workbook.
3. Source `03_Figures_DMFT_Age_BMI.R` using its full script path. It reconstructs and verifies out-of-fold predictions, writes the figures/model CSV and automatically runs variable importance. Set `MAKE_VARIABLE_IMPORTANCE <- FALSE` in script 03 if you want to run script 08 separately later. Change `VI_PERMUTATIONS` in script 08 to adjust the number of shuffles.
4. Source `04_Export_DMFT_Equation_Age_BMI.R` for the full equation, coefficient table and model bundle. Keep the same session; its verified-data check prevents reuse of another run's coefficients.
5. Source `05_Calibration_DMFT_Age_BMI.R` for calibration results and plots. On its own it can read the new `OOF_predictions_k<k>_<tag>.csv` files from `OUT_DIR`; the k=6 file is sufficient to calibrate the selected model alone.
6. The selected age/BMI model is bundled in the Shiny app (`app/app.R` and `06_Shiny_DMFT_Age_BMI.R`). After a refit, compare the new `dmft_age_bmi_model.csv` with the bundled coefficients before updating them; see `PUBLISHING.md`.

Use full script paths when sourcing because preprocessing changes R's working directory. `readxl`, `caret` and `pROC` are needed for figures; equation/summary input uses `readxl`; permutation calculations use base R. Calibration uses base R plus `readxl`. The remaining scripts keep their existing package requirements.

## Available data and verification

The attached workbook contains panel-level metrics only. It has no participant-level predictions, fitted coefficients or calibration estimates. Therefore numerical calibration results, actual permutation importance, a fitted risk equation, ROC curves and participant heatmaps cannot be calculated from this file alone. No previous model values or synthetic test values have been substituted.

The reader was tested against all six actual summary rows. Synthetic fixtures verified raw-equation equivalence and Shiny import for all 111 test participants, held-out baseline metrics, reproducible permutations for eight predictors, calibration/bootstrap calculations, provenance, and rejection of mismatched inputs. All R files parse. The importance plot layout was checked on synthetic data. The full study rerun remains to be done using your matching participant data and saved folds.

From the extracted bundle folder, run `Rscript test_summary_pipeline.R` for the new tests, or `Rscript test_age_bmi.R` for the earlier search/calculator tests. Test data and temporary test equations are synthetic and must not be reported as study results.
