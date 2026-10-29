# Salivary biomarker DMFT calculator

### [Open the calculator](https://qiyan14.github.io/dmft-calculator/)

The selected panel now comprises six salivary biomarkers plus **age and sex**. It was the highest-performing retained panel across the four evaluated configurations (age + sex alone; biomarkers + age + sex; biomarkers + age + sex + BMI; biomarkers alone).

**Model status:** the complete fitted coefficients are bundled for the selected six-biomarker age-and-sex model. Sex is entered as female or male: female (stored code 1) is the reference category, and male (code 2) adds 0.11282243181054701 to the linear predictor. Training on uploaded data and loading complete saved models remain available.

## Selected panel and results

- Bacteroidetes G5 bacterium HMT 511 (B511)
- Prevotella sp. HMT 475 (P475)
- Unclassified Lactobacillus (LAC)
- Unclassified Bacteroidetes G5 (BG5)
- Bacteroides heparinolyticus (BH)
- Peptostreptococcaceae G3 bacterium HMT 950 (P950)

The development cohort contained 111 adults: 56 with DMFT ≤ 13 and 55 with DMFT > 13. Five-fold cross-validation repeated 100 times produced balanced accuracy **88.36%**, sensitivity **91.09%**, specificity **85.63%**, and AUC **0.890**, at a pooled Youden threshold of **0.4511200806**. The 2.5th–97.5th percentile range of repeat-specific balanced accuracy was 84.73%–90.12%; this describes fold-assignment variability, not an independent-validation confidence interval. Panel, adjustment configuration and threshold selection used these predictions, so performance is exploratory and may be optimistic.

Calibration of the same pooled out-of-fold predictions: intercept **−0.01** (95% CI −0.37 to 0.35), slope **0.71** (95% CI 0.39 to 1.40) and Brier score **0.124** (95% CI 0.086 to 0.167), with participant-bootstrap confidence intervals. The slope below 1 suggests that predicted probabilities are too extreme, as expected from overfitting.

## Export the selected model and figures

1. Run the four-model preprocessing script in R with the settings used for the results workbooks.
2. Put the three biomarker-scenario workbooks in the working directory. Keep the `DMFT_Four_Models_Resampling_<registry_key>.rds` file there when available.
3. Run `source("analysis/Figures_DMFT.V3.R")`. It chooses the winning workbook by balanced accuracy, with AUC breaking ties, and verifies that this configuration also wins at each panel size. For these results, that is biomarkers + age + sex. You may instead set `RESULTS_FILE` explicitly.
4. The script reproduces the common cohort (complete age, sex and BMI, even though BMI is not in the selected model), refits the panels, and checks their AUC, balanced accuracy and thresholds against the workbook. It stops on a mismatch.
5. Outputs include the manuscript balanced-accuracy/AUC performance plot, ROC curves (without a covariates-only line or fixed-specificity annotation), radar plot, cyan/pink heatmap with highlighted misclassified rows and black row outlines, individual out-of-fold predictions, and **`dmft_selected_model.csv`**. The CSV includes the unrounded coefficients, log2 offset, threshold and actual sex reference category. The raw-input equation is checked against the fitted model for every participant before export.
6. Run `source("analysis/Calibration_DMFT.R")` in the same session. It computes the calibration intercept, calibration slope and Brier score of each panel and of the age + sex baseline from the pooled out-of-fold predictions, with participant-bootstrap confidence intervals and repeat percentiles, draws the calibration plot (manuscript Figure 1B), and writes paste-ready results text. Run on its own, it reads the out-of-fold prediction files from step 5.

The figure script requires the original preprocessed participant data in the R session; aggregate Excel results alone cannot reconstruct coefficients or ROC curves. No participant data are included in this repository.

## Risk Calculator

Choose **Selected DMFT model**, enter the six raw biomarker measurements, choose **Female** or **Male** and use the **Age (years)** slider. Click **Calculate risk** to display the gauge and calculation. The comparison graph overlays risk curves for females and males, holds biomarkers fixed and marks the selected age. It updates automatically as inputs change; after a calculation, the gauge and equation update too. The adult slider spans 18–100 years for scenario exploration, not a validated prediction range.

You can also choose **Load a saved model** to upload a complete model CSV. The age/sex graph is shown only when the model contains numeric age and categorical sex.

Enter biomarker measurements in the original absolute assay units, age in years, and sex. Models trained on other data show sex in the codes of their own data, with labels if the saved model includes them. Do not transform or standardize the values yourself; enter 0 for a biomarker that was not detected. The app applies log2(raw + 5e-7); LD reference centering and scaling are absorbed into the exported coefficients. The equation is:

```
eta = intercept + sum_j beta_j * log2(raw_j + 5e-7)
                + beta_age * age + beta_male * I(male)      (female is the reference)
p(HD) = 1 / (1 + exp(-eta))
```

Classification is HD when the estimated risk is at least the saved model's threshold. Calculations use the full-precision coefficients supplied from the final selected-model fit.

`tests/test_selected_model.R` checks the implementation: run `Rscript tests/test_selected_model.R app/app.R` from the repository folder. It compares the app's calculation with the published equation for 10,000 random inputs, saves and reloads the model, and checks the sex labels and the example shown in the manuscript figures.

## Train Your Own Data

Upload participant information and biomarker measurements, match IDs, define outcome groups, choose biomarkers and adjustment variables, and run panel search. Age and sex are selected by default when present; sex is treated as a categorical variable even when coded numerically. Train a retained panel and download its complete model for prediction. User-trained models are distinct from the selected study model.

The app runs in your browser through Shinylive. Entered values are not sent to a server. CSV/TSV inputs and CSV outputs are supported; Excel support depends on optional packages. Large searches can consume substantial browser CPU and memory. For large searches or Excel files, run the same app in R on your computer with `shiny::runApp("app")` (requires the shiny package; readxl and writexl are optional).

## Research use

This cross-sectional development model has not been externally validated. It does not replace a dental examination or predict future DMFT progression. Some biomarker estimates are unstable.

## Source

`app/app.R` is the Shiny source; `docs/app.json` contains the browser app and is served by GitHub Pages. `analysis/Figures_DMFT.V3.R` regenerates figures and exports the fitted model; `analysis/Calibration_DMFT.R` adds the calibration analysis. `tests/test_selected_model.R` checks the calculator against the published equation. See [PUBLISHING.md](PUBLISHING.md) for rebuild instructions. Manuscript in preparation.
