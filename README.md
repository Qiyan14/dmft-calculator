# Salivary biomarker DMFT calculator

### [Open the calculator](https://qiyan14.github.io/dmft-calculator/)

The selected panel now comprises six salivary biomarkers plus **age and sex**. It was the highest-performing retained panel across the four evaluated configurations (age + sex alone; biomarkers + age + sex; biomarkers + age + sex + BMI; biomarkers alone).

**Model status:** the uploaded results contain performance metrics and biomarker odds ratios, but not the complete fitted intercept, age coefficient, sex coefficient or reference coding. The old model's coefficients have therefore been removed. The selected model cannot calculate probabilities until its complete export is supplied. Training on uploaded data and loading complete saved models remain available.

## Selected panel and results

- Bacteroidetes G5 bacterium HMT 511 (B511)
- Prevotella sp. HMT 475 (P475)
- Unclassified Lactobacillus (LAC)
- Unclassified Bacteroidetes G5 (BG5)
- Bacteroides heparinolyticus (BH)
- Peptostreptococcaceae G3 bacterium HMT 950 (P950)

The development cohort contained 111 adults: 56 with DMFT ≤ 13 and 55 with DMFT > 13. Five-fold cross-validation repeated 100 times produced balanced accuracy **88.36%**, sensitivity **91.09%**, specificity **85.63%**, and AUC **0.890**, at a pooled Youden threshold of **0.4511200806**. The 2.5th–97.5th percentile range of repeat-specific balanced accuracy was 84.73%–90.12%; this describes fold-assignment variability, not an independent-validation confidence interval. Panel, adjustment configuration and threshold selection used these predictions, so performance is exploratory and may be optimistic.

## Export the selected model and figures

1. Run the four-model preprocessing script in R with the settings used for the results workbooks.
2. Put the three biomarker-scenario workbooks in the working directory. Keep the `DMFT_Four_Models_Resampling_<registry_key>.rds` file there when available.
3. Run `source("analysis/Figures_DMFT.V3.R")`. It chooses the winning workbook by balanced accuracy, with AUC breaking ties, and verifies that this configuration also wins at each panel size. For these results, that is biomarkers + age + sex. You may instead set `RESULTS_FILE` explicitly.
4. The script reproduces the common cohort (complete age, sex and BMI, even though BMI is not in the selected model), refits the panels, and checks their AUC, balanced accuracy and thresholds against the workbook. It stops on a mismatch.
5. Outputs include the manuscript balanced-accuracy/AUC performance plot, ROC curves (without a covariates-only line or fixed-specificity annotation), radar plot, blue/orange heatmap with black outlines around misclassified cells, individual out-of-fold predictions, and **`dmft_selected_model.csv`**. The CSV includes the unrounded coefficients, log2 offset, threshold and actual sex reference category. The raw-input equation is checked against the fitted model for every participant before export.

The figure script requires the original preprocessed participant data in the R session; aggregate Excel results alone cannot reconstruct coefficients or ROC curves. No participant data are included in this repository.

## Risk Calculator

On the website, choose **Load a saved model** and upload `dmft_selected_model.csv`. For desktop use, put the CSV next to `app.R` and restart to load the selected model automatically. To bundle it for all website visitors, place the verified CSV in `app/`, rebuild the site and publish it.

Enter biomarker measurements in the original absolute assay units, age in years, and sex using the categories recorded in the model export. Do not transform or standardize the values yourself. The app applies log2(raw + 5e-7); LD reference centering and scaling are absorbed into the exported coefficients. The equation is:

```
eta = intercept + sum_j beta_j * log2(raw_j + 5e-7)
                + beta_age * age + beta_sex * I(sex = nonreference category)
p(HD) = 1 / (1 + exp(-eta))
```

Classification is HD when the calculated probability is at least the saved model's threshold. Numeric coefficients will be available after the complete export is generated; they are not inferred from the previous model.

## Train Your Own Data

Upload participant information and biomarker measurements, match IDs, define outcome groups, choose biomarkers and adjustment variables, and run panel search. Age and sex are selected by default when present; sex is treated as a categorical variable even when coded numerically. Train a retained panel and download its complete model for prediction. User-trained models are distinct from the selected study model.

The app runs in your browser through Shinylive. Entered values are not sent to a server. CSV/TSV inputs and CSV outputs are supported; Excel support depends on optional packages. Large searches can consume substantial browser CPU and memory.

## Research use

This cross-sectional development model has not been externally validated. It does not replace a dental examination or predict future DMFT progression. Some biomarker estimates are unstable.

## Source

`app/app.R` is the Shiny source; `docs/app.json` contains the browser app and is served by GitHub Pages. `analysis/Figures_DMFT.V3.R` regenerates figures and exports the fitted model. See [PUBLISHING.md](PUBLISHING.md) for rebuild instructions. Manuscript in preparation.
