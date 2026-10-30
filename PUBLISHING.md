# Updating the website

GitHub Pages serves `docs/` on the main branch at https://qiyan14.github.io/dmft-calculator/.

1. Edit `app/app.R` and test locally with Shiny.
2. Run `Rscript tests/test_selected_model.R app/app.R analysis/04_Export_DMFT_Equation_Age_BMI.R` from the repository root. These are synthetic implementation checks, not study validation.
3. Run `source("update_app_source.R")` to copy the app into `docs/app.json` and update the cache version in `docs/index.html`.
4. Commit and push the source and generated files. Wait for the Pages workflow to finish.

For new runtime dependencies or a new Shinylive version, use `build.R`; then restore the versioned app request in `docs/index.html` and run `update_app_source.R`.

## Age/BMI model status

The study uses numeric age and BMI with exhaustive single-marker and pair searches and beam search from size three. The selected six-biomarker age/BMI model is bundled in `PUBLISHED_MODEL` (section 1 of `app/app.R`) with the full-precision coefficients from the final fit (`DMFT_equation.txt`), the workbook threshold 0.4303368621, and the cohort and cross-validated performance from the results workbook. Its biomarker slopes equal the logarithms of the workbook's per-doubling odds ratios to within 10⁻⁹, which the test checks.

If the model is refitted, replace the bundled coefficients, threshold, cohort and performance together, then update the published equation in `tests/test_selected_model.R` and rerun the test.

The desktop app `analysis/06_Shiny_DMFT_Age_BMI.R` is `app/app.R` with one difference: it checks a `dmft_age_bmi_model.csv` in its working directory and uses it only when it matches the bundled model. Other models are loaded with **Load a saved model**.
