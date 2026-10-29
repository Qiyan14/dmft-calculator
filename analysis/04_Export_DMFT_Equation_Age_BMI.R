# Export the best DMFT panel from the ONE workbook specified below.
# Run 01_Preprocessing_DMFT_Age_BMI.R first, then source this file
# in the SAME R session. No panel search or cross-validation is rerun.
# data.dfPH biomarkers must be LD-standardized log2(raw + 5e-7).
# Age_DOS and BMI retain original units.
# Requires readxl for Excel input. All fitting/export operations use base R.
#
# OUTPUTS:
#   dmft_age_bmi_model.csv       complete model for the updated Shiny app
#   dmft_age_bmi_coefficients.csv   identical compatibility copy
#   DMFT_equation.txt             full-precision equation and input definitions
#   DMFT_equation_model.rds       fit, source result, coding, transforms/checks
#   DMFT_panel_selection.csv      best retained panel at each size in this workbook
#
# In Shiny: choose 'Load a saved model' and upload dmft_age_bmi_model.csv.
# Alternatively place that file next to app.R and restart the desktop app.
# A website's built-in model requires bundling the verified CSV and rebuilding.
# Aggregate results alone cannot recover the intercept/covariate coefficients;
# the original preprocessed participant data are required for the final refit.

RESULTS_FILE <- NULL # NULL derives the age/BMI workbook name from registry_key
RESULTS_DIR <- getwd()
OUTPUT_DIR <- getwd()
# Set RESULTS_FILE to its full path if it is not in RESULTS_DIR.
# Only this workbook is read; no other scenario workbooks are required.

read_dmft_workbook <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Install readxl first: install.packages('readxl')")
  pars <- as.data.frame(readxl::read_excel(path, sheet = "Analysis_Parameters", col_types = "text"))
  sheets <- grep("^Results_[0-6]_Biomarkers$", readxl::excel_sheets(path), value = TRUE)
  list(path = path, parameters = pars,
       tables = setNames(lapply(sheets, function(sh)
         as.data.frame(readxl::read_excel(path, sheet = sh))), sheets))
}

dmft_parameter <- function(wb, name) {
  p <- wb$parameters
  if (!all(c("Parameter", "Value") %in% names(p))) stop("Invalid Analysis_Parameters sheet: ", wb$path)
  v <- unique(trimws(as.character(p$Value[which(trimws(p$Parameter) == name)])))
  v <- v[!is.na(v) & nzchar(v)]
  if (length(v) != 1L) stop("Expected one value for '", name, "' in ", wb$path)
  v
}

select_dmft_panel <- function(workbooks, registry_key) {
  candidates <- list()
  for (wb in workbooks) {
    key <- dmft_parameter(wb, "Script 1 registry key")
    if (!identical(trimws(as.character(registry_key)), key))
      stop("Workbook and preprocessing registry keys differ.\nWorkbook: ", key,
           "\nR session: ", registry_key, "\nUse the matching age/BMI preprocessing settings.")
    if (dmft_parameter(wb, "Reference group") != "LD" ||
        dmft_parameter(wb, "Case group (positive)") != "HD") stop("Expected LD reference and HD positive class.")
    config <- dmft_parameter(wb, "Model configuration")
    conf <- dmft_parameter(wb, "Confounders")
    covs <- if (tolower(conf) == "none") character(0) else trimws(strsplit(conf, ",", fixed = TRUE)[[1]])
    allowed <- list(Biomarkers_Age_BMI = c("Age_DOS", "BMI"))
    if (!(config %in% names(allowed)) || !identical(covs, allowed[[config]]) ||
        !grepl("_exhaustive12_beam3$", registry_key) ||
        dmft_parameter(wb, "k = 2 search") != "exhaustive")
      stop("Use rerun age/BMI results with exhaustive two-marker search.")
    expected <- paste0("Results_", 1:6, "_Biomarkers")
    if (!all(expected %in% names(wb$tables))) stop("Missing panel-size sheets in ", wb$path)
    for (sh in expected) {
      k <- as.integer(gsub("[^0-9]", "", sh))
      x <- wb$tables[[sh]]
      need <- c(if (k) paste0("Biomarker", seq_len(k)), "Balanced_Accuracy", "AUC", "Threshold")
      if (!all(need %in% names(x))) stop("Missing result columns in ", wb$path, ": ", sh)
      x$.source_row <- seq_len(nrow(x))
      if (!is.numeric(x$Balanced_Accuracy) || !is.numeric(x$AUC)) stop("Non-numeric performance columns.")
      x <- x[is.finite(x$Balanced_Accuracy) & is.finite(x$AUC), , drop = FALSE]
      if (!nrow(x)) stop("No usable results in ", wb$path, ": ", sh)
      x <- x[order(-x$Balanced_Accuracy, -x$AUC, x$.source_row), , drop = FALSE]
      candidates[[length(candidates) + 1L]] <- list(k = k, row = x[1, , drop = FALSE],
        sheet = sh, workbook = wb, configuration = config, covariates = covs)
    }
  }
  if (!length(candidates)) stop("No candidate panels found.")
  # Rank all configurations/sizes by BA, then AUC. Exact ties: fewer markers,
  # fewer covariates, supplied workbook order, then original workbook row order.
  ranking <- order(-vapply(candidates, function(x) x$row$Balanced_Accuracy, numeric(1)),
                   -vapply(candidates, function(x) x$row$AUC, numeric(1)),
                   vapply(candidates, `[[`, integer(1), "k"),
                   vapply(candidates, function(x) length(x$covariates), integer(1)))
  list(winner = candidates[[ranking[1]]], candidates = candidates[ranking])
}

export_dmft_equation <- function(data, reference_stats, registry_key,
                                output_dir = getwd(), results_file = RESULTS_FILE,
                                results_dir = getwd()) {
  if (length(registry_key) != 1L || is.na(registry_key) || !nzchar(trimws(registry_key)))
    stop("registry_key must be one nonempty string.")
  if (is.null(results_file)) results_file <- sprintf("Results_LogReg_HD_vs_LD_Biomarkers_Age_BMI_%s.xlsx", registry_key)
  if (length(results_file) != 1L || is.na(results_file) || !nzchar(results_file))
    stop("Specify exactly one RESULTS_FILE.")
  results_files <- results_file
  if (!file.exists(results_files)) results_files <- file.path(results_dir, results_file)
  if (!file.exists(results_files)) stop("Results workbook not found: ", results_files)
  selected <- select_dmft_panel(lapply(results_files, read_dmft_workbook), registry_key)
  winner <- selected$winner; wb <- winner$workbook
  markers <- if (winner$k) as.character(unlist(winner$row[paste0("Biomarker", seq_len(winner$k))], use.names = FALSE)) else character(0)
  if (anyNA(markers) || any(!nzchar(markers)) || anyDuplicated(markers)) stop("Invalid selected markers.")
  covariates <- winner$covariates; predictors <- c(markers, covariates)
  if (anyDuplicated(predictors)) stop("Biomarkers and covariates overlap.")
  threshold <- as.numeric(winner$row$Threshold); delta <- 5e-7
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1) stop("Invalid saved threshold.")
  message("Selected ", winner$configuration, ": ", winner$k, " biomarkers; BA = ",
          signif(winner$row$Balanced_Accuracy, 7), "; AUC = ", signif(winner$row$AUC, 7))
  message("Covariates: ", if (length(covariates)) paste(covariates, collapse = ", ") else "none")

  # Replicate the age/BMI complete-case cohort from run_dmft_age_bmi().
  required <- unique(c("condL", "Age_DOS", "BMI", predictors))
  if (anyDuplicated(names(data))) stop("Duplicate data column names.")
  if (length(setdiff(required, names(data)))) stop("Missing data columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  d <- as.data.frame(data)
  if (!all(vapply(d[c("Age_DOS", "BMI", markers)], is.numeric, logical(1)))) stop("Age, BMI and preprocessed biomarkers must be numeric.")
  for (v in c("Age_DOS", "BMI")) if (any(!is.na(d[[v]]) & !is.finite(d[[v]]))) stop(v, " contains infinite values.")
  keep <- complete.cases(d[, c("Age_DOS", "BMI"), drop = FALSE]) & !is.na(d$condL)
  source_rows <- which(keep); d <- droplevels(d[keep, , drop = FALSE])
  if (!setequal(unique(as.character(d$condL)), c("LD", "HD"))) stop("Expected only LD and HD outcomes.")
  d$condL <- factor(as.character(d$condL), levels = c("LD", "HD"))
  if (length(markers) && any(!is.finite(as.matrix(d[markers])))) stop("Non-finite preprocessed biomarkers remain.")
  counts <- table(d$condL)
  # Verify the common-cohort size against the specified workbook.
  for (x in selected$candidates) {
    expected <- as.integer(c(dmft_parameter(x$workbook, "n (LD)"), dmft_parameter(x$workbook, "n (HD)")))
    if (!identical(as.integer(counts), expected) || nrow(d) != as.integer(dmft_parameter(x$workbook, "Total n")))
      stop("Common-cohort counts do not match ", x$workbook$path)
  }
  saved_path <- file.path(results_dir, sprintf("DMFT_Age_BMI_Resampling_%s.rds", registry_key))
  if (file.exists(saved_path)) {
    saved <- readRDS(saved_path)
    if (!identical(saved$registry_key, registry_key) || !identical(saved$source_rows, source_rows) ||
        !identical(as.character(saved$outcome), as.character(d$condL)))
      stop("Saved cohort does not match this preprocessing session.")
  }
  message("Age/BMI cohort: ", nrow(d), "; excluded ", sum(!keep))
  mu <- sd_used <- setNames(numeric(0), character(0))
  if (length(markers)) {
    if (!all(c("Column", "Mean", "SD") %in% names(reference_stats)) || anyDuplicated(reference_stats$Column))
      stop("Reference statistics require unique Column rows and Mean/SD columns.")
    idx <- match(markers, reference_stats$Column)
    if (anyNA(idx)) stop("Missing LD reference statistics for selected biomarkers.")
    mu <- setNames(reference_stats$Mean[idx], markers); sd_used <- setNames(reference_stats$SD[idx], markers)
    if (!is.numeric(mu) || !is.numeric(sd_used) || any(!is.finite(mu)) || any(!is.finite(sd_used)) || any(sd_used <= 0))
      stop("LD means must be finite and SDs positive.")
  }
  fit <- glm(reformulate(predictors, response = "condL"), data = d,
             family = binomial(), na.action = na.fail)
  b <- coef(fit)
  if (!isTRUE(fit$converged) || any(!is.finite(b))) stop("Nonconverged or non-finite fit; no export written.")
  raw_b <- b
  raw_b[markers] <- b[markers] / sd_used
  raw_b["(Intercept)"] <- b["(Intercept)"] - sum(b[markers] * mu / sd_used)
  # Check selected biomarker slopes against the workbook's per-doubling ORs.
  for (j in seq_along(markers)) {
    nm <- paste0("OR_perDoubling_Biomarker", j)
    if (nm %in% names(winner$row)) {
      v <- as.numeric(winner$row[[nm]])
      if (is.finite(v) && v > 0 && abs(raw_b[markers[j]] - log(v)) > 1e-5 * (1 + abs(log(v))))
        stop("Refitted slope differs from the selected workbook for ", markers[j], ". Check data, LD reference statistics and coding.")
    }
  }
  X <- model.matrix(fit)
  for (j in seq_along(markers)) X[, markers[j]] <- X[, markers[j]] * sd_used[j] + mu[j]
  eta_export <- drop(X %*% raw_b); eta_fit <- as.numeric(predict(fit, type = "link"))
  max_error <- max(abs(eta_export - eta_fit))
  if (!is.finite(max_error) || max_error > 1e-8) stop("Raw-input equation failed the fitted-model equivalence check.")
  if (max(abs(b[-1])) > 15) warning("Large coefficients: assess estimation instability before use.")

  terms <- names(b); input <- setNames(rep("value in original units", length(b)), terms)
  input["(Intercept)"] <- ""; input[markers] <- "log2(raw value + 5e-07)"
  means <- sds <- setNames(rep(NA_real_, length(b)), terms); means[markers] <- mu; sds[markers] <- sd_used
  tab <- data.frame(Term = terms, Label = terms, Estimate_model_scale = unname(b),
    Coefficient_raw_input = unname(raw_b), LD_mean_log2 = unname(means), LD_SD_log2 = unname(sds),
    Input = unname(input), stringsAsFactors = FALSE)
  meta <- data.frame(Term = c("(Threshold)", "(Delta)", "(Case group)", "(Reference group)",
                              "(Model configuration)", "(Search design)"),
    Label = c("Classification threshold", "log2 offset", "high DMFT (DMFT > 13)", "low DMFT (DMFT <= 13)",
              "Biomarkers_Age_BMI", "exhaustive12_beam3"),
    Estimate_model_scale = NA_real_, Coefficient_raw_input = c(threshold, delta, NA, NA, NA_real_, NA_real_),
    LD_mean_log2 = NA_real_, LD_SD_log2 = NA_real_, Input = "", stringsAsFactors = FALSE)
  tab <- rbind(tab, meta)
  number <- function(x) sprintf("%.17g", x)
  eqterms <- setNames(terms[-1], terms[-1])
  eqterms[markers] <- paste0("log2(X", seq_along(markers), " + 5e-7)")
  equation <- paste0("eta = ", number(raw_b[1]), paste0(ifelse(raw_b[-1] < 0, " - ", " + "),
    number(abs(raw_b[-1])), " * ", eqterms, collapse = ""))
  lines <- c("Selected DMFT equation: risk of HD (DMFT > 13)",
    paste("Configuration:", winner$configuration), paste("Source workbook:", wb$path),
    paste("Source sheet:", winner$sheet, "| Excel row:", winner$row$.source_row + 1L), "",
    equation, "p_HD = plogis(eta) = 1 / (1 + exp(-eta))",
    paste0("Classify HD when p_HD >= ", number(threshold), "; otherwise LD."), "",
    if (length(markers)) c("Use RAW biomarkers in the original assay units; apply log2(raw + 5e-7) exactly once.",
                           "Do not additionally standardize inputs.", paste0("X", seq_along(markers), " = ", markers)),
    if ("Age_DOS" %in% covariates) "Age_DOS = age in years.",
    if ("BMI" %in% covariates) "BMI = kg/m^2 in original units.",
    "",
    paste0("Common cohort: LD = ", counts[1], "; HD = ", counts[2], "."),
    paste0("Saved cross-validation BA = ", number(winner$row$Balanced_Accuracy), "; AUC = ", number(winner$row$AUC)),
    paste0("Maximum raw-equation versus fitted-model log-odds difference = ", number(max_error)),
    "This full-cohort refit does not supply new validation performance.",
    "Panel, configuration and threshold selection used development predictions; independent validation is needed.")
  selection <- do.call(rbind, lapply(selected$candidates, function(x) data.frame(
    Configuration = x$configuration, Panel_size = x$k, Balanced_Accuracy = x$row$Balanced_Accuracy,
    AUC = x$row$AUC, Threshold = x$row$Threshold, Workbook = x$workbook$path,
    Sheet = x$sheet, Excel_row = x$row$.source_row + 1L)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(digits = 17); on.exit(options(old), add = TRUE)
  for (filename in c("dmft_age_bmi_model.csv", "dmft_age_bmi_coefficients.csv"))
    write.csv(tab, file.path(output_dir, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.csv(selection, file.path(output_dir, "DMFT_panel_selection.csv"), row.names = FALSE)
  writeLines(lines, file.path(output_dir, "DMFT_equation.txt"), useBytes = TRUE)
  bundle <- list(fit = fit, coefficients = tab, markers = markers, covariates = covariates,
    delta = delta, threshold = threshold, reference_mean = mu, reference_sd = sd_used,
    source_rows = source_rows, registry_key = registry_key, group_counts = counts, equation = equation,
    max_link_difference = max_error, selected_result = winner, session_info = sessionInfo())
  saveRDS(bundle, file.path(output_dir, "DMFT_equation_model.rds"))
  cat(paste(lines, collapse = "\n"), "\n\nSaved to: ", normalizePath(output_dir), "\n", sep = "")
  if (!length(markers)) message("Demographics-only model won. Its export is valid; the current app's biomarker-model loader needs adaptation for zero biomarkers.")
  invisible(bundle)
}

if (!isTRUE(getOption("dmft.equation.skip_run"))) {
  if (!exists("data.dfPH") || !exists("registry_key")) stop("Run the matching age/BMI preprocessing script first.")
  .dmft_reference_stats <- if (exists("stats_PH_vif")) stats_PH_vif else if (exists("stats_PH")) stats_PH else NULL
  DMFT_EQUATION <- export_dmft_equation(data.dfPH, .dmft_reference_stats, registry_key,
    output_dir = OUTPUT_DIR, results_file = RESULTS_FILE, results_dir = RESULTS_DIR)
}
