# DMFT BIOMARKER SEARCH ADJUSTED FOR AGE AND BMI
# Run 01_Preprocessing_DMFT_Age_BMI.R first in the same R session.
# All single biomarkers and all pairs are evaluated; beam search starts at size 3.
# Age_DOS and BMI enter every candidate model in original numeric units.
# Separate registry keys prevent accidental reuse of previous results.

# ==========================================================================
# LOGISTIC REGRESSION — FAST COMBINATION SEARCH
# ==========================================================================
# Logistic regression on every combination of min_biomarkers..max_biomarkers
# candidate markers (or a beam search where there are too many), with
# repeated cross-validation, for Script 1's caries, DMFT and perio projects.
# Run Script 1 first, in the same session:
#
#   source("script1_preprocessing.R")   # PROJECT = "caries", "DMFT" or "perio"
#   source("LogReg_Fast.R")             # this file; runs at the bottom
#
#   project  comparison: case vs reference
#   caries   CP vs CF (Cavitated_Caries) or HC vs LC (Clinical_Caries)
#   DMFT     HD vs LD
#   perio    PH: P vs H | GH: G vs H | GP: P vs G   (RUN_PH / RUN_GH / RUN_GP)
#
# Sensitivity always refers to the case group. One workbook per comparison,
# named with Script 1's registry key, e.g.
#   Results_LogReg_PH_PISA_noConfounders_panel_ratios_noVif.xlsx
# The validation project is run in Script 2.
#
# How it works
#   Pass 1, sweep: every combination is cross-validated on the first
#     SWEEP_REPEATS repeats with a bare IRLS on a prebuilt design matrix, and
#     ranked by balanced accuracy. With SWEEP_CV = "onestep" each combination
#     is fitted once on all subjects and each fold's fit is approximated by one
#     Newton step; "exact" refits every fold.
#   Pass 2, exact: the best SWEEP_KEEP per panel size are refitted through
#     caret on all CV_REPEATS repeats. Every exported number comes from pass 2.
#   Both passes use the same folds: caret's own folds after set.seed(123).
#
# Metrics (pass 2), from the out-of-fold predictions of all repeats pooled
#   Youden cut: the cut-off with the highest Sensitivity + Specificity; when
#     several tie, the one closest to 0.5. Accuracy, Balanced_Accuracy,
#     Sensitivity, Specificity, Precision and NPV are at that cut. AUC needs
#     no cut.
#   Balanced_Accuracy 95% CI: the 2.5th and 97.5th percentiles of the
#     balanced accuracies of the CV_REPEATS repeats (each subject predicted
#     once per repeat), at the same cut.
#   Balanced_Accuracy_p_Fisher: Fisher's exact test of true by predicted class
#     at that cut, first repeat (one prediction per subject).
# ==========================================================================


# ==========================================================================
# 1) SETTINGS
# ==========================================================================

# Perio comparisons (caries and DMFT have one comparison each)
RUN_PH <- TRUE
RUN_GH <- TRUE
RUN_GP <- TRUE

min_biomarkers <- 1          # smallest panel exported
max_biomarkers <- 6          # largest panel exported
numCores       <- 10

# Cross-validation for every exported number: CV_FOLDS-fold CV repeated
# CV_REPEATS times. The CI is built from the repeats, so use many.
CV_FOLDS   <- 5
CV_REPEATS <- 100
# Repeats used to rank combinations in the sweep. More gives a steadier
# ranking and is proportionally slower; exported numbers always use CV_REPEATS.
SWEEP_REPEATS <- 5

# Export: at each panel size, the top TOP_N per Balanced_Accuracy, Sensitivity
# and AUC among the refitted models (SWEEP_KEEP, RANK_COLS), sorted by
# Balanced_Accuracy (ties: AUC)
TOP_N <- 50L

# Odds-ratio CIs: profile likelihood (TRUE) or Wald (FALSE). Models flagged for
# separation always get Wald.
USE_PROFILE_CI <- TRUE
# A slope above this (on the standardised scale), or no convergence, flags
# complete or quasi-complete separation; those rows' ORs are not interpretable.
SEPARATION_COEF_LIMIT <- 15

# Checks before the run
RUN_VALIDATION    <- TRUE    # sweep vs caret, and one-step vs exact sweep; results printed
RUN_PREFLIGHT     <- TRUE    # time a sample of models and project the sweep's runtime
MAX_RUNTIME_HOURS <- 12      # stop when the projected sweep is longer (needs RUN_PREFLIGHT); Inf = no limit

# Sweep
# Engine: "fastglm" (C++, if installed), "irls" (base R), "glmfit" (stock glm.fit)
SWEEP_ENGINE <- if (requireNamespace("fastglm", quietly = TRUE)) "fastglm" else "irls"
# Sweep CV: "onestep" fits each combination once on all subjects and gets each
# fold's coefficients from one Newton step that removes that fold's subjects,
# about 3-5x faster than "exact", which refits every fold. It only affects
# which combinations reach pass 2: exported numbers always come from exact
# refits. A combination whose full fit is near separation (not converged, or a
# fitted probability within ~5e-5 of 0 or 1, i.e. |linear predictor| >
# ONESTEP_ETA_MAX) is refitted fold by fold instead.
SWEEP_CV        <- "onestep"
ONESTEP_ETA_MAX <- 10
SWEEP_KEEP   <- TOP_N        # best combinations per ranking metric refitted in pass 2
# Ranking columns of the sweep: 1 Accuracy, 2 Sensitivity, 3 Specificity,
# 4 Precision, 5 NPV, 6 Threshold, 7 AUC, 10 Balanced_Accuracy
RANK_COLS    <- 10L
# Search: "exhaustive" = every combination; "beam" = level k built from the best
# level k-1 models; "auto" = exhaustive while C(p, k) <= EXHAUSTIVE_LIMIT, then
# beam. With all ~419 biomarkers, C(p, 6) is ~7 trillion, so use beam or auto.
SEARCH_MODE      <- "beam"
EXHAUSTIVE_LIMIT <- 2e6
BEAM_WIDTH       <- TOP_N    # parents kept per ranking metric

# This file's name for the parameters sheet; found automatically when the file
# is run with source() or Rscript, otherwise this value is used.
SCRIPT_FILE <- "LogReg_Fast.R"


# ==========================================================================
# 2) PACKAGES AND PREREQUISITES
# ==========================================================================
# Run this file with source(): it stops at the first unmet requirement. Line by
# line (Ctrl+Enter), RStudio carries on after an error, and the message left on
# screen is then the last failure, not the first.

.need_pkgs <- c("caret", "parallel", "doSNOW", "foreach", "writexl", "dplyr", "pROC")
.miss_pkgs <- .need_pkgs[!vapply(.need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkgs) > 0)
  stop("Missing package(s): ", paste(.miss_pkgs, collapse = ", "),
       "\n  install.packages(c(", paste0('"', .miss_pkgs, '"', collapse = ", "), "))", call. = FALSE)
suppressMessages({
  library(caret); library(parallel); library(doSNOW); library(foreach)
  library(writexl); library(dplyr); library(pROC)
})

script_file <- local({
  f <- NULL
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { f <- fr$ofile; break }   # source()
  if (is.null(f)) {                                                                   # Rscript
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) f <- gsub("~+~", " ", sub("^--file=", "", a[1]), fixed = TRUE)
  }
  if (is.null(f) || !nzchar(f)) SCRIPT_FILE else basename(f)
})

PERIO_DEFS <- list(
  PH = list(ref = "H", case = "P", label = "Periodontitis vs Healthy"),
  GH = list(ref = "H", case = "G", label = "Gingivitis vs Healthy"),
  GP = list(ref = "G", case = "P", label = "Periodontitis vs Gingivitis")
)
PERIO_RUN <- names(PERIO_DEFS)[c(isTRUE(RUN_PH), isTRUE(RUN_GH), isTRUE(RUN_GP))]

.msgs <- character(0)
if (!exists("PROJECT") || !exists("registry_key")) {
  .msgs <- "Script 1 has not run in this session: source(\"script1_preprocessing.R\") first."
} else if (identical(PROJECT, "validation")) {
  .msgs <- "PROJECT = \"validation\": the validation project is run in Script 2, not here."
} else if (!PROJECT %in% c("caries", "DMFT", "perio")) {
  .msgs <- sprintf("PROJECT = \"%s\"; this script runs \"caries\", \"DMFT\" and \"perio\".", PROJECT)
} else {
  .cmps <- if (PROJECT == "perio") PERIO_RUN else "PH"
  if (length(.cmps) == 0)
    .msgs <- c(.msgs, "RUN_PH, RUN_GH and RUN_GP are all FALSE: nothing to run.")
  .need <- c("confounders", "biomarker_registry",
             unlist(lapply(.cmps, function(x) paste0(c("data.df", "biomarkers_", "stats_"), x,
                                                     c("", "_vif", "_vif")))))
  .miss <- .need[!vapply(.need, exists, logical(1))]
  if (length(.miss) > 0) {
    .msgs <- c(.msgs, paste("Script 1 object(s) missing:", paste(.miss, collapse = ", "),
                            "- run Script 1 to the end in this session."))
  } else {
    for (.x in .cmps) {
      .d <- get(paste0("data.df", .x))
      if (is.null(.d) || nrow(.d) == 0)
        .msgs <- c(.msgs, sprintf("data.df%s has no rows; check Script 1's settings for %s.", .x, PROJECT))
      if (length(get(paste0("biomarkers_", .x, "_vif"))) == 0)
        .msgs <- c(.msgs, sprintf("biomarkers_%s_vif is empty; check Script 1's candidate markers.", .x))
    }
  }
}
if (!isTRUE(SWEEP_CV %in% c("onestep", "exact")))
  .msgs <- c(.msgs, sprintf("SWEEP_CV is \"%s\"; use \"onestep\" or \"exact\".", paste(SWEEP_CV, collapse = ", ")))
if (length(.msgs) > 0)
  stop("\n\n===== PREREQUISITES NOT MET =====\n  ", paste(.msgs, collapse = "\n  "), "\n", call. = FALSE)


# ==========================================================================
# 3) COMPARISONS FROM SCRIPT 1
# ==========================================================================
# One spec per comparison: Script 1's data, candidate markers and z-score stats,
# the reference and case labels, and the rows describing the outcome on the
# parameters sheet. condL becomes a factor with the reference first: caret
# models P(second level), so the case is always the second level.

make_spec <- function(tag, data, biomarkers, stats, ref, case, info, selected_size = NA_integer_) {
  found <- unique(as.character(data$condL))
  bad   <- setdiff(found, c(ref, case))
  if (length(bad) > 0)
    stop(sprintf("%s: condL contains %s; expected only %s and %s. Re-run Script 1.",
                 tag, paste(bad, collapse = ", "), ref, case), call. = FALSE)
  if (!all(c(ref, case) %in% found))
    stop(sprintf("%s: condL needs both %s and %s; found only %s.",
                 tag, ref, case, paste(found, collapse = ", ")), call. = FALSE)
  data$condL <- factor(as.character(data$condL), levels = c(ref, case))
  info$reference <- ref
  info$case      <- case
  list(tag = tag, data = data, biomarkers = biomarkers, stats = stats, info = info,
       selected_size = selected_size)
}

# Size of Script 1's selected-marker list (for the runtime hint only)
.selected_size <- function(markers) if (is.null(markers)) NA_integer_ else length(markers)
.sel <- get0("selected_markers", ifnotfound = list())

if (PROJECT == "caries") {

  clinical <- CARIES_SOURCE_COL == "Clinical_Caries"
  ref  <- if (clinical) "LC" else "CF"
  case <- if (clinical) "HC" else "CP"
  rule <- if (clinical) {
    sprintf("<= %g -> LC (low caries burden) | > %g -> HC (high caries burden)",
            CARIES_BINARY_CUTOFF, CARIES_BINARY_CUTOFF)
  } else {
    sprintf("<= %g -> CF (no caries) | >= %g -> CP (high caries experience) | in between dropped",
            CARIES_NOCARIES_MAX, CARIES_HIGH_MIN)
  }
  info <- list(project = PROJECT, condl_source = CONDL_SOURCE, rows = c(
    "Caries source column" = CARIES_SOURCE_COL,
    "Caries cut-off"       = rule,
    "Periodontal filter"   = if (is.null(CARIES_PERIODONTAL_FILTER)) "none"
                             else paste(CARIES_PERIODONTAL_FILTER, collapse = ", ")))
  SPECS <- list(make_spec(paste0(case, "_vs_", ref), data.dfPH, biomarkers_PH_vif, stats_PH_vif,
                          ref, case, info, .selected_size(.sel$Caries[[CARIES_SOURCE_COL]])))
  # for the ROC and heatmap scripts, which look for data.dfCaries
  data.dfCaries <- SPECS[[1]]$data

} else if (PROJECT == "DMFT") {

  info <- list(project = PROJECT, condl_source = CONDL_SOURCE, rows = c(
    "DMFT"               = paste(dmft_cols, collapse = " + "),
    "DMFT cut-off"       = sprintf("<= %g -> LD (low DMFT) | > %g -> HD (high DMFT) (%s)",
                                   DMFT_CUTOFF_USED, DMFT_CUTOFF_USED,
                                   if (is.null(DMFT_CUTOFF)) "median of the analysed subjects" else "fixed"),
    "Periodontal filter" = if (is.null(DMFT_PERIODONTAL_FILTER)) "none"
                           else paste(DMFT_PERIODONTAL_FILTER, collapse = ", ")))
  SPECS <- list(make_spec("HD_vs_LD", data.dfPH, biomarkers_PH_vif, stats_PH_vif,
                          "LD", "HD", info, .selected_size(.sel$DMFT$PH)))

} else {   # perio

  .score_rule <- function(score, lo, hi)
    sprintf("%s < %g -> Healthy | %g to %g -> Gingivitis | > %g -> Periodontitis", score, lo, lo, hi, hi)
  rule <- switch(PERIO_DEFINITION,
                 original  = "Original clinical diagnosis",
                 PISA      = .score_rule("PISA", PISA_HEALTHY_MAX, PISA_MOD_MAX),
                 AdjustPRA = .score_rule("AdjustPRA", ADJUSTPRA_HEALTHY_MAX, ADJUSTPRA_MOD_MAX))
  SPECS <- lapply(PERIO_RUN, function(x) {
    d <- PERIO_DEFS[[x]]
    make_spec(x, get(paste0("data.df", x)), get(paste0("biomarkers_", x, "_vif")),
              get(paste0("stats_", x, "_vif")), d$ref, d$case,
              list(project = PROJECT, condl_source = CONDL_SOURCE,
                   rows = c("condL definition" = rule, "Contrast" = d$label)),
              .selected_size(.sel[[PERIO_DEFINITION]][[x]]))
  })
  names(SPECS) <- PERIO_RUN
}

cat(sprintf("\n===== Logistic regression: %s (condL from %s) =====\n", PROJECT, CONDL_SOURCE))
for (.sp in SPECS) {
  .tab <- table(.sp$data$condL)
  cat(sprintf("  %s: %s vs %s (reference) | n = %d (%s)\n", .sp$tag, .sp$info$case, .sp$info$reference,
              sum(.tab), paste(names(.tab), .tab, sep = " ", collapse = ", ")))
}


# ==========================================================================
# 4) CHECKS AGAINST SCRIPT 1
# ==========================================================================
cat("\n===== Checks =====\n")

# The models use `confounders` from this session. Stop if it no longer matches
# what Script 1 ran with, since the file names (Script 1's registry key) would
# then describe the wrong models.
.conf_script1 <- as.character(biomarker_registry[[registry_key]]$confounders)
if (!setequal(.conf_script1, as.character(confounders)))
  stop(sprintf("confounders in this session (%s) differ from the ones Script 1 ran with (%s). Re-run Script 1.",
               if (length(confounders)) paste(confounders, collapse = ", ") else "none",
               if (length(.conf_script1)) paste(.conf_script1, collapse = ", ") else "none"), call. = FALSE)
cat(sprintf("  confounders: %s\n", if (length(confounders)) paste(confounders, collapse = ", ") else "none"))
cat(sprintf("  CV: %d folds x %d repeats (sweep ranks on the first %d)\n",
            CV_FOLDS, CV_REPEATS, min(SWEEP_REPEATS, CV_REPEATS)))
if (CV_REPEATS < 20)
  cat("  NOTE: the CI is built from the repeats; with fewer than 20 its percentiles are rough.\n")

# Script 1 imputes the biomarkers but not the confounders; rows with a missing
# confounder are left out once, before the folds are built.
.max_combos <- 0
for (.sp in SPECS) {
  .bm    <- intersect(.sp$biomarkers, colnames(.sp$data))
  .conf  <- intersect(confounders, colnames(.sp$data))
  .n_cc  <- if (length(.conf)) sum(!complete.cases(as.data.frame(.sp$data)[, .conf, drop = FALSE])) else 0L
  .n_cmb <- sum(choose(length(.bm), min_biomarkers:max_biomarkers))
  .max_combos <- max(.max_combos, .n_cmb)
  cat(sprintf("  %s: %d candidate markers | combinations (k = %d..%d): %s\n", .sp$tag, length(.bm),
              min_biomarkers, max_biomarkers, format(.n_cmb, big.mark = ",", scientific = FALSE)))
  if (anyNA(.sp$data[, .bm, drop = FALSE]))
    warning(sprintf("%s: NAs in the candidate markers despite Script 1's imputation.", .sp$tag))
  if (.n_cc > 0)
    cat(sprintf("    %d row(s) with a missing confounder will be left out.\n", .n_cc))
}
if (.max_combos > 1e8)
  cat("  Exhaustive search is infeasible at this size; SEARCH_MODE \"beam\" or \"auto\" is needed.\n")


# ==========================================================================
# 5) SHARED HELPERS
# ==========================================================================

format_duration <- function(seconds) {
  hours   <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs    <- round(seconds %% 60)
  if (hours > 0)   return(sprintf("%dh %dm %ds", hours, minutes, secs))
  if (minutes > 0) return(sprintf("%dm %ds", minutes, secs))
  sprintf("%ds", secs)
}

# caret's per-resample summary. twoClassSummary treats the first level as
# positive, so its Sens and Spec are swapped to make the case (second level)
# positive. glm has nothing to tune, so caret selects nothing with it.
customSummary <- function(data, lev = NULL, model = NULL) {
  ts  <- twoClassSummary(data, lev, model)
  ds  <- defaultSummary(data, lev, model)
  out <- c(ts[["ROC"]], ts[["Spec"]], ts[["Sens"]], ds[["Accuracy"]], ds[["Kappa"]],
           mean(c(ts[["Spec"]], ts[["Sens"]])))
  names(out) <- c("ROC", "Sens", "Spec", "Accuracy", "Kappa", "BalancedAccuracy")
  out
}

control <- trainControl(
  method          = "repeatedcv",
  number          = CV_FOLDS,
  repeats         = CV_REPEATS,
  summaryFunction = customSummary,
  classProbs      = TRUE,
  allowParallel   = FALSE,
  savePredictions = "final"
)

# Ratio features (Script 1's bacterial ratios) are log2 ratios, not z-scored
is_ratio_feature <- function(name) startsWith(name, "ratio__")

# Share as a percentage, "<0.1%" for small non-zero shares
.pct <- function(x) ifelse(x > 0 & x < 0.001, "<0.1%", sprintf("%.1f%%", 100 * x))

# Highest Balanced_Accuracy first, ties by AUC, NAs last
sort_by_balanced_accuracy <- function(df) {
  out <- df[order(-df$Balanced_Accuracy, -df$AUC, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Rows in the top n for at least one metric, sorted by Balanced_Accuracy
get_top_n_per_metric <- function(df, n = TOP_N,
                                 metrics = c("Balanced_Accuracy", "Sensitivity", "AUC")) {
  top_rows <- unique(unlist(lapply(metrics, function(m) {
    if (m %in% colnames(df)) order(df[[m]], decreasing = TRUE, na.last = TRUE)[seq_len(min(n, nrow(df)))]
  })))
  sort_by_balanced_accuracy(df[top_rows, , drop = FALSE])
}


# ==========================================================================
# 6) PASS 2: EXACT REFIT THROUGH CARET
# ==========================================================================

# Youden cut from the pooled out-of-fold predictions; among tied cut-offs the
# one closest to 0.5 (equally close: the higher). Ties are found with a
# tolerance: distinct values of Sensitivity + Specificity differ by at least
# 1 / (n_case x n_control), but equal ones can differ in the last bit, which
# pROC's "best" (exact equality) would not count as a tie.
youden_best <- function(roc_obj) {
  m <- coords(roc_obj, "all",
              ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv"),
              transpose = FALSE)
  J    <- m$sensitivity + m$specificity
  tol  <- 0.5 / (length(roc_obj$cases) * length(roc_obj$controls))
  tied <- which(J >= max(J, na.rm = TRUE) - tol)
  m[tied[order(abs(m$threshold[tied] - 0.5), -m$threshold[tied])[1]], , drop = FALSE]
}

# Balanced accuracy of each repeat at `cut` (case if P(case) >= cut), the
# 2.5th / 97.5th percentiles across repeats, and Fisher's exact test on the
# first repeat
repeat_ci <- function(pred, lev, cut) {
  rep_id <- sub("^.*\\.(Rep[0-9]+)$", "\\1", pred$Resample)
  ba_r <- vapply(split(pred, rep_id), function(d) {
    p <- d[[lev[2]]]
    (mean(p[d$obs == lev[2]] >= cut) + mean(p[d$obs == lev[1]] < cut)) / 2
  }, numeric(1))
  ci <- if (length(ba_r) >= 2) quantile(ba_r, c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
  first <- pred[rep_id == sort(unique(rep_id))[1], , drop = FALSE]
  p_fis <- if (!is.finite(cut)) NA_real_ else
    tryCatch(fisher.test(table(factor(first$obs, levels = lev),
                               factor(first[[lev[2]]] >= cut, levels = c(FALSE, TRUE))))$p.value,
             error = function(e) NA_real_)
  c(ci, p_fis)
}

# Fits every combination in chunk_combinations through caret on ctrl's folds.
# Returns, per combination, the metrics, the CI and Fisher p, the odds ratios
# and the separation diagnostics.
run_foreach_chunk <- function(chunk_combinations, cl, data_df, stats_df, confounders, ctrl = control) {
  foreach(
    biomarker_set  = chunk_combinations,
    .packages      = c("caret", "pROC"),
    .export        = c("customSummary", "is_ratio_feature", "youden_best", "repeat_ci",
                       "USE_PROFILE_CI", "SEPARATION_COEF_LIMIT"),
    .errorhandling = "pass"
  ) %dopar% {

    k <- length(biomarker_set)
    tryCatch({
      predictors  <- c(biomarker_set, confounders)
      subset_data <- data_df[, c(predictors, "condL")]
      set.seed(123)
      logit_model <- train(as.formula(paste("condL ~", paste(predictors, collapse = " + "))),
                           data = subset_data, method = "glm", family = "binomial",
                           metric = "BalancedAccuracy", trControl = ctrl)

      lev  <- levels(subset_data$condL)
      pred <- logit_model$pred

      if (length(unique(pred$obs)) < 2 || !all(lev %in% colnames(pred))) {
        results    <- rep(NA_real_, 10)
        or_results <- rep(NA_real_, k * 6)
        diagnostics <- rep(NA_real_, 3)
      } else {
        # direction fixed: a higher P(case) means case
        roc_obj <- suppressMessages(roc(pred$obs, pred[[lev[2]]], levels = lev,
                                        direction = "<", quiet = TRUE))
        m <- youden_best(roc_obj)
        results <- c(m$accuracy, m$sensitivity, m$specificity, m$precision, m$npv,
                     m$threshold, as.numeric(roc_obj$auc),
                     repeat_ci(pred, lev, m$threshold))

        final_model <- logit_model$finalModel
        coefs       <- coef(final_model)
        slopes      <- coefs[names(coefs) != "(Intercept)"]
        max_abs     <- if (length(slopes)) max(abs(slopes), na.rm = TRUE) else NA_real_
        converged   <- isTRUE(final_model$converged)
        separated   <- (!is.na(max_abs) && max_abs > SEPARATION_COEF_LIMIT) || !converged
        conf_ints   <- suppressMessages(tryCatch(
          if (USE_PROFILE_CI && !separated) confint(final_model, level = 0.95)
          else confint.default(final_model, level = 0.95),
          error = function(e) confint.default(final_model, level = 0.95)))

        # The model is fitted on Z = (log2 X - Mean) / SD, so the OR per
        # doubling of the raw value is exp(beta / SD) and the OR per SD of log2
        # abundance is exp(beta). Ratio features are already log2 and not
        # z-scored (SD = 1).
        or_results <- rep(NA_real_, k * 6)
        for (i in seq_len(k)) {
          bm <- biomarker_set[i]
          if (!bm %in% names(coefs)) next
          b <- coefs[[bm]]; lo <- conf_ints[bm, 1]; hi <- conf_ints[bm, 2]
          base <- (i - 1) * 6
          or_results[base + 4:6] <- exp(c(b, lo, hi))
          sd_bm <- if (is_ratio_feature(bm)) 1 else {
            s <- stats_df$SD[stats_df$Column == bm]
            if (length(s) == 1 && !is.na(s) && s > 0) s else NA_real_
          }
          if (!is.na(sd_bm)) or_results[base + 1:3] <- exp(c(b, lo, hi) / sd_bm)
        }
        diagnostics <- c(max_abs, separated, converged)
      }
      list(biomarkers = biomarker_set, metrics = as.numeric(results),
           or_metrics = as.numeric(or_results), diagnostics = as.numeric(diagnostics), success = TRUE)

    }, error = function(e) {
      list(biomarkers = biomarker_set, metrics = rep(NA_real_, 10),
           or_metrics = rep(NA_real_, k * 6), diagnostics = rep(NA_real_, 3),
           success = FALSE, error_msg = e$message)
    })
  }
}

assemble_results <- function(results_list, num_biomarkers) {
  biomarker_cols <- paste0("Biomarker", seq_len(num_biomarkers))
  metric_cols    <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "Threshold", "AUC",
                      "Balanced_Accuracy_CI_Lower", "Balanced_Accuracy_CI_Upper", "Balanced_Accuracy_p_Fisher")
  or_cols <- unlist(lapply(seq_len(num_biomarkers), function(i)
    paste0(c("OR_perDoubling_", "OR_perDoubling_Lower_CI_", "OR_perDoubling_Upper_CI_",
             "OR_perSD_", "OR_perSD_Lower_CI_", "OR_perSD_Upper_CI_"), "Biomarker", i)))
  diag_cols <- c("Max_Abs_Coef", "Separation_Flag", "Converged")
  all_cols  <- c(biomarker_cols, metric_cols, or_cols, diag_cols)

  assembled <- do.call(rbind, lapply(results_list, function(r) {
    if (is.null(r) || inherits(r, "error")) return(matrix(NA, 1, length(all_cols)))
    matrix(c(r$biomarkers, r$metrics, r$or_metrics, r$diagnostics), nrow = 1)
  }))
  df <- as.data.frame(assembled, stringsAsFactors = FALSE)
  colnames(df) <- all_cols
  for (col in c(metric_cols, or_cols, diag_cols)) df[[col]] <- as.numeric(df[[col]])
  df$Balanced_Accuracy <- (df$Sensitivity + df$Specificity) / 2
  df[, c(biomarker_cols, "Balanced_Accuracy", "Balanced_Accuracy_CI_Lower", "Balanced_Accuracy_CI_Upper",
         "Balanced_Accuracy_p_Fisher", "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV",
         "Threshold", "AUC", or_cols, diag_cols)]
}

build_params_sheet <- function(spec, biomarkers, confounders, start_time, end_time,
                               group_counts, n_dropped = 0L) {
  info <- spec$info
  rows <- list(
    c("=== OUTCOME ===", ""),
    c("Project", info$project),
    c("condL source", info$condl_source))
  for (nm in names(info$rows)) rows <- c(rows, list(c(nm, unname(info$rows[[nm]]))))
  rows <- c(rows, list(
    c("Comparison", spec$tag),
    c("Case group (positive)", info$case),
    c("Reference group", info$reference),
    c("Script 1 registry key", registry_key),
    c("Script file", script_file),
    c("", ""),
    c("=== SAMPLE ===", ""),
    c("Total n", as.character(sum(unlist(group_counts))))))
  for (g in names(group_counts)) rows <- c(rows, list(c(paste0("  n (", g, ")"), as.character(group_counts[[g]]))))
  if (n_dropped > 0) rows <- c(rows, list(c("Rows left out (missing confounder)", as.character(n_dropped))))
  rows <- c(rows, list(
    c("", ""),
    c("=== MODEL ===", ""),
    c("Model", "Multivariable logistic regression (glm, binomial, caret)"),
    c("Cross-validation", sprintf("%d-fold CV x %d repeats, seed 123 (caret's own folds)", CV_FOLDS, CV_REPEATS)),
    c("Sweep", sprintf("ranked on the first %d repeat(s) of the same folds; best %d per panel size refitted",
                       min(SWEEP_REPEATS, CV_REPEATS), SWEEP_KEEP)),
    c("Sweep CV", if (SWEEP_CV == "onestep") {
      sprintf(paste("one-step: each combination fitted once on all subjects; each fold's coefficients from one",
                    "Newton step removing that fold's subjects. Combinations whose full fit did not converge or",
                    "had |linear predictor| > %g refitted fold by fold. Affects only which combinations are",
                    "refitted; every exported number comes from the exact refits"), ONESTEP_ETA_MAX)
    } else "exact: every fold refitted"),
    c("Candidate markers", as.character(length(biomarkers))),
    c("Markers per model", sprintf("%d to %d", min_biomarkers, max_biomarkers)),
    c("Confounders", if (length(confounders)) paste(confounders, collapse = ", ") else "none"),
    c("", ""),
    c("=== METRICS ===", ""),
    c("Youden cut", paste("Cut-off with the highest Sensitivity + Specificity on the pooled out-of-fold",
                          "predictions of all repeats; among tied cut-offs, the one closest to 0.5")),
    c("Accuracy ... NPV", "At the Youden cut; Sensitivity and Precision refer to the case group"),
    c("Balanced_Accuracy", "(Sensitivity + Specificity) / 2 at the Youden cut; sheets sorted by it (ties: AUC)"),
    c("Balanced_Accuracy 95% CI", sprintf(paste("2.5th and 97.5th percentiles of the balanced accuracies of the %d repeats",
                                                "(each subject predicted once per repeat) at the Youden cut; they reflect",
                                                "the variation from the random split into folds"), CV_REPEATS)),
    c("Balanced_Accuracy_p_Fisher", paste("Fisher's exact test of true class x predicted class at the Youden cut, first repeat",
                                          "(one prediction per subject); leans optimistic, as the cut is chosen on the same predictions")),
    c("AUC", "Pooled out-of-fold predictions of all repeats; direction fixed (higher P(case) = case)"),
    c("OR_perDoubling", paste("Odds ratio per doubling of the raw value, adjusted for the rest of the model: exp(beta / SD),",
                              "SD of log2 abundance in the reference group. Ratio features are log2 ratios (SD = 1)")),
    c("OR_perSD", "Odds ratio per 1 SD of log2 abundance (reference group's SD)"),
    c("OR CI", if (USE_PROFILE_CI) "95%, profile likelihood (Wald on error or separation)" else "95%, Wald"),
    c("Separation_Flag", sprintf("1 = max |slope| > %g or glm did not converge; ORs on flagged rows are not interpretable",
                                 SEPARATION_COEF_LIMIT)),
    c("", ""),
    c("=== TIMING ===", ""),
    c("Start", as.character(start_time)),
    c("End", as.character(end_time)),
    c("Duration", format_duration(as.numeric(difftime(end_time, start_time, units = "secs"))))))
  do.call(rbind, lapply(rows, function(r) data.frame(Parameter = r[1], Value = r[2], stringsAsFactors = FALSE)))
}


# ==========================================================================
# 7) PASS 1: THE SWEEP
# ==========================================================================

# --- Folds ----------------------------------------------------------------
# caret draws a private seed inside train() before building its resamples, so
# set.seed(123) + createMultiFolds() does not reproduce them. The folds are
# read off a throwaway train() instead and passed to both passes. Folds depend
# only on the outcome and the seed, not on the predictors.
caret_fold_index <- function(y, ctrl = control, seed = 123) {
  d <- data.frame(condL = y, z = sin(seq_along(y)))
  set.seed(seed)
  fit <- suppressWarnings(caret::train(condL ~ z, data = d, method = "glm", family = "binomial",
                                       metric = "BalancedAccuracy", trControl = ctrl))
  fit$control$index
}

# --- Confounders -----------------------------------------------------------
resolve_confounders <- function(data_df, requested = confounders) {
  confounders <- requested
  conf <- intersect(confounders, colnames(data_df))
  if (length(conf) < length(confounders))
    warning("Confounders not in the data, dropped: ", paste(setdiff(confounders, conf), collapse = ", "))
  if (length(conf) == 0) NULL else conf
}

# A row with a missing confounder would be dropped by glm from every model, so
# it is dropped once here and the folds are built on the rows actually used.
drop_incomplete_confounders <- function(data_df, conf) {
  if (is.null(conf)) return(list(data = data_df, n_dropped = 0L))
  ok <- stats::complete.cases(as.data.frame(data_df)[, conf, drop = FALSE])
  if (!all(ok)) {
    has_na <- conf[vapply(conf, function(v) anyNA(data_df[[v]]), logical(1))]
    cat(sprintf("  Leaving out %d of %d row(s) with a missing confounder (%s).\n",
                sum(!ok), length(ok), paste(has_na, collapse = ", ")))
  }
  data_df <- data_df[ok, , drop = FALSE]
  for (v in conf) if (is.factor(data_df[[v]])) data_df[[v]] <- droplevels(data_df[[v]])
  list(data = data_df, n_dropped = sum(!ok))
}

# --- Design ------------------------------------------------------------------
# folds_idx: training rows per fold, as trainControl(index) takes them.
# id identifies the design, so workers rebuild their copy when it changes.
build_design <- function(data_df, biomarkers, conf = NULL, folds_idx, id = "") {
  data_df$condL <- factor(data_df$condL)
  levels(data_df$condL) <- make.names(levels(data_df$condL))
  lev <- levels(data_df$condL)
  stopifnot(length(lev) == 2L)

  X <- as.matrix(data_df[, biomarkers, drop = FALSE])
  storage.mode(X) <- "double"
  dimnames(X) <- NULL
  if (anyNA(X))
    warning("NAs in the marker matrix: caret would drop those rows per combination, the sweep cannot.")

  Z <- if (length(conf)) {
    mm <- model.matrix(~ ., data = data_df[, conf, drop = FALSE])[, -1, drop = FALSE]
    dimnames(mm) <- NULL
    mm
  } else NULL

  y <- as.numeric(data_df$condL == lev[2L])   # caret models P(second level)
  n <- nrow(X)
  folds <- lapply(folds_idx, function(tr) { tr <- as.integer(tr); list(tr = tr, te = setdiff(seq_len(n), tr)) })
  list(X = X, Z = Z, y = y, lev = lev, folds = folds, n = n, p = length(biomarkers),
       bm_names = biomarkers, id = id)
}

# --- Worker-side precomputation (once per worker, not per combination) -----
# Rows are subset per fold once, so each combination only gathers columns.
prepare_worker <- function() {
  n    <- .DSN$n
  base <- if (is.null(.DSN$Z)) matrix(1, n, 1L) else cbind(1, .DSN$Z)
  Xall <- cbind(base, .DSN$X)
  trs  <- lapply(.DSN$folds, `[[`, "tr")
  tes  <- lapply(.DSN$folds, `[[`, "te")
  pooled_y <- unlist(lapply(tes, function(i) .DSN$y[i]), use.names = FALSE)
  W <- list(
    id = .DSN$id, nb = ncol(base),
    Xall = Xall, y = .DSN$y, tes = tes,          # for the one-step approximation
    Xtr = lapply(trs, function(i) Xall[i, , drop = FALSE]),
    Xte = lapply(tes, function(i) Xall[i, , drop = FALSE]),
    ytr = lapply(trs, function(i) .DSN$y[i]),
    nte = vapply(tes, length, integer(1)),
    nf  = length(trs),
    pooled_y = pooled_y, np = sum(pooled_y), nn = length(pooled_y) - sum(pooled_y))
  assign("W", W, envir = globalenv())
  invisible(TRUE)
}

# --- Fit engines -----------------------------------------------------------
# A fit an engine cannot finish cleanly (warm start overshooting, no
# convergence, separation) is redone with glm.fit from glm's own start, the fit
# caret makes.
.fit_glmfit <- function(X, y, start = NULL) {
  f <- tryCatch(suppressWarnings(glm.fit(X, y, family = binomial(), start = start,
                                         control = list(epsilon = 1e-8, maxit = 25L))),
                error = function(e) NULL)
  if (!is.null(start) &&
      (is.null(f) || !isTRUE(f$converged) || anyNA(f$coefficients) || max(abs(X %*% f$coefficients)) > 30))
    return(.fit_glmfit(X, y))
  if (is.null(f) || anyNA(f$coefficients)) return(NULL)
  list(coef = f$coefficients, conv = isTRUE(f$converged))
}

.fit_fastglm <- function(X, y, start = NULL) {
  f <- tryCatch(fastglm::fastglm(X, y, family = binomial(), start = start, method = 2L),
                error = function(e) NULL)
  if (is.null(f) || anyNA(f$coefficients) || !(is.null(f$converged) || isTRUE(f$converged)) ||
      max(abs(X %*% f$coefficients)) > 30)
    return(.fit_glmfit(X, y))
  list(coef = f$coefficients, conv = TRUE)
}

# IRLS with normal equations and a Cholesky solve, cheaper than glm.fit's QR for
# a handful of columns. |eta| > 30 (probability within 1e-13 of 0 or 1) hands
# over to glm.fit.
.fit_irls <- function(X, y, start = NULL, maxit = 30L, eps = 1e-8) {
  p <- ncol(X)
  b <- if (is.null(start) || length(start) != p) numeric(p) else start
  dev_old <- Inf
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% b)
    if (max(abs(eta)) > 30) return(.fit_glmfit(X, y))
    mu  <- 1 / (1 + exp(-eta))
    dev <- -2 * sum(y * log(mu) + (1 - y) * log1p(-mu))
    if (!is.finite(dev)) return(.fit_glmfit(X, y))
    if (abs(dev - dev_old) < eps * (abs(dev) + 0.1)) return(list(coef = b, conv = TRUE))
    dev_old <- dev
    w <- mu * (1 - mu)
    w[w < 1e-10] <- 1e-10
    A  <- crossprod(X, X * w)
    bb <- crossprod(X, w * eta + (y - mu))
    R  <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(R)) return(.fit_glmfit(X, y))
    b <- as.vector(backsolve(R, backsolve(R, bb, transpose = TRUE)))
  }
  .fit_glmfit(X, y)
}

.fit <- switch(SWEEP_ENGINE, fastglm = .fit_fastglm, irls = .fit_irls, glmfit = .fit_glmfit,
               stop("SWEEP_ENGINE must be \"fastglm\", \"irls\" or \"glmfit\"."))

# --- Pooled ROC metrics with one sort ----------------------------------------
# AUC is the trapezoid under the tie-collapsed ROC (pROC's value). Cut-offs are
# pROC's: midpoints between adjacent distinct scores, +Inf (nobody positive) and
# -Inf (everyone). Youden ties go to the cut-off closest to 0.5, as in pass 2.
.roc_stats <- function(p, y, np, nn) {
  if (anyNA(p)) return(rep(NA_real_, 7L))
  o  <- sort.list(p, method = "radix", decreasing = TRUE)
  ps <- p[o]; ys <- y[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  m    <- length(ps)
  keep <- c(ps[-1L] != ps[-m], TRUE)          # last index of each tied run
  tp   <- tp[keep]; fp <- fp[keep]; s <- ps[keep]
  nk   <- length(tp)
  auc  <- sum((c(0, tp[-nk]) + tp) * (fp - c(0, fp[-nk]))) / (2 * np * nn)

  tp  <- c(0, tp); fp <- c(0, fp)
  thr <- c(Inf, if (nk > 1) (s[-nk] + s[-1L]) / 2, -Inf)   # cut for "score >= s_j"
  sens <- tp / np
  spec <- 1 - fp / nn
  J    <- sens + spec
  tied <- which(J >= max(J) - 0.5 / (np * nn))          # same tie rule as pass 2
  j    <- tied[order(abs(thr[tied] - 0.5), -thr[tied])[1]]

  TP <- tp[j]; FP <- fp[j]; FN <- np - TP; TN <- nn - FP
  c((TP + TN) / (np + nn), sens[j], spec[j],
    if (TP + FP > 0) TP / (TP + FP) else NA_real_,
    if (TN + FN > 0) TN / (TN + FN) else NA_real_,
    thr[j], auc)
}

# --- Combination tasks (never materialise C(p, k) on the master) -------------
combo_tasks <- function(p, k, max_block = 5e3) {
  tasks <- list()
  rec <- function(prefix, start) {
    need  <- k - length(prefix)
    avail <- p - start + 1L
    if (avail < need) return(invisible(NULL))
    if (need == 0L || choose(avail, need) <= max_block) {
      tasks[[length(tasks) + 1L]] <<- list(prefix = prefix, start = start, need = need)
      return(invisible(NULL))
    }
    for (i in start:(p - need + 1L)) rec(c(prefix, i), i + 1L)
  }
  rec(integer(0), 1L)
  tasks
}

expand_task <- function(t, p) {
  if (t$need == 0L) return(matrix(as.integer(t$prefix), ncol = 1L))
  m <- combn(t$start:p, t$need)
  if (length(t$prefix))
    m <- rbind(matrix(as.integer(t$prefix), nrow = length(t$prefix), ncol = ncol(m)), m)
  storage.mode(m) <- "integer"
  m
}

.task_size <- function(t, p) if (t$need == 0L) 1 else choose(p - t$start + 1L, t$need)

# Top N per ranking column; the top N of a union of top-N sets is the global top N
keep_top <- function(cmb, met, N = SWEEP_RETAIN) {
  if (ncol(cmb) <= N) return(list(cmb = cmb, met = met))
  idx <- unique(unlist(lapply(RANK_COLS, function(cc) {
    o <- order(met[, cc], decreasing = TRUE, na.last = NA)
    o[seq_len(min(N, length(o)))]
  })))
  list(cmb = cmb[, idx, drop = FALSE], met = met[idx, , drop = FALSE])
}

# --- Sweep worker ----------------------------------------------------------
# Pooled out-of-fold predictions of one combination (design columns `sel`).
#
# Exact: nf fits, each fold warm-started from the previous fold's converged
# coefficients (training sets overlap by ~80%).
.exact_pred <- function(sel, linkinv) {
  out <- numeric(length(W$pooled_y)); b <- NULL; pos <- 0L
  for (f in seq_len(W$nf)) {
    fit <- .fit(W$Xtr[[f]][, sel, drop = FALSE], W$ytr[[f]], start = b)
    if (is.null(fit)) return(NULL)
    b  <- if (fit$conv) fit$coef else NULL
    nt <- W$nte[f]
    out[pos + seq_len(nt)] <- linkinv(W$Xte[[f]][, sel, drop = FALSE] %*% fit$coef)
    pos <- pos + nt
  }
  out
}

# One-step: one fit on all subjects (coefficients b, fitted mu, weights
# w = mu(1 - mu)); then for each fold f, with its subjects left out,
#   b_f = b - (X'WX - X_f'W_f X_f)^-1 X_f'(y_f - mu_f)
# This is one IRLS iteration on the training rows started from b: at b the
# score over all subjects is zero, so the training score is minus the fold's
# score, and the training information is the full information minus the
# fold's share. Only the fold's own rows are touched. Returns NULL near
# separation, where one step is not enough; the caller then refits per fold.
.onestep_pred <- function(sel, linkinv) {
  X    <- W$Xall[, sel, drop = FALSE]
  full <- .fit(X, W$y)
  if (is.null(full) || !full$conv) return(NULL)
  b   <- full$coef
  eta <- as.vector(X %*% b)
  if (max(abs(eta)) > ONESTEP_ETA_MAX) return(NULL)
  mu <- 1 / (1 + exp(-eta)); w <- mu * (1 - mu); r <- W$y - mu
  H  <- crossprod(X, X * w)
  out <- numeric(length(W$pooled_y)); pos <- 0L
  for (f in seq_len(W$nf)) {
    S  <- W$tes[[f]]; XS <- X[S, , drop = FALSE]
    R  <- tryCatch(chol(H - crossprod(XS, XS * w[S])), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    bf <- b - backsolve(R, backsolve(R, crossprod(XS, r[S]), transpose = TRUE))
    out[pos + seq_along(S)] <- linkinv(XS %*% bf)
    pos <- pos + length(S)
  }
  out
}

# Sweep metrics for the combinations in cmb (one per column). Column 9 is 1
# when the combination was refitted per fold, 0 when the one-step was used;
# column 8 is unused (pass 2 supplies the diagnostics).
sweep_cmb <- function(cmb, onestep = SWEEP_CV == "onestep") {
  m <- ncol(cmb); k <- nrow(cmb); nb <- W$nb
  met     <- matrix(NA_real_, m, 10L)
  linkinv <- stats::binomial()$linkinv                # as predict.glm, so saturated predictions tie as in caret
  sel     <- c(1L, integer(k), if (nb > 1L) 2:nb)     # caret's column order: intercept, markers, confounders
  bm_slot <- 1L + seq_len(k)
  for (j in seq_len(m)) {
    sel[bm_slot] <- nb + cmb[, j]
    pred <- if (onestep) .onestep_pred(sel, linkinv) else NULL
    met[j, 9L] <- is.null(pred)
    if (is.null(pred)) pred <- .exact_pred(sel, linkinv)
    if (!is.null(pred)) met[j, 1:7] <- .roc_stats(pred, W$pooled_y, W$np, W$nn)
  }
  met[, 10L] <- (met[, 2L] + met[, 3L]) / 2
  met
}

# Sweep a block and keep its top, counting the combinations refitted per fold
.sweep_keep <- function(cmb) {
  met <- sweep_cmb(cmb)
  out <- keep_top(cmb, met)
  out$n_exact <- as.numeric(sum(met[, 9L] == 1))
  out
}
sweep_task        <- function(t)   .sweep_keep(expand_task(t, .DSN$p))
sweep_matrix_task <- function(cmb) .sweep_keep(cmb)

.setup_workers <- function(cl, extra) {
  clusterExport(cl, c(".DSN", ".fit", ".fit_glmfit", ".roc_stats", "sweep_cmb", "keep_top",
                      ".exact_pred", ".onestep_pred", ".sweep_keep", "prepare_worker",
                      "SWEEP_RETAIN", "RANK_COLS", "SWEEP_ENGINE", "SWEEP_CV", "ONESTEP_ETA_MAX", extra),
                envir = globalenv())
  clusterEvalQ(cl, {
    if (SWEEP_ENGINE == "fastglm") suppressMessages(library(fastglm))
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)   # one BLAS thread per worker
    if (!exists("W", envir = globalenv()) || !identical(get("W", envir = globalenv())$id, .DSN$id))
      prepare_worker()
    NULL
  })
}

# Runs `tasks` in groups across the workers, keeping the top as it goes, with
# progress counted in combinations. Returns the top, the number of combinations
# swept (n_total) and how many of them were refitted per fold (n_exact).
.run_tasks <- function(cl, tasks, sizes, task_fun, k, group_mult = 4L) {
  groups <- split(seq_along(tasks), ceiling(seq_along(tasks) / max(1L, length(cl) * group_mult)))
  total  <- sum(sizes)
  t0 <- Sys.time(); done <- 0; n_exact <- 0
  buf_cmb <- NULL; buf_met <- NULL
  for (g in groups) {
    res <- parLapplyLB(cl, tasks[g], task_fun)
    n_exact <- n_exact + sum(vapply(res, `[[`, numeric(1), "n_exact"))
    buf_cmb <- cbind(buf_cmb, do.call(cbind, lapply(res, `[[`, "cmb")))
    buf_met <- rbind(buf_met, do.call(rbind, lapply(res, `[[`, "met")))
    if (ncol(buf_cmb) > 20000L) { red <- keep_top(buf_cmb, buf_met); buf_cmb <- red$cmb; buf_met <- red$met }
    done <- done + sum(sizes[g])
    el   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- done / max(el, 1e-6)
    cat(sprintf("\r  [k=%d] %s/%s (%5.1f%%) | %s/sec | elapsed %s | ETA %s      ",
                k, format(done, big.mark = ","), format(total, big.mark = ","), 100 * done / total,
                format(round(rate), big.mark = ","), format_duration(el), format_duration((total - done) / rate)))
    flush.console()
  }
  top <- keep_top(buf_cmb, buf_met)
  top$n_total <- total
  top$n_exact <- n_exact
  cat(sprintf("\n  k = %d done in %s | %d retained%s\n", k,
              format_duration(as.numeric(difftime(Sys.time(), t0, units = "secs"))), ncol(top$cmb),
              if (SWEEP_CV == "onestep")
                sprintf(" | refitted per fold (near separation): %s (%s)",
                        format(n_exact, big.mark = ","), .pct(n_exact / total)) else ""))
  top
}

# Exhaustive level k
run_sweep <- function(k, cl) {
  tasks <- combo_tasks(.DSN$p, k)
  sizes <- vapply(tasks, .task_size, numeric(1), p = .DSN$p)
  cat(sprintf("  k = %d | %s combinations | engine: %s\n", k, format(sum(sizes), big.mark = ","), SWEEP_ENGINE))
  .setup_workers(cl, c("expand_task", "sweep_task"))
  .run_tasks(cl, tasks, sizes, sweep_task, k)
}

# Beam level k: every retained (k-1)-marker model extended by each marker it
# lacks. It cannot find a k-marker model none of whose (k-1)-marker sub-models
# made the parents; a wider BEAM_WIDTH lowers that risk.
extend_parents <- function(parents, p) {
  km1 <- nrow(parents)
  pieces <- lapply(seq_len(ncol(parents)), function(i) {
    par <- parents[, i]
    js  <- setdiff(seq_len(p), par)
    if (!length(js)) return(NULL)
    M <- apply(rbind(matrix(par, nrow = km1, ncol = length(js)), js), 2L, sort.int)
    if (!is.matrix(M)) M <- matrix(M, nrow = km1 + 1L)
    M
  })
  cand <- do.call(cbind, pieces)
  dimnames(cand) <- NULL
  storage.mode(cand) <- "integer"
  cand[, !duplicated(t(cand)), drop = FALSE]
}

run_beam_level <- function(parents, k, cl, chunk = 2000L) {
  cand  <- extend_parents(parents, .DSN$p)
  total <- ncol(cand)
  cat(sprintf("  k = %d | beam: %d parents -> %s candidates (exhaustive: %s)\n", k, ncol(parents),
              format(total, big.mark = ","), format(choose(.DSN$p, k), big.mark = ",", scientific = FALSE)))
  .setup_workers(cl, "sweep_matrix_task")
  starts <- seq(1L, total, by = chunk)
  tasks  <- lapply(starts, function(st) cand[, st:min(st + chunk - 1L, total), drop = FALSE])
  .run_tasks(cl, tasks, vapply(tasks, ncol, integer(1)), sweep_matrix_task, k)
}

# --- Preflight: time a sample at the largest k and project the runtime ------
preflight <- function(ks, counts, modes, selected_size = NA_integer_) {
  if (!exists("W", envir = globalenv()) || !identical(W$id, .DSN$id)) prepare_worker()
  k_probe <- max(ks)
  set.seed(1)
  cmb <- vapply(seq_len(100), function(i) sort(sample.int(.DSN$p, k_probe)), integer(k_probe))
  if (!is.matrix(cmb)) cmb <- matrix(cmb, nrow = k_probe)
  per  <- system.time(sweep_cmb(cmb))[["elapsed"]] / 100
  secs <- counts * per / numCores

  cat("\n===== Preflight =====\n")
  cat(sprintf("  Markers: %d | n: %d | sweep folds: %d | engine: %s | %.3f ms per model (1 core, k = %d)\n",
              .DSN$p, .DSN$n, length(.DSN$folds), SWEEP_ENGINE, per * 1000, k_probe))
  for (i in seq_along(ks))
    cat(sprintf("  k = %-2d %-10s %18s models  %s\n", ks[i], modes[i],
                format(counts[i], big.mark = ",", scientific = FALSE), format_duration(secs[i])))
  cat(sprintf("  all                %18s models  %s (sweep only)\n",
              format(sum(counts), big.mark = ",", scientific = FALSE), format_duration(sum(secs))))
  if (any(modes == "beam")) cat("  (beam counts are upper bounds)\n")

  if (is.finite(MAX_RUNTIME_HOURS) && sum(secs) / 3600 > MAX_RUNTIME_HOURS)
    stop(sprintf(paste0("Projected sweep time %.1f h exceeds MAX_RUNTIME_HOURS (%.1f). Options:\n",
                        "  - Script 1: %s_MARKERS <- \"selected\" (%s instead of all %d markers)\n",
                        "  - SEARCH_MODE <- \"beam\" or \"auto\"\n",
                        "  - fewer SWEEP_REPEATS or a lower max_biomarkers (now %d)\n",
                        "  - a higher MAX_RUNTIME_HOURS"),
                 sum(secs) / 3600, MAX_RUNTIME_HOURS, toupper(PROJECT),
                 if (is.na(selected_size)) "the selected list" else sprintf("%d selected markers", selected_size),
                 .DSN$p, max(ks)), call. = FALSE)
  invisible(per)
}


# ==========================================================================
# 8) PIPELINE (one comparison)
# ==========================================================================
run_pipeline_fast <- function(spec) {

  cat("\n", strrep("=", 60), "\n", sprintf("  %s | condL from %s\n", spec$tag, spec$info$condl_source),
      strrep("=", 60), "\n", sep = "")

  data_df    <- spec$data
  biomarkers <- spec$biomarkers
  conf       <- resolve_confounders(data_df, if (!is.null(spec$model_name)) spec$covariates else confounders)
  cat(sprintf("  Confounders: %s\n", if (is.null(conf)) "none" else paste(conf, collapse = ", ")))
  cc <- drop_incomplete_confounders(data_df, conf)
  data_df   <- cc$data
  n_dropped <- cc$n_dropped + if (!is.null(spec$common_dropped)) spec$common_dropped else 0L

  missing_bm <- setdiff(biomarkers, colnames(data_df))
  if (length(missing_bm) > 0) {
    warning(sprintf("%d marker(s) not in the data and dropped: %s", length(missing_bm),
                    paste(missing_bm, collapse = ", ")))
    biomarkers <- intersect(biomarkers, colnames(data_df))
  }
  if (length(biomarkers) < max_biomarkers)
    stop(sprintf("max_biomarkers (%d) exceeds the number of candidate markers (%d).",
                 max_biomarkers, length(biomarkers)))
  cat(sprintf("  Candidate markers: %d\n", length(biomarkers)))

  group_counts <- as.list(table(data_df$condL))
  for (g in names(group_counts))
    cat(sprintf("    %s: %d%s\n", g, group_counts[[g]], if (g == spec$info$case) " (case)" else " (reference)"))

  # factor() keeps the reference-first order set by make_spec()
  data_df$condL <- factor(data_df$condL)
  levels(data_df$condL) <- make.names(levels(data_df$condL))

  # one set of folds for both passes; the sweep uses the first SWEEP_REPEATS repeats
  idx <- if (!is.null(spec$folds)) spec$folds else caret_fold_index(data_df$condL)
  ctrl_exact       <- control
  ctrl_exact$index <- idx
  n_sweep_folds    <- CV_FOLDS * min(SWEEP_REPEATS, CV_REPEATS)
  assign(".DSN", build_design(data_df, biomarkers, conf, folds_idx = idx[seq_len(n_sweep_folds)],
                              id = paste(spec$tag, format(Sys.time(), "%Y%m%d %H:%M:%OS6"))),
         envir = globalenv())
  prepare_worker()
  assign("SWEEP_RETAIN", max(SWEEP_KEEP, BEAM_WIDTH), envir = globalenv())

  p  <- .DSN$p
  ks <- min_biomarkers:max_biomarkers
  all_levels <- seq_len(max_biomarkers)
  modes <- vapply(all_levels, function(k) {
    if (k <= 2L) return("exhaustive")
    switch(SEARCH_MODE, exhaustive = "exhaustive", beam = "beam",
           auto = if (choose(p, k) <= EXHAUSTIVE_LIMIT) "exhaustive" else "beam",
           stop("SEARCH_MODE must be \"auto\", \"exhaustive\" or \"beam\"."))
  }, character(1))
  # beam at level k needs level k-1, so then every level from 1 is swept
  levels_run  <- if (any(modes[ks] == "beam")) all_levels else ks
  max_parents <- length(RANK_COLS) * BEAM_WIDTH
  counts <- vapply(levels_run, function(k) {
    if (modes[k] == "exhaustive") choose(p, k)
    else min(choose(p, k), min(max_parents, choose(p, k - 1)) * (p - k + 1))
  }, numeric(1))
  if (RUN_PREFLIGHT) preflight(levels_run, counts, modes[levels_run], spec$selected_size)

  t_start <- Sys.time()
  cl <- makeCluster(numCores)
  on.exit({ try(stopCluster(cl), silent = TRUE); registerDoSEQ() }, add = TRUE)

  all_results <- list(); n_evaluated <- 0; best_by_k <- c(); exact_share <- c(); parents <- NULL
  for (k in levels_run) {
    top <- if (modes[k] == "beam") run_beam_level(parents, k, cl) else run_sweep(k, cl)
    n_evaluated <- n_evaluated + top$n_total
    exact_share[as.character(k)] <- top$n_exact / top$n_total
    best_by_k[as.character(k)] <- suppressWarnings(max(top$met[, 10], na.rm = TRUE))
    cat(sprintf("  [k=%d] best sweep balanced accuracy: %.4f\n", k, best_by_k[[as.character(k)]]))
    parents <- keep_top(top$cmb, top$met, BEAM_WIDTH)$cmb
    if (!(k %in% ks)) next                       # swept only to seed the beam

    exp_top <- keep_top(top$cmb, top$met, SWEEP_KEEP)
    combos  <- lapply(seq_len(ncol(exp_top$cmb)), function(j) biomarkers[exp_top$cmb[, j]])
    cat(sprintf("  [k=%d] refitting %d through caret (%d x %d CV)...\n", k, length(combos), CV_FOLDS, CV_REPEATS))
    registerDoSNOW(cl)
    exact <- run_foreach_chunk(combos, cl, data_df, spec$stats, conf, ctrl_exact)
    all_results[[as.character(k)]] <- assemble_results(exact, k) %>%
      mutate(across(where(is.numeric), ~ round(., 10)))
  }
  t_end <- Sys.time()

  cat(sprintf("\n  Total: %s | models swept: %s (exhaustive equivalent: %s)\n",
              format_duration(as.numeric(difftime(t_end, t_start, units = "secs"))),
              format(n_evaluated, big.mark = ",", scientific = FALSE),
              format(sum(choose(p, ks)), big.mark = ",", scientific = FALSE)))
  # if this flattens, extra markers are fitting noise
  cat("  Best sweep balanced accuracy by number of markers:\n")
  for (nm in names(best_by_k))
    cat(sprintf("    k = %s: %.4f%s\n", nm, best_by_k[[nm]],
                if (as.integer(nm) > 1) sprintf("  (%+.4f)", best_by_k[[nm]] - best_by_k[[as.character(as.integer(nm) - 1)]]) else ""))

  params <- build_params_sheet(spec, biomarkers, conf, t_start, t_end, group_counts, n_dropped)
  onestep <- SWEEP_CV == "onestep"
  params <- rbind(params, data.frame(
    Parameter = c("", "=== SEARCH ===", "Search mode", "Exhaustive limit", "Beam width (per metric)",
                  "Models swept", "Exhaustive equivalent", "Sweep ranking", "Export",
                  paste0("  k = ", levels_run, " search"), paste0("  k = ", levels_run, " best sweep balanced accuracy"),
                  if (onestep) paste0("  k = ", levels_run, " refitted per fold in the sweep (near separation)")),
    Value = c("", "", SEARCH_MODE, format(EXHAUSTIVE_LIMIT, big.mark = ",", scientific = FALSE),
              as.character(BEAM_WIDTH), format(n_evaluated, scientific = FALSE),
              format(sum(choose(p, ks)), scientific = FALSE),
              paste(c("Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "Threshold", "AUC",
                      "Max_Abs_Coef", "Converged", "Balanced_Accuracy")[RANK_COLS], collapse = ", "),
              sprintf("Top %d per Balanced_Accuracy, Sensitivity and AUC among the refitted models", TOP_N),
              modes[levels_run], sprintf("%.4f", best_by_k[as.character(levels_run)]),
              if (onestep) .pct(exact_share[as.character(levels_run)])),
    stringsAsFactors = FALSE))

  if (!is.null(spec$model_name)) params <- rbind(params, data.frame(
    Parameter = c("Model configuration", "Cohort policy"),
    Value = c(spec$model_name, "Participants complete for Age_DOS and BMI; no Sex-based exclusion")))
  export_list <- list(Analysis_Parameters = params)
  for (k in ks) {
    df_k <- get_top_n_per_metric(all_results[[as.character(k)]])
    cat(sprintf("  %d marker(s): %d refitted -> %d exported\n", k, nrow(all_results[[as.character(k)]]), nrow(df_k)))
    export_list[[paste0("Results_", k, "_Biomarkers")]] <- df_k
  }
  out_file <- sprintf("Results_LogReg_%s_%s.xlsx", spec$tag, registry_key)
  write_xlsx(export_list, path = out_file)
  cat("  Saved:", out_file, "\n")
  invisible(all_results)
}


# ==========================================================================
# 9) CHECK THE SWEEP
# ==========================================================================
# 1) The sweep's exact path against caret: ~50 random combinations through
#    both, on the same folds (the sweep's). Accuracy, Sensitivity,
#    Specificity, Precision, NPV and Threshold should agree exactly; AUC to
#    rounding, or in the 4th-5th decimal where predictions nearly tie (both
#    stop at glm's 1e-8 tolerance from different starts). Rows flagged for
#    separation are reported apart: glm stops at 25 iterations, the IRLS at 30.
# 2) SWEEP_CV = "onestep": the one-step against the exact path on up to
#    n_rank combinations of k markers (all of them when there are fewer): how
#    close the metrics are, how similar the ranking is, and whether the exact
#    top 10 are inside the one-step's top SWEEP_KEEP, the set pass 2 refits.

# n combinations of k of p markers, one per column; all of them when there are
# no more than n
.sample_combos <- function(p, k, n, seed) {
  if (choose(p, k) <= n) return(combn(p, k))
  set.seed(seed)
  pick <- matrix(integer(0), nrow = k, ncol = 0)
  while (ncol(pick) < n) {
    pick <- cbind(pick, matrix(replicate(n, sort(sample.int(p, k))), nrow = k))
    pick <- pick[, !duplicated(t(pick)), drop = FALSE]
  }
  pick[, seq_len(n), drop = FALSE]
}

validate_fast_vs_caret <- function(spec, k = 3L, n_sample = 50L, n_rank = 2000L) {
  conf       <- resolve_confounders(spec$data, if (!is.null(spec$model_name)) spec$covariates else confounders)
  data_df    <- drop_incomplete_confounders(spec$data, conf)$data
  biomarkers <- intersect(spec$biomarkers, colnames(data_df))
  k          <- min(as.integer(k), length(biomarkers))
  data_df$condL <- factor(data_df$condL)
  levels(data_df$condL) <- make.names(levels(data_df$condL))

  idx  <- caret_fold_index(data_df$condL)
  idx  <- idx[seq_len(CV_FOLDS * min(SWEEP_REPEATS, CV_REPEATS))]
  ctrl <- control
  ctrl$index   <- idx
  ctrl$repeats <- min(SWEEP_REPEATS, CV_REPEATS)
  assign(".DSN", build_design(data_df, biomarkers, conf, folds_idx = idx,
                              id = paste("validate", format(Sys.time(), "%H:%M:%OS6"))), envir = globalenv())
  prepare_worker()

  p    <- length(biomarkers)
  pick <- .sample_combos(p, k, n_sample, seed = 99)

  cat(sprintf("  1) Exact sweep vs caret: %d combinations of %d markers\n", ncol(pick), k))
  tf   <- system.time(fast <- sweep_cmb(pick, onestep = FALSE))[["elapsed"]]
  cl   <- makeCluster(2L); registerDoSNOW(cl)
  ts   <- system.time(slow <- run_foreach_chunk(lapply(seq_len(ncol(pick)), function(j) biomarkers[pick[, j]]),
                                                cl, data_df, spec$stats, conf, ctrl))[["elapsed"]]
  stopCluster(cl); registerDoSEQ()

  slow_m <- t(vapply(slow, function(r)
    if (is.list(r) && length(r$metrics) == 10) as.numeric(r$metrics[1:7]) else rep(NA_real_, 7), numeric(7)))
  sep    <- vapply(slow, function(r) is.list(r) && length(r$diagnostics) == 3 && isTRUE(r$diagnostics[2] == 1),
                   logical(1))
  failed <- is.na(slow_m[, 7])
  use    <- !sep & !failed
  nm  <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "Threshold", "AUC")
  cmp <- data.frame(metric = nm, max_abs_diff = vapply(seq_len(7), function(i)
    if (any(use)) max(abs(fast[use, i] - slow_m[use, i]), na.rm = TRUE) else NA_real_, numeric(1)))
  print(cmp, row.names = FALSE)
  cat(sprintf("\n  %d combinations compared | %d flagged for separation (left out) | %d caret failures\n",
              sum(use), sum(sep), sum(failed)))
  if (any(sep))
    cat(sprintf("  Separation-flagged rows, max |AUC difference|: %.4f\n",
                max(abs(fast[sep, 7] - slow_m[sep, 7]), na.rm = TRUE)))
  cat(sprintf("  sweep (1 core): %.2fs | caret (2 cores): %.2fs\n", tf, ts))

  if (SWEEP_CV == "onestep") {
    pick2 <- .sample_combos(p, k, n_rank, seed = 98)
    cat(sprintf("\n  2) One-step vs exact sweep: %s combinations of %d markers%s\n",
                format(ncol(pick2), big.mark = ","), k,
                if (ncol(pick2) == choose(p, k)) " (all of them)" else " (random)"))
    t_ex <- system.time(ex <- sweep_cmb(pick2, onestep = FALSE))[["elapsed"]]
    t_os <- system.time(os <- sweep_cmb(pick2, onestep = TRUE))[["elapsed"]]
    used   <- os[, 9] == 0 & !is.na(ex[, 10]) & !is.na(os[, 10])
    keep_n <- min(SWEEP_KEEP, ncol(pick2))
    top10  <- order(-ex[, 10], -ex[, 7])[seq_len(min(10, ncol(pick2)))]
    rk     <- rank(-os[, 10], ties.method = "min", na.last = "keep")
    cat(sprintf("    refitted per fold (near separation): %d (%s)\n", sum(os[, 9] == 1), .pct(mean(os[, 9] == 1))))
    if (sum(used) >= 2) {
      d_auc <- abs(ex[used, 7] - os[used, 7]); d_ba <- abs(ex[used, 10] - os[used, 10])
      cat(sprintf("    |AUC difference|: median %.4f, 95th percentile %.4f\n", median(d_auc), quantile(d_auc, 0.95)))
      cat(sprintf("    balanced accuracy: identical in %.0f%%, |difference| 95th percentile %.4f\n",
                  100 * mean(d_ba < 1e-12), quantile(d_ba, 0.95)))
      cat(sprintf("    rank correlation (Spearman) of balanced accuracy: %.3f\n",
                  cor(ex[used, 10], os[used, 10], method = "spearman")))
    }
    cat(sprintf("    exact top %d inside the one-step top %d (the set pass 2 refits): %d of %d\n",
                length(top10), keep_n, sum(rk[top10] <= keep_n, na.rm = TRUE), length(top10)))
    missed <- top10[is.na(rk[top10]) | rk[top10] > keep_n]
    if (length(missed))
      cat(sprintf("    left out: exact balanced accuracy %s; lowest kept: %.4f (close values = near-ties)\n",
                  paste(sprintf("%.4f", ex[missed, 10]), collapse = ", "),
                  min(ex[which(rk <= keep_n), 10], na.rm = TRUE)))
    cat(sprintf("    time (1 core): exact %.1fs, one-step %.1fs (%.1fx faster)\n", t_ex, t_os, t_ex / max(t_os, 1e-6)))
  }
  invisible(list(fast = fast, slow = slow_m, cmp = cmp, separated = sep))
}


# ==========================================================================
# 10) RUN AGE-AND-BMI ADJUSTED DMFT SEARCH
# ==========================================================================
# One age/BMI-adjusted search uses shared participants and folds at every panel size.
# Selection and Youden cutoff optimization still use development predictions;
# this comparison is exploratory, not independent validation.
if (PROJECT != "DMFT") stop("This age/BMI runner requires PROJECT = 'DMFT'.")
if (!exists("DMFT_MODEL_SETS"))
  stop("Run 01_Preprocessing_DMFT_Age_BMI.R first in this session.")

run_dmft_age_bmi <- function(base_spec) {
  if (SEARCH_MODE != "beam") stop("This study version requires SEARCH_MODE = beam (from size 3).")
  sets <- DMFT_MODEL_SETS
  if (!identical(names(sets), "Biomarkers_Age_BMI") ||
      !identical(sets[[1]]$covariates, c("Age_DOS", "BMI")) ||
      !grepl("_exhaustive12_beam3$", registry_key))
    stop("Run 01_Preprocessing_DMFT_Age_BMI.R before this revised analysis.")
  required <- unique(unlist(lapply(sets, `[[`, "covariates")))
  missing <- setdiff(required, names(base_spec$data))
  if (length(missing)) stop("Required columns missing: ", paste(missing, collapse = ", "))
  d <- as.data.frame(base_spec$data)
  for (v in c("Age_DOS", "BMI")) {
    if (!is.numeric(d[[v]])) stop(v, " must be numeric in original units.")
    if (any(!is.na(d[[v]]) & !is.finite(d[[v]]))) stop(v, " contains infinite values.")
  }
  ok <- complete.cases(d[, required, drop = FALSE]) & !is.na(d$condL)
  common_dropped <- sum(!ok)
  source_rows <- which(ok)
  d <- droplevels(d[ok, , drop = FALSE])
  if (!setequal(as.character(unique(d$condL)), c("LD", "HD")))
    stop("Expected LD and HD outcome labels.")
  d$condL <- factor(d$condL, levels = c("LD", "HD"))
  if (min(table(d$condL)) < CV_FOLDS) stop("Too few participants per group for CV_FOLDS.")
  if (any(!is.finite(as.matrix(d[, base_spec$biomarkers, drop = FALSE]))))
    stop("Missing or non-finite preprocessed biomarkers remain.")
  folds <- caret_fold_index(d$condL)
  ctl <- control; ctl$index <- folds
  # Persist exact fold assignment and original row mapping for reproducibility.
  saveRDS(list(index = folds, source_rows = source_rows, outcome = d$condL,
               registry_key = registry_key,
               model_sets = sets, common_dropped = common_dropped),
          sprintf("DMFT_Age_BMI_Resampling_%s.rds", registry_key))
  message("Age/BMI cohort: ", nrow(d), " participants; excluded ", common_dropped)
  all_summary <- list(); status <- character(0); all_runs <- list()
  metric_cols <- c("Balanced_Accuracy", "Accuracy", "Sensitivity", "Specificity",
                   "Precision", "NPV", "Threshold", "AUC",
                   "Balanced_Accuracy_CI_Lower", "Balanced_Accuracy_CI_Upper")
  for (nm in names(sets)) {
    spec <- base_spec; spec$data <- d; spec$model_name <- nm
    spec$covariates <- sets[[nm]]$covariates
    spec$common_dropped <- common_dropped; spec$folds <- folds
    spec$tag <- paste0(base_spec$tag, "_", nm)
    status[nm] <- tryCatch({
      if (!isTRUE(sets[[nm]]$biomarkers)) {
        started <- Sys.time(); set.seed(123)
        fit <- caret::train(reformulate(spec$covariates, response = "condL"),
                            data = d, method = "glm", family = "binomial",
                            metric = "BalancedAccuracy", trControl = ctl)
        pred <- fit$pred; lev <- levels(d$condL)
        if (any(!is.finite(pred[[lev[2]]]))) stop("Non-finite baseline predictions.")
        roc_obj <- pROC::roc(pred$obs, pred[[lev[2]]], levels = lev, direction = "<", quiet = TRUE)
        m <- youden_best(roc_obj); ci <- repeat_ci(pred, lev, m$threshold)
        row <- data.frame(Balanced_Accuracy = (m$sensitivity + m$specificity) / 2,
                          Accuracy = m$accuracy, Sensitivity = m$sensitivity,
                          Specificity = m$specificity, Precision = m$precision,
                          NPV = m$npv, Threshold = m$threshold, AUC = as.numeric(pROC::auc(roc_obj)),
                          Balanced_Accuracy_CI_Lower = ci[1], Balanced_Accuracy_CI_Upper = ci[2],
                          Balanced_Accuracy_p_Fisher = ci[3], Converged = fit$finalModel$converged)
        pars <- build_params_sheet(spec, character(0), spec$covariates, started, Sys.time(),
                                   as.list(table(d$condL)), common_dropped)
        pars$Value[which(pars$Parameter == "Markers per model")] <- "0 (demographic baseline)"
        pars$Value[which(pars$Parameter %in% c("Sweep", "Sweep CV"))] <- "Not applicable; fixed age + BMI model"
        pars <- rbind(pars, data.frame(Parameter = c("Model configuration", "Cohort policy"),
                Value = c(nm, "Participants complete for Age_DOS and BMI; no Sex-based exclusion")))
        co <- coef(summary(fit$finalModel))
        co <- data.frame(Term = rownames(co), co, row.names = NULL, check.names = FALSE)
        write_xlsx(list(Analysis_Parameters = pars, Results_0_Biomarkers = row,
                        Coefficients = co, CV_Predictions = pred),
                   sprintf("Results_LogReg_%s_%s.xlsx", spec$tag, registry_key))
        all_runs[[nm]] <- list(`0` = row)
      } else {
        if (RUN_VALIDATION) validate_fast_vs_caret(spec)
        all_runs[[nm]] <- run_pipeline_fast(spec)
      }
      summaries <- lapply(names(all_runs[[nm]]), function(k) {
        res <- all_runs[[nm]][[k]]
        valid <- is.finite(res$Balanced_Accuracy) & is.finite(res$AUC)
        if (!any(valid)) stop("No valid exact refits for ", nm, " at size ", k)
        best <- sort_by_balanced_accuracy(res[valid, , drop = FALSE])[1, , drop = FALSE]
        bcols <- grep("^Biomarker[0-9]+$", names(best), value = TRUE)
        data.frame(Model = nm, N = nrow(d), LD = sum(d$condL == "LD"), HD = sum(d$condL == "HD"),
                   N_Biomarkers = as.integer(k),
                   Covariates = if (length(spec$covariates)) paste(spec$covariates, collapse = ", ") else "none",
                   Panel = if (length(bcols)) paste(unlist(best[bcols]), collapse = "; ") else "none",
                   best[, metric_cols, drop = FALSE], row.names = NULL)
      })
      all_summary[[nm]] <- do.call(rbind, summaries)
      "done"
    }, error = function(e) {
      message("FAILED ", nm, ": ", conditionMessage(e))
      paste("FAILED:", conditionMessage(e))
    })
  }
  successful <- all_summary[status[names(all_summary)] == "done"]
  if (length(successful)) {
    by_size <- do.call(rbind, successful); rownames(by_size) <- NULL
    best_overall <- do.call(rbind, lapply(successful, function(x) {
      x[order(-x$Balanced_Accuracy, -x$AUC, x$N_Biomarkers), , drop = FALSE][1, , drop = FALSE]
    }))
    rownames(best_overall) <- NULL
  } else { by_size <- data.frame(); best_overall <- data.frame() }
  notes <- data.frame(Note = c(
    "Age and BMI adjustment; exhaustive one- and two-marker search, beam from three markers.",
    "Only biomarkers are transformed/standardized; Age_DOS and BMI retain original units.",
    "Each biomarker configuration has its own 1-6 marker search. Best means highest balanced accuracy, then AUC.",
    "Each model uses its own Youden threshold selected from pooled out-of-fold predictions.",
    "Preprocessing and model/cutoff selection use development data; results are exploratory and may be optimistic.",
    "Existing preprocessing defines DMFT labels before common-cohort exclusion; labels/cutoff are not recalculated.",
    "Repeat percentile intervals describe variability over fold assignments, not independent-validation confidence intervals."))
  write_xlsx(list(Best_Overall = best_overall, Best_By_Size = by_size,
                  Run_Status = data.frame(Model = names(status), Status = unname(status)), Notes = notes),
             sprintf("DMFT_Age_BMI_Summary_%s.xlsx", registry_key))
  print(data.frame(Model = names(status), Status = unname(status)))
  if (any(startsWith(status, "FAILED"))) warning("Some models failed; review Run_Status before comparing results.")
  invisible(list(results = all_runs, summary = best_overall, by_size = by_size, status = status))
}
DMFT_AGE_BMI_RESULTS <- run_dmft_age_bmi(SPECS[[1]])
