# ==========================================================================
# 04b_Export_DMFT_Equation_Single_Biomarker.R
# ==========================================================================
# Exports the calculator equation for the best SINGLE-biomarker model of the
# 10-repeat run (one biomarker with age and BMI). Under nested cross-validation
# this model generalized better than the larger panels, so the manuscript
# (Section 3.4) and the calculator now use it instead of the six-marker panel.
#
# HOW TO RUN
#   1. In a fresh R session, run the preprocessing script with the settings of
#      the 10-repeat run (it creates data.dfPH, stats_PH and registry_key).
#   2. Put the run's Biomarkers_Age_BMI results workbook in the working
#      directory, or give its full path in RESULTS_FILE below.
#   3. source("04b_Export_DMFT_Equation_Single_Biomarker.R")
#
# WHAT IT DOES
#   - takes the top row of the Results_1_Biomarkers sheet (the best single
#     biomarker, Bacteroides heparinolyticus in this run) and its threshold;
#   - refits the logistic regression of HD on that biomarker, Age_DOS and BMI
#     in the run's cohort (glm, like the search's final fits) and checks the
#     refit against the odds ratios in the workbook;
#   - converts the coefficients to raw-input form (the LD-group centring and
#     scaling are absorbed into the intercept and slope) and checks that the
#     raw-input equation reproduces the fitted log-odds;
#   - computes odds ratios with 95% profile-likelihood confidence intervals,
#     the calculator's worked example, and the values for the manuscript.
#
# OUTPUTS (in OUT_DIR)
#   DMFT_equation_single_biomarker.txt   equation, odds ratios, threshold, checks, example,
#                                        and the values to paste into Sections 3.4 and 3.7
#   dmft_age_bmi_model.csv               the model as a calculator file (Patient prediction ->
#                                        "Saved model"); same format as script 04's export
#   PUBLISHED_MODEL_single_biomarker.R   replacement PUBLISHED_MODEL and EXAMPLE blocks for app/app.R
#   DMFT_equation_single_biomarker.rds   the fit and everything above
#
# The refit uses all participants and adds no validation: the threshold and the
# panel were chosen on the development data, and the model needs testing in
# independent cohorts before use.
# ==========================================================================


# ==========================================================================
# 1) SETTINGS
# ==========================================================================
RESULTS_FILE  <- "Results_LogReg_HD_vs_LD_Biomarkers_Age_BMI_DMFT_cutoff13median_conf_AgeDOS_BMI_Sex_allBiomarkers_noRatios_noVif_exhaustive12_beam3.xlsx"
OUT_DIR       <- file.path(getwd(), "DMFT_single_biomarker_equation")
PANEL_SIZE    <- 1L                     # the best panel of this size is exported (1 = best single biomarker)
EXAMPLE_INPUT <- list(biomarkers = 0,   # the calculator's example values (one per biomarker), not a participant
                      Age_DOS = 40, BMI = 25)
DELTA         <- 5e-7                   # offset of log2(raw + DELTA), as in the preprocessing
CHECK_TOL     <- 1e-4                   # largest relative difference from the workbook's odds ratios
REQUIRE_SAME_REGISTRY <- TRUE           # stop when registry_key in this session is not the run's


# ==========================================================================
# 2) THE RUN: WORKBOOK, PANEL AND COHORT
# ==========================================================================
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package readxl is needed: install.packages(\"readxl\")", call. = FALSE)
if (!exists("data.dfPH"))
  stop("data.dfPH is not in this R session. Run the preprocessing script (PROJECT = \"DMFT\", ",
       "with the settings of the 10-repeat run) in this session first.", call. = FALSE)
if (!file.exists(RESULTS_FILE))
  stop(sprintf("Workbook not found: %s\n  Working directory: %s\n  Put it there or give its full path in RESULTS_FILE.",
               RESULTS_FILE, getwd()), call. = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.P <- as.data.frame(readxl::read_excel(RESULTS_FILE, sheet = "Analysis_Parameters", col_types = "text"))
.P <- data.frame(Parameter = trimws(.P[[1]]), Value = trimws(.P[[2]]), stringsAsFactors = FALSE)
.pget <- function(name) {
  v <- .P$Value[which(.P$Parameter == name)]
  if (length(v) && !is.na(v[1])) v[1] else NA_character_
}
.stopm <- character(0)
if (!identical(.pget("Model configuration"), "Biomarkers_Age_BMI"))
  .stopm <- c(.stopm, sprintf("RESULTS_FILE should be the Biomarkers_Age_BMI workbook (its model configuration: %s).",
                              .pget("Model configuration")))
REGISTRY_KEY <- .pget("Script 1 registry key")
CONF <- trimws(strsplit(.pget("Confounders"), ",", fixed = TRUE)[[1]])
if (!setequal(CONF, c("Age_DOS", "BMI")))
  .stopm <- c(.stopm, sprintf("Expected a model adjusted for Age_DOS and BMI; the workbook's confounders: %s.",
                              .pget("Confounders")))
CONF <- c("Age_DOS", "BMI")
REF  <- .pget("Reference group"); CASE <- .pget("Case group (positive)")
if (is.na(REF) || is.na(CASE)) { REF <- "LD"; CASE <- "HD" }
if (exists("registry_key") && !identical(get("registry_key"), REGISTRY_KEY)) {
  .msg <- sprintf("registry_key in this session (%s) is not the run's (%s).", get("registry_key"), REGISTRY_KEY)
  if (isTRUE(REQUIRE_SAME_REGISTRY)) {
    .stopm <- c(.stopm, paste(.msg, "Re-run the preprocessing with the run's settings, or set REQUIRE_SAME_REGISTRY <- FALSE."))
  } else warning(.msg, call. = FALSE)
}

# The best panel of size PANEL_SIZE: top row of its sheet (balanced accuracy, then AUC)
.sheet <- sprintf("Results_%d_Biomarkers", as.integer(PANEL_SIZE))
if (!.sheet %in% readxl::excel_sheets(RESULTS_FILE)) stop("RESULTS_FILE has no sheet ", .sheet, call. = FALSE)
.R <- suppressWarnings(as.data.frame(readxl::read_excel(RESULTS_FILE, sheet = .sheet)))
.num <- function(v) suppressWarnings(as.numeric(v))
for (.v in intersect(names(.R), c("Balanced_Accuracy", "AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV",
                                   "Threshold", "Balanced_Accuracy_Repeat_Min", "Balanced_Accuracy_Repeat_Max")))
  .R[[.v]] <- .num(.R[[.v]])
.need <- c(paste0("Biomarker", seq_len(PANEL_SIZE)), "Balanced_Accuracy", "AUC", "Threshold")
if (length(setdiff(.need, names(.R)))) stop("Sheet ", .sheet, " lacks: ", paste(setdiff(.need, names(.R)), collapse = ", "), call. = FALSE)
.R <- .R[order(-.R$Balanced_Accuracy, -.R$AUC, na.last = TRUE), , drop = FALSE]
ROW <- .R[1, , drop = FALSE]
MARKERS <- as.character(unlist(ROW[paste0("Biomarker", seq_len(PANEL_SIZE))], use.names = FALSE))
THRESHOLD <- ROW$Threshold
if (anyNA(MARKERS) || any(!nzchar(MARKERS))) stop("The top row of ", .sheet, " has no biomarker names.", call. = FALSE)
if (!is.finite(THRESHOLD) || THRESHOLD <= 0 || THRESHOLD >= 1) stop("Invalid threshold in ", .sheet, call. = FALSE)
if (length(EXAMPLE_INPUT$biomarkers) != length(MARKERS))
  stop(sprintf("EXAMPLE_INPUT$biomarkers needs %d value(s), one per biomarker.", length(MARKERS)), call. = FALSE)

# The cohort: participants complete for every covariate named in the run's cohort policy
.policy <- .pget("Cohort policy")
COHORT_VARS <- if (!is.na(.policy) && grepl("^Participants complete for ", .policy)) {
  trimws(strsplit(sub("\\s*\\(.*$", "", sub("^Participants complete for ", "", .policy)), ",", fixed = TRUE)[[1]])
} else character(0)
COHORT_VARS <- unique(c(COHORT_VARS, CONF))
DAT <- as.data.frame(data.dfPH)
.miss <- setdiff(c("condL", COHORT_VARS, MARKERS), names(DAT))
if (length(.miss))
  stop("Not in data.dfPH: ", paste(.miss, collapse = ", "),
       ".\n  Run the preprocessing with the settings of the run (registry key ", REGISTRY_KEY, ").", call. = FALSE)
.ok <- stats::complete.cases(DAT[, COHORT_VARS, drop = FALSE]) & !is.na(DAT$condL)
DAT <- DAT[.ok, , drop = FALSE]; rownames(DAT) <- NULL
DAT$condL <- factor(as.character(DAT$condL), levels = c(REF, CASE))
if (anyNA(DAT$condL)) stop(sprintf("condL in data.dfPH has values other than %s and %s.", REF, CASE), call. = FALSE)
N <- nrow(DAT); N_REF <- sum(DAT$condL == REF); N_CASE <- sum(DAT$condL == CASE)
.n_run <- suppressWarnings(as.integer(c(.pget("Total n"), .pget("n (LD)"), .pget("n (HD)"))))
if (!identical(c(N, N_REF, N_CASE), .n_run))
  .stopm <- c(.stopm, sprintf("Cohort in this session: n = %d (%s %d, %s %d); the run: n = %s (%s %s, %s %s), complete for %s.",
                              N, REF, N_REF, CASE, N_CASE, .n_run[1], REF, .n_run[2], CASE, .n_run[3],
                              paste(COHORT_VARS, collapse = ", ")))
.bad <- c(MARKERS, CONF)[!vapply(DAT[c(MARKERS, CONF)], function(z) is.numeric(z) && all(is.finite(z)), logical(1))]
if (length(.bad)) .stopm <- c(.stopm, paste("Not numeric, or missing values, in the cohort:", paste(.bad, collapse = ", ")))

# LD reference statistics of log2(raw + DELTA): the centring and scaling used for the biomarkers
.stats_name <- NULL
for (.s in c("stats_PH", "stats_PH_vif")) {
  if (exists(.s)) {
    .st <- get(.s)
    if (is.data.frame(.st) && all(c("Column", "Mean", "SD") %in% names(.st)) && all(MARKERS %in% .st$Column)) { .stats_name <- .s; break }
  }
}
if (is.null(.stats_name))
  .stopm <- c(.stopm, "No LD reference statistics (stats_PH with columns Column, Mean and SD) for the biomarker in this session.")
if (length(.stopm))
  stop("\n\n===== THIS SESSION DOES NOT MATCH THE RUN =====\n  ", paste(.stopm, collapse = "\n  "), "\n", call. = FALSE)
.st <- get(.stats_name)
MU <- setNames(.st$Mean[match(MARKERS, .st$Column)], MARKERS)
SDV <- setNames(.st$SD[match(MARKERS, .st$Column)], MARKERS)
if (any(!is.finite(MU)) || any(!is.finite(SDV)) || any(SDV <= 0))
  stop("The LD reference mean and SD of the biomarker must be finite, with a positive SD.", call. = FALSE)
# The biomarker columns should hold the LD-standardized values. The preprocessing standardizes
# over the LD group of data.dfPH (before the cohort is restricted to complete covariates), so
# check there; the SD is also checked against the workbook below.
for (.m in MARKERS) {
  .z <- as.numeric(data.dfPH[[.m]][as.character(data.dfPH$condL) == REF])
  .z <- .z[is.finite(.z)]
  if (abs(mean(.z)) > 1e-6 || abs(stats::sd(.z) - 1) > 1e-6)
    warning(sprintf(paste("%s in data.dfPH has LD mean %.3g and SD %.3g, not 0 and 1: check that it holds the values",
                          "standardized with %s (the conversion to raw-input coefficients assumes it does)."),
                    .m, mean(.z), stats::sd(.z), .stats_name), call. = FALSE)
}


# ==========================================================================
# 3) THE REFIT, ITS CHECKS AND THE RAW-INPUT EQUATION
# ==========================================================================
.d <- data.frame(y = DAT$condL, setNames(DAT[MARKERS], paste0("B", seq_along(MARKERS))),
                 Age_DOS = DAT$Age_DOS, BMI = DAT$BMI)
FIT <- stats::glm(y ~ ., data = .d, family = stats::binomial())
B <- stats::coef(FIT)
if (!isTRUE(FIT$converged) || any(!is.finite(B))) stop("The refit did not converge; nothing was written.", call. = FALSE)
.bn <- paste0("B", seq_along(MARKERS))
B_Z <- setNames(B[.bn], MARKERS)                       # per LD-SD of log2(raw + DELTA)

# Against the workbook: odds ratio per SD and per doubling, and the SD they imply
CHECKS <- data.frame(Check = character(0), Workbook = numeric(0), This_session = numeric(0), Relative_difference = numeric(0))
.add_check <- function(what, wb, here) {
  rd <- if (is.finite(wb)) abs(here - wb) / max(abs(wb), 1e-12) else NA_real_
  CHECKS[nrow(CHECKS) + 1L, ] <<- list(what, wb, here, rd)
}
for (j in seq_along(MARKERS)) {
  sd_wb <- .num(ROW[[paste0("OR_perSD_Biomarker", j)]]); dbl_wb <- .num(ROW[[paste0("OR_perDoubling_Biomarker", j)]])
  if (length(sd_wb)) .add_check(sprintf("OR per LD-SD, %s", MARKERS[j]), sd_wb, exp(B_Z[j]))
  if (length(dbl_wb)) .add_check(sprintf("OR per doubling, %s", MARKERS[j]), dbl_wb, exp(B_Z[j] / SDV[j]))
  if (length(sd_wb) && length(dbl_wb) && is.finite(sd_wb) && is.finite(dbl_wb) && sd_wb > 0 && dbl_wb > 0 && dbl_wb != 1)
    .add_check(sprintf("LD SD of log2(raw + delta), %s", MARKERS[j]), log(sd_wb) / log(dbl_wb), SDV[j])
}
if (!nrow(CHECKS)) warning("The workbook has no odds-ratio columns to check the refit against.", call. = FALSE)
.fail <- CHECKS[is.finite(CHECKS$Relative_difference) & CHECKS$Relative_difference > CHECK_TOL, , drop = FALSE]
if (nrow(.fail))
  stop("\n\nThe refit does not reproduce the workbook:\n  ",
       paste(sprintf("%s: workbook %.10g, this session %.10g", .fail$Check, .fail$Workbook, .fail$This_session), collapse = "\n  "),
       "\n  Check that the preprocessing used the run's settings (cohort, imputation, LD standardization).", call. = FALSE)

# Raw-input form: logit p = b0 + sum_j b_j * log2(raw_j + DELTA) + b_age * Age + b_bmi * BMI
B_RAW <- setNames(B_Z / SDV, MARKERS)
B0_RAW <- unname(B["(Intercept)"] - sum(B_Z * MU / SDV))
B_COV <- setNames(unname(B[CONF]), CONF)
.L <- sweep(sweep(as.matrix(DAT[MARKERS]), 2, SDV, `*`), 2, MU, `+`)         # log2(raw + DELTA) of each participant
.eta_raw <- drop(B0_RAW + .L %*% B_RAW + as.matrix(DAT[CONF]) %*% B_COV)
MAX_LINK_DIFF <- max(abs(.eta_raw - stats::predict(FIT, type = "link")))
if (!is.finite(MAX_LINK_DIFF) || MAX_LINK_DIFF > 1e-8)
  stop("The raw-input equation does not reproduce the fitted log-odds (largest difference ", signif(MAX_LINK_DIFF, 3), ").", call. = FALSE)
if (max(abs(B[-1])) > 15) warning("Very large coefficient: the estimates may be unstable (separation).", call. = FALSE)

# Detection (raw value above 0, after imputation) by group: with a biomarker detected in
# one group only, the fit shows quasi-complete separation and the size of its slope is
# set by where glm stopped iterating rather than by the data
DETECT <- do.call(rbind, lapply(MARKERS, function(m) {
  L <- DAT[[m]] * SDV[[m]] + MU[[m]]; det <- L > log2(DELTA) + 1e-9
  data.frame(Biomarker = m, Detected_ref = sum(det[DAT$condL == REF]), N_ref = N_REF,
             Detected_case = sum(det[DAT$condL == CASE]), N_case = N_CASE)
}))
DETECT$Separation <- with(DETECT, (Detected_ref > 0 & Detected_case == 0) | (Detected_case > 0 & Detected_ref == 0) |
                                  (Detected_ref == N_ref & Detected_case < N_case) | (Detected_case == N_case & Detected_ref < N_ref))
.sep_txt <- vapply(seq_len(nrow(DETECT)), function(i) with(DETECT[i, ], sprintf(
  "  %s above 0 in %d of %d %s and %d of %d %s participants%s", Biomarker, Detected_ref, N_ref, REF, Detected_case, N_case, CASE,
  if (Separation) " -> detected in one group only: quasi-complete separation (the slope's size is not estimable)" else "")), "")

# Likelihood-ratio test of the biomarker(s) given age and BMI (in-sample; it does not allow for
# the biomarker having been chosen among all candidates)
FIT0 <- stats::glm(y ~ Age_DOS + BMI, data = .d, family = stats::binomial())
LRT <- c(chisq = FIT0$deviance - FIT$deviance, df = length(MARKERS))
LRT["p"] <- stats::pchisq(LRT[["chisq"]], LRT[["df"]], lower.tail = FALSE)

# Odds ratios with 95% confidence intervals (profile likelihood; Wald where profiling fails).
# Under quasi-complete separation neither profile nor Wald limits are meaningful for that
# biomarker, so none are reported for it.
.se <- sqrt(diag(stats::vcov(FIT)))
.ci <- tryCatch(suppressMessages(suppressWarnings(stats::confint(FIT))), error = function(e) NULL)
.ci_one <- function(term) {
  lo <- hi <- NA_real_; how <- "profile likelihood"
  if (!is.null(.ci) && term %in% rownames(.ci)) { lo <- .ci[term, 1]; hi <- .ci[term, 2] }
  if (!is.finite(lo) || !is.finite(hi)) {
    how <- if (is.finite(lo) || is.finite(hi)) "profile likelihood, one limit not found (shown as NA)" else "Wald (profiling failed)"
    if (!is.finite(lo) && !is.finite(hi)) { lo <- B[term] - 1.96 * .se[term]; hi <- B[term] + 1.96 * .se[term] }
  }
  c(lo = unname(lo), hi = unname(hi), how = how)
}
OR <- do.call(rbind, c(
  lapply(seq_along(MARKERS), function(j) {
    ci <- .ci_one(.bn[j]); s <- SDV[j]
    if (isTRUE(DETECT$Separation[j])) ci <- c(lo = NA, hi = NA, how = "not estimable: quasi-complete separation")
    data.frame(Term = MARKERS[j], Unit = "per doubling of the raw value", OR = exp(B_RAW[j]),
               Lower = exp(as.numeric(ci["lo"]) / s), Upper = exp(as.numeric(ci["hi"]) / s), Method = ci["how"])
  }),
  lapply(CONF, function(v) {
    ci <- .ci_one(v)
    data.frame(Term = v, Unit = if (v == "BMI") "per kg/m2" else "per year", OR = exp(B_COV[v]),
               Lower = exp(as.numeric(ci["lo"])), Upper = exp(as.numeric(ci["hi"])), Method = ci["how"])
  })))
rownames(OR) <- NULL

# The calculator's worked example (illustrative values, not a participant)
EX_L <- log2(as.numeric(EXAMPLE_INPUT$biomarkers) + DELTA)
EX_TERMS <- c(Intercept = B0_RAW, setNames(B_RAW * EX_L, MARKERS),
              Age_DOS = B_COV[["Age_DOS"]] * EXAMPLE_INPUT$Age_DOS, BMI = B_COV[["BMI"]] * EXAMPLE_INPUT$BMI)
EX_ETA <- sum(EX_TERMS); EX_P <- stats::plogis(EX_ETA)


# ==========================================================================
# 4) OUTPUTS
# ==========================================================================
.f17 <- function(x) sprintf("%.17g", x)
.f6  <- function(x) formatC(signif(x, 6), digits = 6, format = "fg", flag = "#")    # six significant digits
.sgn <- function(x) if (x < 0) " \u2212 " else " + "
.label <- function(k) {                      # display name: italic genus and species, spaces for underscores
  if (grepl("^[A-Z][a-z]+_[a-z]+$", k)) sprintf("*%s*", sub("_", " ", k)) else gsub("_", " ", k)
}
.lsym <- if (length(MARKERS) == 1L) "L" else paste0("L", seq_along(MARKERS))
EQ_FULL <- paste0("eta = ", .f17(B0_RAW),
                  paste0(ifelse(B_RAW < 0, " - ", " + "), .f17(abs(B_RAW)), " * log2(X", if (length(MARKERS) > 1) seq_along(MARKERS) else "", " + 5e-7)", collapse = ""),
                  ifelse(B_COV[["Age_DOS"]] < 0, " - ", " + "), .f17(abs(B_COV[["Age_DOS"]])), " * Age",
                  ifelse(B_COV[["BMI"]] < 0, " - ", " + "), .f17(abs(B_COV[["BMI"]])), " * BMI")
EQ_MS <- paste0("\u03b7 = ", sub("^-", "\u2212", .f6(B0_RAW)),
                paste0(vapply(B_RAW, .sgn, ""), .f6(abs(B_RAW)), .lsym, collapse = ""),
                .sgn(B_COV[["Age_DOS"]]), .f6(abs(B_COV[["Age_DOS"]])), "Age",
                .sgn(B_COV[["BMI"]]), .f6(abs(B_COV[["BMI"]])), "BMI")
.or_txt <- function(i, d = 2) sprintf("%s (95%% CI %s to %s; %s)", formatC(OR$OR[i], digits = d, format = "f"),
                                      ifelse(is.finite(OR$Lower[i]), formatC(OR$Lower[i], digits = d, format = "f"), "NA"),
                                      ifelse(is.finite(OR$Upper[i]), formatC(OR$Upper[i], digits = d, format = "f"), "NA"), OR$Method[i])
.perf <- function(v) if (v %in% names(ROW) && is.finite(ROW[[v]])) sprintf("%.2f%%", 100 * ROW[[v]]) else "NA"
.above <- if (EX_P >= THRESHOLD) "above" else "below"

LINES <- c(
  "DMFT calculator equation: best single-biomarker model (risk of HD, DMFT > 13)",
  paste("Biomarker:", paste(MARKERS, collapse = ", ")),
  paste("Covariates: Age_DOS (years), BMI (kg/m2)"),
  paste("Source:", basename(RESULTS_FILE), "| sheet", .sheet, "| top row (balanced accuracy, then AUC)"),
  paste("Registry key:", REGISTRY_KEY),
  sprintf("Cohort: n = %d (%s %d, %s %d); LD reference statistics from %s", N, REF, N_REF, CASE, N_CASE, .stats_name),
  "",
  "EQUATION (full precision; X = raw biomarker value in the original assay units)",
  EQ_FULL,
  "p_HD = plogis(eta) = 1 / (1 + exp(-eta))",
  paste0("Classify HD when p_HD >= ", .f17(THRESHOLD), " (the panel's pooled Youden threshold in the run); otherwise LD."),
  "Apply log2(raw + 5e-7) once; do not standardize the inputs (the LD centring and scaling are in the coefficients).",
  "",
  "COEFFICIENTS",
  sprintf("  %-34s model scale %12.6g (SE %.4g) | raw-input %s", "(Intercept)", B[["(Intercept)"]], .se[["(Intercept)"]], .f17(B0_RAW)),
  vapply(seq_along(MARKERS), function(j) sprintf("  %-34s model scale %12.6g (SE %.4g, per LD-SD) | raw-input %s per log2 unit | LD mean %.6g, SD %.6g",
                                                  MARKERS[j], B_Z[j], .se[[.bn[j]]], .f17(B_RAW[j]), MU[j], SDV[j]), ""),
  vapply(CONF, function(v) sprintf("  %-34s model scale %12.6g (SE %.4g) | raw-input %s", v, B[[v]], .se[[v]], .f17(B_COV[[v]])), ""),
  "",
  "ODDS RATIOS (95% CI)",
  vapply(seq_len(nrow(OR)), function(i) sprintf("  %-34s %s: %s", OR$Term[i], OR$Unit[i], .or_txt(i, if (OR$Term[i] == "BMI") 3 else 2)), ""),
  "",
  "DETECTION BY GROUP",
  .sep_txt,
  "",
  sprintf("  Likelihood-ratio test of %s given age and BMI: chi-square %.2f on %d df, p = %.3g (in-sample; ignores the selection among candidates)",
          paste(MARKERS, collapse = " + "), LRT[["chisq"]], as.integer(LRT[["df"]]), LRT[["p"]]),
  "",
  "CHECKS",
  if (nrow(CHECKS)) sprintf("  %-55s workbook %.10g | this session %.10g | relative difference %.2g",
                            CHECKS$Check, CHECKS$Workbook, CHECKS$This_session, CHECKS$Relative_difference) else "  (no odds ratios in the workbook)",
  sprintf("  Raw-input equation vs fitted log-odds, largest difference: %.2g", MAX_LINK_DIFF),
  "",
  "APPARENT CROSS-VALIDATED PERFORMANCE (from the workbook; panel and threshold chosen on these predictions)",
  sprintf("  Balanced accuracy %s (range across repeats %s to %s), accuracy %s, sensitivity %s, specificity %s, PPV %s, NPV %s, AUC %.4f",
          .perf("Balanced_Accuracy"), .perf("Balanced_Accuracy_Repeat_Min"), .perf("Balanced_Accuracy_Repeat_Max"),
          .perf("Accuracy"), .perf("Sensitivity"), .perf("Specificity"), .perf("Precision"), .perf("NPV"), ROW$AUC),
  "",
  sprintf("CALCULATOR EXAMPLE (not a participant): %s = %s, age %g years, BMI %g kg/m2",
          paste(MARKERS, collapse = ", "), paste(EXAMPLE_INPUT$biomarkers, collapse = ", "), EXAMPLE_INPUT$Age_DOS, EXAMPLE_INPUT$BMI),
  sprintf("  log2(x + 5e-7) = %s", paste(sprintf("%.4f", EX_L), collapse = ", ")),
  sprintf("  terms: %s", paste(sprintf("%s %.4f", names(EX_TERMS), EX_TERMS), collapse = "; ")),
  sprintf("  eta = %.4f; p_HD = %.4f (%.1f%%), %s the %.1f%% threshold -> %s", EX_ETA, EX_P, 100 * EX_P, .above, 100 * THRESHOLD,
          if (EX_P >= THRESHOLD) CASE else REF),
  "",
  "FOR THE MANUSCRIPT",
  paste0("  Section 3.4 equation:  ", EQ_MS),
  sprintf("  Section 3.4 odds ratios: age %s per year; BMI %s per kg/m2; %s %s per doubling",
          formatC(OR$OR[OR$Term == "Age_DOS"], digits = 2, format = "f"), formatC(OR$OR[OR$Term == "BMI"], digits = 2, format = "f"),
          paste(MARKERS, collapse = ", "), paste(formatC(OR$OR[OR$Term %in% MARKERS], digits = 2, format = "f"), collapse = ", ")),
  sprintf("  Section 3.4 threshold: p(HD) >= %.4f (%.1f%%)", THRESHOLD, 100 * THRESHOLD),
  sprintf("  Section 3.4 detection: %s", paste(sprintf("%s above 0 in %d of %d %s and %d of %d %s participants", DETECT$Biomarker,
          DETECT$Detected_ref, DETECT$N_ref, REF, DETECT$Detected_case, DETECT$N_case, CASE), collapse = "; ")),
  sprintf("  Section 3.7 example: biomarker term %.1f, age term %.1f, BMI term %.1f, intercept %.1f, linear predictor %.3f, risk %.1f%% (%s the %.1f%% threshold)",
          sum(EX_TERMS[MARKERS]), EX_TERMS[["Age_DOS"]], EX_TERMS[["BMI"]], B0_RAW, EX_ETA, 100 * EX_P, .above, 100 * THRESHOLD),
  "",
  "This refit uses all participants and adds no validation; the panel and threshold were chosen on the development data.")
writeLines(LINES, file.path(OUT_DIR, "DMFT_equation_single_biomarker.txt"), useBytes = TRUE)

# Calculator model file (same columns and metadata rows as script 04's dmft_age_bmi_model.csv)
.terms <- c("(Intercept)", MARKERS, CONF)
TAB <- data.frame(Term = .terms,
                  Label = c("Intercept", vapply(MARKERS, .label, ""), "Age", "Body mass index"),
                  Estimate_model_scale = c(B[["(Intercept)"]], unname(B_Z), B[["Age_DOS"]], B[["BMI"]]),
                  Coefficient_raw_input = c(B0_RAW, unname(B_RAW), B_COV[["Age_DOS"]], B_COV[["BMI"]]),
                  LD_mean_log2 = c(NA, unname(MU), NA, NA), LD_SD_log2 = c(NA, unname(SDV), NA, NA),
                  Input = c("", rep(sprintf("log2(raw value + %s)", format(DELTA)), length(MARKERS)), rep("value in original units", 2)),
                  stringsAsFactors = FALSE)
TAB <- rbind(TAB, data.frame(
  Term = c("(Threshold)", "(Delta)", "(Case group)", "(Reference group)", "(Model configuration)", "(Search design)"),
  Label = c("Classification threshold", "log2 offset", "high DMFT (DMFT > 13)", "low DMFT (DMFT <= 13)", "Biomarkers_Age_BMI",
            sub("^.*_(exhaustive[^_]*_beam[^_]*)$", "\\1", REGISTRY_KEY)),
  Estimate_model_scale = NA_real_, Coefficient_raw_input = c(THRESHOLD, DELTA, NA, NA, NA, NA),
  LD_mean_log2 = NA_real_, LD_SD_log2 = NA_real_, Input = "", stringsAsFactors = FALSE))
.old <- options(digits = 17)
utils::write.csv(TAB, file.path(OUT_DIR, "dmft_age_bmi_model.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
options(.old)

# PUBLISHED_MODEL and EXAMPLE blocks for app/app.R (replace the existing ones)
.cvm <- regmatches(.pget("Cross-validation"), regexec("([0-9]+)-fold CV x ([0-9]+) repeat", .pget("Cross-validation")))[[1]]
.cv_text <- if (length(.cvm) == 3L) sprintf("%s-fold cross-validation repeated %s times", .cvm[2], .cvm[3]) else "Repeated cross-validation"
.abbr <- vapply(MARKERS, function(k) paste(substr(strsplit(k, "_")[[1]][1:2], 1, 1), collapse = ""), "")
.q <- function(x) paste0("\"", x, "\"")
.vec <- function(x, f = .f17) paste0("c(", paste(f(x), collapse = ", "), ")")
APP_BLOCK <- c(
  "# ---- generated by 04b_Export_DMFT_Equation_Single_Biomarker.R: replace PUBLISHED_MODEL and EXAMPLE in app/app.R ----",
  "PUBLISHED_MODEL <- list(",
  "  source = \"published\", name = \"Study DMFT model, trained on the 111 study participants\",",
  sprintf("  description = \"%s with age and BMI (logistic regression).\",",
          if (length(MARKERS) == 1L) "One salivary biomarker" else sprintf("%d salivary biomarkers", length(MARKERS))),
  "  case_label = \"high DMFT (DMFT > 13)\", ref_label = \"low DMFT (DMFT \\u2264 13)\",",
  sprintf("  case_short = %s, ref_short = %s, transform = \"log2\", delta = %s,", .q(CASE), .q(REF), format(DELTA)),
  sprintf("  analysis_configuration = \"Biomarkers_Age_BMI\", search_design = %s,", .q(TAB$Label[TAB$Term == "(Search design)"])),
  sprintf("  intercept = %s,", .f17(B0_RAW)),
  "  markers = data.frame(",
  sprintf("    key = c(%s),", paste(.q(MARKERS), collapse = ", ")),
  sprintf("    abbr = c(%s),", paste(.q(toupper(.abbr)), collapse = ", ")),
  sprintf("    label = c(%s),", paste(.q(vapply(MARKERS, .label, "")), collapse = ", ")),
  sprintf("    coef = %s, stringsAsFactors = FALSE),", .vec(unname(B_RAW))),
  "  covariates = data.frame(var = c(\"Age_DOS\", \"BMI\"), term = c(\"Age_DOS\", \"BMI\"),",
  "    abbr = c(\"Age\", \"BMI\"), label = c(\"Age\", \"Body mass index\"), unit = c(\"years\", \"kg/m\\u00b2\"),",
  "    type = c(\"numeric\", \"numeric\"), level = c(NA_character_, NA_character_),",
  sprintf("    coef = %s, stringsAsFactors = FALSE),", .vec(unname(B_COV))),
  sprintf("  cov_levels = list(), level_labels = list(), threshold = %s,", .f17(THRESHOLD)),
  sprintf("  cohort = \"%d adults: %d with DMFT \\u2264 13 and %d with DMFT > 13.\",", N, N_REF, N_CASE),
  sprintf("  performance = c(\"Balanced accuracy\" = %s, \"Sensitivity\" = %s,", .f17(ROW$Balanced_Accuracy), .f17(ROW$Sensitivity)),
  sprintf("                  \"Specificity\" = %s, \"AUC\" = %s),", .f17(ROW$Specificity), .f17(ROW$AUC)),
  sprintf("  performance_note = paste(\"%s. Panel and threshold were chosen on these\",", .cv_text),
  "    \"predictions, so development performance may be optimistic; not externally validated.\"),",
  "  caution = \"Development model; the biomarker coefficient is imprecise and independent validation is required. Age is in years and BMI in kg/m\\u00b2.\"",
  ")",
  "",
  sprintf("EXAMPLE <- list(markers = %s, covariates = c(%s, %s))   # illustrative, not a participant",
          .vec(as.numeric(EXAMPLE_INPUT$biomarkers), function(x) format(x)), format(EXAMPLE_INPUT$Age_DOS), format(EXAMPLE_INPUT$BMI)))
writeLines(APP_BLOCK, file.path(OUT_DIR, "PUBLISHED_MODEL_single_biomarker.R"), useBytes = TRUE)

DMFT_SINGLE_EQUATION <- list(fit = FIT, markers = MARKERS, covariates = CONF, coefficients = TAB, odds_ratios = OR,
                             checks = CHECKS, detection = DETECT, lrt = LRT, threshold = THRESHOLD, delta = DELTA, reference_mean = MU, reference_sd = SDV,
                             equation_full = EQ_FULL, equation_manuscript = EQ_MS, example = list(input = EXAMPLE_INPUT,
                             terms = EX_TERMS, eta = EX_ETA, p = EX_P), max_link_difference = MAX_LINK_DIFF,
                             workbook_row = ROW, registry_key = REGISTRY_KEY, session_info = utils::sessionInfo())
saveRDS(DMFT_SINGLE_EQUATION, file.path(OUT_DIR, "DMFT_equation_single_biomarker.rds"))
cat(paste(LINES, collapse = "\n"), "\n\nSaved in: ", normalizePath(OUT_DIR), "\n", sep = "")
