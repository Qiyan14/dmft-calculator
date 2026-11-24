# ==========================================================================
# 04b_Export_DMFT_Equation_Single_Biomarker.V9.R
# ==========================================================================
# Exports the calculator equation for the best SINGLE-biomarker model of the
# run, for one covariate set (MODEL_SET): "Biomarkers_Age" (one biomarker with
# age; the calculator's model since manuscript 11.20 V5) or
# "Biomarkers_Age_BMI" (one biomarker with age and BMI). Under nested
# cross-validation the single biomarker generalized better than the larger
# panels, so the manuscript (Section 3.4) and the calculator use it.
#
# V9: MODEL_SET (section 1). The covariates come from the workbook (age, or
# age and BMI), and the equation, odds ratios, worked example, model file and
# app block follow them. Outputs for Biomarkers_Age_BMI go to
# DMFT_single_biomarker_equation, as before (model file
# dmft_age_bmi_model.csv); for another set to
# DMFT_single_biomarker_equation_<MODEL_SET> (for Biomarkers_Age the model
# file is dmft_age_model.csv). With MODEL_SET <- "Biomarkers_Age_BMI" the
# outputs are those of V8.
#
# HOW TO RUN
#   1. In a fresh R session, run the preprocessing script with the settings of
#      the run (it creates data.dfPH, stats_PH and registry_key).
#   2. Put the run's results workbook for MODEL_SET in the working directory
#      (RESULTS_FILE = NULL finds it from registry_key), or give its full path
#      in RESULTS_FILE below.
#   3. source("04b_Export_DMFT_Equation_Single_Biomarker.V9.R")
#
# WHAT IT DOES
#   - takes the top row of the Results_1_Biomarkers sheet (the best single
#     biomarker, Bacteroides heparinolyticus in this run) and its threshold;
#   - refits the logistic regression of HD on that biomarker and the set's
#     covariates (Age_DOS, or Age_DOS and BMI) in the run's cohort (glm, like
#     the search's final fits) and checks the refit against the odds ratios in
#     the workbook;
#   - converts the coefficients to raw-input form (the LD-group centring and
#     scaling are absorbed into the intercept and slope) and checks that the
#     raw-input equation reproduces the fitted log-odds;
#   - computes odds ratios with 95% profile-likelihood confidence intervals,
#     the calculator's worked example, and the values for the manuscript.
#
# OUTPUTS (in OUT_DIR)
#   DMFT_equation_single_biomarker.txt   equation, odds ratios, threshold, checks, example,
#                                        and the values to paste into Sections 3.4 and 3.7
#   dmft_age_model.csv                   the model as a calculator file (Patient prediction ->
#     (dmft_age_bmi_model.csv with BMI)  "Saved model"); same format as script 04's export
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
MODEL_SET     <- "Biomarkers_Age"       # covariate set: "Biomarkers_Age" (age) or "Biomarkers_Age_BMI" (age and BMI)
# NULL = the workbook script 2 wrote for this session's settings:
# Results_LogReg_HD_vs_LD_<MODEL_SET>_<registry_key>.xlsx
RESULTS_FILE  <- NULL
OUT_DIR       <- NULL                   # NULL = DMFT_single_biomarker_equation for Biomarkers_Age_BMI (as V8),
                                        # DMFT_single_biomarker_equation_<MODEL_SET> for any other set
PANEL_SIZE    <- 1L                     # the best panel of this size is exported (1 = best single biomarker)
EXAMPLE_INPUT <- list(biomarkers = 0,   # the calculator's example values (one per biomarker), not a participant;
                      Age_DOS = 40, BMI = 25)   # a covariate that is not in the model is not used
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
if (!is.character(MODEL_SET) || length(MODEL_SET) != 1 || is.na(MODEL_SET) || !nzchar(MODEL_SET))
  stop("MODEL_SET: the name of one covariate set, e.g. \"Biomarkers_Age\".", call. = FALSE)
if (is.null(RESULTS_FILE)) {
  if (!exists("registry_key"))
    stop("registry_key is not in this R session: run the preprocessing first, or give RESULTS_FILE.", call. = FALSE)
  RESULTS_FILE <- sprintf("Results_LogReg_HD_vs_LD_%s_%s.xlsx", MODEL_SET, registry_key)
}
if (!file.exists(RESULTS_FILE))
  stop(sprintf("Workbook not found: %s\n  Working directory: %s\n  Put it there or give its full path in RESULTS_FILE.",
               RESULTS_FILE, getwd()), call. = FALSE)
if (is.null(OUT_DIR))
  OUT_DIR <- file.path(getwd(), if (identical(MODEL_SET, "Biomarkers_Age_BMI")) "DMFT_single_biomarker_equation" else
    sprintf("DMFT_single_biomarker_equation_%s", MODEL_SET))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.P <- as.data.frame(readxl::read_excel(RESULTS_FILE, sheet = "Analysis_Parameters", col_types = "text"))
.P <- data.frame(Parameter = trimws(.P[[1]]), Value = trimws(.P[[2]]), stringsAsFactors = FALSE)
.pget <- function(name) {
  v <- .P$Value[which(.P$Parameter == name)]
  if (length(v) && !is.na(v[1])) v[1] else NA_character_
}
.stopm <- character(0)
if (!identical(.pget("Model configuration"), MODEL_SET))
  .stopm <- c(.stopm, sprintf("RESULTS_FILE should be the %s workbook (its model configuration: %s).",
                              MODEL_SET, .pget("Model configuration")))
REGISTRY_KEY <- .pget("Script 1 registry key")
CONF <- trimws(strsplit(.pget("Confounders"), ",", fixed = TRUE)[[1]])
CONF <- CONF[!is.na(CONF) & nzchar(CONF) & CONF != "none"]
.COV_KNOWN <- c("Age_DOS", "BMI")       # the covariates the calculator takes, in its order
if (!length(CONF) || length(setdiff(CONF, .COV_KNOWN)))
  .stopm <- c(.stopm, sprintf("Expected a model adjusted for Age_DOS and/or BMI; the workbook's confounders: %s.",
                              .pget("Confounders")))
CONF <- intersect(.COV_KNOWN, CONF)
# Each covariate's symbol in the equation, wording, label, unit and odds-ratio unit
.COV_SYM   <- c(Age_DOS = "Age", BMI = "BMI")
.COV_WORD  <- c(Age_DOS = "age", BMI = "BMI")
.COV_LABEL <- c(Age_DOS = "Age", BMI = "Body mass index")
.COV_UNIT  <- c(Age_DOS = "years", BMI = "kg/m2")
.COV_PER   <- c(Age_DOS = "per year", BMI = "per kg/m2")
.COV_APPU  <- c(Age_DOS = "years", BMI = "kg/m\\u00b2")   # app.R code: the escape for the superscript 2
COV_AND <- if (length(CONF) <= 1) unname(.COV_WORD[CONF]) else
  paste(paste(.COV_WORD[CONF][-length(CONF)], collapse = ", "), "and", .COV_WORD[CONF][length(CONF)])
REF  <- .pget("Reference group"); CASE <- .pget("Case group (positive)")
if (is.na(REF) || is.na(CASE)) { REF <- "LD"; CASE <- "HD" }
# LD / HD labels from the run's cut-off: the workbook's "DMFT cut-off" row, e.g.
# "<= 13.9 -> LD (low DMFT) | > 13.9 -> HD (high DMFT) (fixed)"
.cut <- suppressWarnings(as.numeric(sub("^<=\\s*([-+0-9.eE]+)\\s*->.*$", "\\1", .pget("DMFT cut-off"))))
if (!is.finite(.cut)) .cut <- get0("DMFT_CUTOFF_USED", ifnotfound = NA_real_)
if (length(.cut) != 1 || !is.numeric(.cut) || !is.finite(.cut))
  stop("The run's DMFT cut-off is neither in RESULTS_FILE (Analysis_Parameters) nor in this session.", call. = FALSE)
RULE_HD     <- sprintf("DMFT > %g", .cut)
RULE_LD     <- sprintf("DMFT <= %g", .cut)
RULE_LD_APP <- sprintf("DMFT \\u2264 %g", .cut)   # app.R code: the escape for "less than or equal"
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
.ex_bad <- CONF[!vapply(CONF, function(v) is.numeric(EXAMPLE_INPUT[[v]]) && length(EXAMPLE_INPUT[[v]]) == 1 &&
                          is.finite(EXAMPLE_INPUT[[v]]), logical(1))]
if (length(.ex_bad))
  stop("EXAMPLE_INPUT needs one number for each covariate of the model: ", paste(.ex_bad, collapse = ", "), call. = FALSE)

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
.d <- data.frame(y = DAT$condL, setNames(DAT[MARKERS], paste0("B", seq_along(MARKERS))), DAT[CONF])
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

# Raw-input form: logit p = b0 + sum_j b_j * log2(raw_j + DELTA) + b_age * Age (+ b_bmi * BMI)
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

# Likelihood-ratio test of the biomarker(s) given the covariates (in-sample; it does not allow for
# the biomarker having been chosen among all candidates)
FIT0 <- stats::glm(stats::reformulate(CONF, response = "y"), data = .d, family = stats::binomial())
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
    data.frame(Term = v, Unit = .COV_PER[[v]], OR = exp(B_COV[v]),
               Lower = exp(as.numeric(ci["lo"])), Upper = exp(as.numeric(ci["hi"])), Method = ci["how"])
  })))
rownames(OR) <- NULL

# The calculator's worked example (illustrative values, not a participant)
EX_L <- log2(as.numeric(EXAMPLE_INPUT$biomarkers) + DELTA)
EX_TERMS <- c(Intercept = B0_RAW, setNames(B_RAW * EX_L, MARKERS),
              setNames(vapply(CONF, function(v) B_COV[[v]] * EXAMPLE_INPUT[[v]], numeric(1)), CONF))
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
                  paste0(ifelse(B_COV[CONF] < 0, " - ", " + "), .f17(abs(B_COV[CONF])), " * ", .COV_SYM[CONF], collapse = ""))
EQ_MS <- paste0("\u03b7 = ", sub("^-", "\u2212", .f6(B0_RAW)),
                paste0(vapply(B_RAW, .sgn, ""), .f6(abs(B_RAW)), .lsym, collapse = ""),
                paste0(vapply(B_COV[CONF], .sgn, ""), .f6(abs(B_COV[CONF])), .COV_SYM[CONF], collapse = ""))
.or_txt <- function(i, d = 2) sprintf("%s (95%% CI %s to %s; %s)", formatC(OR$OR[i], digits = d, format = "f"),
                                      ifelse(is.finite(OR$Lower[i]), formatC(OR$Lower[i], digits = d, format = "f"), "NA"),
                                      ifelse(is.finite(OR$Upper[i]), formatC(OR$Upper[i], digits = d, format = "f"), "NA"), OR$Method[i])
.perf <- function(v) if (v %in% names(ROW) && is.finite(ROW[[v]])) sprintf("%.2f%%", 100 * ROW[[v]]) else "NA"
.above <- if (EX_P >= THRESHOLD) "above" else "below"

LINES <- c(
  sprintf("DMFT calculator equation: best single-biomarker model (risk of HD, %s)", RULE_HD),
  paste("Biomarker:", paste(MARKERS, collapse = ", ")),
  paste("Covariates:", paste(sprintf("%s (%s)", CONF, .COV_UNIT[CONF]), collapse = ", ")),
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
  sprintf("  Likelihood-ratio test of %s given %s: chi-square %.2f on %d df, p = %.3g (in-sample; ignores the selection among candidates)",
          paste(MARKERS, collapse = " + "), COV_AND, LRT[["chisq"]], as.integer(LRT[["df"]]), LRT[["p"]]),
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
  sprintf("CALCULATOR EXAMPLE (not a participant): %s = %s%s",
          paste(MARKERS, collapse = ", "), paste(EXAMPLE_INPUT$biomarkers, collapse = ", "),
          paste0(", ", .COV_WORD[CONF], " ", vapply(CONF, function(v) sprintf("%g", EXAMPLE_INPUT[[v]]), ""), " ",
                 .COV_UNIT[CONF], collapse = "")),
  sprintf("  log2(x + 5e-7) = %s", paste(sprintf("%.4f", EX_L), collapse = ", ")),
  sprintf("  terms: %s", paste(sprintf("%s %.4f", names(EX_TERMS), EX_TERMS), collapse = "; ")),
  sprintf("  eta = %.4f; p_HD = %.4f (%.1f%%), %s the %.1f%% threshold -> %s", EX_ETA, EX_P, 100 * EX_P, .above, 100 * THRESHOLD,
          if (EX_P >= THRESHOLD) CASE else REF),
  "",
  "FOR THE MANUSCRIPT",
  paste0("  Section 3.4 equation:  ", EQ_MS),
  sprintf("  Section 3.4 odds ratios: %s; %s %s per doubling",
          paste(sprintf("%s %s %s", .COV_WORD[CONF], formatC(OR$OR[match(CONF, OR$Term)], digits = 2, format = "f"),
                        .COV_PER[CONF]), collapse = "; "),
          paste(MARKERS, collapse = ", "), paste(formatC(OR$OR[OR$Term %in% MARKERS], digits = 2, format = "f"), collapse = ", ")),
  sprintf("  Section 3.4 threshold: p(HD) >= %.4f (%.1f%%)", THRESHOLD, 100 * THRESHOLD),
  sprintf("  Section 3.4 detection: %s", paste(sprintf("%s above 0 in %d of %d %s and %d of %d %s participants", DETECT$Biomarker,
          DETECT$Detected_ref, DETECT$N_ref, REF, DETECT$Detected_case, DETECT$N_case, CASE), collapse = "; ")),
  sprintf("  Section 3.7 example: biomarker term %.1f%s, intercept %.1f, linear predictor %.3f, risk %.1f%% (%s the %.1f%% threshold)",
          sum(EX_TERMS[MARKERS]), paste0(", ", .COV_WORD[CONF], " term ", sprintf("%.1f", EX_TERMS[CONF]), collapse = ""),
          B0_RAW, EX_ETA, 100 * EX_P, .above, 100 * THRESHOLD),
  "",
  "This refit uses all participants and adds no validation; the panel and threshold were chosen on the development data.")
writeLines(LINES, file.path(OUT_DIR, "DMFT_equation_single_biomarker.txt"), useBytes = TRUE)

# Calculator model file (same columns and metadata rows as script 04's dmft_age_bmi_model.csv)
.terms <- c("(Intercept)", MARKERS, CONF)
.na_cov <- rep(NA_real_, length(CONF))
TAB <- data.frame(Term = .terms,
                  Label = unname(c("Intercept", vapply(MARKERS, .label, ""), .COV_LABEL[CONF])),
                  Estimate_model_scale = unname(c(B[["(Intercept)"]], B_Z, B[CONF])),
                  Coefficient_raw_input = unname(c(B0_RAW, B_RAW, B_COV[CONF])),
                  LD_mean_log2 = c(NA, unname(MU), .na_cov), LD_SD_log2 = c(NA, unname(SDV), .na_cov),
                  Input = c("", rep(sprintf("log2(raw value + %s)", format(DELTA)), length(MARKERS)),
                            rep("value in original units", length(CONF))),
                  stringsAsFactors = FALSE)
TAB <- rbind(TAB, data.frame(
  Term = c("(Threshold)", "(Delta)", "(Case group)", "(Reference group)", "(Model configuration)", "(Search design)"),
  Label = c("Classification threshold", "log2 offset", sprintf("high DMFT (%s)", RULE_HD), sprintf("low DMFT (%s)", RULE_LD),
            MODEL_SET,
            sub("^.*_(exhaustive[^_]*_beam[^_]*)$", "\\1", REGISTRY_KEY)),
  Estimate_model_scale = NA_real_, Coefficient_raw_input = c(THRESHOLD, DELTA, NA, NA, NA, NA),
  LD_mean_log2 = NA_real_, LD_SD_log2 = NA_real_, Input = "", stringsAsFactors = FALSE))
.old <- options(digits = 17)
MODEL_CSV <- sprintf("dmft_%s_model.csv", paste(tolower(.COV_SYM[CONF]), collapse = "_"))   # dmft_age_model.csv, dmft_age_bmi_model.csv
utils::write.csv(TAB, file.path(OUT_DIR, MODEL_CSV), row.names = FALSE, na = "", fileEncoding = "UTF-8")
options(.old)

# PUBLISHED_MODEL and EXAMPLE blocks for app/app.R (replace the existing ones)
.cvm <- regmatches(.pget("Cross-validation"), regexec("([0-9]+)-fold CV x ([0-9]+) repeat", .pget("Cross-validation")))[[1]]
.cv_text <- if (length(.cvm) == 3L) sprintf("%s-fold cross-validation repeated %s times", .cvm[2], .cvm[3]) else "Repeated cross-validation"
.abbr <- vapply(MARKERS, function(k) paste(substr(strsplit(k, "_")[[1]][1:2], 1, 1), collapse = ""), "")
.q <- function(x) paste0("\"", x, "\"")
.cvec <- function(x) paste0("c(", paste(.q(unname(x)), collapse = ", "), ")")
.units_sentence <- paste0(.COV_SYM[CONF][1], " is in ", .COV_APPU[CONF][1],
                          if (length(CONF) > 1) paste0(" and ", paste(sprintf("%s in %s", .COV_SYM[CONF][-1], .COV_APPU[CONF][-1]),
                                                                      collapse = " and ")) else "")
.vec <- function(x, f = .f17) paste0("c(", paste(f(x), collapse = ", "), ")")
APP_BLOCK <- c(
  "# ---- generated by 04b_Export_DMFT_Equation_Single_Biomarker.R: replace PUBLISHED_MODEL and EXAMPLE in app/app.R ----",
  "PUBLISHED_MODEL <- list(",
  "  source = \"published\", name = \"Study DMFT model, trained on the 111 study participants\",",
  sprintf("  description = \"%s with %s (logistic regression).\",",
          if (length(MARKERS) == 1L) "One salivary biomarker" else sprintf("%d salivary biomarkers", length(MARKERS)), COV_AND),
  sprintf("  case_label = \"high DMFT (%s)\", ref_label = \"low DMFT (%s)\",", RULE_HD, RULE_LD_APP),
  sprintf("  case_short = %s, ref_short = %s, transform = \"log2\", delta = %s,", .q(CASE), .q(REF), format(DELTA)),
  sprintf("  analysis_configuration = %s, search_design = %s,", .q(MODEL_SET), .q(TAB$Label[TAB$Term == "(Search design)"])),
  sprintf("  intercept = %s,", .f17(B0_RAW)),
  "  markers = data.frame(",
  sprintf("    key = c(%s),", paste(.q(MARKERS), collapse = ", ")),
  sprintf("    abbr = c(%s),", paste(.q(toupper(.abbr)), collapse = ", ")),
  sprintf("    label = c(%s),", paste(.q(vapply(MARKERS, .label, "")), collapse = ", ")),
  sprintf("    coef = %s, stringsAsFactors = FALSE),", .vec(unname(B_RAW))),
  sprintf("  covariates = data.frame(var = %s, term = %s,", .cvec(CONF), .cvec(CONF)),
  sprintf("    abbr = %s, label = %s, unit = %s,", .cvec(.COV_SYM[CONF]), .cvec(.COV_LABEL[CONF]), .cvec(.COV_APPU[CONF])),
  sprintf("    type = %s, level = c(%s),", .cvec(rep("numeric", length(CONF))),
          paste(rep("NA_character_", length(CONF)), collapse = ", ")),
  sprintf("    coef = %s, stringsAsFactors = FALSE),", .vec(unname(B_COV[CONF]))),
  sprintf("  cov_levels = list(), level_labels = list(), threshold = %s,", .f17(THRESHOLD)),
  sprintf("  cohort = \"%d adults: %d with %s and %d with %s.\",", N, N_REF, RULE_LD_APP, N_CASE, RULE_HD),
  sprintf("  performance = c(\"Balanced accuracy\" = %s, \"Sensitivity\" = %s,", .f17(ROW$Balanced_Accuracy), .f17(ROW$Sensitivity)),
  sprintf("                  \"Specificity\" = %s, \"AUC\" = %s),", .f17(ROW$Specificity), .f17(ROW$AUC)),
  sprintf("  performance_note = paste(\"%s. Panel and threshold were chosen on these\",", .cv_text),
  "    \"predictions, so development performance may be optimistic; not externally validated.\"),",
  sprintf("  caution = \"Development model; the biomarker coefficient is imprecise and independent validation is required. %s.\"",
          .units_sentence),
  ")",
  "",
  sprintf("EXAMPLE <- list(markers = %s, covariates = c(%s))   # illustrative, not a participant",
          .vec(as.numeric(EXAMPLE_INPUT$biomarkers), function(x) format(x)),
          paste(vapply(CONF, function(v) format(EXAMPLE_INPUT[[v]]), ""), collapse = ", ")))
writeLines(APP_BLOCK, file.path(OUT_DIR, "PUBLISHED_MODEL_single_biomarker.R"), useBytes = TRUE)

DMFT_SINGLE_EQUATION <- list(fit = FIT, model_set = MODEL_SET, markers = MARKERS, covariates = CONF, coefficients = TAB, odds_ratios = OR,
                             checks = CHECKS, detection = DETECT, lrt = LRT, threshold = THRESHOLD, delta = DELTA, reference_mean = MU, reference_sd = SDV,
                             equation_full = EQ_FULL, equation_manuscript = EQ_MS, example = list(input = EXAMPLE_INPUT,
                             terms = EX_TERMS, eta = EX_ETA, p = EX_P), max_link_difference = MAX_LINK_DIFF,
                             workbook_row = ROW, registry_key = REGISTRY_KEY, session_info = utils::sessionInfo())
saveRDS(DMFT_SINGLE_EQUATION, file.path(OUT_DIR, "DMFT_equation_single_biomarker.rds"))
cat(paste(LINES, collapse = "\n"), "\n\nSaved in: ", normalizePath(OUT_DIR), "\n", sep = "")
