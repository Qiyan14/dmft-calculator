# ==========================================================================
# CALIBRATION AND FIGURE 1A-B FOR THE DMFT RUN (one covariate set: MODEL_SET)
# 05b_Calibration_10rep_DMFT.V11.R
# ==========================================================================
# Run in the same R session as the preprocessing (script 01, with the run's
# settings), after it:
#
#   source("05b_Calibration_10rep_DMFT.V11.R")
#
# V10: MODEL_SET (section 1) names the covariate set: "Biomarkers_Age" (the
# panel search with age only, the calculator's model since manuscript 11.20
# V5) or "Biomarkers_Age_BMI" (with age and BMI, as V9). Its covariates-only
# workbook (Age, or Age_BMI_Only) is found among the run's workbooks: the one
# with the same covariates. Outputs for Biomarkers_Age_BMI go to
# DMFT_10rep_outputs, as before; for another set to
# DMFT_10rep_outputs_<MODEL_SET>, so the earlier files stay as they are.
# With MODEL_SET <- "Biomarkers_Age_BMI" the results are those of V9.
# V11: the calibration intercept is the root of its score equation, which
# always exists and is found reliably, instead of glm.fit's iterations: with
# some predictions at 0 or 1 these diverged in some bootstrap resamples and
# gave intervals such as -9e14 to 0.30 (the four-marker panel of the age
# set); the calibration slope falls back to a guarded Newton fit when glm.fit
# does not converge. Estimates are unchanged (to 1e-8); intervals change
# where resamples had failed. New setting CAL_SEED: the bootstrap seed
# (default 123, the seed of the manuscript's intervals; NULL = script 01's
# SEED, as V10).
#
# Why a separate script: the 10-repeat run (LogReg_Fast.R) reports the
# predictions of its own fitting routine, refitted in every fold (fastglm or
# IRLS, glm.fit as the fallback), not caret's. 03_Figures and 05_Calibration
# refit through caret and demand an exact match, so they stop on this run's
# workbook. This script rebuilds the run's out-of-fold predictions with the
# run's own routine, on the run's own folds, and checks them against the
# workbooks before using them.
#
# Inputs (settings below)
#   data.dfPH       from the preprocessing (same settings as the run)
#   RESULTS_FILE    the run's workbook for MODEL_SET (Results_1..6 sheets)
#   BASELINE_FILE   the run's covariates-only workbook with the same covariates
#                   (Age for Biomarkers_Age, Age_BMI_Only for
#                   Biomarkers_Age_BMI); its CV_Predictions sheet holds the
#                   run's folds and that model's predictions
# Steps
#   1. Cohort: participants complete for the run's covariates (Age_DOS, BMI,
#      Sex); n, LD and HD must equal the workbook's.
#   2. Folds: read from CV_Predictions; the covariates-only model is
#      refitted on them and must reproduce the run's predictions.
#   3. Best panel at each size (top row of each Results sheet): pooled
#      out-of-fold predictions with each fitting routine; the routine that
#      reproduces the workbook (balanced accuracy, accuracy, sensitivity,
#      specificity, PPV, NPV, threshold, AUC) is used.
#   4. Calibration intercept and slope, Brier and scaled Brier scores, with
#      AUC, balanced accuracy, sensitivity, specificity, accuracy, PPV and NPV
#      at the panel's pooled Youden threshold: 95% CIs from 2,000
#      participant-level bootstrap resamples (within HD and LD, predictions
#      held fixed), and the lowest and highest repeat-specific values.
#   5. Figure 1A (ROC curves, 1-6 biomarkers) and 1B (calibration of the
#      selected panel), the prediction file and ready-to-paste text.
#   6. Nested cross-validation (section 7), when script 09's workbook for this
#      run is in the working directory: the same measures and figures from
#      09's held-out predictions (each outer fold's panel, fit and threshold
#      chosen without its held-out participants), checked against 09's own
#      numbers, and figures that set the repeated and nested results side by side.
# Outputs (OUT_DIR; <tag> = the number of repeats, e.g. 5rep)
#   Figure_1AB_ROC_Calibration_<tag>.png/.pdf   Figure 1, panels A and B
#   Figure_ROC_<tag>, Figure_Calibration_<tag>  the two panels on their own
#   S2_OOF_predictions_<tag>.csv      out-of-fold predictions, all models
#   S3_Calibration_metrics_<tag>.csv  all models x measures, CIs, repeat range
#   Performance_<tag>.csv             Table 2 values with PPV, NPV, repeat range
#   Calibration_text_<tag>.txt        ready-to-paste sentences
#   Refit_check_<tag>.csv, Fold_diagnostics_<tag>.csv, Run_info_<tag>.txt
#   Nested cross-validation (with 09's workbook):
#   Figure_1_Repeated_Nested_<tag>    ROC curves (A repeated, B nested) and
#                                     calibration (C best repeated, D best nested)
#   Figure_ROC_nested_<tag>, Figure_Calibration_nested_<tag>   nested panels on their own
#   Figure_BA_by_size_<tag>           balanced accuracy by panel size, repeated and nested
#   S4_Nested_calibration_metrics_<tag>.csv   nested models x measures, bootstrap CIs
#   Nested_performance_<tag>.csv      one row per nested model (Table 2 and 3 values)
#   Nested_text_<tag>.txt             ready-to-paste sentences
# Packages: readxl. fastglm is used when installed (as by LogReg_Fast);
# caret only if the folds have to be regenerated. About 1-2 minutes.
# ==========================================================================

# ==========================================================================
# 1) SETTINGS
# ==========================================================================

# Covariate set (an add_set() name of script 01): "Biomarkers_Age" (panels
# with age; the calculator's model) or "Biomarkers_Age_BMI" (with age and BMI)
MODEL_SET     <- "Biomarkers_Age"

# The run's two workbooks. NULL = the ones script 2 wrote for this session's
# settings (script 1's registry_key), in the working directory:
#   Results_LogReg_HD_vs_LD_<MODEL_SET>_<registry_key>.xlsx
#   the covariates-only workbook of the run with the same covariates, e.g.
#   Results_LogReg_HD_vs_LD_Age_<registry_key>.xlsx for Biomarkers_Age or
#   Results_LogReg_HD_vs_LD_Age_BMI_Only_<registry_key>.xlsx for Biomarkers_Age_BMI
# Or give the file names (full paths if they are not in the working directory).
RESULTS_FILE  <- NULL
BASELINE_FILE <- NULL
FOLD_FILE     <- NULL     # optional: the run's DMFT_Resampling_<registry key>.rds (cross-checked)
OUT_DIR       <- NULL     # NULL = DMFT_10rep_outputs for Biomarkers_Age_BMI (as V9),
                          # DMFT_10rep_outputs_<MODEL_SET> for any other set
FIG_FORMATS   <- c("png", "pdf")   # any of "png", "pdf", "tiff"
PNG_DPI       <- 300

# Fitting routine. "auto" tries LogReg_Fast's routines ("fastglm" when that
# package is installed, "irls", "glmfit") and caret's ("glm"), and keeps the
# one that reproduces the workbook most closely. Or name one of them.
ENGINE        <- "auto"
MATCH_TOL     <- 1e-6     # largest difference that counts as reproducing the workbook
ACCEPT_TOL    <- 0.005    # stop when no routine comes closer than this
REQUIRE_SAME_REGISTRY <- TRUE   # stop when registry_key in this session is not the run's

# Marker names on figures and in tables
ITALIC_TAXA    <- TRUE    # genus and species in italics (figures)
NAME_OVERRIDES <- c()     # e.g. c(Moraxella_Moraxella_unclassified = "*Moraxella* sp. (unclassified)")

# Calibration
CAL_BOOT      <- 2000     # participant-level bootstrap resamples
CAL_SEED      <- 123      # bootstrap seed (123: the manuscript's intervals); NULL = script 01's SEED
CAL_SPAN      <- 0.75     # loess span of the calibration curves
CAL_BINS      <- 10       # risk groups for the observed proportions
CAL_EPS       <- 1e-12    # predictions clipped to [CAL_EPS, 1 - CAL_EPS] before the logit
CAL_HIST_BINS <- 40       # bins of the histogram under the calibration plot
CAL_BAND      <- "range"  # shaded band: "range" (lowest to highest repeat curve) or "p95" (2.5th-97.5th percentiles)
CAL_K         <- NULL     # panel size in Figure 1B; NULL = best overall panel
S2_INCLUDE_STUDYID <- FALSE  # TRUE adds StudyID to the prediction file

# Nested cross-validation (section 7): script 09's workbook for this run.
# NULL = DMFT_Nested_CV_<registry_key>.xlsx in the working directory; when it
# is not there, section 7 is skipped and everything above is unchanged.
NESTED_FILE   <- NULL
NESTED_SET    <- NULL     # covariate set in 09's workbook; NULL = MODEL_SET
NESTED_CAL_K  <- NULL     # panel size of the nested calibration figure; NULL = best size under nested CV
NESTED_MATCH_TOL <- 1e-6  # largest difference from 09's own numbers that counts as a match


# ==========================================================================
# 2) SHARED SET-UP: THE RUN'S COHORT, FOLDS AND PANELS, REFITS CHECKED
#    (as in 08b_Variable_Importance_10rep_DMFT, which still expects age and
#    BMI; V10 adds MODEL_SET)
# ==========================================================================

# --- 2a) Prerequisites ------------------------------------------------------
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package readxl is needed: install.packages(\"readxl\")", call. = FALSE)
if (!exists("data.dfPH"))
  stop("data.dfPH is not in this R session. Run the preprocessing script (PROJECT = \"DMFT\", ",
       "with the settings of the 10-repeat run) in this session first.", call. = FALSE)
if (!is.null(CAL_SEED) && (length(CAL_SEED) != 1 || is.na(CAL_SEED) || CAL_SEED != round(CAL_SEED)))
  stop("CAL_SEED: a whole number, or NULL for script 01's SEED.", call. = FALSE)
if (!exists("SEED"))
  stop("SEED is not in this R session. Run the preprocessing script (script 01, which sets SEED) ",
       "in this session first.", call. = FALSE)
if (!is.character(MODEL_SET) || length(MODEL_SET) != 1 || is.na(MODEL_SET) || !nzchar(MODEL_SET))
  stop("MODEL_SET: the name of one covariate set, e.g. \"Biomarkers_Age\".", call. = FALSE)
# Workbook not named in the settings: the one of this session's settings
if (is.null(RESULTS_FILE)) {
  if (!exists("registry_key"))
    stop("registry_key is not in this R session: run the preprocessing first, or give RESULTS_FILE and BASELINE_FILE.",
         call. = FALSE)
  RESULTS_FILE <- sprintf("Results_LogReg_HD_vs_LD_%s_%s.xlsx", MODEL_SET, registry_key)
}
.no_file <- function(f) stop(sprintf("Workbook not found: %s\n  Working directory: %s\n  Put it there or give its full path in the settings.",
                                     f, getwd()), call. = FALSE)
if (!file.exists(RESULTS_FILE)) .no_file(RESULTS_FILE)
if (is.null(OUT_DIR))
  OUT_DIR <- file.path(getwd(), if (identical(MODEL_SET, "Biomarkers_Age_BMI")) "DMFT_10rep_outputs" else
    sprintf("DMFT_10rep_outputs_%s", MODEL_SET))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- 2b) The workbooks --------------------------------------------------------
.read_params <- function(path) {
  p <- as.data.frame(readxl::read_excel(path, sheet = "Analysis_Parameters", col_types = "text"))
  data.frame(Parameter = trimws(p[[1]]), Value = trimws(p[[2]]), stringsAsFactors = FALSE)
}
.pget <- function(p, name) {
  v <- p$Value[which(p$Parameter == name)]
  if (length(v) && !is.na(v[1])) v[1] else NA_character_
}
.P  <- .read_params(RESULTS_FILE)
.covs_of <- function(p) {
  v <- .pget(p, "Confounders")
  if (is.na(v) || !nzchar(v) || v == "none") character(0) else trimws(strsplit(v, ",", fixed = TRUE)[[1]])
}
if (!length(.covs_of(.P)))
  stop(sprintf(paste("%s has no covariates; this script needs a set with covariates (their covariates-only",
                     "model holds the run's folds)."), .pget(.P, "Model configuration")), call. = FALSE)
# Covariates-only workbook not named in the settings: the run's workbook without
# biomarkers whose covariates are MODEL_SET's (script 02 pairs the sets the same
# way), among the workbooks next to RESULTS_FILE
if (is.null(BASELINE_FILE)) {
  .key  <- .pget(.P, "Script 1 registry key")
  .cand <- list.files(dirname(RESULTS_FILE), pattern = "^Results_LogReg_HD_vs_LD_.*\\.xlsx$", full.names = TRUE)
  .cand <- .cand[endsWith(basename(.cand), sprintf("_%s.xlsx", .key)) &
                 normalizePath(.cand) != normalizePath(RESULTS_FILE)]
  .cand <- .cand[vapply(.cand, function(f) {
    sh <- tryCatch(readxl::excel_sheets(f), error = function(e) character(0))
    if (!"Results_0_Biomarkers" %in% sh || any(grepl("^Results_[1-9][0-9]*_Biomarkers$", sh))) return(FALSE)
    p <- tryCatch(.read_params(f), error = function(e) NULL)
    !is.null(p) && setequal(.covs_of(p), .covs_of(.P))
  }, logical(1))]
  if (length(.cand) != 1)
    stop(sprintf("%s covariates-only workbook of this run with the covariates of %s (%s) in %s%s.\n  Give it in BASELINE_FILE%s.",
                 if (length(.cand)) "More than one" else "No", MODEL_SET, paste(.covs_of(.P), collapse = ", "),
                 normalizePath(dirname(RESULTS_FILE)),
                 if (length(.cand)) paste0(": ", paste(basename(.cand), collapse = ", ")) else "",
                 if (length(.cand)) "" else
                   " (script 02 writes it when script 01 has a covariates-only add_set() line with those covariates)"),
         call. = FALSE)
  BASELINE_FILE <- .cand
}
if (!file.exists(BASELINE_FILE)) .no_file(BASELINE_FILE)
.P0 <- .read_params(BASELINE_FILE)
.stopm <- character(0)
if (!identical(.pget(.P, "Model configuration"), MODEL_SET))
  .stopm <- c(.stopm, sprintf("RESULTS_FILE should be the %s workbook (its model configuration: %s).",
                              MODEL_SET, .pget(.P, "Model configuration")))
for (.nm in c("Script 1 registry key", "Total n", "n (LD)", "n (HD)", "Cross-validation", "Cohort policy"))
  if (!identical(.pget(.P, .nm), .pget(.P0, .nm)))
    .stopm <- c(.stopm, sprintf("The two workbooks are from different runs (%s: %s vs %s).",
                                .nm, .pget(.P, .nm), .pget(.P0, .nm)))
if (!setequal(.covs_of(.P), .covs_of(.P0)))
  .stopm <- c(.stopm, sprintf("BASELINE_FILE's covariates (%s) are not those of RESULTS_FILE (%s).",
                              .pget(.P0, "Confounders"), .pget(.P, "Confounders")))
.cvs <- .pget(.P, "Cross-validation")
.cvm <- regmatches(.cvs, regexec("([0-9]+)-fold CV x ([0-9]+) repeats?, seed ([0-9]+)", .cvs))[[1]]
if (length(.cvm) != 4L) stop("Cannot read the cross-validation settings in RESULTS_FILE: ", .cvs, call. = FALSE)
CV_FOLDS   <- as.integer(.cvm[2])
CV_REPEATS <- as.integer(.cvm[3])
CV_SEED    <- as.integer(.cvm[4])
BOOT_SEED  <- if (is.null(CAL_SEED)) as.integer(SEED) else as.integer(CAL_SEED)
if (!identical(CV_SEED, as.integer(SEED)))
  warning(sprintf(paste("The run's folds were made with seed %d, but script 01's SEED is now %d. Run script 02",
                        "again to use one seed throughout (the bootstrap uses CAL_SEED)."), CV_SEED, as.integer(SEED)),
          call. = FALSE, immediate. = TRUE)
TAG <- sprintf("%drep", CV_REPEATS)
REGISTRY_KEY <- .pget(.P, "Script 1 registry key")
CONF <- .covs_of(.P)
REF  <- .pget(.P, "Reference group")
CASE <- .pget(.P, "Case group (positive)")
if (is.na(REF) || is.na(CASE)) { REF <- "LD"; CASE <- "HD" }
# The covariates in words, for labels and text: "age + BMI" / "age and BMI", or "age"
.COV_WORD <- c(Age_DOS = "age", BMI = "BMI", Sex = "sex")
COV_WORDS <- unname(ifelse(CONF %in% names(.COV_WORD), .COV_WORD[CONF], CONF))
COV_PLUS  <- paste(COV_WORDS, collapse = " + ")
COV_AND   <- if (length(COV_WORDS) <= 1) paste(COV_WORDS, collapse = "") else
  paste(paste(COV_WORDS[-length(COV_WORDS)], collapse = ", "), "and", COV_WORDS[length(COV_WORDS)])
COV_PLUS_CAP <- paste0(toupper(substring(COV_PLUS, 1, 1)), substring(COV_PLUS, 2))
COV_ONLY  <- paste(COV_PLUS_CAP, "only")   # "Age + BMI only", "Age only"
if (exists("registry_key") && !identical(get("registry_key"), REGISTRY_KEY)) {
  .msg <- sprintf("registry_key in this session (%s) is not the run's (%s).", get("registry_key"), REGISTRY_KEY)
  if (isTRUE(REQUIRE_SAME_REGISTRY)) {
    .stopm <- c(.stopm, paste(.msg, "Re-run the preprocessing with the run's settings, or set REQUIRE_SAME_REGISTRY <- FALSE."))
  } else warning(.msg, call. = FALSE)
}

# Best panel at each size: the top row of each Results_<k>_Biomarkers sheet
# (sorted again by balanced accuracy, then AUC)
.METRICS  <- c("Balanced_Accuracy", "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "Threshold", "AUC")
.NUM_COLS <- c(.METRICS, "Balanced_Accuracy_Repeat_Min", "Balanced_Accuracy_Repeat_Max")
.read_sheet <- function(path, sheet) {
  d <- suppressWarnings(as.data.frame(readxl::read_excel(path, sheet = sheet)))
  for (v in intersect(.NUM_COLS, names(d))) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d
}
.sheets <- grep("^Results_[1-9][0-9]*_Biomarkers$", readxl::excel_sheets(RESULTS_FILE), value = TRUE)
if (!length(.sheets)) stop("No Results_<k>_Biomarkers sheets in RESULTS_FILE.", call. = FALSE)
PANELS <- lapply(.sheets, function(s) {
  k <- as.integer(sub("^Results_([0-9]+)_Biomarkers$", "\\1", s))
  d <- .read_sheet(RESULTS_FILE, s)
  miss <- setdiff(c(paste0("Biomarker", seq_len(k)), .NUM_COLS), names(d))
  if (length(miss)) stop(sprintf("Sheet %s lacks column(s): %s", s, paste(miss, collapse = ", ")), call. = FALSE)
  d <- d[order(-d$Balanced_Accuracy, -d$AUC, na.last = TRUE), , drop = FALSE]
  list(k = k, markers = as.character(unlist(d[1, paste0("Biomarker", seq_len(k))], use.names = FALSE)),
       row = d[1, .NUM_COLS, drop = FALSE])
})
names(PANELS) <- vapply(PANELS, function(x) as.character(x$k), "")
PANELS <- PANELS[order(as.integer(names(PANELS)))]
.ks <- as.integer(names(PANELS))
BEST_K <- .ks[order(-vapply(PANELS, function(x) x$row$Balanced_Accuracy, numeric(1)),
                    -vapply(PANELS, function(x) x$row$AUC, numeric(1)), .ks)[1]]

.sh0 <- readxl::excel_sheets(BASELINE_FILE)
if (!"Results_0_Biomarkers" %in% .sh0) stop("BASELINE_FILE has no Results_0_Biomarkers sheet.", call. = FALSE)
BASE_ROW <- .read_sheet(BASELINE_FILE, "Results_0_Biomarkers")
if (length(setdiff(.NUM_COLS, names(BASE_ROW))))
  stop("Results_0_Biomarkers lacks: ", paste(setdiff(.NUM_COLS, names(BASE_ROW)), collapse = ", "), call. = FALSE)
BASE_ROW <- BASE_ROW[1, .NUM_COLS, drop = FALSE]
CVP <- if ("CV_Predictions" %in% .sh0) as.data.frame(readxl::read_excel(BASELINE_FILE, sheet = "CV_Predictions")) else NULL
if (!is.null(CVP) && !all(c("rowIndex", "Resample", "obs", CASE) %in% names(CVP))) CVP <- NULL

# --- 2c) The cohort -------------------------------------------------------------
# Participants complete for every covariate named in the run's cohort policy,
# in data.dfPH's row order (the row numbers of the run's folds).
.policy <- .pget(.P, "Cohort policy")
COHORT_VARS <- if (!is.na(.policy) && grepl("^Participants complete for ", .policy)) {
  trimws(strsplit(sub("\\s*\\(.*$", "", sub("^Participants complete for ", "", .policy)), ",", fixed = TRUE)[[1]])
} else character(0)
COHORT_VARS <- unique(c(COHORT_VARS, CONF))
DAT <- as.data.frame(data.dfPH)
.all_markers <- unique(unlist(lapply(PANELS, `[[`, "markers"), use.names = FALSE))
.miss <- setdiff(c("condL", COHORT_VARS, .all_markers), names(DAT))
if (length(.miss))
  stop("Not in data.dfPH: ", paste(.miss, collapse = ", "),
       ".\n  Run the preprocessing with the settings of the run (registry key ", REGISTRY_KEY, ").", call. = FALSE)
.ok <- stats::complete.cases(DAT[, COHORT_VARS, drop = FALSE]) & !is.na(DAT$condL)
SOURCE_ROWS <- which(.ok)
DAT <- DAT[.ok, , drop = FALSE]
rownames(DAT) <- NULL
DAT$condL <- factor(as.character(DAT$condL), levels = c(REF, CASE))
if (anyNA(DAT$condL)) stop(sprintf("condL in data.dfPH has values other than %s and %s.", REF, CASE), call. = FALSE)
N <- nrow(DAT)
Y <- as.integer(DAT$condL == CASE)
.n_run <- suppressWarnings(as.integer(c(.pget(.P, "Total n"), .pget(.P, "n (LD)"), .pget(.P, "n (HD)"))))
if (!identical(c(N, sum(Y == 0L), sum(Y == 1L)), .n_run))
  .stopm <- c(.stopm, sprintf("Cohort in this session: n = %d (%s %d, %s %d); the run: n = %s (%s %s, %s %s), complete for %s.",
                              N, REF, sum(Y == 0L), CASE, sum(Y == 1L), .n_run[1], REF, .n_run[2], CASE, .n_run[3],
                              paste(COHORT_VARS, collapse = ", ")))
.bad <- c(.all_markers, CONF)[!vapply(DAT[c(.all_markers, CONF)], function(z) is.numeric(z) && all(is.finite(z)), logical(1))]
if (length(.bad)) .stopm <- c(.stopm, paste("Not numeric, or missing values, in the cohort:", paste(.bad, collapse = ", ")))
if (length(.stopm))
  stop("\n\n===== THIS SESSION DOES NOT MATCH THE RUN =====\n  ", paste(.stopm, collapse = "\n  "), "\n", call. = FALSE)

# --- 2d) The run's folds ------------------------------------------------------------
# From the CV_Predictions sheet (caret's folds from the run's seed: the held-out rows
# of each FoldX.RepYY), else FOLD_FILE, else regenerated with caret. Ordered
# as caret lists them (Fold1.Rep01, ..., Fold5.Rep01, Fold1.Rep02, ...).
.fold_order <- function(nm) order(as.integer(sub("^.*\\.Rep", "", nm)), as.integer(sub("^Fold([0-9]+)\\..*$", "\\1", nm)))
.check_index <- function(index) {
  nm <- names(index)
  if (is.null(nm) || any(!grepl("^Fold[0-9]+\\.Rep[0-9]+$", nm)) || anyDuplicated(nm))
    stop("Fold names must look like Fold1.Rep01.", call. = FALSE)
  if (length(index) != CV_FOLDS * CV_REPEATS)
    stop(sprintf("%d folds found; the run used %d-fold CV x %d repeats.", length(index), CV_FOLDS, CV_REPEATS), call. = FALSE)
  reps <- sub("^.*\\.", "", nm)
  for (r in unique(reps)) {
    te <- unlist(lapply(index[reps == r], function(tr) setdiff(seq_len(N), tr)), use.names = FALSE)
    if (!identical(sort(as.integer(te)), seq_len(N)))
      stop("Each participant must be held out exactly once per repeat (", r, ").", call. = FALSE)
  }
  invisible(TRUE)
}
INDEX <- NULL
FOLD_SOURCE <- NA_character_
if (!is.null(CVP)) {
  CVP$rowIndex <- as.integer(CVP$rowIndex)
  CVP$Resample <- as.character(CVP$Resample)
  if (anyNA(CVP$rowIndex) || any(CVP$rowIndex < 1L | CVP$rowIndex > N))
    stop("CV_Predictions has row indices outside 1..", N, ".", call. = FALSE)
  if (any(as.character(CVP$obs) != as.character(DAT$condL[CVP$rowIndex])))
    stop("The outcomes in CV_Predictions do not match condL in this session: the cohort or its row order differs from the run.",
         call. = FALSE)
  .res <- unique(CVP$Resample)
  .res <- .res[.fold_order(.res)]
  INDEX <- setNames(lapply(.res, function(r) setdiff(seq_len(N), CVP$rowIndex[CVP$Resample == r])), .res)
  FOLD_SOURCE <- sprintf("CV_Predictions sheet of %s", basename(BASELINE_FILE))
}
if (!is.null(FOLD_FILE)) {
  if (!file.exists(FOLD_FILE)) stop("FOLD_FILE not found: ", FOLD_FILE, call. = FALSE)
  .sv <- readRDS(FOLD_FILE)
  .ix <- if (is.list(.sv) && !is.null(.sv$index)) .sv$index else .sv
  if (!is.list(.ix) || is.null(names(.ix))) stop("FOLD_FILE does not hold a named list of training rows.", call. = FALSE)
  .ix <- lapply(.ix[.fold_order(names(.ix))], function(v) sort(as.integer(v)))
  if (!is.null(.sv$source_rows) && !identical(as.integer(.sv$source_rows), SOURCE_ROWS))
    stop("FOLD_FILE was made for other rows of data.dfPH.", call. = FALSE)
  if (is.null(INDEX)) {
    INDEX <- .ix
    FOLD_SOURCE <- basename(FOLD_FILE)
  } else if (!identical(names(.ix), names(INDEX)) || !all(mapply(identical, .ix, INDEX))) {
    stop("FOLD_FILE and the CV_Predictions sheet hold different folds; they are not from the same run.", call. = FALSE)
  } else FOLD_SOURCE <- sprintf("%s (identical in %s)", FOLD_SOURCE, basename(FOLD_FILE))
}
if (is.null(INDEX)) {
  if (!requireNamespace("caret", quietly = TRUE))
    stop("No folds: BASELINE_FILE has no CV_Predictions sheet and FOLD_FILE is not set. ",
         "Install caret to regenerate them.", call. = FALSE)
  .summ <- function(data, lev = NULL, model = NULL) {
    ts <- caret::twoClassSummary(data, lev, model)
    ds <- caret::defaultSummary(data, lev, model)
    c(ROC = ts[["ROC"]], Sens = ts[["Sens"]], Spec = ts[["Spec"]], Accuracy = ds[["Accuracy"]],
      Kappa = ds[["Kappa"]], BalancedAccuracy = mean(c(ts[["Sens"]], ts[["Spec"]])))
  }
  .ctrl <- caret::trainControl(method = "repeatedcv", number = CV_FOLDS, repeats = CV_REPEATS,
                               summaryFunction = .summ, classProbs = TRUE, allowParallel = FALSE,
                               savePredictions = "final")
  .d <- data.frame(condL = DAT$condL, z = sin(seq_len(N)))
  levels(.d$condL) <- make.names(levels(.d$condL))
  set.seed(CV_SEED)
  .fit <- suppressWarnings(caret::train(condL ~ z, data = .d, method = "glm", family = "binomial",
                                        metric = "BalancedAccuracy", trControl = .ctrl))
  INDEX <- lapply(.fit$control$index[.fold_order(names(.fit$control$index))], function(v) sort(as.integer(v)))
  FOLD_SOURCE <- sprintf("regenerated with caret::train after set.seed(%d)", CV_SEED)
}
.check_index(INDEX)
FOLDS <- lapply(INDEX, function(tr) list(tr = as.integer(tr), te = setdiff(seq_len(N), tr)))

# --- 2e) Fitting routines -------------------------------------------------------------
# .fit_glmfit, .fit_fastglm and .fit_irls are LogReg_Fast.R's own (its
# "sweep"): each fold is refitted, warm-started from the previous fold's
# converged coefficients; a fit that does not finish cleanly (no
# convergence, |linear predictor| > 30, separation) is redone with glm.fit.
# .fit_glm is caret's fit: glm() with its defaults, from glm's own start.
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
.fit_glm <- function(X, y) {
  f <- tryCatch(suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(), control = stats::glm.control())),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  b <- f$coefficients
  b[is.na(b)] <- 0        # predict.glm leaves aliased columns out
  list(coef = b, conv = isTRUE(f$converged))
}
.linkinv <- stats::binomial()$linkinv   # as predict.glm: saturated predictions tie as in caret

# Design: intercept, markers, then the covariates (caret's column order)
.Z <- stats::model.matrix(stats::reformulate(CONF), data = DAT)[, -1, drop = FALSE]
.design <- function(markers) {
  X <- cbind(1, as.matrix(DAT[, markers, drop = FALSE]), .Z)
  storage.mode(X) <- "double"
  colnames(X) <- c("(Intercept)", markers, colnames(.Z))
  X
}

# Pooled out-of-fold predictions over all folds of all repeats, in fold order
.oof <- function(X, y, folds, engine) {
  nf <- length(folds)
  nte <- vapply(folds, function(f) length(f$te), integer(1))
  p <- numeric(sum(nte)); row <- integer(sum(nte)); fold <- integer(sum(nte))
  coef <- matrix(NA_real_, nf, ncol(X), dimnames = list(names(folds), colnames(X)))
  conv <- logical(nf); eta_tr <- numeric(nf)
  b <- NULL; pos <- 0L
  for (f in seq_len(nf)) {
    tr <- folds[[f]]$tr; te <- folds[[f]]$te
    Xtr <- X[tr, , drop = FALSE]; ytr <- y[tr]
    fit <- switch(engine,
                  fastglm = .fit_fastglm(Xtr, ytr, start = b),
                  irls    = .fit_irls(Xtr, ytr, start = b),
                  glmfit  = .fit_glmfit(Xtr, ytr, start = b),
                  glm     = .fit_glm(Xtr, ytr),
                  stop("Unknown fitting routine: ", engine))
    if (is.null(fit)) return(NULL)
    b <- if (fit$conv) fit$coef else NULL
    ix <- pos + seq_along(te); pos <- pos + length(te)
    p[ix] <- .linkinv(drop(X[te, , drop = FALSE] %*% fit$coef))
    row[ix] <- te; fold[ix] <- f
    coef[f, ] <- fit$coef; conv[f] <- fit$conv
    eta_tr[f] <- max(abs(Xtr %*% fit$coef))
  }
  list(p = p, row = row, fold = fold, resample = names(folds)[fold], coef = coef,
       converged = conv, max_abs_eta_train = eta_tr)
}

# Pooled metrics as LogReg_Fast computes them: AUC under the tie-collapsed
# ROC; cut-offs at midpoints between distinct predictions (HD when p >= cut);
# the cut with the highest sensitivity + specificity, ties to the one closest
# to 0.5 (equally close: the higher)
.roc_stats <- function(p, y, np, nn) {
  if (anyNA(p)) return(rep(NA_real_, 7L))
  o  <- sort.list(p, method = "radix", decreasing = TRUE)
  ps <- p[o]; ys <- y[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  m    <- length(ps)
  keep <- c(ps[-1L] != ps[-m], TRUE)
  tp   <- tp[keep]; fp <- fp[keep]; s <- ps[keep]
  nk   <- length(tp)
  auc  <- sum((c(0, tp[-nk]) + tp) * (fp - c(0, fp[-nk]))) / (2 * np * nn)
  tp  <- c(0, tp); fp <- c(0, fp)
  thr <- c(Inf, if (nk > 1) (s[-nk] + s[-1L]) / 2, -Inf)
  sens <- tp / np
  spec <- 1 - fp / nn
  J    <- sens + spec
  tied <- which(J >= max(J) - 0.5 / (np * nn))
  j    <- tied[order(abs(thr[tied] - 0.5), -thr[tied])[1]]
  TP <- tp[j]; FP <- fp[j]; FN <- np - TP; TN <- nn - FP
  c((TP + TN) / (np + nn), sens[j], spec[j],
    if (TP + FP > 0) TP / (TP + FP) else NA_real_,
    if (TN + FN > 0) TN / (TN + FN) else NA_real_,
    thr[j], auc)
}
.perf <- function(p, y) {
  s <- .roc_stats(p, y, sum(y == 1L), sum(y == 0L))
  setNames(c((s[2] + s[3]) / 2, s), .METRICS)
}

# --- 2f) Labels -------------------------------------------------------------------------
# As in 03_Figures: Actinomyces_sp._HMT_169 -> *Actinomyces* sp. HMT 169;
# *...* marks italics (genus and species), higher ranks stay upright.
.CYTOKINE_LABELS <- c(IL6 = "IL-6", MMP9 = "MMP-9", S100A8 = "S100A8", MMP8 = "MMP-8",
                      IL1B_IL1F2 = "IL-1\u03b2", CCL3_MIP1a = "MIP-1\u03b1", BAFF = "BAFF",
                      Resistin = "Resistin", IFNa = "IFN-\u03b1", TIMP1 = "TIMP-1", PGE2 = "PGE2",
                      Adiponectin = "Adiponectin")
.COVARIATE_LABELS <- c(Age_DOS = "Age", Teeth = "Teeth", BMI = "BMI", Sex = "Sex",
                       Cigsmoke = "Smoking", Diabetes = "Diabetes")
.HIGHER_TAXA <- c("Bacteria", "Bacteroidetes", "Bacteroidota", "Firmicutes", "Bacillota",
                  "Proteobacteria", "Pseudomonadota", "Actinobacteria", "Actinomycetota",
                  "Fusobacteria", "Fusobacteriota", "Spirochaetes", "Spirochaetota", "Synergistetes",
                  "Synergistota", "Saccharibacteria", "Absconditabacteria", "Gracilibacteria",
                  "Chloroflexi", "Tenericutes", "Mollicutes", "Clostridia", "Bacilli", "Negativicutes",
                  "Erysipelotrichia", "Bacteroidia", "Flavobacteriia", "Coriobacteriia",
                  "Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria",
                  "Epsilonproteobacteria", "Deltaproteobacteria", "Campylobacterota", "TM7", "SR1", "GN02")
.is_higher_taxon <- function(tok) grepl("(aceae|ales|ineae|idae|etes|ota)$", tok) || tok %in% .HIGHER_TAXA
marker_label <- function(name) {
  if (length(NAME_OVERRIDES) && name %in% names(NAME_OVERRIDES)) return(NAME_OVERRIDES[[name]])
  if (name %in% names(.CYTOKINE_LABELS))  return(.CYTOKINE_LABELS[[name]])
  if (name %in% names(.COVARIATE_LABELS)) return(.COVARIATE_LABELS[[name]])
  if (startsWith(name, "ratio__")) {
    parts <- strsplit(sub("^ratio__", "", name), "__", fixed = TRUE)[[1]]
    if (length(parts) == 2) return(paste(marker_label(parts[1]), "/", marker_label(parts[2])))
  }
  if (grepl("_div_TIMP1$", name)) return(paste0(marker_label(sub("_div_TIMP1$", "", name)), "/TIMP-1"))
  tok <- strsplit(name, "_", fixed = TRUE)[[1]]
  tok <- tok[nzchar(tok)]
  if (!length(tok)) return(name)
  head    <- tok[1]
  it_head <- ITALIC_TAXA && grepl("^[A-Z][a-z]+$", head) && !.is_higher_taxon(head)
  fhead   <- if (it_head) paste0("*", head, "*") else head
  if ("unclassified" %in% tolower(tok)) return(paste("unclassified", fhead))
  out <- fhead; i <- 2L; epithet_ok <- TRUE
  while (i <= length(tok)) {
    t <- tok[i]
    if (grepl("^G[0-9]+$", t)) {
      out <- c(out, sprintf("[G-%s]", sub("^G", "", t))); epithet_ok <- FALSE
    } else if (t == "HMT" && i < length(tok) && grepl("^[0-9]+$", tok[i + 1])) {
      out <- c(out, paste("HMT", tok[i + 1])); i <- i + 1L; epithet_ok <- FALSE
    } else if (t %in% c("sp.", "sp", "spp.", "bacterium", "clone", "oral", "taxon", "genomosp.", "cf.")) {
      out <- c(out, if (t == "sp") "sp." else t); epithet_ok <- FALSE
    } else if (t %in% c("subsp.", "subsp")) {
      out <- c(out, "subsp."); epithet_ok <- TRUE
    } else if (it_head && epithet_ok && grepl("^[a-z]+$", t)) {
      out <- c(out, paste0("*", t, "*"))
    } else {
      out <- c(out, t)
    }
    i <- i + 1L
  }
  paste(out, collapse = " ")
}
label_plain <- function(s) gsub("*", "", s, fixed = TRUE)
.label_k <- function(k) if (k == 0) COV_ONLY else sprintf("%d biomarker%s + %s", k, if (k == 1) "" else "s", COV_PLUS)

# --- 2g) Reproduce the run ---------------------------------------------------------------
cat(sprintf("\n===== DMFT %d-repeat run: %s =====\n", CV_REPEATS, basename(RESULTS_FILE)))
cat(sprintf("  Cohort: n = %d (%s %d, %s %d); %d-fold CV x %d repeats; folds: %s\n",
            N, REF, sum(Y == 0L), CASE, sum(Y == 1L), CV_FOLDS, CV_REPEATS, FOLD_SOURCE))

# Covariates only (age, or age + BMI): refitted on the folds, against the run's own predictions
.b0 <- .oof(.design(character(0)), Y, FOLDS, "glm")
if (!is.null(CVP)) {
  .m0 <- match(paste(CVP$rowIndex, CVP$Resample), paste(.b0$row, .b0$resample))
  BASE_PRED_DIFF <- if (anyNA(.m0)) Inf else max(abs(.b0$p[.m0] - as.numeric(CVP[[CASE]])))
  if (!is.finite(BASE_PRED_DIFF) || BASE_PRED_DIFF > 1e-6)
    stop(sprintf(paste0("The refitted %s predictions differ from the run's (largest difference %.2g): ",
                        "data.dfPH (rows, row order, or %s values) is not the run's data."),
                 COV_PLUS, BASE_PRED_DIFF, paste(CONF, collapse = "/")), call. = FALSE)
  BASE_DATA <- data.frame(id = CVP$rowIndex, resample = CVP$Resample, p = as.numeric(CVP[[CASE]]),
                          stringsAsFactors = FALSE)
  cat(sprintf("  %s predictions reproduced (largest difference %.1e)\n", COV_PLUS_CAP, BASE_PRED_DIFF))
} else {
  BASE_PRED_DIFF <- NA_real_
  BASE_DATA <- data.frame(id = .b0$row, resample = .b0$resample, p = .b0$p, stringsAsFactors = FALSE)
}
.bp <- .perf(BASE_DATA$p, Y[BASE_DATA$id])
.bd <- abs(.bp - unlist(BASE_ROW[1, .METRICS]))
BASE_METRIC_DIFF <- max(.bd, na.rm = TRUE)
if (BASE_METRIC_DIFF > MATCH_TOL)
  stop(sprintf("%s metrics differ from Results_0_Biomarkers (largest difference %.2g).", COV_PLUS_CAP, BASE_METRIC_DIFF),
       call. = FALSE)

# Panels: every available routine; keep the one that reproduces the workbook
.avail <- c(if (requireNamespace("fastglm", quietly = TRUE)) "fastglm", "irls", "glmfit", "glm")
.ENGINES <- if (identical(ENGINE, "auto")) .avail else {
  if (!ENGINE %in% c("fastglm", "irls", "glmfit", "glm"))
    stop("ENGINE must be \"auto\", \"fastglm\", \"irls\", \"glmfit\" or \"glm\".", call. = FALSE)
  if (ENGINE == "fastglm" && !"fastglm" %in% .avail)
    stop("ENGINE = \"fastglm\", but the fastglm package is not installed.", call. = FALSE)
  ENGINE
}
.FITS <- list(); .rows <- list()
for (.e in .ENGINES) {
  .FITS[[.e]] <- list()
  for (.kk in names(PANELS)) {
    .pn <- PANELS[[.kk]]
    .X  <- .design(.pn$markers)
    .o  <- tryCatch(.oof(.X, Y, FOLDS, .e), error = function(err) NULL)
    .pf <- if (is.null(.o)) setNames(rep(NA_real_, length(.METRICS)), .METRICS) else .perf(.o$p, Y[.o$row])
    .wb <- unlist(.pn$row[1, .METRICS])
    .d  <- abs(.pf - .wb)
    .d[is.na(.pf) & is.na(.wb)] <- 0
    .rows[[length(.rows) + 1L]] <- data.frame(
      Engine = .e, Panel_size = .pn$k, Panel = paste(.pn$markers, collapse = "; "),
      as.list(setNames(.pf, paste0(.METRICS, "_refit"))), as.list(setNames(.wb, paste0(.METRICS, "_workbook"))),
      Max_abs_difference = if (anyNA(.d)) Inf else max(.d), check.names = FALSE, stringsAsFactors = FALSE)
    .FITS[[.e]][[.kk]] <- list(X = .X, oof = .o, perf = .pf)
  }
}
REFIT_CHECK <- do.call(rbind, .rows)
.maxd <- vapply(.ENGINES, function(e) max(REFIT_CHECK$Max_abs_difference[REFIT_CHECK$Engine == e]), numeric(1))
ENGINE_USED <- .ENGINES[which.min(.maxd)]   # closest; ties go to the first (LogReg_Fast's default order)
REFIT_MAXDIFF <- .maxd[[ENGINE_USED]]
REFIT_CHECK$Used <- REFIT_CHECK$Engine == ENGINE_USED
utils::write.csv(REFIT_CHECK, file.path(OUT_DIR, sprintf("Refit_check_%s.csv", TAG)), row.names = FALSE, na = "")
cat("  Largest difference from the workbook over all panels and metrics, by fitting routine:\n")
for (.e in .ENGINES) cat(sprintf("    %-8s %.2g%s\n", .e, .maxd[[.e]], if (.e == ENGINE_USED) "   <- used" else ""))
if (REFIT_MAXDIFF > ACCEPT_TOL)
  stop(sprintf(paste0("No fitting routine reproduces the workbook (closest: %s, largest difference %.3g). ",
                      "See Refit_check_%s.csv. The data in this session are probably not those of the run."),
               ENGINE_USED, REFIT_MAXDIFF, TAG), call. = FALSE)
if (REFIT_MAXDIFF > MATCH_TOL)
  warning(sprintf(paste0("The refits come close to the workbook but do not reproduce it exactly (routine %s, largest ",
                         "difference %.2g; see Refit_check_%s.csv). Outputs use these refits."),
                  ENGINE_USED, REFIT_MAXDIFF, TAG), call. = FALSE, immediate. = TRUE)

# One entry per model: the out-of-fold predictions (id = row of the cohort)
.frame <- function(id, resample, p) {
  d <- data.frame(id = as.integer(id), resample = as.character(resample), y = Y[as.integer(id)], p = as.numeric(p),
                  stringsAsFactors = FALSE)
  d$rep <- sub("^.*\\.", "", d$resample)
  d
}
MODELS <- list(`0` = list(k = 0L, markers = character(0), label = .label_k(0),
                          panel = paste(unname(ifelse(CONF %in% names(.COVARIATE_LABELS), .COVARIATE_LABELS[CONF], CONF)),
                                        collapse = "; "),
                          data = .frame(BASE_DATA$id, BASE_DATA$resample, BASE_DATA$p), workbook = BASE_ROW,
                          X = .design(character(0)), coef = .b0$coef, converged = .b0$converged,
                          max_abs_eta = .b0$max_abs_eta_train, engine = "glm"))
for (.kk in names(PANELS)) {
  .f <- .FITS[[ENGINE_USED]][[.kk]]
  MODELS[[.kk]] <- list(k = PANELS[[.kk]]$k, markers = PANELS[[.kk]]$markers, label = .label_k(PANELS[[.kk]]$k),
                        panel = paste(label_plain(vapply(PANELS[[.kk]]$markers, marker_label, "")), collapse = "; "),
                        data = .frame(.f$oof$row, .f$oof$resample, .f$oof$p), workbook = PANELS[[.kk]]$row,
                        X = .f$X, coef = .f$oof$coef, converged = .f$oof$converged,
                        max_abs_eta = .f$oof$max_abs_eta_train, engine = ENGINE_USED)
}
for (.kk in names(MODELS)) MODELS[[.kk]]$threshold <- .perf(MODELS[[.kk]]$data$p, MODELS[[.kk]]$data$y)[["Threshold"]]

FOLD_DIAGNOSTICS <- do.call(rbind, lapply(MODELS, function(m) data.frame(
  Model = m$label, Panel_size = m$k, Routine = m$engine, Folds = length(m$converged),
  Fits_not_converged = sum(!m$converged),
  Max_abs_slope = if (ncol(m$coef) > 1) max(abs(m$coef[, -1, drop = FALSE])) else NA_real_,
  Max_abs_linear_predictor_training = max(m$max_abs_eta),
  Folds_with_training_linear_predictor_over_15 = sum(m$max_abs_eta > 15), stringsAsFactors = FALSE)))
utils::write.csv(FOLD_DIAGNOSTICS, file.path(OUT_DIR, sprintf("Fold_diagnostics_%s.csv", TAG)), row.names = FALSE)

.pkgv <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"
RUN_INFO <- c(
  sprintf("Run (registry key): %s", REGISTRY_KEY),
  sprintf("Covariate set: %s (covariates %s); covariates-only model: %s", MODEL_SET, paste(CONF, collapse = ", "),
          .pget(.P0, "Model configuration")),
  sprintf("LD / HD: DMFT %s", .pget(.P, "DMFT cut-off")),
  sprintf("Workbooks: %s (md5 %s); %s (md5 %s)", basename(RESULTS_FILE), unname(tools::md5sum(RESULTS_FILE)),
          basename(BASELINE_FILE), unname(tools::md5sum(BASELINE_FILE))),
  sprintf("Cross-validation: %d-fold x %d repeats, seed %d; folds from the %s", CV_FOLDS, CV_REPEATS, CV_SEED, FOLD_SOURCE),
  sprintf("Bootstrap seed (%s): %d", if (is.null(CAL_SEED)) "script 01's SEED" else "CAL_SEED", BOOT_SEED),
  sprintf("Cohort: n = %d (%s %d, %s %d), complete for %s", N, REF, sum(Y == 0L), CASE, sum(Y == 1L),
          paste(COHORT_VARS, collapse = ", ")),
  sprintf("%s: refitted predictions vs the run's, largest difference %.2g; metrics vs workbook %.2g",
          COV_ONLY, BASE_PRED_DIFF, BASE_METRIC_DIFF),
  sprintf("Panels: fitting routine %s reproduces the workbook (largest difference %.2g over %d panels x %d metrics)",
          ENGINE_USED, REFIT_MAXDIFF, length(PANELS), length(.METRICS)),
  sprintf("Routines tried: %s", paste(sprintf("%s %.2g", .ENGINES, .maxd), collapse = "; ")),
  sprintf("Best overall panel: %d biomarkers (%s)", BEST_K, paste(PANELS[[as.character(BEST_K)]]$markers, collapse = ", ")),
  sprintf("%s; readxl %s; fastglm %s; caret %s", R.version.string, .pkgv("readxl"), .pkgv("fastglm"), .pkgv("caret")),
  sprintf("Written %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines(RUN_INFO, file.path(OUT_DIR, sprintf("Run_info_%s.txt", TAG)))

.save_figure <- function(stem, width, height, draw, pointsize = 12) {
  for (fmt in FIG_FORMATS) {
    f <- file.path(OUT_DIR, paste0(stem, ".", fmt))
    switch(fmt,
           png  = grDevices::png(f, width = width, height = height, units = "in", res = PNG_DPI, pointsize = pointsize),
           pdf  = grDevices::pdf(f, width = width, height = height, pointsize = pointsize, useDingbats = FALSE),
           tiff = grDevices::tiff(f, width = width, height = height, units = "in", res = PNG_DPI,
                                  pointsize = pointsize, compression = "lzw"),
           stop("FIG_FORMATS: use png, pdf or tiff."))
    tryCatch(draw(), finally = grDevices::dev.off())
    cat("  Saved:", f, "\n")
  }
}
.write_csv <- function(df, stem) {
  f <- file.path(OUT_DIR, paste0(stem, ".csv"))
  utils::write.csv(df, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  cat("  Saved:", f, "\n")
}


# ==========================================================================
# 3) CALIBRATION MEASURES (every model)
# ==========================================================================
# Calibration intercept: logistic regression of HD on an intercept with
#   logit(p) as offset; 0 is ideal, > 0 means predicted risks too low on average
# Calibration slope: logistic regression of HD on logit(p); 1 is ideal, < 1
#   means predictions too extreme
# Brier score: mean (p - y)^2; scaled Brier: 1 - Brier / [ybar (1 - ybar)]
# Classification measures at the model's pooled Youden threshold.
# 95% CI: participant-level bootstrap within HD and LD, each participant
# keeping all of its out-of-fold predictions (predictions held fixed, so the
# interval leaves out panel and threshold selection). Repeat_Min/Max: lowest
# and highest repeat-specific values (fold-assignment variability).

# Calibration intercept: the root of its score, sum(y - plogis(a + logit p)), which falls steadily from
# sum(y) to sum(y) - n as a grows, so it has exactly one root whenever both groups are present. (glm.fit's
# iterations, used before, could diverge when some predictions were extreme, 0 or 1 clipped to CAL_EPS,
# and then returned values such as -4e13 in some bootstrap resamples.)
.cal_intercept <- function(y, lp) {
  if (!any(y == 1) || !any(y == 0)) return(NA_real_)
  f <- function(a) sum(y - stats::plogis(a + lp))
  lo <- -1; hi <- 1
  while (f(lo) <= 0 && lo > -1e8) lo <- 2 * lo
  while (f(hi) >= 0 && hi < 1e8) hi <- 2 * hi
  if (f(lo) <= 0 || f(hi) >= 0) return(NA_real_)
  stats::uniroot(f, c(lo, hi), tol = 1e-12, maxiter = 1000)$root
}
# Calibration slope: glm.fit as before or, when it does not converge, Newton-Raphson with step-halving on
# the (concave) log-likelihood, from the ideal intercept 0 and slope 1
.cal_slope <- function(y, lp) {
  X <- cbind(1, lp)
  f <- suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(), control = list(maxit = 50)))
  if (isTRUE(f$converged) && all(is.finite(f$coefficients))) return(unname(f$coefficients[2]))
  ll <- function(b) { eta <- drop(X %*% b); sum(y * stats::plogis(eta, log.p = TRUE) + (1 - y) * stats::plogis(-eta, log.p = TRUE)) }
  b <- c(0, 1); cur <- ll(b)
  for (it in 1:200) {
    mu <- stats::plogis(drop(X %*% b))
    st <- tryCatch(drop(solve(crossprod(X, X * (mu * (1 - mu))), crossprod(X, y - mu))), error = function(e) NULL)
    if (is.null(st)) break
    t <- 1
    repeat { nb <- b + t * st; nl <- ll(nb); if (is.finite(nl) && nl >= cur) break; t <- t / 2; if (t < 1e-10) break }
    if (t < 1e-10) break
    done <- abs(nl - cur) < 1e-12; b <- nb; cur <- nl
    if (done) break
  }
  b[2]
}
.auc <- function(p, y) {
  r <- rank(p); n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cal_metrics <- function(p, y, thr) {
  pc <- pmin(pmax(p, CAL_EPS), 1 - CAL_EPS); lp <- stats::qlogis(pc)
  a <- .cal_intercept(y, lp)
  b <- .cal_slope(y, lp)
  brier <- mean((p - y)^2); ybar <- mean(y)
  pos <- p >= thr
  sens <- mean(pos[y == 1L]); spec <- mean(!pos[y == 0L])
  c(Calibration_intercept = a, Calibration_slope = b, Brier = brier,
    Scaled_Brier = 1 - brier / (ybar * (1 - ybar)), AUC = .auc(p, y),
    Balanced_accuracy = (sens + spec) / 2, Sensitivity = sens, Specificity = spec,
    Accuracy = mean(pos == (y == 1L)),
    PPV = if (any(pos)) mean(y[pos] == 1L) else NA_real_,
    NPV = if (any(!pos)) mean(y[!pos] == 0L) else NA_real_,
    Mean_predicted = mean(p), Observed = ybar)
}
cal_evaluate <- function(m, B = CAL_BOOT, seed = BOOT_SEED) {
  d <- m$data; thr <- m$threshold
  est <- cal_metrics(d$p, d$y, thr)
  reps <- split(seq_len(nrow(d)), d$rep)
  rv <- t(vapply(reps, function(ix) cal_metrics(d$p[ix], d$y[ix], thr), numeric(length(est))))
  rq <- apply(rv, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  rr <- apply(rv, 2, function(v) suppressWarnings(range(v, na.rm = TRUE)))
  by_id <- split(seq_len(nrow(d)), d$id); ids <- names(by_id)
  y_id <- vapply(by_id, function(ix) d$y[ix[1]], integer(1))
  set.seed(seed)
  boot <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {
    s <- c(sample(ids[y_id == 1L], replace = TRUE), sample(ids[y_id == 0L], replace = TRUE))
    ix <- unlist(by_id[s], use.names = FALSE)
    boot[b, ] <- cal_metrics(d$p[ix], d$y[ix], thr)
  }
  bq <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  data.frame(Model = m$label, Panel_size = m$k, Panel = m$panel, Participants = length(by_id),
             HD = sum(y_id == 1L), LD = sum(y_id == 0L), Predictions = nrow(d), Repeats = length(reps),
             Threshold = thr, Metric = names(est), Estimate = unname(est),
             Bootstrap_CI_Lower = bq[1, ], Bootstrap_CI_Upper = bq[2, ],
             Repeat_Min = rr[1, ], Repeat_Max = rr[2, ], Repeat_P2.5 = rq[1, ], Repeat_P97.5 = rq[2, ],
             Share_predictions_below_0.001_or_above_0.999 = mean(d$p < 0.001 | d$p > 0.999),
             Clipped_predictions = sum(d$p < CAL_EPS | d$p > 1 - CAL_EPS),
             Bootstrap_resamples = B, stringsAsFactors = FALSE, row.names = NULL)
}

cat(sprintf("\n  Calibration of %d models, %d bootstrap resamples each\n", length(MODELS), CAL_BOOT))
CAL_RESULTS <- do.call(rbind, lapply(MODELS, function(m) {
  t0 <- Sys.time(); r <- cal_evaluate(m)
  g <- function(x) r$Estimate[r$Metric == x]
  cat(sprintf("    %-28s intercept %6.2f  slope %5.2f  Brier %.3f  AUC %.3f  (%.0fs)\n", m$label,
              g("Calibration_intercept"), g("Calibration_slope"), g("Brier"), g("AUC"),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  r
}))
rownames(CAL_RESULTS) <- NULL
.write_csv(CAL_RESULTS, sprintf("S3_Calibration_metrics_%s", TAG))


# ==========================================================================
# 4) PERFORMANCE (TABLE 2) AND THE PREDICTION FILE (S2)
# ==========================================================================
.rep_ba <- function(d, thr) vapply(split(seq_len(nrow(d)), d$rep), function(ix) {
  p <- d$p[ix]; y <- d$y[ix]
  (mean(p[y == 1L] >= thr) + mean(p[y == 0L] < thr)) / 2
}, numeric(1))
PERFORMANCE <- do.call(rbind, lapply(MODELS, function(m) {
  pf <- .perf(m$data$p, m$data$y); rb <- .rep_ba(m$data, m$threshold)
  s5 <- mean(m$data$p[m$data$y == 1L] >= 0.5); c5 <- mean(m$data$p[m$data$y == 0L] < 0.5)
  data.frame(Model = m$label, Panel_size = m$k, Panel = m$panel, Routine = m$engine,
             Threshold = pf[["Threshold"]], Balanced_accuracy = pf[["Balanced_Accuracy"]],
             Sensitivity = pf[["Sensitivity"]], Specificity = pf[["Specificity"]], Accuracy = pf[["Accuracy"]],
             PPV = pf[["Precision"]], NPV = pf[["NPV"]], AUC = pf[["AUC"]],
             BA_repeat_min = min(rb), BA_repeat_max = max(rb),
             BA_repeat_P2.5 = stats::quantile(rb, 0.025, names = FALSE),
             BA_repeat_P97.5 = stats::quantile(rb, 0.975, names = FALSE),
             Workbook_BA_repeat_min = m$workbook$Balanced_Accuracy_Repeat_Min,
             Workbook_BA_repeat_max = m$workbook$Balanced_Accuracy_Repeat_Max,
             BA_at_0.5 = (s5 + c5) / 2, Sensitivity_at_0.5 = s5, Specificity_at_0.5 = c5,
             Repeats = length(rb), Predictions = nrow(m$data), stringsAsFactors = FALSE)
}))
rownames(PERFORMANCE) <- NULL
.write_csv(PERFORMANCE, sprintf("Performance_%s", TAG))

# One row per participant and repeat; one prediction column per model
S2 <- do.call(rbind, lapply(seq_along(FOLDS), function(f)
  data.frame(Participant = FOLDS[[f]]$te, Resample = names(FOLDS)[f], stringsAsFactors = FALSE)))
S2$Repeat   <- as.integer(sub("^.*\\.Rep", "", S2$Resample))
S2$Fold     <- as.integer(sub("^Fold([0-9]+)\\..*$", "\\1", S2$Resample))
S2$Observed <- c(REF, CASE)[Y[S2$Participant] + 1L]
.pcols <- character(0)
for (.kk in names(MODELS)) {
  .m <- MODELS[[.kk]]
  .col <- if (.m$k == 0) sprintf("P_%s_%s_only", CASE, paste(COV_WORDS, collapse = "_")) else
    sprintf("P_%s_%d_biomarker%s", CASE, .m$k, if (.m$k == 1) "" else "s")
  S2[[.col]] <- .m$data$p[match(paste(S2$Participant, S2$Resample), paste(.m$data$id, .m$data$resample))]
  .pcols <- c(.pcols, .col)
}
S2 <- S2[order(S2$Repeat, S2$Participant), c("Participant", "Repeat", "Fold", "Resample", "Observed", .pcols)]
if (isTRUE(S2_INCLUDE_STUDYID) && "StudyID" %in% names(DAT)) S2 <- cbind(StudyID = DAT$StudyID[S2$Participant], S2)
if (anyNA(S2[.pcols])) stop("Internal check: a prediction is missing from the prediction file.", call. = FALSE)
.write_csv(S2, sprintf("S2_OOF_predictions_%s", TAG))


# ==========================================================================
# 5) FIGURE 1A-B: ROC CURVES AND CALIBRATION
# ==========================================================================
SIZE_COLS <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")   # Okabe-Ito
SIZE_PCH  <- c(16, 17, 15, 18, 4, 8)
.col_k <- function(k) SIZE_COLS[(k - 1L) %% length(SIZE_COLS) + 1L]
.pch_k <- function(k) SIZE_PCH[(k - 1L) %% length(SIZE_PCH) + 1L]
.n_markers <- function(k) sprintf("%d biomarker%s", k, if (k == 1) "" else "s")
MODEL_TITLE <- sprintf("Logistic regression, adjusted for %s", COV_AND)

.roc_xy <- function(p, y) {
  o <- order(p, decreasing = TRUE); ps <- p[o]; ys <- y[o]
  keep <- c(ps[-1] != ps[-length(ps)], TRUE)
  list(fpr = c(0, cumsum(1 - ys)[keep] / sum(y == 0L)), tpr = c(0, cumsum(ys)[keep] / sum(y == 1L)))
}
# models: list keyed by panel size, each with $data (p, y; thr = each
# prediction's own threshold when it has one) and $threshold
# letter: panel letter drawn above the plot's top-left corner (combined figures)
.letter <- function(letter) if (!is.null(letter)) {
  op <- par(xpd = NA); on.exit(par(op))
  x <- grconvertX(grconvertX(par("usr")[1], "user", "ndc") - 0.045, "ndc", "user")
  mtext(letter, side = 3, line = 1.9, at = x, adj = 0, font = 2, cex = 1.5)
}
draw_roc <- function(models = MODELS, ks = .ks, main = MODEL_TITLE, note = "Symbols: Youden cut of each panel",
                     letter = NULL) {
  op <- par(mar = c(4.6, 4.6, 3.2, 1.6), pty = "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
       main = main, cex.main = 0.98)
  .letter(letter)
  abline(0, 1, col = "grey65")
  for (k in ks) {
    d <- models[[as.character(k)]]$data; r <- .roc_xy(d$p, d$y)
    lines(r$fpr, r$tpr, col = .col_k(k), lwd = 2)
  }
  for (k in ks) {   # each panel's threshold
    m <- models[[as.character(k)]]; d <- m$data
    thr <- if (!is.null(d$thr)) d$thr else m$threshold
    se <- mean((d$p >= thr)[d$y == 1L]); sp <- mean((d$p < thr)[d$y == 0L])
    points(1 - sp, se, pch = .pch_k(k), col = "white", cex = 1.9, lwd = 3)
    points(1 - sp, se, pch = .pch_k(k), col = .col_k(k), cex = 1.35, lwd = 1.5)
  }
  auc <- vapply(ks, function(k) { d <- models[[as.character(k)]]$data; .perf(d$p, d$y)[["AUC"]] }, numeric(1))
  legend("bottomright", legend = sprintf("%s (AUC = %.3f)", vapply(ks, .n_markers, ""), auc),
         col = vapply(ks, .col_k, ""), lty = 1, lwd = 2, pch = vapply(ks, .pch_k, numeric(1)), pt.cex = 1.2,
         bty = "n", inset = c(0.01, 0.06), cex = 0.9, seg.len = 2.6)
  mtext(note, side = 1, line = -1.3, adj = 0.98, cex = 0.7, col = "grey35")
}

.k_fig <- if (is.null(CAL_K)) BEST_K else as.integer(CAL_K)
if (!as.character(.k_fig) %in% names(MODELS)) stop("CAL_K = ", .k_fig, " is not among the models.", call. = FALSE)
CAL_FIG_MODEL <- MODELS[[as.character(.k_fig)]]

.cal_loess <- function(p, y, grid) {
  fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = CAL_SPAN, degree = 2,
                                                control = stats::loess.control(surface = "direct"))),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, length(grid)))
  pmin(pmax(suppressWarnings(stats::predict(fit, newdata = data.frame(p = grid))), 0), 1)
}
prepare_calibration <- function(m, res = NULL) {   # res: the model's rows of CAL_RESULTS (default) or of the nested results
  d <- m$data
  rng  <- range(d$p)
  grid <- seq(rng[1], rng[2], length.out = 200)
  curve <- .cal_loess(d$p, d$y, grid)
  reps <- split(seq_len(nrow(d)), d$rep)
  rc <- vapply(reps, function(ix) {
    inr <- grid >= min(d$p[ix]) & grid <= max(d$p[ix])
    out <- rep(NA_real_, length(grid)); out[inr] <- .cal_loess(d$p[ix], d$y[ix], grid[inr]); out
  }, numeric(length(grid)))
  band <- t(apply(rc, 1, function(v) {
    if (sum(is.finite(v)) < length(reps)) return(c(NA_real_, NA_real_))   # every repeat must cover the point
    if (identical(CAL_BAND, "p95")) stats::quantile(v, c(0.025, 0.975), names = FALSE) else range(v)
  }))
  br  <- unique(stats::quantile(d$p, seq(0, 1, length.out = CAL_BINS + 1), names = FALSE))
  grp <- cut(d$p, br, include.lowest = TRUE)
  pts <- data.frame(pred = tapply(d$p, grp, mean), obs = tapply(d$y, grp, mean))
  list(d = d, grid = grid, curve = curve, band = band, pts = pts[stats::complete.cases(pts), ],
       res = if (is.null(res)) CAL_RESULTS[CAL_RESULTS$Model == m$label, ] else res, label = m$label,
       repeats = length(reps))
}
.fmtn <- function(x, digits) {   # fixed decimals, without "-0.00"
  s <- formatC(x, format = "f", digits = digits)
  z <- grepl("^-0\\.?0*$", s); s[z] <- sub("^-", "", s[z]); s
}
.fmt_ci <- function(res, metric, digits) {
  r <- res[res$Metric == metric, ]
  f <- function(x) .fmtn(x, digits)
  sprintf("%s (%s to %s)", f(r$Estimate), f(r$Bootstrap_CI_Lower), f(r$Bootstrap_CI_Upper))
}
draw_calibration <- function(cal, main = NULL, hist = TRUE, letter = NULL) {
  d <- cal$d
  op <- par(mar = c(if (hist) 0.6 else 4.6, 4.6, 3.2, 1.6), pty = if (hist) "m" else "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = if (hist) "n" else "s",
       xlab = if (hist) "" else sprintf("Predicted risk of %s", CASE), ylab = sprintf("Observed proportion %s", CASE),
       main = main %||% "", cex.main = 0.98)
  .letter(letter)
  abline(h = seq(0.2, 0.8, 0.2), v = seq(0.2, 0.8, 0.2), col = "grey92", lwd = 0.8)
  abline(0, 1, col = "grey55", lty = 2, lwd = 1.5)
  has_band <- cal$repeats > 1   # one repeat (e.g. nested CV with one outer split): no band
  ok <- is.finite(cal$band[, 1]) & is.finite(cal$band[, 2])
  if (has_band && any(ok)) polygon(c(cal$grid[ok], rev(cal$grid[ok])), c(cal$band[ok, 1], rev(cal$band[ok, 2])),
                                   col = grDevices::adjustcolor("#0072B2", 0.18), border = NA)
  lines(cal$grid, cal$curve, col = "#0072B2", lwd = 2.5)
  points(cal$pts$pred, cal$pts$obs, pch = 21, bg = "white", col = "#0b0b0b", cex = 1.2, lwd = 1.4)
  txt <- c(sprintf("Calibration intercept %s", .fmt_ci(cal$res, "Calibration_intercept", 2)),
           sprintf("Calibration slope %s", .fmt_ci(cal$res, "Calibration_slope", 2)),
           sprintf("Brier score %s", .fmt_ci(cal$res, "Brier", 3)))
  legend("bottomright", legend = txt, bty = "n", cex = 0.78, inset = c(0.01, 0.02),
         title = "Estimate (95% bootstrap CI)", title.adj = 0)
  band_lab <- if (identical(CAL_BAND, "p95")) "Middle 95% across repeats" else sprintf("Range across the %d repeats", cal$repeats)
  keep_leg <- if (has_band) 1:4 else c(1, 2, 4)
  legend("topleft", bty = "n", cex = 0.78, inset = c(0.01, 0.01),
         legend = c("Ideal", "Loess calibration curve", band_lab, sprintf("Observed in %d risk groups", nrow(cal$pts)))[keep_leg],
         lty = c(2, 1, NA, NA)[keep_leg], lwd = c(1.5, 2.5, NA, 1.4)[keep_leg], pch = c(NA, NA, 15, 21)[keep_leg],
         col = c("grey55", "#0072B2", grDevices::adjustcolor("#0072B2", 0.25), "#0b0b0b")[keep_leg],
         pt.bg = c(NA, NA, NA, "white")[keep_leg], pt.cex = c(1, 1, 2, 1.2)[keep_leg])
  if (hist) {   # prediction distributions, per repeat: HD up, LD down
    par(mar = c(4.6, 4.6, 0.3, 1.6))
    brk <- seq(0, 1, length.out = CAL_HIST_BINS + 1)
    h1 <- graphics::hist(d$p[d$y == 1L], breaks = brk, plot = FALSE)$counts / cal$repeats
    h0 <- graphics::hist(d$p[d$y == 0L], breaks = brk, plot = FALSE)$counts / cal$repeats
    ym <- max(c(h1, h0), 1)
    plot(NA, xlim = c(0, 1), ylim = c(-ym, ym), xaxs = "i", yaxt = "n",
         xlab = sprintf("Predicted risk of %s", CASE), ylab = "")
    rect(brk[-length(brk)], 0, brk[-1], h1, col = "#D55E00", border = "white", lwd = 0.4)
    rect(brk[-length(brk)], -h0, brk[-1], 0, col = "#56B4E9", border = "white", lwd = 0.4)
    abline(h = 0, col = "grey40", lwd = 0.8)
    mtext(CASE, side = 2, line = 0.6, at = ym / 2, las = 1, cex = 0.7, col = "#D55E00")
    mtext(REF, side = 2, line = 0.6, at = -ym / 2, las = 1, cex = 0.7, col = "#0072B2")
  }
}

cat("\n  Figures\n")
CAL_FIG <- prepare_calibration(CAL_FIG_MODEL)
.title <- sprintf("Calibration: %s", CAL_FIG_MODEL$label)
.save_figure(sprintf("Figure_ROC_%s", TAG), 6.5, 6.5, draw_roc)
.save_figure(sprintf("Figure_Calibration_%s", TAG), 6.5, 7.6, function() {
  layout(matrix(1:2, 2), heights = c(4.6, 1.4))
  draw_calibration(CAL_FIG, main = .title)
})
.save_figure(sprintf("Figure_1AB_ROC_Calibration_%s", TAG), 12.5, 6.9, function() {
  layout(matrix(c(1, 1, 2, 3), 2, 2), widths = c(1, 1), heights = c(4.6, 1.4))
  par(oma = c(0, 0, 1.6, 0))
  draw_roc()
  draw_calibration(CAL_FIG, main = .title)
  mtext("A", outer = TRUE, side = 3, at = 0.012, adj = 0, font = 2, cex = 1.5, line = 0)
  mtext("B", outer = TRUE, side = 3, at = 0.512, adj = 0, font = 2, cex = 1.5, line = 0)
})


# ==========================================================================
# 6) READY-TO-PASTE TEXT
# ==========================================================================
.r  <- CAL_RESULTS[CAL_RESULTS$Model == CAL_FIG_MODEL$label, ]
.g  <- function(metric, col = "Estimate") .r[[col]][.r$Metric == metric]
.f2 <- function(x) .fmtn(x, 2)
.f3 <- function(x) .fmtn(x, 3)
.f4 <- function(x) .fmtn(x, 4)
.pc <- function(x, d = 1) .fmtn(100 * x, d)
.ci <- function(metric, f) sprintf("%s (95%% CI %s to %s)", f(.g(metric)), f(.g(metric, "Bootstrap_CI_Lower")),
                                   f(.g(metric, "Bootstrap_CI_Upper")))
.rg <- function(metric, f) sprintf("%s to %s", f(.g(metric, "Repeat_Min")), f(.g(metric, "Repeat_Max")))
.pf <- PERFORMANCE[PERFORMANCE$Panel_size == CAL_FIG_MODEL$k, ]
.slope_lo <- .g("Calibration_slope", "Bootstrap_CI_Lower"); .slope_hi <- .g("Calibration_slope", "Bootstrap_CI_Upper")
.int_lo <- .g("Calibration_intercept", "Bootstrap_CI_Lower"); .int_hi <- .g("Calibration_intercept", "Bootstrap_CI_Upper")
.interp_slope <- if (.slope_hi < 1) {
  "The calibration slope below 1 indicates predictions that were too extreme, with high risks overestimated and low risks underestimated, consistent with overfitting."
} else if (.slope_lo > 1) {
  "The calibration slope above 1 indicates predictions that were too moderate, with high risks underestimated and low risks overestimated."
} else "The slope's confidence interval included 1."
.interp_int <- if (.int_lo > 0) {
  "The positive calibration intercept indicates that predicted risks were, on average, lower than the observed HD proportion."
} else if (.int_hi < 0) {
  "The negative calibration intercept indicates that predicted risks were, on average, higher than the observed HD proportion."
} else "The intercept's confidence interval included 0, so systematic over- or underestimation was not detected."
.base <- CAL_RESULTS[CAL_RESULTS$Panel_size == 0, ]
.bg <- function(metric) .base$Estimate[.base$Metric == metric]
.table2 <- vapply(.ks, function(k) {
  p <- PERFORMANCE[PERFORMANCE$Panel_size == k, ]
  sprintf("  %d | %s | BA %s | Sens %s | Spec %s | Acc %s | AUC %s | PPV %s | NPV %s | threshold %s | BA range across repeats %s-%s",
          k, p$Panel, .pc(p$Balanced_accuracy), .pc(p$Sensitivity), .pc(p$Specificity), .pc(p$Accuracy), .f3(p$AUC),
          .pc(p$PPV), .pc(p$NPV), .f4(p$Threshold), .pc(p$BA_repeat_min), .pc(p$BA_repeat_max))
}, "")
.b0p <- PERFORMANCE[PERFORMANCE$Panel_size == 0, ]

.text <- c(
  "RUN AND CHECKS",
  paste0("  ", RUN_INFO),
  "",
  "TABLE 2 (best panel at each size; pooled out-of-fold predictions; %, AUC, threshold):",
  .table2,
  sprintf("  0 | %s only | BA %s | Sens %s | Spec %s | Acc %s | AUC %s | PPV %s | NPV %s | threshold %s | BA range across repeats %s-%s",
          gsub("; ", ", ", MODELS[["0"]]$panel, fixed = TRUE), .pc(.b0p$Balanced_accuracy), .pc(.b0p$Sensitivity), .pc(.b0p$Specificity), .pc(.b0p$Accuracy), .f3(.b0p$AUC),
          .pc(.b0p$PPV), .pc(.b0p$NPV), .f4(.b0p$Threshold), .pc(.b0p$BA_repeat_min), .pc(.b0p$BA_repeat_max)),
  "",
  "RESULTS 3.2 (selected panel):",
  paste0("The best overall panel comprised ", .pf$Panel, ". At its probability threshold of ", .f4(.pf$Threshold),
         ", positive and negative predictive values were ", .pc(.pf$PPV, 2), "% and ", .pc(.pf$NPV, 2),
         "%, respectively. Balanced accuracy ranged from ", .pc(.pf$BA_repeat_min, 2), "% to ", .pc(.pf$BA_repeat_max, 2),
         "% across the ", .pf$Repeats, " cross-validation repeats."),
  "",
  "RESULTS 3.2 (calibration, Figure 1B):",
  paste0("Figure 1B shows the calibration of the selected ", CAL_FIG_MODEL$k, "-marker model. The calibration intercept was ",
         .ci("Calibration_intercept", .f2), ", the calibration slope ", .ci("Calibration_slope", .f2),
         ", and the Brier score ", .ci("Brier", .f3), "; the scaled Brier score was ", .f2(.g("Scaled_Brier")),
         ". Across the ", .r$Repeats[1], " repeats, the slope ranged from ", .rg("Calibration_slope", .f2),
         " and the Brier score from ", .rg("Brier", .f3), ". ", .interp_slope, " ", .interp_int,
         sprintf(" Of the pooled predictions, %s%% were below 0.001 or above 0.999.",
                 formatC(100 * .r$Share_predictions_below_0.001_or_above_0.999[1], format = "f", digits = 1)),
         if (nrow(.base)) sprintf(" For comparison, %s alone had a calibration intercept of %s, slope of %s, and Brier score of %s (AUC %s).",
                                  COV_AND, .f2(.bg("Calibration_intercept")), .f2(.bg("Calibration_slope")), .f3(.bg("Brier")), .f3(.bg("AUC"))) else ""),
  "",
  "ABSTRACT (Results):",
  paste0("Calibration intercept and slope were ", .f2(.g("Calibration_intercept")), " and ", .f2(.g("Calibration_slope")),
         ", and the Brier score was ", .f3(.g("Brier")), "."),
  "",
  "BOOTSTRAP INTERVALS FOR DISCRIMINATION (selected panel; optional):",
  paste0("AUC ", .ci("AUC", .f3), "; balanced accuracy ", .ci("Balanced_accuracy", .f3), "; PPV ", .ci("PPV", .f3),
         "; NPV ", .ci("NPV", .f3), " at the ", .f4(.r$Threshold[1]), " threshold."),
  "",
  sprintf("Participants %d (%s %d, %s %d); pooled predictions %d from %d repeats; %d bootstrap resamples; loess span %.2f; band: %s; predictions clipped to [%g, 1 - %g] before the logit (%d clipped).",
          .r$Participants[1], CASE, .r$HD[1], REF, .r$LD[1], .r$Predictions[1], .r$Repeats[1], CAL_BOOT, CAL_SPAN,
          if (identical(CAL_BAND, "p95")) "2.5th-97.5th percentiles of repeat curves" else "lowest to highest repeat curve",
          CAL_EPS, CAL_EPS, .r$Clipped_predictions[1])
)
writeLines(.text, file.path(OUT_DIR, sprintf("Calibration_text_%s.txt", TAG)), useBytes = TRUE)
cat("\n", paste(.text, collapse = "\n"), "\n", sep = "")
# ==========================================================================
# 7) NESTED CROSS-VALIDATION (script 09's workbook)
# ==========================================================================
# Script 09's held-out predictions for the same covariate set: in each outer
# fold, the best panel of each size, its fit and its Youden threshold were
# chosen in the outer training participants only, then applied to the
# held-out participants. The measures of section 3 are computed from those
# predictions; each prediction is classified at the threshold chosen in its own
# outer training set. 95% CIs: participant-level bootstrap within HD and LD,
# predictions held fixed (as in section 3). With more than one outer repeat,
# each measure is the mean over the repeats, as in script 09. The figures set
# the repeated and nested results side by side: the best repeated-CV panel
# (CAL_K) and the best nested-CV panel (NESTED_CAL_K).

.cov_words <- c(Age_DOS = "age", BMI = "BMI", Sex = "sex")
.cov_text_of <- function(covs) unname(ifelse(covs %in% names(.cov_words), .cov_words[covs], covs))
.nested_label <- function(k, covs) {
  w <- .cov_text_of(covs)
  if (k == 0) { s <- paste(w, collapse = " + "); return(paste0(toupper(substring(s, 1, 1)), substring(s, 2), " only")) }
  paste(c(.n_markers(k), w), collapse = " + ")
}
.and_list <- function(w) if (length(w) <= 1) paste(w, collapse = "") else
  paste(paste(w[-length(w)], collapse = ", "), "and", w[length(w)])

# Measures of one nested model (d: id, y, p, thr, rep); mean over outer repeats
.nested_est <- function(d, ix = seq_len(nrow(d))) {
  reps <- split(ix, d$rep[ix])
  colMeans(do.call(rbind, lapply(reps, function(j) cal_metrics(d$p[j], d$y[j], d$thr[j]))))
}
nested_evaluate <- function(m, B = CAL_BOOT, seed = BOOT_SEED) {
  d <- m$data
  est <- .nested_est(d)
  by_id <- split(seq_len(nrow(d)), d$id); ids <- names(by_id)
  y_id <- vapply(by_id, function(ix) d$y[ix[1]], integer(1))
  set.seed(seed)
  boot <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {
    s <- c(sample(ids[y_id == 1L], replace = TRUE), sample(ids[y_id == 0L], replace = TRUE))
    boot[b, ] <- .nested_est(d, unlist(by_id[s], use.names = FALSE))
  }
  bq <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  data.frame(Set = m$set, Method = m$method, Model = m$label, Panel_size = m$k, Participants = length(by_id),
             HD = sum(y_id == 1L), LD = sum(y_id == 0L), Predictions = nrow(d), Outer_repeats = length(unique(d$rep)),
             Metric = names(est), Estimate = unname(est), Bootstrap_CI_Lower = bq[1, ], Bootstrap_CI_Upper = bq[2, ],
             Share_predictions_below_0.05_or_above_0.95 = mean(d$p < 0.05 | d$p > 0.95),
             Clipped_predictions = sum(d$p < CAL_EPS | d$p > 1 - CAL_EPS),
             Bootstrap_resamples = B, stringsAsFactors = FALSE, row.names = NULL)
}

NESTED_FILE_USED <- if (is.null(NESTED_FILE)) sprintf("DMFT_Nested_CV_%s.xlsx", REGISTRY_KEY) else NESTED_FILE
NESTED_DONE <- FALSE
if (!file.exists(NESTED_FILE_USED)) {
  cat(sprintf(paste0("\n  Nested cross-validation: %s is not in %s, so section 7 was skipped.\n",
                     "  Run script 09 for this run first, or give its workbook in NESTED_FILE.\n"),
              basename(NESTED_FILE_USED), getwd()))
} else {
  # --- 7a) 09's workbook, checked against this run ---------------------------
  .nsheets <- readxl::excel_sheets(NESTED_FILE_USED)
  .miss_sh <- setdiff(c("Predictions", "Overview", "Performance_by_repeat", "Settings"), .nsheets)
  if (length(.miss_sh))
    stop(NESTED_FILE_USED, " lacks sheet(s) ", paste(.miss_sh, collapse = ", "), "; it is not script 09's workbook.", call. = FALSE)
  .nread <- function(s) as.data.frame(readxl::read_excel(NESTED_FILE_USED, sheet = s))
  N_SETTINGS <- .nread("Settings"); N_PRED <- .nread("Predictions"); N_OV <- .nread("Overview"); N_PBR <- .nread("Performance_by_repeat")
  .nkey <- as.character(N_SETTINGS$Value[N_SETTINGS$Parameter == "Registry key"])
  if (!identical(.nkey, REGISTRY_KEY))
    stop(sprintf("%s is from another run (registry key %s; this run: %s).", basename(NESTED_FILE_USED), .nkey, REGISTRY_KEY),
         call. = FALSE)
  NESTED_SET_USED <- if (is.null(NESTED_SET)) .pget(.P, "Model configuration") else NESTED_SET
  if (!NESTED_SET_USED %in% N_OV$Set)
    stop(sprintf("Covariate set %s is not among the sets with biomarkers in %s (%s).", NESTED_SET_USED,
                 basename(NESTED_FILE_USED), paste(N_OV$Set, collapse = ", ")), call. = FALSE)
  .ov <- N_OV[N_OV$Set == NESTED_SET_USED, , drop = FALSE][1, ]
  .ncovs <- if (is.na(.ov$Covariates) || .ov$Covariates == "none") character(0) else
    trimws(strsplit(as.character(.ov$Covariates), ",", fixed = TRUE)[[1]])
  .ncov_model <- if ("Covariates_only_model" %in% names(.ov)) as.character(.ov$Covariates_only_model) else NA_character_
  # 09's Source_row is the row of data.dfPH: the same participants and outcomes as this run
  N_PRED$id <- match(N_PRED$Source_row, SOURCE_ROWS)
  if (anyNA(N_PRED$id) || length(unique(N_PRED$id)) != N)
    stop("The participants in 09's workbook are not this run's cohort.", call. = FALSE)
  if (any(as.character(N_PRED$Outcome) != as.character(DAT$condL[N_PRED$id])) || any(N_PRED$y != Y[N_PRED$id]))
    stop("The outcomes in 09's workbook do not match condL in this session.", call. = FALSE)
  .nmodel <- function(set, method, k, label) {
    d <- N_PRED[N_PRED$Set == set & N_PRED$Method == method, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    dd <- data.frame(id = as.integer(d$id), y = as.integer(d$y), p = as.numeric(d$p), thr = as.numeric(d$Threshold),
                     rep = as.integer(d$Outer_repeat), fold = as.integer(d$Outer_fold))
    if (anyDuplicated(dd[c("id", "rep")]) || any(table(dd$rep) != N) || anyNA(dd$p) || anyNA(dd$thr))
      stop(sprintf("09's predictions for %s / %s do not cover every participant once per outer repeat.", set, method),
           call. = FALSE)
    list(k = k, label = label, set = set, method = method, data = dd, threshold = dd$thr)
  }
  NESTED_MODELS <- list()
  for (.k in 1:6) {
    .m <- .nmodel(NESTED_SET_USED, sprintf("Best %d-marker panel", .k), .k, .nested_label(.k, .ncovs))
    if (!is.null(.m)) NESTED_MODELS[[as.character(.k)]] <- .m
  }
  if (!length(NESTED_MODELS)) stop("No best-panel predictions for ", NESTED_SET_USED, " in 09's workbook.", call. = FALSE)
  .nks <- as.integer(names(NESTED_MODELS))
  NESTED_OTHER <- Filter(Negate(is.null), list(
    .nmodel(NESTED_SET_USED, "Panel search", NA_integer_, "Panel search (best panel of any size)"),
    .nmodel(NESTED_SET_USED, "Lasso", NA_integer_, "Lasso"),
    .nmodel(NESTED_SET_USED, "Elastic net (alpha 0.5)", NA_integer_, "Elastic net (alpha 0.5)"),
    .nmodel(NESTED_SET_USED, "Random forest", NA_integer_, "Random forest"),
    if (!is.na(.ncov_model)) .nmodel(.ncov_model, "Covariates only", 0L, .nested_label(0, .ncovs))))
  .n_outer <- as.character(N_SETTINGS$Value[N_SETTINGS$Parameter == "Outer folds x repeats"])

  # --- 7b) Measures --------------------------------------------------------------
  cat(sprintf("\n  Nested cross-validation (%s; set %s; outer folds x repeats %s): %d models, %d bootstrap resamples each\n",
              basename(NESTED_FILE_USED), NESTED_SET_USED, .n_outer, length(NESTED_MODELS) + length(NESTED_OTHER), CAL_BOOT))
  NESTED_RESULTS <- do.call(rbind, lapply(c(NESTED_MODELS, NESTED_OTHER), function(m) {
    r <- nested_evaluate(m)
    g <- function(x) r$Estimate[r$Metric == x]
    cat(sprintf("    %-42s balanced accuracy %.3f  AUC %.3f  slope %5.2f  Brier %.3f\n", m$label,
                g("Balanced_accuracy"), g("AUC"), g("Calibration_slope"), g("Brier")))
    r
  }))
  rownames(NESTED_RESULTS) <- NULL
  .write_csv(NESTED_RESULTS, sprintf("S4_Nested_calibration_metrics_%s", TAG))

  # Check: 09's own numbers (Performance_by_repeat; mean over outer repeats)
  .pairs <- c(AUC = "AUC", Balanced_accuracy = "Balanced_Accuracy", Brier = "Brier",
              Calibration_intercept = "Calibration_intercept", Calibration_slope = "Calibration_slope")
  NESTED_CHECK <- do.call(rbind, lapply(c(NESTED_MODELS, NESTED_OTHER), function(m) {
    w <- N_PBR[N_PBR$Set == m$set & N_PBR$Method == m$method, , drop = FALSE]
    r <- NESTED_RESULTS[NESTED_RESULTS$Set == m$set & NESTED_RESULTS$Method == m$method, ]
    data.frame(Model = m$label, Measure = names(.pairs),
               Script_09 = vapply(.pairs, function(cn) if (nrow(w) && cn %in% names(w)) mean(as.numeric(w[[cn]])) else NA_real_,
                                  numeric(1)),
               This_script = vapply(names(.pairs), function(mn) r$Estimate[r$Metric == mn], numeric(1)),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  NESTED_CHECK$Difference <- abs(NESTED_CHECK$Script_09 - NESTED_CHECK$This_script)
  NESTED_MAXDIFF <- suppressWarnings(max(NESTED_CHECK$Difference, na.rm = TRUE))
  .write_csv(NESTED_CHECK, sprintf("Nested_check_%s", TAG))
  if (!is.finite(NESTED_MAXDIFF) || NESTED_MAXDIFF > NESTED_MATCH_TOL)
    warning(sprintf("The nested measures differ from 09's own (largest difference %.2g; see Nested_check_%s.csv).",
                    NESTED_MAXDIFF, TAG), call. = FALSE, immediate. = TRUE)
  cat(sprintf("  09's own numbers reproduced: largest difference %.1e\n", NESTED_MAXDIFF))

  # One row per model, with the repeated-CV values of the same panel size (Tables 2 and 3)
  NESTED_PERFORMANCE <- do.call(rbind, lapply(c(NESTED_MODELS, NESTED_OTHER), function(m) {
    r <- NESTED_RESULTS[NESTED_RESULTS$Set == m$set & NESTED_RESULTS$Method == m$method, ]
    g <- function(x, col = "Estimate") r[[col]][r$Metric == x]
    rk <- if (!is.na(m$k) && as.character(m$k) %in% names(MODELS)) MODELS[[as.character(m$k)]] else NULL
    rp <- if (!is.null(rk)) PERFORMANCE[PERFORMANCE$Panel_size == m$k, ] else NULL
    rc <- if (!is.null(rk)) CAL_RESULTS[CAL_RESULTS$Model == rk$label, ] else NULL
    gc <- function(x) if (is.null(rc)) NA_real_ else rc$Estimate[rc$Metric == x]
    data.frame(Model = m$label, Set = m$set, Method = m$method, Panel_size = m$k,
               Nested_balanced_accuracy = g("Balanced_accuracy"),
               Nested_balanced_accuracy_CI_Lower = g("Balanced_accuracy", "Bootstrap_CI_Lower"),
               Nested_balanced_accuracy_CI_Upper = g("Balanced_accuracy", "Bootstrap_CI_Upper"),
               Nested_AUC = g("AUC"), Nested_calibration_intercept = g("Calibration_intercept"),
               Nested_calibration_slope = g("Calibration_slope"),
               Nested_calibration_slope_CI_Lower = g("Calibration_slope", "Bootstrap_CI_Lower"),
               Nested_calibration_slope_CI_Upper = g("Calibration_slope", "Bootstrap_CI_Upper"),
               Nested_Brier = g("Brier"), Nested_Brier_CI_Lower = g("Brier", "Bootstrap_CI_Lower"),
               Nested_Brier_CI_Upper = g("Brier", "Bootstrap_CI_Upper"),
               Nested_share_below_0.05_or_above_0.95 = r$Share_predictions_below_0.05_or_above_0.95[1],
               Repeated_balanced_accuracy = if (is.null(rp)) NA_real_ else rp$Balanced_accuracy,
               Repeated_AUC = if (is.null(rp)) NA_real_ else rp$AUC,
               Repeated_calibration_slope = gc("Calibration_slope"), Repeated_Brier = gc("Brier"),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  .write_csv(NESTED_PERFORMANCE, sprintf("Nested_performance_%s", TAG))

  # Best size under nested CV: highest balanced accuracy, then AUC, then fewer markers (as script 09)
  .nsc <- NESTED_PERFORMANCE[NESTED_PERFORMANCE$Method %in% sprintf("Best %d-marker panel", .nks), ]
  NESTED_BEST_K <- .nsc$Panel_size[order(-.nsc$Nested_balanced_accuracy, -.nsc$Nested_AUC, .nsc$Panel_size)][1]
  .nk_fig <- if (is.null(NESTED_CAL_K)) NESTED_BEST_K else as.integer(NESTED_CAL_K)
  if (!as.character(.nk_fig) %in% names(NESTED_MODELS)) stop("NESTED_CAL_K = ", .nk_fig, " is not among the nested models.", call. = FALSE)
  NESTED_CAL_MODEL <- NESTED_MODELS[[as.character(.nk_fig)]]
  NESTED_CAL_FIG <- prepare_calibration(NESTED_CAL_MODEL, res = NESTED_RESULTS[NESTED_RESULTS$Set == NESTED_CAL_MODEL$set &
                                                                                 NESTED_RESULTS$Method == NESTED_CAL_MODEL$method, ])

  # --- 7c) Figures ---------------------------------------------------------------
  .adj <- if (length(.ncovs)) paste("adjusted for", .and_list(.cov_text_of(.ncovs))) else "biomarkers only"
  .ttl_rep_roc <- sprintf("Repeated cross-validation, adjusted for %s", COV_AND)
  .ttl_nes_roc <- sprintf("Nested cross-validation, %s", .adj)
  .ttl_rep_cal <- sprintf("Calibration, repeated cross-validation: %s", CAL_FIG_MODEL$label)
  .ttl_nes_cal <- sprintf("Calibration, nested cross-validation: %s", NESTED_CAL_MODEL$label)
  .nes_note <- "Symbols: thresholds chosen in each outer training set"
  draw_ba_by_size <- function(letter = NULL) {
    kk  <- intersect(.ks, .nks)
    rba <- 100 * vapply(kk, function(k) PERFORMANCE$Balanced_accuracy[PERFORMANCE$Panel_size == k], numeric(1))
    nba <- 100 * vapply(kk, function(k) NESTED_PERFORMANCE$Nested_balanced_accuracy[
      NESTED_PERFORMANCE$Method == sprintf("Best %d-marker panel", k)], numeric(1))
    yl <- range(c(rba, nba)); yl <- c(floor(yl[1] / 10) * 10 - 5, ceiling(yl[2] / 10) * 10 + 5)
    op <- par(mar = c(4.6, 4.6, 3.2, 1.6)); on.exit(par(op))
    plot(NA, xlim = range(kk) + c(-0.35, 0.35), ylim = yl, xaxt = "n",
         xlab = if (length(.ncovs)) sprintf("Number of biomarkers (with %s)", .and_list(.cov_text_of(.ncovs))) else "Number of biomarkers",
         ylab = "Balanced accuracy (%)", main = "Balanced accuracy by panel size", cex.main = 0.98)
    .letter(letter)
    axis(1, at = kk)
    abline(h = pretty(yl), col = "grey92", lwd = 0.8)
    lines(kk, rba, col = "#0072B2", lwd = 2); points(kk, rba, pch = 16, col = "#0072B2", cex = 1.3)
    lines(kk, nba, col = "#D55E00", lwd = 2); points(kk, nba, pch = 17, col = "#D55E00", cex = 1.3)
    legend("bottomleft", legend = c("Repeated cross-validation", "Nested cross-validation"), col = c("#0072B2", "#D55E00"),
           pch = c(16, 17), lwd = 2, bty = "n", inset = c(0.01, 0.02), cex = 0.9)
    box()
  }
  cat("\n  Nested figures\n")
  .save_figure(sprintf("Figure_ROC_nested_%s", TAG), 6.5, 6.5, function()
    draw_roc(NESTED_MODELS, .nks, main = .ttl_nes_roc, note = .nes_note))
  .save_figure(sprintf("Figure_Calibration_nested_%s", TAG), 6.5, 7.6, function() {
    layout(matrix(1:2, 2), heights = c(4.6, 1.4))
    draw_calibration(NESTED_CAL_FIG, main = .ttl_nes_cal)
  })
  .save_figure(sprintf("Figure_BA_by_size_%s", TAG), 6.5, 5.6, function() draw_ba_by_size())
  .save_figure(sprintf("Figure_1_Repeated_Nested_%s", TAG), 12.5, 12.8, function() {
    layout(matrix(c(1, 2, 3, 5, 4, 6), 3, 2, byrow = TRUE), heights = c(5.2, 4.6, 1.4))
    par(oma = c(0, 0, 0.6, 0))
    draw_roc(main = .ttl_rep_roc, letter = "A")
    draw_roc(NESTED_MODELS, .nks, main = .ttl_nes_roc, note = .nes_note, letter = "B")
    draw_calibration(CAL_FIG, main = .ttl_rep_cal, letter = "C")
    draw_calibration(NESTED_CAL_FIG, main = .ttl_nes_cal, letter = "D")
  })

  # --- 7d) Ready-to-paste text and run information ---------------------------------
  .np <- function(meth) NESTED_PERFORMANCE[NESTED_PERFORMANCE$Method == meth & NESTED_PERFORMANCE$Set == NESTED_SET_USED, ]
  .nr <- NESTED_RESULTS[NESTED_RESULTS$Set == NESTED_CAL_MODEL$set & NESTED_RESULTS$Method == NESTED_CAL_MODEL$method, ]
  .nci <- function(metric, f) { r <- .nr[.nr$Metric == metric, ]
    sprintf("%s (95%% CI %s to %s)", f(r$Estimate), f(r$Bootstrap_CI_Lower), f(r$Bootstrap_CI_Upper)) }
  .row_txt <- function(p) sprintf("balanced accuracy %s%% | AUC %s | calibration slope %s (%s to %s) | intercept %s | Brier %s (%s to %s) | predictions < 0.05 or > 0.95: %s%%",
                                  .pc(p$Nested_balanced_accuracy), .f3(p$Nested_AUC), .f2(p$Nested_calibration_slope),
                                  .f2(p$Nested_calibration_slope_CI_Lower), .f2(p$Nested_calibration_slope_CI_Upper),
                                  .f2(p$Nested_calibration_intercept), .f3(p$Nested_Brier), .f3(p$Nested_Brier_CI_Lower),
                                  .f3(p$Nested_Brier_CI_Upper), .pc(p$Nested_share_below_0.05_or_above_0.95))
  .nmeth_other <- vapply(NESTED_OTHER, `[[`, "", "method")
  NESTED_TEXT <- c(
    "NESTED CROSS-VALIDATION",
    sprintf("  Workbook %s; set %s; outer folds x repeats %s; %d bootstrap resamples; 09's own numbers reproduced (largest difference %.1e)",
            basename(NESTED_FILE_USED), NESTED_SET_USED, .n_outer, CAL_BOOT, NESTED_MAXDIFF),
    "",
    "TABLE 2 (best panel of each size): repeated CV | nested CV",
    vapply(.nks, function(k) {
      p <- .np(sprintf("Best %d-marker panel", k))
      sprintf("  %d | repeated: balanced accuracy %s%%, AUC %s, calibration slope %s | nested: %s", k,
              .pc(p$Repeated_balanced_accuracy), .f3(p$Repeated_AUC), .f2(p$Repeated_calibration_slope), .row_txt(p))
    }, ""),
    "",
    "TABLE 3 (comparators, nested):",
    vapply(seq_along(NESTED_OTHER), function(i) sprintf("  %s | %s", NESTED_OTHER[[i]]$label,
                                                         .row_txt(NESTED_PERFORMANCE[NESTED_PERFORMANCE$Model == NESTED_OTHER[[i]]$label &
                                                                                       NESTED_PERFORMANCE$Method == .nmeth_other[i], ])), ""),
    "",
    sprintf("CALIBRATION, nested (best size under nested CV: %d; figure: %s):", NESTED_BEST_K, NESTED_CAL_MODEL$label),
    paste0("The calibration intercept was ", .nci("Calibration_intercept", .f2), ", the calibration slope ",
           .nci("Calibration_slope", .f2), ", and the Brier score ", .nci("Brier", .f3), "; ",
           .pc(.nr$Share_predictions_below_0.05_or_above_0.95[1]), "% of the held-out predictions were below 0.05 or above 0.95."),
    "",
    sprintf("FIGURE 1 (repeated and nested): A, ROC curves under repeated CV; B, under nested CV; C, calibration of the best repeated-CV panel (%s); D, of the best nested-CV panel (%s).",
            CAL_FIG_MODEL$label, NESTED_CAL_MODEL$label))
  writeLines(NESTED_TEXT, file.path(OUT_DIR, sprintf("Nested_text_%s.txt", TAG)), useBytes = TRUE)
  cat("\n", paste(NESTED_TEXT, collapse = "\n"), "\n", sep = "")
  RUN_INFO <- c(RUN_INFO[!startsWith(RUN_INFO, "Written ")],
                sprintf("Nested CV: %s (md5 %s); set %s; outer folds x repeats %s; 09's numbers reproduced, largest difference %.2g",
                        basename(NESTED_FILE_USED), unname(tools::md5sum(NESTED_FILE_USED)), NESTED_SET_USED, .n_outer, NESTED_MAXDIFF),
                RUN_INFO[startsWith(RUN_INFO, "Written ")])
  writeLines(RUN_INFO, file.path(OUT_DIR, sprintf("Run_info_%s.txt", TAG)))
  NESTED_DONE <- TRUE
}
cat(sprintf("\nDone. Upload the files in %s (the .pdf copies are optional).\n", OUT_DIR))
