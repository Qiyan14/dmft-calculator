# ==========================================================================
# ROC CURVES AND CALIBRATION UNDER REPEATED AND NESTED CROSS-VALIDATION,
# EACH WITH ITS OWN BEST PANEL (one covariate set: NESTED_SET)
# 05c_Calibration_Repeated_Nested_DMFT.V16.R
# ==========================================================================
# Run after scripts 05b and 09, in the same working directory (and, ideally,
# the same R session):
#
#   source("05c_Calibration_Repeated_Nested_DMFT.V16.R")
#
# V10: new setting MAX_PANEL_SIZE (section 1). Set to 4, the figures and the
# choice of the best panel use panels of 1 to 4 biomarkers only, as in the
# manuscript; nothing is refitted (05b's and 09's results for those sizes do
# not depend on the larger panels). Output files then end in _max4, so the
# files of an earlier run stay as they are. NULL = every size, as V9.
# V11: 05b's folder may hold the files of more than one run (for example
# ..._10rep.csv and ..._5rep.csv). The run whose balanced accuracy and AUC by
# panel size equal the repeated cross-validation values stored in 09's
# workbook is used, and named in the console; REPEATED_TAG (section 1) names
# it by hand.
# V12: NESTED_SET defaults to "Biomarkers_Age" (the panel search with age
# only, the calculator's model since manuscript 11.20 V5). REPEATED_DIR and
# OUT_DIR default to that set's folders: 05b V10 writes DMFT_10rep_outputs
# for Biomarkers_Age_BMI and DMFT_10rep_outputs_<set> for another set, and
# this script writes DMFT_calibration_repeated_nested (Biomarkers_Age_BMI) or
# DMFT_calibration_repeated_nested_<set>, so earlier figures stay. It stops
# when 05b's run is for another covariate set (05b V10 records the set in its
# Run_info file; 05b V9 ran Biomarkers_Age_BMI only).
# V13: calibration intercept and slope fitted as in 05b V11 (the intercept as
# the root of its score equation; glm.fit's iterations diverged in some
# bootstrap resamples when predictions were 0 or 1 and drew intervals such as
# -9e14 to 0.30); CAL_SEED is 123 by default, the seed of the manuscript's
# intervals (NULL = script 01's SEED, as V12).
# V14: SHOW_CI (default FALSE): the calibration panels print the estimates
# without their bootstrap intervals or the box title, as in the manuscript;
# the metrics file still gives the intervals. TRUE draws V13's panels.
# V15: MAX_PANEL_SIZE defaults to 5 (the manuscript reports panels of one to
# five biomarkers with age); outputs end in _max5.
# V16: one row per analysis: Figure_ROC_Calibration_nested and
# Figure_ROC_Calibration_repeated put the ROC curves (A) beside the calibration
# of that analysis's best panel (B), with its prediction histogram below it
# (the manuscript's Figure 1 and Supplementary Figure S1). Other outputs as V15.
#
# ROC curves: the best panel of each size (1 to MAX_PANEL_SIZE biomarkers) under
# repeated cross-validation (05b's out-of-fold predictions, pooled over the repeats;
# symbols at each panel's Youden cut) and under nested cross-validation (09's
# held-out predictions; symbols at the thresholds chosen in the outer
# training sets). The best panel of each analysis is drawn thicker and
# marked "best" in the legend.
# Calibration: two calibrations, side by side, each of the panel that its own
# analysis ranks best:
#   A. Repeated cross-validation: the best panel under repeated
#      cross-validation (highest balanced accuracy, then AUC, then fewer
#      markers, as scripts 02 and 05b), from its out-of-fold predictions over
#      the repeats (05b's S2_OOF_predictions file). Shaded band: range of
#      the repeat-specific calibration curves.
#   B. Nested cross-validation: the best panel size under nested
#      cross-validation (as script 09), from 09's held-out predictions. In
#      each outer fold the panel, its fit and its threshold were chosen in
#      the outer training set, so the panel can differ between folds (they
#      are listed). No band: one outer split has no repeats to span.
# For each: calibration intercept (logistic regression of HD on an intercept
# with logit(p) as offset; 0 is ideal), calibration slope (on logit(p); 1 is
# ideal, below 1 means predictions too extreme), Brier and scaled Brier
# scores and AUC, with 95% CIs from CAL_BOOT bootstrap resamples of
# participants drawn within HD and LD (each participant keeping all of its
# predictions; predictions held fixed). Measures and figures as in 05b; the
# numbers are checked against 05b's S3 file and 09's own.
#
# Inputs (settings below)
#   REPEATED_DIR  05b's output folder (S2_OOF_predictions, Performance,
#                 Run_info and, for the check, S3_Calibration_metrics)
#   NESTED_FILE   09's workbook DMFT_Nested_CV_<registry_key>.xlsx
# Outputs (OUT_DIR)
#   Figure_ROC_Calibration_nested.png/.pdf    A ROC curves, B calibration of the best panel, nested CV
#   Figure_ROC_Calibration_repeated.png/.pdf  the same for repeated CV
#   Figure_ROC_Calibration_repeated_vs_nested.png/.pdf   A, B ROC curves (repeated, nested);
#                                     C, D calibration of each one's best panel
#   Figure_ROC_repeated_vs_nested.png/.pdf           A repeated, B nested ROC curves
#   Figure_Calibration_repeated_vs_nested.png/.pdf   A repeated, B nested calibration
#   Figure_Calibration_repeated_<k>, Figure_Calibration_nested_<k>   each alone
#   Calibration_repeated_nested_metrics.csv   estimates, 95% CIs, 05b's / 09's values
#   Calibration_repeated_nested_text.txt      ready-to-paste sentences (AUCs by size, calibration)
#   (with MAX_PANEL_SIZE set, every file name ends in _max<size>, e.g. ..._max4.png)
# Packages: readxl. About 1 minute.
# ==========================================================================

# ==========================================================================
# 1) SETTINGS
# ==========================================================================
REPEATED_DIR  <- NULL     # 05b's OUT_DIR; NULL = the folder 05b V10 writes for NESTED_SET
                          # (DMFT_10rep_outputs for Biomarkers_Age_BMI, DMFT_10rep_outputs_<set> otherwise)
REPEATED_TAG  <- NULL     # which run in that folder, e.g. "5rep" (the end of 05b's file names);
                          # NULL = the only one, or the one that matches 09's workbook
NESTED_FILE   <- NULL     # NULL = DMFT_Nested_CV_<registry_key>.xlsx in the working directory
NESTED_SET    <- "Biomarkers_Age"       # covariate set in 09's workbook ("Biomarkers_Age_BMI": age and BMI)
MAX_PANEL_SIZE <- 5       # largest panel size drawn and considered (the manuscript reports 1 to 5);
                          # NULL = every size in 05b's and 09's results
REPEATED_K    <- NULL     # panel size under repeated CV; NULL = the best of the sizes considered
NESTED_K      <- NULL     # panel size under nested CV; NULL = the best of the sizes considered
CAL_BOOT      <- 2000     # bootstrap resamples
CAL_SEED      <- 123      # bootstrap seed (123: the manuscript's intervals); NULL = script 01's SEED,
                          # or the seed recorded in 09's workbook
CAL_SPAN      <- 0.75     # loess span of the calibration curves
CAL_BINS      <- 10       # risk groups for the observed proportions
CAL_EPS       <- 1e-12    # predictions clipped to [CAL_EPS, 1 - CAL_EPS] before the logit
CAL_HIST_BINS <- 40       # bins of the histograms under the calibration plots
CAL_BAND      <- "range"  # repeated CV band: "range" (lowest to highest repeat curve) or "p95"
SHOW_CI       <- FALSE    # print the 95% bootstrap intervals in the calibration panels (the metrics file always has them)
PANEL_LETTERS <- c("A", "B")   # letters of the side-by-side figure; NULL = none
OUT_DIR       <- NULL     # NULL = DMFT_calibration_repeated_nested for Biomarkers_Age_BMI,
                          # DMFT_calibration_repeated_nested_<NESTED_SET> otherwise
FIG_FORMATS   <- c("png", "pdf")   # any of "png", "pdf", "tiff"
PNG_DPI       <- 300
MATCH_TOL     <- 1e-6     # largest difference from 05b's / 09's numbers that counts as a match
CASE <- "HD"; REF <- "LD"


# ==========================================================================
# 2) INPUTS: 05b'S OUT-OF-FOLD PREDICTIONS AND 09'S WORKBOOK
# ==========================================================================
if (!requireNamespace("readxl", quietly = TRUE)) stop('Install readxl first: install.packages("readxl")', call. = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
.must <- function(ok, ...) if (!ok) stop(..., call. = FALSE)
.must(is.character(NESTED_SET) && length(NESTED_SET) == 1 && !is.na(NESTED_SET) && nzchar(NESTED_SET),
      "NESTED_SET: the name of one covariate set, e.g. \"Biomarkers_Age\".")
.set_dir <- function(stem) file.path(getwd(), if (identical(NESTED_SET, "Biomarkers_Age_BMI")) stem else
  sprintf("%s_%s", stem, NESTED_SET))
if (is.null(REPEATED_DIR)) REPEATED_DIR <- .set_dir("DMFT_10rep_outputs")
if (is.null(OUT_DIR)) OUT_DIR <- .set_dir("DMFT_calibration_repeated_nested")
.must(is.null(MAX_PANEL_SIZE) || (is.numeric(MAX_PANEL_SIZE) && length(MAX_PANEL_SIZE) == 1 && !is.na(MAX_PANEL_SIZE) &&
                                  MAX_PANEL_SIZE >= 1 && MAX_PANEL_SIZE == round(MAX_PANEL_SIZE)),
      "MAX_PANEL_SIZE: a whole number of biomarkers (1 or more), or NULL for every size.")
KMAX <- if (is.null(MAX_PANEL_SIZE)) 12L else as.integer(MAX_PANEL_SIZE)   # 12: more sizes than scripts 02 and 09 search
LIMITED <- !is.null(MAX_PANEL_SIZE)
FILE_SUFFIX <- if (LIMITED) sprintf("_max%d", KMAX) else ""

# --- 05b's outputs ------------------------------------------------------------------
.must(dir.exists(REPEATED_DIR), sprintf("%s does not exist: run script 05b first, or set REPEATED_DIR to its output folder.",
                                        REPEATED_DIR))
.s2 <- list.files(REPEATED_DIR, pattern = "^S2_OOF_predictions_.*\\.csv$")
.must(length(.s2) >= 1, sprintf("%s holds no S2_OOF_predictions_<tag>.csv: run script 05b first, or set REPEATED_DIR.", REPEATED_DIR))
TAGS_05B <- sub("^S2_OOF_predictions_(.*)\\.csv$", "\\1", .s2)
if (!is.null(REPEATED_TAG)) {
  .must(length(REPEATED_TAG) == 1 && REPEATED_TAG %in% TAGS_05B,
        sprintf("REPEATED_TAG \"%s\" is not in %s (runs there: %s).", paste(REPEATED_TAG, collapse = " "), REPEATED_DIR,
                paste(TAGS_05B, collapse = ", ")))
  TAGS_05B <- as.character(REPEATED_TAG)
}
.f05t <- function(stem, tag) file.path(REPEATED_DIR, sprintf("%s_%s.%s", stem, tag, if (stem == "Run_info") "txt" else "csv"))
.key05 <- function(tag) {   # registry key recorded in a run's Run_info file (NA when there is none)
  f <- .f05t("Run_info", tag)
  if (!file.exists(f)) return(NA_character_)
  sub("^Run \\(registry key\\): ", "", grep("^Run \\(registry key\\): ", readLines(f, warn = FALSE), value = TRUE)[1])
}
KEYS_05B <- vapply(TAGS_05B, .key05, "")
.ukey <- unique(KEYS_05B[!is.na(KEYS_05B)])
KEY_05B <- if (length(.ukey) == 1) .ukey else NA_character_   # set again below, once the run is chosen

# --- 09's workbook --------------------------------------------------------------------
NESTED_FILE_USED <- NESTED_FILE
if (is.null(NESTED_FILE_USED)) {
  .k <- if (exists("registry_key")) get("registry_key") else KEY_05B
  .must(!is.na(.k) && length(.k) == 1,
        "registry_key is not in this R session and 05b's Run_info file does not settle the run: set NESTED_FILE to 09's workbook.")
  NESTED_FILE_USED <- sprintf("DMFT_Nested_CV_%s.xlsx", .k)
}
.must(file.exists(NESTED_FILE_USED), sprintf("%s is not in %s. Run script 09 first, or set NESTED_FILE.", NESTED_FILE_USED, getwd()))
.miss <- setdiff(c("Predictions", "Performance_by_repeat", "Overview", "Settings", "Selected_panels"),
                 readxl::excel_sheets(NESTED_FILE_USED))
.must(!length(.miss), NESTED_FILE_USED, " lacks sheet(s) ", paste(.miss, collapse = ", "), "; it is not script 09's workbook.")
.rd <- function(s) as.data.frame(withCallingHandlers(readxl::read_excel(NESTED_FILE_USED, sheet = s),
  warning = function(w) if (grepl("^Coercing text to numeric", conditionMessage(w))) invokeRestart("muffleWarning")))
# (09 writes the threshold of a model with constant predictions as the text "Inf"; read as Inf)
N_PRED <- .rd("Predictions"); N_PBR <- .rd("Performance_by_repeat"); N_OV <- .rd("Overview"); N_SET <- .rd("Settings")
N_SEL <- .rd("Selected_panels")
.setting <- function(prefix) {
  v <- N_SET$Value[startsWith(as.character(N_SET$Parameter), prefix)]
  if (length(v)) as.character(v[1]) else NA_character_
}
KEY_09 <- .setting("Registry key")

# --- which of 05b's runs (when its folder holds several) ----------------------------------
if (length(TAGS_05B) > 1) {
  .bs <- if ("By_size" %in% readxl::excel_sheets(NESTED_FILE_USED)) .rd("By_size") else NULL
  .bs <- if (!is.null(.bs) && all(c("Set", "Panel_size", "Apparent_Balanced_Accuracy", "Apparent_AUC") %in% names(.bs)))
    .bs[.bs$Set == NESTED_SET, , drop = FALSE] else NULL
  .fits <- vapply(TAGS_05B, function(tag) {   # same repeated cross-validation as the one 09's workbook records?
    if (!is.na(KEYS_05B[[tag]]) && !is.na(KEY_09) && !identical(KEYS_05B[[tag]], KEY_09)) return(FALSE)
    f <- .f05t("Performance", tag)
    if (is.null(.bs) || !nrow(.bs) || !file.exists(f)) return(NA)
    p <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    m <- merge(data.frame(k = as.integer(.bs$Panel_size), ba = as.numeric(.bs$Apparent_Balanced_Accuracy),
                          auc = as.numeric(.bs$Apparent_AUC)),
               data.frame(k = as.integer(p$Panel_size), ba5 = as.numeric(p$Balanced_accuracy), auc5 = as.numeric(p$AUC)), by = "k")
    nrow(m) > 0 && max(abs(m$ba - m$ba5), abs(m$auc - m$auc5)) < 1e-6
  }, logical(1))
  .hit <- TAGS_05B[!is.na(.fits) & .fits]
  .must(length(.hit) == 1, sprintf(
    "%s holds %d runs of script 05b (%s) and %s. Set REPEATED_TAG in section 1 to the run to use, for example REPEATED_TAG <- \"%s\".",
    REPEATED_DIR, length(TAGS_05B), paste(TAGS_05B, collapse = ", "),
    if (length(.hit) > 1) "more than one matches 09's workbook" else "none could be matched to 09's workbook", TAGS_05B[1]))
  cat(sprintf("  05b's folder holds %d runs (%s): using \"%s\", whose balanced accuracy and AUC by panel size equal those in 09's workbook.\n",
              length(TAGS_05B), paste(TAGS_05B, collapse = ", "), .hit))
  TAGS_05B <- .hit
}
TAG_05B <- TAGS_05B
.f05 <- function(stem) .f05t(stem, TAG_05B)
.must(file.exists(.f05("Performance")), "05b's Performance_", TAG_05B, ".csv is missing from ", REPEATED_DIR, ".")
S2  <- utils::read.csv(.f05("S2_OOF_predictions"), check.names = FALSE, stringsAsFactors = FALSE)
PERF05 <- utils::read.csv(.f05("Performance"), check.names = FALSE, stringsAsFactors = FALSE)
S3_05B <- if (file.exists(.f05("S3_Calibration_metrics")))
  utils::read.csv(.f05("S3_Calibration_metrics"), check.names = FALSE, stringsAsFactors = FALSE) else NULL
.rinfo <- if (file.exists(.f05("Run_info"))) readLines(.f05("Run_info"), warn = FALSE) else character(0)
KEY_05B <- sub("^Run \\(registry key\\): ", "", grep("^Run \\(registry key\\): ", .rinfo, value = TRUE)[1])
SEED_05B <- suppressWarnings(as.integer(sub("^.*: *", "", grep("^Bootstrap seed", .rinfo, value = TRUE)[1])))
# Covariate set of 05b's run: recorded by 05b V10; 05b V9 and earlier ran Biomarkers_Age_BMI only
SET_05B <- sub("^Covariate set: ([^ ]+).*$", "\\1", grep("^Covariate set: ", .rinfo, value = TRUE)[1])
if (is.na(SET_05B)) SET_05B <- "Biomarkers_Age_BMI"
.must(identical(SET_05B, NESTED_SET),
      sprintf(paste0("05b's outputs in %s are for covariate set %s, not %s. Run 05b V10 with MODEL_SET <- \"%s\", ",
                     "or set REPEATED_DIR to the folder of that run."), REPEATED_DIR, SET_05B, NESTED_SET, NESTED_SET))

# --- same run? -------------------------------------------------------------------------
.must(is.na(KEY_05B) || is.na(KEY_09) || identical(KEY_05B, KEY_09),
      sprintf("05b's outputs are from run %s but 09's workbook from run %s.", KEY_05B, KEY_09))
if (exists("registry_key") && !is.na(KEY_09) && !identical(KEY_09, get("registry_key")))
  warning(sprintf("09's workbook is from run %s; this session's registry_key is %s.", KEY_09, get("registry_key")),
          call. = FALSE, immediate. = TRUE)
.wb_seed <- suppressWarnings(as.integer(.setting("Seed")))
SEED_USED <- if (!is.null(CAL_SEED)) as.integer(CAL_SEED) else if (exists("SEED")) as.integer(get("SEED")) else .wb_seed
.must(length(SEED_USED) == 1 && !is.na(SEED_USED), "No seed: set CAL_SEED, or run script 01 first (its SEED).")
if (is.null(CAL_SEED) && !is.na(.wb_seed) && .wb_seed != SEED_USED)
  warning(sprintf("The bootstrap seed (%d) is not the seed of 09's run (%d).", SEED_USED, .wb_seed),
          call. = FALSE, immediate. = TRUE)
.must(NESTED_SET %in% N_OV$Set, sprintf("Covariate set %s is not in %s (sets: %s).", NESTED_SET,
                                        basename(NESTED_FILE_USED), paste(N_OV$Set, collapse = ", ")))
.covs <- as.character(N_OV$Covariates[N_OV$Set == NESTED_SET][1])
.covs <- if (is.na(.covs) || .covs == "none") character(0) else trimws(strsplit(.covs, ",", fixed = TRUE)[[1]])
.cov_words <- c(Age_DOS = "age", BMI = "BMI", Sex = "sex")
.cov_text <- unname(ifelse(.covs %in% names(.cov_words), .cov_words[.covs], .covs))
.label <- function(k) paste(c(sprintf("%d biomarker%s", k, if (k == 1) "" else "s"), .cov_text), collapse = " + ")
.s2col <- function(k) sprintf("P_%s_%d_biomarker%s", CASE, k, if (k == 1) "" else "s")

# --- the best panel of each analysis -----------------------------------------------------
.rp <- PERF05[PERF05$Panel_size >= 1 & PERF05$Panel_size <= KMAX, , drop = FALSE]
.must(nrow(.rp) > 0, "05b's Performance file has no biomarker panels", if (LIMITED) sprintf(" of 1 to %d markers", KMAX), ".")
BEST_REPEATED_K <- .rp$Panel_size[order(-.rp$Balanced_accuracy, -.rp$AUC, .rp$Panel_size)][1]
.p02 <- as.character(N_OV$Script02_panel[N_OV$Set == NESTED_SET][1])
.k02 <- if (is.na(.p02) || !nzchar(.p02)) NA_integer_ else length(strsplit(.p02, ";", fixed = TRUE)[[1]])
if (!is.na(.k02) && .k02 != BEST_REPEATED_K && !(LIMITED && .k02 > KMAX))   # (a larger script-02 panel is expected under a limit)
  warning(sprintf("05b's best panel under repeated CV has %d markers; script 02's panel in 09's workbook has %d.",
                  BEST_REPEATED_K, .k02), call. = FALSE, immediate. = TRUE)
.nsz <- do.call(rbind, lapply(seq_len(KMAX), function(k) {
  w <- N_PBR[N_PBR$Set == NESTED_SET & N_PBR$Method == sprintf("Best %d-marker panel", k), , drop = FALSE]
  if (!nrow(w)) return(NULL)
  data.frame(k = k, BA = mean(as.numeric(w$Balanced_Accuracy)), AUC = mean(as.numeric(w$AUC)))
}))
.must(!is.null(.nsz), "09's workbook has no best-panel results for ", NESTED_SET, ".")
SIZES_REP <- sort(unique(.rp$Panel_size)); SIZES_NES <- .nsz$k
if (LIMITED && (max(SIZES_REP) < KMAX || max(SIZES_NES) < KMAX))
  warning(sprintf("MAX_PANEL_SIZE is %d, but the results stop at %d (repeated) and %d (nested) markers.", KMAX, max(SIZES_REP),
                  max(SIZES_NES)), call. = FALSE, immediate. = TRUE)
BEST_NESTED_K <- .nsz$k[order(-.nsz$BA, -.nsz$AUC, .nsz$k)][1]
K_REP <- if (is.null(REPEATED_K)) BEST_REPEATED_K else as.integer(REPEATED_K)
K_NES <- if (is.null(NESTED_K)) BEST_NESTED_K else as.integer(NESTED_K)
.must(!LIMITED || (K_REP <= KMAX && K_NES <= KMAX),
      sprintf("REPEATED_K / NESTED_K ask for a panel larger than MAX_PANEL_SIZE (%d): raise it, or set it to NULL.", KMAX))
.best_words <- if (LIMITED && KMAX > 1) sprintf("best panel of 1 to %d biomarkers", KMAX) else "best panel"
.must(.s2col(K_REP) %in% names(S2), "05b's S2 file has no column ", .s2col(K_REP), ".")
.must(K_REP %in% PERF05$Panel_size, "05b's Performance file has no ", K_REP, "-marker panel.")

# Panels chosen in 09's outer training sets
.chosen <- function(k) {
  w <- N_SEL[N_SEL$Set == NESTED_SET & N_SEL$Choice == sprintf("Best %d-marker panel", k), , drop = FALSE]
  if (!nrow(w)) return(NA_character_)
  w <- w[order(w$Outer_repeat, w$Outer_fold), ]
  tag <- if (length(unique(w$Outer_repeat)) > 1) sprintf("%d.%d", w$Outer_repeat, w$Outer_fold) else as.character(w$Outer_fold)
  g <- split(tag, factor(as.character(w$Panel), levels = unique(as.character(w$Panel))))   # by first fold
  paste(sprintf("fold%s %s: %s", ifelse(lengths(g) > 1, "s", ""), vapply(g, paste, "", collapse = ", "),
                gsub("_", " ", names(g))), collapse = " | ")
}

# --- the two models --------------------------------------------------------------------------
.r05 <- PERF05[PERF05$Panel_size == K_REP, ][1, ]
REP_MODEL <- list(
  analysis = "Repeated", k = K_REP, label = .label(K_REP), panel = gsub("_", " ", as.character(.r05$Panel)),
  role = if (K_REP == BEST_REPEATED_K) paste(.best_words, "under repeated cross-validation") else "",
  data = data.frame(id = as.integer(S2$Participant), y = as.integer(S2$Observed == CASE), p = as.numeric(S2[[.s2col(K_REP)]]),
                    thr = as.numeric(.r05$Threshold), rep = as.integer(S2$Repeat)))
.nd <- N_PRED[N_PRED$Set == NESTED_SET & N_PRED$Method == sprintf("Best %d-marker panel", K_NES), , drop = FALSE]
.must(nrow(.nd) > 0, sprintf("09's workbook has no predictions for the best %d-marker panel (%s).", K_NES, NESTED_SET))
.nd <- .nd[order(.nd$Outer_repeat, .nd$Source_row), ]
NES_MODEL <- list(
  analysis = "Nested", k = K_NES, label = .label(K_NES), panel = .chosen(K_NES),
  role = if (K_NES == BEST_NESTED_K) "best panel size under nested cross-validation; panel chosen in each outer fold" else
    "panel chosen in each outer fold",
  data = data.frame(id = as.integer(.nd$Source_row), y = as.integer(.nd$y), p = as.numeric(.nd$p),
                    thr = as.numeric(.nd$Threshold), rep = as.integer(.nd$Outer_repeat)))
for (m in list(REP_MODEL, NES_MODEL)) {
  d <- m$data; n_id <- length(unique(d$id))
  .must(!anyNA(d$p) && !anyNA(d$y) && !anyDuplicated(d[c("id", "rep")]) && all(table(d$rep) == n_id),
        sprintf("The %s predictions do not cover every participant once per repeat.", tolower(m$analysis)))
}
.yr <- REP_MODEL$data[REP_MODEL$data$rep == REP_MODEL$data$rep[1], ]; .yn <- NES_MODEL$data[NES_MODEL$data$rep == 1, ]
.must(nrow(.yr) == nrow(.yn) && sum(.yr$y) == sum(.yn$y),
      sprintf("05b's cohort (%d participants, %d %s) is not 09's (%d, %d).", nrow(.yr), sum(.yr$y), CASE, nrow(.yn), sum(.yn$y)))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\n  Run %s | %d participants (%s %d, %s %d) | %d bootstrap resamples, seed %d\n", KEY_09 %||% "?", nrow(.yn),
            CASE, sum(.yn$y), REF, nrow(.yn) - sum(.yn$y), CAL_BOOT, SEED_USED))
if (LIMITED) cat(sprintf("  MAX_PANEL_SIZE = %d: panels of 1 to %d biomarkers only (files end in %s)\n", KMAX, KMAX, FILE_SUFFIX))
cat(sprintf("  Repeated CV (%s, %d repeats): best panel %d markers -> %s\n", basename(.f05("S2_OOF_predictions")),
            length(unique(REP_MODEL$data$rep)), BEST_REPEATED_K, REP_MODEL$panel))
cat(sprintf("  Nested CV (%s): best size %d markers -> %s\n", basename(NESTED_FILE_USED), BEST_NESTED_K, NES_MODEL$panel))

# ROC curves: every panel size of each analysis
REP_ROC <- list(); NES_ROC <- list()
for (k in seq_len(KMAX)) {
  if (.s2col(k) %in% names(S2) && k %in% PERF05$Panel_size)
    REP_ROC[[as.character(k)]] <- list(k = k, data = data.frame(p = as.numeric(S2[[.s2col(k)]]),
      y = as.integer(S2$Observed == CASE), thr = as.numeric(PERF05$Threshold[PERF05$Panel_size == k][1])))
  .d <- N_PRED[N_PRED$Set == NESTED_SET & N_PRED$Method == sprintf("Best %d-marker panel", k), , drop = FALSE]
  if (nrow(.d))
    NES_ROC[[as.character(k)]] <- list(k = k, data = data.frame(p = as.numeric(.d$p), y = as.integer(.d$y),
                                                               thr = as.numeric(.d$Threshold)))
}
.must(length(REP_ROC) > 0 && length(NES_ROC) > 0, "No panels for the ROC curves.")


# ==========================================================================
# 3) CALIBRATION MEASURES (as script 05b)
# ==========================================================================
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
    Balanced_accuracy = (sens + spec) / 2, Sensitivity = sens, Specificity = spec)
}
# repeated CV: pooled over the repeats (as 05b's section 3); nested CV: mean over outer repeats (as 05b's section 7 and 09)
.estimate <- function(m, ix = seq_len(nrow(m$data))) {
  d <- m$data
  if (m$analysis == "Repeated") return(cal_metrics(d$p[ix], d$y[ix], d$thr[ix]))
  reps <- split(ix, d$rep[ix])
  colMeans(do.call(rbind, lapply(reps, function(j) cal_metrics(d$p[j], d$y[j], d$thr[j]))))
}
evaluate <- function(m, B = CAL_BOOT, seed = SEED_USED) {
  d <- m$data
  est <- .estimate(m)
  by_id <- split(seq_len(nrow(d)), d$id); ids <- names(by_id)
  y_id <- vapply(by_id, function(ix) d$y[ix[1]], integer(1))
  set.seed(seed)
  boot <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {   # participants resampled within HD and LD, with all their predictions
    s <- c(sample(ids[y_id == 1L], replace = TRUE), sample(ids[y_id == 0L], replace = TRUE))
    boot[b, ] <- .estimate(m, unlist(by_id[s], use.names = FALSE))
  }
  bq <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  # reference values: 05b's S3 file (repeated) or 09's Performance_by_repeat (nested)
  ref <- setNames(rep(NA_real_, length(est)), names(est)); ref_lo <- ref; ref_hi <- ref
  if (m$analysis == "Repeated" && !is.null(S3_05B)) {
    s3 <- S3_05B[S3_05B$Panel_size == m$k, , drop = FALSE]
    same_boot <- nrow(s3) > 0 && all(s3$Bootstrap_resamples == B) && !is.na(SEED_05B) && SEED_05B == seed
    for (x in intersect(names(est), s3$Metric)) {
      ref[x] <- s3$Estimate[s3$Metric == x]
      if (same_boot) { ref_lo[x] <- s3$Bootstrap_CI_Lower[s3$Metric == x]; ref_hi[x] <- s3$Bootstrap_CI_Upper[s3$Metric == x] }
    }
    if (!same_boot) cat(sprintf("  (05b's intervals not compared: it used %s resamples and seed %s; here %d and %d)\n",
                                paste(unique(s3$Bootstrap_resamples), collapse = "/"), SEED_05B, B, seed))
  }
  if (m$analysis == "Nested") {
    w <- N_PBR[N_PBR$Set == NESTED_SET & N_PBR$Method == sprintf("Best %d-marker panel", m$k), , drop = FALSE]
    s09 <- c(Calibration_intercept = "Calibration_intercept", Calibration_slope = "Calibration_slope", Brier = "Brier",
             Scaled_Brier = "Scaled_Brier", AUC = "AUC", Balanced_accuracy = "Balanced_Accuracy",
             Sensitivity = "Sensitivity", Specificity = "Specificity")
    for (x in names(est)) if (nrow(w) && s09[[x]] %in% names(w)) ref[x] <- mean(as.numeric(w[[s09[[x]]]]))
  }
  data.frame(Analysis = m$analysis, Model = m$label, Panel_size = m$k, Role = m$role, Panel = m$panel,
             Participants = length(by_id), HD = sum(y_id == 1L), LD = sum(y_id == 0L), Predictions = nrow(d),
             Repeats = length(unique(d$rep)), Metric = names(est), Estimate = unname(est), CI_Lower = bq[1, ],
             CI_Upper = bq[2, ], Reference = if (m$analysis == "Repeated") "05b S3" else "09",
             Reference_value = unname(ref), Reference_CI_Lower = unname(ref_lo), Reference_CI_Upper = unname(ref_hi),
             Share_predictions_below_0.05_or_above_0.95 = mean(d$p < 0.05 | d$p > 0.95),
             Confidently_wrong = sum((d$p > 0.95 & d$y == 0L) | (d$p < 0.05 & d$y == 1L)),
             Bootstrap_resamples = B, Seed = seed, stringsAsFactors = FALSE, row.names = NULL)
}
RESULTS <- do.call(rbind, lapply(list(REP_MODEL, NES_MODEL), function(m) {
  r <- evaluate(m)
  g <- function(x, col = "Estimate") r[[col]][r$Metric == x]
  cat(sprintf("    %-8s %-26s intercept %5.2f (%5.2f to %5.2f)  slope %4.2f (%4.2f to %4.2f)  Brier %.3f (%.3f to %.3f)\n",
              m$analysis, m$label, g("Calibration_intercept"), g("Calibration_intercept", "CI_Lower"),
              g("Calibration_intercept", "CI_Upper"), g("Calibration_slope"), g("Calibration_slope", "CI_Lower"),
              g("Calibration_slope", "CI_Upper"), g("Brier"), g("Brier", "CI_Lower"), g("Brier", "CI_Upper")))
  r
}))
rownames(RESULTS) <- NULL
.chk_metrics <- c("Calibration_intercept", "Calibration_slope", "Brier", "Scaled_Brier", "AUC")
.chk <- with(RESULTS[RESULTS$Metric %in% .chk_metrics, ],
             c(abs(Estimate - Reference_value), abs(CI_Lower - Reference_CI_Lower), abs(CI_Upper - Reference_CI_Upper)))
MAXDIFF <- suppressWarnings(max(.chk, na.rm = TRUE))
if (!is.finite(MAXDIFF)) {
  cat("  (no reference values to check against: 05b's S3 file is missing)\n")
} else {
  if (MAXDIFF > MATCH_TOL)
    warning(sprintf("The measures differ from 05b's / 09's (largest difference %.2g); see the metrics file.", MAXDIFF),
            call. = FALSE, immediate. = TRUE)
  cat(sprintf("  05b's and 09's numbers %s: largest difference %.1e\n",
              if (MAXDIFF > MATCH_TOL) "NOT reproduced (an older 05b run? see the metrics file)" else "reproduced", MAXDIFF))
}
.f <- file.path(OUT_DIR, sprintf("Calibration_repeated_nested_metrics%s.csv", FILE_SUFFIX))
utils::write.csv(RESULTS, .f, row.names = FALSE, na = "", fileEncoding = "UTF-8"); cat("  Saved:", .f, "\n")


# ==========================================================================
# 4) FIGURES (as Figure 1B of script 05b)
# ==========================================================================
.cal_loess <- function(p, y, grid) {
  fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = CAL_SPAN, degree = 2,
                                                control = stats::loess.control(surface = "direct"))),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, length(grid)))
  pmin(pmax(suppressWarnings(stats::predict(fit, newdata = data.frame(p = grid))), 0), 1)
}
prepare_calibration <- function(m) {
  d <- m$data
  grid <- seq(min(d$p), max(d$p), length.out = 200)
  reps <- split(seq_len(nrow(d)), d$rep)
  band <- if (length(reps) > 1) {
    rc <- vapply(reps, function(ix) {
      inr <- grid >= min(d$p[ix]) & grid <= max(d$p[ix])
      out <- rep(NA_real_, length(grid)); out[inr] <- .cal_loess(d$p[ix], d$y[ix], grid[inr]); out
    }, numeric(length(grid)))
    t(apply(rc, 1, function(v) {
      if (sum(is.finite(v)) < length(reps)) return(c(NA_real_, NA_real_))   # every repeat must cover the point
      if (identical(CAL_BAND, "p95")) stats::quantile(v, c(0.025, 0.975), names = FALSE) else range(v)
    }))
  } else matrix(NA_real_, length(grid), 2)
  br  <- unique(stats::quantile(d$p, seq(0, 1, length.out = CAL_BINS + 1), names = FALSE))
  grp <- cut(d$p, br, include.lowest = TRUE)
  pts <- data.frame(pred = tapply(d$p, grp, mean), obs = tapply(d$y, grp, mean))
  list(d = d, grid = grid, curve = .cal_loess(d$p, d$y, grid), band = band, pts = pts[stats::complete.cases(pts), ],
       res = RESULTS[RESULTS$Analysis == m$analysis, ], label = m$label, role = m$role, analysis = m$analysis,
       repeats = length(reps))
}
.fmtn <- function(x, digits) {   # fixed decimals, without "-0.00"
  s <- formatC(x, format = "f", digits = digits)
  z <- grepl("^-0\\.?0*$", s); s[z] <- sub("^-", "", s[z]); s
}
.fmt_ci <- function(res, metric, digits) {
  r <- res[res$Metric == metric, ]
  sprintf("%s (%s to %s)", .fmtn(r$Estimate, digits), .fmtn(r$CI_Lower, digits), .fmtn(r$CI_Upper, digits))
}
# letter: panel letter drawn above the plot's top-left corner (combined figures)
.letter <- function(letter) if (!is.null(letter)) {
  op <- par(xpd = NA); on.exit(par(op))
  x <- grconvertX(grconvertX(par("usr")[1], "user", "ndc") - 0.045, "ndc", "user")
  mtext(letter, side = 3, line = 1.9, at = x, adj = 0, font = 2, cex = 1.5)
}
draw_calibration <- function(cal, main = NULL, hist = TRUE, role = TRUE, letter = NULL) {
  d <- cal$d
  op <- par(mar = c(if (hist) 0.6 else 4.6, 4.6, 3.2, 1.6), pty = if (hist) "m" else "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = if (hist) "n" else "s",
       xlab = if (hist) "" else sprintf("Predicted risk of %s", CASE), ylab = sprintf("Observed proportion %s", CASE),
       main = main %||% "", cex.main = 0.98)
  if (role && nzchar(cal$role %||% "")) mtext(sprintf("(%s)", cal$role), side = 3, line = 0.3, cex = 0.87 * par("cex"), col = "grey25")   # scales with the layout
  .letter(letter)
  abline(h = seq(0.2, 0.8, 0.2), v = seq(0.2, 0.8, 0.2), col = "grey92", lwd = 0.8)
  abline(0, 1, col = "grey55", lty = 2, lwd = 1.5)
  has_band <- cal$repeats > 1
  ok <- is.finite(cal$band[, 1]) & is.finite(cal$band[, 2])
  if (has_band && any(ok)) polygon(c(cal$grid[ok], rev(cal$grid[ok])), c(cal$band[ok, 1], rev(cal$band[ok, 2])),
                                   col = grDevices::adjustcolor("#0072B2", 0.18), border = NA)
  lines(cal$grid, cal$curve, col = "#0072B2", lwd = 2.5)
  points(cal$pts$pred, cal$pts$obs, pch = 21, bg = "white", col = "#0b0b0b", cex = 1.2, lwd = 1.4)
  .fmt_stat <- if (isTRUE(SHOW_CI)) .fmt_ci else function(res, metric, digits) .fmtn(res$Estimate[res$Metric == metric], digits)
  txt <- c(sprintf("Calibration intercept %s", .fmt_stat(cal$res, "Calibration_intercept", 2)),
           sprintf("Calibration slope %s", .fmt_stat(cal$res, "Calibration_slope", 2)),
           sprintf("Brier score %s", .fmt_stat(cal$res, "Brier", 3)))
  legend("bottomright", legend = txt, bty = "n", cex = 0.78, inset = c(0.01, 0.02),
         title = if (isTRUE(SHOW_CI)) "Estimate (95% bootstrap CI)" else NULL, title.adj = 0)
  band_lab <- if (identical(CAL_BAND, "p95")) "Middle 95% across repeats" else sprintf("Range across the %d repeats", cal$repeats)
  keep_leg <- if (has_band) 1:4 else c(1, 2, 4)
  legend("topleft", bty = "n", cex = 0.78, inset = c(0.01, 0.01),
         legend = c("Ideal", "Loess calibration curve", band_lab, sprintf("Observed in %d risk groups", nrow(cal$pts)))[keep_leg],
         lty = c(2, 1, NA, NA)[keep_leg], lwd = c(1.5, 2.5, NA, 1.4)[keep_leg], pch = c(NA, NA, 15, 21)[keep_leg],
         col = c("grey55", "#0072B2", grDevices::adjustcolor("#0072B2", 0.25), "#0b0b0b")[keep_leg],
         pt.bg = c(NA, NA, NA, "white")[keep_leg], pt.cex = c(1, 1, 2, 1.2)[keep_leg])
  if (hist) {   # prediction distributions (per repeat): HD up, LD down
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
.save_figure <- function(stem, width, height, draw, pointsize = 12) {
  for (fmt in FIG_FORMATS) {
    f <- file.path(OUT_DIR, paste0(stem, FILE_SUFFIX, ".", fmt))
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
# ROC curves, as Figure 1A of script 05b
SIZE_COLS <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")   # Okabe-Ito
SIZE_PCH  <- c(16, 17, 15, 18, 4, 8)
.col_k <- function(k) SIZE_COLS[(k - 1L) %% length(SIZE_COLS) + 1L]
.pch_k <- function(k) SIZE_PCH[(k - 1L) %% length(SIZE_PCH) + 1L]
.n_markers <- function(k) sprintf("%d biomarker%s", k, if (k == 1) "" else "s")
.roc_xy <- function(p, y) {
  o <- order(p, decreasing = TRUE); ps <- p[o]; ys <- y[o]
  keep <- c(ps[-1] != ps[-length(ps)], TRUE)
  list(fpr = c(0, cumsum(1 - ys)[keep] / sum(y == 0L)), tpr = c(0, cumsum(ys)[keep] / sum(y == 1L)))
}
draw_roc <- function(models, main, note, best = NA, letter = NULL) {
  op <- par(mar = c(4.6, 4.6, 3.2, 1.6), pty = "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity", main = main, cex.main = 0.98)
  .letter(letter)
  abline(0, 1, col = "grey65")
  ks <- as.integer(names(models))
  lw <- ifelse(ks == best, 3.2, 2)
  for (i in seq_along(ks)) {
    d <- models[[i]]$data; r <- .roc_xy(d$p, d$y)
    lines(r$fpr, r$tpr, col = .col_k(ks[i]), lwd = lw[i])
  }
  for (i in seq_along(ks)) {   # each panel's threshold (nested: each prediction's own fold threshold)
    d <- models[[i]]$data
    se <- mean((d$p >= d$thr)[d$y == 1L]); sp <- mean((d$p < d$thr)[d$y == 0L])
    points(1 - sp, se, pch = .pch_k(ks[i]), col = "white", cex = 1.9, lwd = 3)
    points(1 - sp, se, pch = .pch_k(ks[i]), col = .col_k(ks[i]), cex = 1.35, lwd = 1.5)
  }
  auc <- vapply(models, function(m) .auc(m$data$p, m$data$y), numeric(1))
  legend("bottomright", legend = sprintf("%s (AUC = %.3f)%s", vapply(ks, .n_markers, ""), auc, ifelse(ks == best, ", best", "")),
         col = vapply(ks, .col_k, ""), lty = 1, lwd = lw, pch = vapply(ks, .pch_k, numeric(1)), pt.cex = 1.2,
         bty = "n", inset = c(0.01, 0.06), cex = 0.9, seg.len = 2.6)
  mtext(note, side = 1, line = -1.3, adj = 0.98, cex = 0.7, col = "grey35")
}
.adj <- if (length(.cov_text)) paste("adjusted for", if (length(.cov_text) == 1) .cov_text else
  paste(paste(.cov_text[-length(.cov_text)], collapse = ", "), "and", .cov_text[length(.cov_text)])) else "biomarkers only"
TITLE_ROC_REP <- sprintf("Repeated cross-validation, %s", .adj)
TITLE_ROC_NES <- sprintf("Nested cross-validation, %s", .adj)
NOTE_ROC_REP <- "Symbols: Youden cut of each panel"
NOTE_ROC_NES <- "Symbols: thresholds chosen in each outer training set"
.save_figure("Figure_ROC_repeated_vs_nested", 12.5, 6.5, function() {
  layout(matrix(1:2, 1))
  par(oma = c(0, 0, if (is.null(PANEL_LETTERS)) 0 else 1.6, 0))
  draw_roc(REP_ROC, TITLE_ROC_REP, NOTE_ROC_REP, best = BEST_REPEATED_K)
  draw_roc(NES_ROC, TITLE_ROC_NES, NOTE_ROC_NES, best = BEST_NESTED_K)
  if (!is.null(PANEL_LETTERS))
    for (i in 1:2) mtext(PANEL_LETTERS[i], outer = TRUE, side = 3, at = (i - 1) / 2 + 0.012, adj = 0, font = 2, cex = 1.5,
                         line = 0)
})

CAL_REP <- prepare_calibration(REP_MODEL)
CAL_NES <- prepare_calibration(NES_MODEL)
TITLE_REP <- sprintf("Calibration, repeated cross-validation: %s", REP_MODEL$label)
TITLE_NES <- sprintf("Calibration, nested cross-validation: %s", NES_MODEL$label)
.save_figure(sprintf("Figure_Calibration_repeated_%d", K_REP), 6.5, 7.6, function() {
  layout(matrix(1:2, 2), heights = c(4.6, 1.4)); draw_calibration(CAL_REP, main = TITLE_REP)
})
.save_figure(sprintf("Figure_Calibration_nested_%d", K_NES), 6.5, 7.6, function() {
  layout(matrix(1:2, 2), heights = c(4.6, 1.4)); draw_calibration(CAL_NES, main = TITLE_NES)
})
.save_figure("Figure_Calibration_repeated_vs_nested", 12.5, 6.9, function() {
  layout(matrix(1:4, 2), heights = c(4.6, 1.4))
  par(oma = c(0, 0, if (is.null(PANEL_LETTERS)) 0 else 1.6, 0))
  draw_calibration(CAL_REP, main = TITLE_REP)
  draw_calibration(CAL_NES, main = TITLE_NES)
  if (!is.null(PANEL_LETTERS))
    for (i in 1:2) mtext(PANEL_LETTERS[i], outer = TRUE, side = 3, at = (i - 1) / 2 + 0.012, adj = 0, font = 2, cex = 1.5,
                         line = 0)
})
.save_figure("Figure_ROC_Calibration_repeated_vs_nested", 12.5, 12.8, function() {
  layout(matrix(c(1, 2, 3, 5, 4, 6), 3, 2, byrow = TRUE), heights = c(5.2, 4.6, 1.4))
  par(oma = c(0, 0, 0.6, 0))
  draw_roc(REP_ROC, TITLE_ROC_REP, NOTE_ROC_REP, best = BEST_REPEATED_K, letter = "A")
  draw_roc(NES_ROC, TITLE_ROC_NES, NOTE_ROC_NES, best = BEST_NESTED_K, letter = "B")
  draw_calibration(CAL_REP, main = TITLE_REP, letter = "C")
  draw_calibration(CAL_NES, main = TITLE_NES, letter = "D")
})
# One row per analysis: ROC curves (A) beside the calibration of its best panel (B), with the
# prediction histogram under the calibration plot (Figure 1: nested; Supplementary Figure S1: repeated)
.one_row <- function(roc, roc_title, roc_note, best, cal, cal_title) function() {
  layout(matrix(c(1, 1, 2, 3), 2, 2), heights = c(4.6, 1.4))
  par(oma = c(0, 0, 0.6, 0))
  draw_roc(roc, roc_title, roc_note, best = best, letter = "A")
  draw_calibration(cal, main = cal_title, letter = "B")
}
.save_figure("Figure_ROC_Calibration_nested", 12.5, 6.9,
             .one_row(NES_ROC, TITLE_ROC_NES, NOTE_ROC_NES, BEST_NESTED_K, CAL_NES, TITLE_NES))
.save_figure("Figure_ROC_Calibration_repeated", 12.5, 6.9,
             .one_row(REP_ROC, TITLE_ROC_REP, NOTE_ROC_REP, BEST_REPEATED_K, CAL_REP, TITLE_REP))


# ==========================================================================
# 5) READY-TO-PASTE TEXT
# ==========================================================================
.ci <- function(r, metric, digits) {
  x <- r[r$Metric == metric, ]
  sprintf("%s (95%% CI %s to %s)", .fmtn(x$Estimate, digits), .fmtn(x$CI_Lower, digits), .fmtn(x$CI_Upper, digits))
}
.sent <- function(m) {
  r <- RESULTS[RESULTS$Analysis == m$analysis, ]
  sprintf("%s cross-validation, %s (%s): calibration intercept %s, calibration slope %s, Brier score %s; %s%% of the predictions were below 0.05 or above 0.95. Panel: %s.",
          m$analysis, m$label, m$role, .ci(r, "Calibration_intercept", 2), .ci(r, "Calibration_slope", 2), .ci(r, "Brier", 3),
          formatC(100 * r$Share_predictions_below_0.05_or_above_0.95[1], format = "f", digits = 1), m$panel)
}
TEXT <- c(
  "CALIBRATION, REPEATED AND NESTED CROSS-VALIDATION (each with its own best panel)",
  sprintf("  Run %s; repeated: %s (%d repeats); nested: %s; %d bootstrap resamples (seed %d); %s",
          KEY_09 %||% "?", basename(.f05("S2_OOF_predictions")), length(unique(REP_MODEL$data$rep)), basename(NESTED_FILE_USED),
          CAL_BOOT, SEED_USED, if (is.finite(MAXDIFF)) sprintf("05b's and 09's numbers %s (largest difference %.1e)",
                                                    if (MAXDIFF > MATCH_TOL) "NOT reproduced" else "reproduced", MAXDIFF)
          else "not checked (05b's S3 file missing)"),
  if (LIMITED) sprintf("  Panel sizes considered: 1 to %d biomarkers (MAX_PANEL_SIZE); larger panels are left out of the figures and of the choice of the best panel.", KMAX),
  "",
  sprintf("ROC, repeated cross-validation (pooled over the repeats): %s; best panel: %d marker%s.",
          paste(sprintf("%s AUC %.3f", vapply(REP_ROC, function(m) .n_markers(m$k), ""),
                        vapply(REP_ROC, function(m) .auc(m$data$p, m$data$y), numeric(1))), collapse = ", "),
          BEST_REPEATED_K, if (BEST_REPEATED_K == 1) "" else "s"),
  sprintf("ROC, nested cross-validation (held-out predictions): %s; best size: %d marker%s.",
          paste(sprintf("%s AUC %.3f", vapply(NES_ROC, function(m) .n_markers(m$k), ""),
                        vapply(NES_ROC, function(m) .auc(m$data$p, m$data$y), numeric(1))), collapse = ", "),
          BEST_NESTED_K, if (BEST_NESTED_K == 1) "" else "s"),
  "", .sent(REP_MODEL), "", .sent(NES_MODEL), "",
  "FIGURE (draft caption):",
  paste0("Calibration of the best panel of each analysis", if (LIMITED) sprintf(" (panels of 1 to %d biomarkers)", KMAX),
         ". (", PANEL_LETTERS[1] %||% "A", ") Repeated cross-validation: the ",
         REP_MODEL$label, " panel ranked highest under repeated cross-validation, from its out-of-fold predictions over the ",
         length(unique(REP_MODEL$data$rep)), " repeats; the band spans the repeat-specific curves. (", PANEL_LETTERS[2] %||% "B",
         ") Nested cross-validation: the panel size ranked highest under nested cross-validation (", NES_MODEL$label,
         "), with the panel, its fit, and its threshold chosen within each outer training set; held-out predictions pooled over ",
         "the outer folds. Loess calibration curve (solid), observed HD proportions by predicted-risk decile (circles), and ideal ",
         "calibration (dashed); bars show the HD (upward) and LD (downward) prediction distributions.",
         if (isTRUE(SHOW_CI)) " Statistics with 95% bootstrap confidence intervals." else ""))
writeLines(TEXT, file.path(OUT_DIR, sprintf("Calibration_repeated_nested_text%s.txt", FILE_SUFFIX)), useBytes = TRUE)
cat("\n", paste(TEXT, collapse = "\n"), "\n", sep = "")
cat(sprintf("\nDone. The files are in %s.\n", OUT_DIR))
