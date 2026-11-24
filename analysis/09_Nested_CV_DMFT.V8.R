# ==========================================================================
# NESTED CROSS-VALIDATION FOR EVERY COVARIATE SET (stand-alone)
# ==========================================================================
# Honest performance estimates (reviewer comment 7) and comparison with
# standard methods (comment 6). This script does not source script 01 or 02. It
# needs
#   - script 01 run in this R session: the preprocessed data (data.dfPH) and
#     the covariate sets (the add_set() lines), and
#   - script 02's Excel output DMFT_Covariate_Sets_Summary_<registry_key>.xlsx,
#     with the apparent results of those sets.
#
#   source("01_Preprocessing_DMFT.R")   # (file names without their version numbers)
#   source("02_LogReg_DMFT.R")
#   source("09_Nested_CV_DMFT.R")       # reads the Excel file script 02 just wrote
#
# The covariate sets come from script 01 (its add_set() lines), together with
# the covariates-only model script 02 paired with each biomarker set. Script
# 02's file must hold every one of them; sets in the file that script 01 no
# longer has are left out (listed when the script runs). So after adding a set
# to script 01, run script 02 again before this script; after removing one, the
# earlier file of script 02 can still be used.
#
# To use an Excel file from an earlier run of script 02 instead, import it with
# NESTED_RESULTS_FILE (section 1): "choose" opens a file browser, or give the
# file's path. Script 01 must still be run first, with the settings of that run.
#
# For every covariate set, on the same folds:
#   1. The participants (script 02's common cohort) are split into
#      NESTED_OUTER_FOLDS groups, with LD and HD balanced across groups.
#   2. One group is set aside. Script 02's procedure is run again from the
#      start on the other groups only:
#        - sets with biomarkers: every single marker and every pair, then beam
#          search to NESTED_MAX_MARKERS markers keeping NESTED_BEAM_WIDTH panels
#          per size, each panel scored by NESTED_INNER_FOLDS-fold CV repeated
#          NESTED_INNER_REPEATS times. The best panel is chosen as in script 02
#          (balanced accuracy, then AUC, then fewer markers), with its Youden
#          threshold from the same predictions (ties: closest to 0.5). The best
#          panel of each size is kept as well;
#        - sets without biomarkers: the covariates are fitted, with the Youden
#          threshold of their cross-validated predictions;
#        - benchmarks, for sets with biomarkers: lasso, elastic net (alpha 0.5)
#          and random forest, each with the set's covariates (unpenalized in
#          glmnet).
#   3. The chosen models predict the set-aside group, using their thresholds;
#      then the next group is set aside.
#   4. The predictions of the set-aside groups are scored. Those participants
#      took no part in choosing the panel, the threshold or any tuning.
# The panels chosen in different groups can differ from each other and from
# script 02's panel. Script 02's panel (from all participants) is the panel to
# report; these scores estimate how well a panel found by script 02's
# procedure predicts participants who were not used to find it.
#
# Run time: with the defaults (script 02's settings), each covariate set with
# biomarkers costs about NESTED_OUTER_FOLDS x NESTED_OUTER_REPEATS times script
# 02's sweep for that set. A projection is printed after the first set. For a
# quicker, rougher run, lower NESTED_INNER_REPEATS (e.g. 2) and
# NESTED_BEAM_WIDTH (e.g. 10); the Settings sheet records the values used.
#
# Outputs (working directory)
#   DMFT_Nested_CV_<registry_key>.xlsx      Overview, Performance, By_size,
#     Selected_panels, Marker_frequency, Benchmark_markers,
#     Performance_by_repeat, Predictions, Settings, Notes
#   DMFT_Nested_CV_text_<registry_key>.txt  ready-to-paste sentences
# ==========================================================================


# ==========================================================================
# 1) SETTINGS
# ==========================================================================

# Script 02's results (DMFT_Covariate_Sets_Summary_<registry_key>.xlsx)
#   NULL      the file script 02 wrote in the working directory (run 09 right after 02)
#   "choose"  import an Excel file chosen in a file browser
#   a path    import that Excel file, e.g. "C:/Users/me/Results/DMFT_Covariate_Sets_Summary_<registry_key>.xlsx"
NESTED_RESULTS_FILE  <- NULL
# NULL = every covariate set of script 01 (its add_set() lines) and the
# covariates-only models script 02 paired with them; or some of those, by name,
# e.g. c("Biomarkers_Age_BMI", "Age_BMI_Only")
NESTED_SETS          <- NULL
NESTED_OUTER_FOLDS   <- 5L
NESTED_OUTER_REPEATS <- 1L     # > 1 repeats everything on new outer splits (steadier, slower)
# Inside each outer training set; the defaults are script 02's values
NESTED_INNER_FOLDS   <- 5L
NESTED_INNER_REPEATS <- 5L
NESTED_BEAM_WIDTH    <- 50L    # panels kept per size
NESTED_MAX_MARKERS   <- 5L     # panels of one to five biomarkers
# Seed: script 01's SEED (its section 1) sets the outer and inner folds and the random forests
NESTED_CORES         <- 10L
NESTED_BENCHMARKS    <- TRUE   # lasso, elastic net, random forest
RF_TREES             <- 500L
NESTED_CAL_EPS       <- 1e-12  # predictions clipped to [eps, 1 - eps] before the logit (calibration)


# ==========================================================================
# 2) PREREQUISITES AND SCRIPT 02'S EXCEL FILE
# ==========================================================================

.need_pkgs <- c("readxl", "writexl", "parallel", if (NESTED_BENCHMARKS) c("glmnet", "ranger"))
.miss_pkgs <- .need_pkgs[!vapply(.need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkgs) > 0)
  stop("Missing package(s): ", paste(.miss_pkgs, collapse = ", "),
       "\n  install.packages(c(", paste0('"', .miss_pkgs, '"', collapse = ", "), "))", call. = FALSE)
if (!exists("PROJECT") || !identical(PROJECT, "DMFT") || !exists("registry_key") ||
    !exists("data.dfPH") || !exists("biomarkers_PH_vif") || !exists("DMFT_MODEL_SETS") || !exists("SEED"))
  stop("Run script 01 (preprocessing, PROJECT = \"DMFT\") in this session first.", call. = FALSE)
.SEED <- as.integer(SEED)   # script 01's seed

# File browser: RStudio's when available, otherwise R's own
.choose_results_file <- function() {
  if (!interactive())
    stop("No file browser outside an interactive session: set NESTED_RESULTS_FILE to the path of script 02's Excel file.",
         call. = FALSE)
  cap <- "Choose script 02's Excel file (DMFT_Covariate_Sets_Summary_....xlsx)"
  message(cap)
  f <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudioapi::selectFile(caption = cap, path = getwd(), filter = "Excel files (*.xlsx)", existing = TRUE)
  } else if (.Platform$OS.type == "windows") {
    utils::choose.files(default = gsub("/", "\\\\", file.path(getwd(), "*.xlsx")), caption = cap, multi = FALSE,
                        filters = matrix(c("Excel files (*.xlsx)", "*.xlsx"), 1))
  } else tryCatch(file.choose(), error = function(e) NULL)
  if (!length(f) || is.na(f[1]) || !nzchar(f[1])) stop("No file chosen.", call. = FALSE)
  f[1]
}
.res_file <- if (is.null(NESTED_RESULTS_FILE)) sprintf("DMFT_Covariate_Sets_Summary_%s.xlsx", registry_key) else
  if (identical(NESTED_RESULTS_FILE, "choose")) .choose_results_file() else NESTED_RESULTS_FILE
if (!file.exists(.res_file))
  stop(if (is.null(NESTED_RESULTS_FILE))
         paste0("Script 02's Excel file is not in the working directory (", getwd(), "):\n  ", .res_file,
                "\n  Run 02_LogReg_DMFT_Age_BMI.R first, or import the file with NESTED_RESULTS_FILE",
                " (\"choose\" for a file browser, or its path).")
       else paste0("File not found: ", .res_file), call. = FALSE)
cat("Script 02's results:", normalizePath(.res_file), "\n")
.file_key <- sub("^DMFT_Covariate_Sets_Summary_(.*)\\.xlsx$", "\\1", basename(.res_file))
if (.file_key != basename(.res_file) && .file_key != registry_key)
  warning("The file name's settings (", .file_key, ") differ from this session's (", registry_key,
          "). Script 01 should be run with the settings of script 02's run.", call. = FALSE, immediate. = TRUE)
.sheets <- tryCatch(readxl::excel_sheets(.res_file), error = function(e)
  stop(.res_file, " could not be read as an Excel file: ", conditionMessage(e), call. = FALSE))
if (!all(c("Run_Status", "Overview") %in% .sheets))
  stop(.res_file, " is not script 02's covariate-set summary (needs the Run_Status and Overview sheets).", call. = FALSE)
.read_sheet <- function(s) if (s %in% .sheets) as.data.frame(readxl::read_excel(.res_file, sheet = s)) else data.frame()
.status   <- .read_sheet("Run_Status")
.apparent <- .read_sheet("Overview")         # best panel of each set with biomarkers
.app_size <- .read_sheet("Best_By_Size")     # best panel of each size
.app_cov  <- .read_sheet("Covariates_Only")  # covariates-only models

.parse_covs <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(trimws(x)) || trimws(x) == "none") character(0) else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}
.all_sets <- lapply(seq_len(nrow(.status)), function(i)
  list(name = as.character(.status$Model[i]), biomarkers = isTRUE(as.logical(.status$Biomarkers[i])),
       covariates = .parse_covs(.status$Covariates[i]), done = identical(as.character(.status$Status[i]), "done")))
names(.all_sets) <- vapply(.all_sets, `[[`, character(1), "name")

# The covariate sets to run are script 01's (its add_set() lines, as script 02
# reads them), each biomarker set with the covariates-only model script 02
# paired with it. Script 02's file only supplies their apparent results, so a
# file from a run with more sets can be used; one without a set cannot.
.sets01 <- get0("biomarker_registry", ifnotfound = list())[[registry_key]]$model_sets
if (is.null(.sets01)) .sets01 <- DMFT_MODEL_SETS
if (!is.list(.sets01) || !length(.sets01) || is.null(names(.sets01)))
  stop("Script 01's covariate sets (DMFT_MODEL_SETS) are missing: run script 01 to the end in this session.",
       call. = FALSE)
.covtxt <- function(v) if (length(v)) paste(v, collapse = ", ") else "none"
.covkey <- function(v) paste(sort(as.character(v)), collapse = "+")
.in01 <- names(.sets01)
.miss <- setdiff(.in01, names(.all_sets))
if (length(.miss))
  stop("Script 02's file has no results for these covariate sets of script 01: ", paste(.miss, collapse = ", "),
       ".\n  ", basename(.res_file), " is from a run of script 02 with other add_set() lines in script 01.",
       "\n  Run script 02 again (after script 01), then this script.", call. = FALSE)
.changed <- .in01[vapply(.in01, function(nm)
  !identical(.covkey(.sets01[[nm]]$covariates), .covkey(.all_sets[[nm]]$covariates)) ||
    !identical(isTRUE(.sets01[[nm]]$biomarkers), .all_sets[[nm]]$biomarkers), logical(1))]
if (length(.changed))
  stop("Covariate set(s) that differ between script 01 and script 02's file: ",
       paste(sprintf("%s (script 01: %s%s; file: %s%s)", .changed,
                     vapply(.changed, function(nm) if (isTRUE(.sets01[[nm]]$biomarkers)) "biomarkers + " else "", ""),
                     vapply(.changed, function(nm) .covtxt(.sets01[[nm]]$covariates), ""),
                     vapply(.changed, function(nm) if (.all_sets[[nm]]$biomarkers) "biomarkers + " else "", ""),
                     vapply(.changed, function(nm) .covtxt(.all_sets[[nm]]$covariates), "")), collapse = "; "),
       ".\n  Run script 02 again (after script 01), then this script.", call. = FALSE)
.paired <- if (nrow(.apparent) && all(c("Set", "Covariates_only_model") %in% names(.apparent)))
  as.character(.apparent$Covariates_only_model[as.character(.apparent$Set) %in% .in01]) else character(0)
.paired <- intersect(setdiff(unique(.paired[!is.na(.paired) & nzchar(.paired)]), .in01), names(.all_sets))
.run <- names(.all_sets)[names(.all_sets) %in% c(.in01, .paired)]   # in script 02's order
.left_out <- setdiff(names(.all_sets), .run)
.failed <- .run[!vapply(.all_sets[.run], `[[`, logical(1), "done")]
if (length(.failed))
  warning("Not completed in script 02's run, so left out: ", paste(.failed, collapse = ", "),
          " (see the Run_Status sheet of ", basename(.res_file), ").", call. = FALSE, immediate. = TRUE)
NESTED_SET_LIST <- .all_sets[setdiff(.run, .failed)]
cat("Covariate sets from script 01:", paste(.in01, collapse = ", "), "\n")
if (length(.paired)) cat("  with script 02's covariates-only model(s):", paste(.paired, collapse = ", "), "\n")
if (length(.left_out))
  cat("  In script 02's file but not among script 01's sets, so left out:", paste(.left_out, collapse = ", "), "\n")
if (!is.null(NESTED_SETS)) {
  if (length(setdiff(NESTED_SETS, names(NESTED_SET_LIST))))
    stop("NESTED_SETS not among the sets above: ", paste(setdiff(NESTED_SETS, names(NESTED_SET_LIST)), collapse = ", "),
         call. = FALSE)
  NESTED_SET_LIST <- NESTED_SET_LIST[NESTED_SETS]
}
if (!length(NESTED_SET_LIST)) stop("No completed covariate sets to run.", call. = FALSE)
.bm_sets <- names(Filter(function(s) s$biomarkers, NESTED_SET_LIST))


# ==========================================================================
# 3) DATA: SCRIPT 02'S COMMON COHORT
# ==========================================================================
# Participants complete for every covariate of every set that script 02 ran
# (as in script 02), with LD / HD outcome; checked against the Excel counts.
# These are script 01's sets' participants too when the file is the one script
# 02 wrote for this session's settings (the registry key lists every covariate).

.cov_all <- unique(unlist(lapply(.all_sets, `[[`, "covariates")))
.cov01   <- unique(unlist(lapply(.sets01, function(s) as.character(s$covariates))))
if (!setequal(.cov_all, .cov01))
  warning(sprintf(paste("Script 02's run used the participants complete for %s (every covariate of its sets);",
                        "script 01's sets need %s. The nested results use script 02's participants, to match its",
                        "apparent results; run script 02 again for script 01's sets alone."),
                  .covtxt(.cov_all), .covtxt(.cov01)), call. = FALSE, immediate. = TRUE)
ND <- as.data.frame(data.dfPH)
if (length(setdiff(c("condL", .cov_all), names(ND))))
  stop("Column(s) missing from data.dfPH: ", paste(setdiff(c("condL", .cov_all), names(ND)), collapse = ", "), call. = FALSE)
for (.v in .cov_all) if (is.character(ND[[.v]]) || is.logical(ND[[.v]])) ND[[.v]] <- factor(ND[[.v]])
.ok <- (if (length(.cov_all)) stats::complete.cases(ND[, .cov_all, drop = FALSE]) else rep(TRUE, nrow(ND))) & !is.na(ND$condL)
ND_SOURCE_ROWS <- which(.ok)
ND <- droplevels(ND[.ok, , drop = FALSE])
if (!setequal(as.character(unique(ND$condL)), c("LD", "HD"))) stop("Expected LD and HD outcome labels.", call. = FALSE)
ND$condL <- factor(as.character(ND$condL), levels = c("LD", "HD"))
ND_Y <- as.integer(ND$condL == "HD")
.n <- nrow(ND)
.chk <- rbind(if (nrow(.apparent)) .apparent[, c("N", "LD", "HD")], if (nrow(.app_cov)) .app_cov[, c("N", "LD", "HD")])
if (is.null(.chk) || !nrow(.chk)) stop(.res_file, " has no results to compare with.", call. = FALSE)
if (any(.chk$N != .n) || any(.chk$LD != sum(ND_Y == 0)) || any(.chk$HD != sum(ND_Y == 1)))
  stop(sprintf(paste("This session's cohort (N %d, LD %d, HD %d) differs from the Excel file's (N %s, LD %s, HD %s).",
                     "Run script 01 with the settings of script 02's run, or import the matching Excel file."),
               .n, sum(ND_Y == 0), sum(ND_Y == 1), .chk$N[1], .chk$LD[1], .chk$HD[1]), call. = FALSE)
ND_MARKERS <- intersect(biomarkers_PH_vif, colnames(ND))
.panel_markers <- trimws(unlist(strsplit(as.character(c(.apparent$Best_panel[.apparent$Set %in% names(NESTED_SET_LIST)],
                                                        .app_size$Panel[.app_size$Model %in% names(NESTED_SET_LIST)])),
                                         ";", fixed = TRUE)))
.panel_markers <- setdiff(unique(.panel_markers[!is.na(.panel_markers) & nzchar(.panel_markers)]), ND_MARKERS)
if (length(.panel_markers))
  stop("Script 02's panels contain marker(s) that are not candidate markers in this session (",
       paste(head(.panel_markers, 5), collapse = ", "), if (length(.panel_markers) > 5) ", ...",
       "). Run script 01 with the settings of script 02's run, or import the matching Excel file.", call. = FALSE)
XM <- as.matrix(ND[, ND_MARKERS, drop = FALSE]); storage.mode(XM) <- "double"
if (any(!is.finite(XM))) stop("Missing or non-finite preprocessed biomarkers.", call. = FALSE)
.max_k <- min(NESTED_MAX_MARKERS, ncol(XM))

# Intercept + covariates (categorical ones as indicators) for one set, all rows
.base_matrix <- function(covs) {
  if (!length(covs)) return(matrix(1, .n, 1, dimnames = list(NULL, "(Intercept)")))
  mm <- stats::model.matrix(~ ., data = ND[, covs, drop = FALSE])
  if (nrow(mm) != .n) stop("Covariates have missing values in the cohort.")
  mm
}

cat(sprintf("\n===== Nested cross-validation: %d participants (LD %d, HD %d), %d candidate markers, %d covariate set(s) =====\n",
            .n, sum(ND_Y == 0), sum(ND_Y == 1), length(ND_MARKERS), length(NESTED_SET_LIST)))
for (.s in NESTED_SET_LIST)
  cat(sprintf("  %-28s %s %s\n", .s$name, if (.s$biomarkers) "panel search +" else "covariates only:",
              if (length(.s$covariates)) paste(.s$covariates, collapse = ", ") else "none"))
cat(sprintf("  Outer: %d folds x %d repeat(s) | inside: %d-fold CV x %d repeats, beam %d, up to %d markers\n",
            NESTED_OUTER_FOLDS, NESTED_OUTER_REPEATS, NESTED_INNER_FOLDS, NESTED_INNER_REPEATS, NESTED_BEAM_WIDTH, .max_k))


# ==========================================================================
# 4) THE SEARCH ENGINE (the same calculations as script 02's sweep)
# ==========================================================================

# Logistic fit with glm.fit, used when IRLS cannot finish cleanly
.nest_glmfit <- function(X, y, start = NULL) {
  f <- tryCatch(suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(), start = start,
                                                control = list(epsilon = 1e-8, maxit = 25L))),
                error = function(e) NULL)
  if (!is.null(start) &&
      (is.null(f) || !isTRUE(f$converged) || anyNA(f$coefficients) || max(abs(X %*% f$coefficients)) > 30))
    return(.nest_glmfit(X, y))
  if (is.null(f) || anyNA(f$coefficients)) return(NULL)
  list(coef = f$coefficients, conv = isTRUE(f$converged))
}

# Logistic fit by IRLS (normal equations, Cholesky); |eta| > 30 hands over to glm.fit
.nest_fit <- function(X, y, start = NULL, maxit = 30L, eps = 1e-8) {
  p <- ncol(X)
  b <- if (is.null(start) || length(start) != p) numeric(p) else start
  dev_old <- Inf
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% b)
    if (max(abs(eta)) > 30) return(.nest_glmfit(X, y))
    mu  <- 1 / (1 + exp(-eta))
    dev <- -2 * sum(y * log(mu) + (1 - y) * log1p(-mu))
    if (!is.finite(dev)) return(.nest_glmfit(X, y))
    if (abs(dev - dev_old) < eps * (abs(dev) + 0.1)) return(list(coef = b, conv = TRUE))
    dev_old <- dev
    w <- mu * (1 - mu); w[w < 1e-10] <- 1e-10
    A <- crossprod(X, X * w); bb <- crossprod(X, w * eta + (y - mu))
    R <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(R)) return(.nest_glmfit(X, y))
    b <- as.vector(backsolve(R, backsolve(R, bb, transpose = TRUE)))
  }
  .nest_glmfit(X, y)
}

# Pooled ROC metrics: balanced accuracy, AUC, Youden threshold (ties: closest
# to 0.5, then the higher), sensitivity, specificity. Case when p >= threshold.
.nest_roc <- function(p, y) {
  if (is.null(p) || anyNA(p)) return(rep(NA_real_, 5L))
  np <- sum(y); nn <- length(y) - np
  o  <- sort.list(p, method = "radix", decreasing = TRUE)
  ps <- p[o]; ys <- y[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  m <- length(ps); keep <- c(ps[-1L] != ps[-m], TRUE)
  tp <- tp[keep]; fp <- fp[keep]; s <- ps[keep]; nk <- length(tp)
  auc <- sum((c(0, tp[-nk]) + tp) * (fp - c(0, fp[-nk]))) / (2 * np * nn)
  tp <- c(0, tp); fp <- c(0, fp)
  thr  <- c(Inf, if (nk > 1) (s[-nk] + s[-1L]) / 2, -Inf)
  sens <- tp / np; spec <- 1 - fp / nn; J <- sens + spec
  tied <- which(J >= max(J) - 0.5 / (np * nn))
  j <- tied[order(abs(thr[tied] - 0.5), -thr[tied])[1]]
  unname(c((sens[j] + spec[j]) / 2, auc, thr[j], sens[j], spec[j]))
}

# Training and test matrices for every fold of one training set
# B: intercept + covariates; X: markers; folds: list(tr, te) within these rows
.nest_prepare <- function(B, X, y, folds) {
  Xall <- cbind(B, X)
  list(nb = ncol(B),
       Xtr = lapply(folds, function(f) Xall[f$tr, , drop = FALSE]),
       Xte = lapply(folds, function(f) Xall[f$te, , drop = FALSE]),
       ytr = lapply(folds, function(f) y[f$tr]),
       nte = vapply(folds, function(f) length(f$te), integer(1)),
       ypool = unlist(lapply(folds, function(f) y[f$te]), use.names = FALSE))
}

# Pooled out-of-fold predictions for the marker columns `cols` (each fold
# warm-started from the previous fold's fit). binomial()$linkinv, as in script
# 02, so saturated predictions tie the same way.
.nest_cv_pred <- function(cols, D) {
  sel <- c(seq_len(D$nb), D$nb + cols)
  linkinv <- stats::binomial()$linkinv
  out <- numeric(length(D$ypool)); pos <- 0L; b <- NULL
  for (f in seq_along(D$Xtr)) {
    fit <- .nest_fit(D$Xtr[[f]][, sel, drop = FALSE], D$ytr[[f]], start = b)
    if (is.null(fit)) return(NULL)
    b  <- if (fit$conv) fit$coef else NULL
    nt <- D$nte[f]
    out[pos + seq_len(nt)] <- linkinv(as.vector(D$Xte[[f]][, sel, drop = FALSE] %*% fit$coef))
    pos <- pos + nt
  }
  out
}
.nest_score_one <- function(cols, D) .nest_roc(.nest_cv_pred(cols, D), D$ypool)

# Scores of the combinations in `cmb` (one per column) on a worker, using the
# design .NW loaded there
.nest_score_cmb <- function(cmb) {
  out <- matrix(NA_real_, ncol(cmb), 5L)
  for (j in seq_len(ncol(cmb))) out[j, ] <- .nest_score_one(cmb[, j], .NW)
  out
}
.nest_score_par <- function(cl, cmb, chunk = 1000L) {
  if (!ncol(cmb)) return(matrix(NA_real_, 0L, 5L))
  starts <- seq(1L, ncol(cmb), by = chunk)
  pieces <- lapply(starts, function(s) cmb[, s:min(s + chunk - 1L, ncol(cmb)), drop = FALSE])
  do.call(rbind, parallel::parLapplyLB(cl, pieces, .nest_score_cmb))
}

# Best n by balanced accuracy, then AUC, with their scores; failed fits dropped
.nest_keep <- function(cmb, sc, n) {
  o <- order(-sc[, 1], -sc[, 2], na.last = NA)
  o <- o[seq_len(min(n, length(o)))]
  list(cmb = cmb[, o, drop = FALSE], sc = sc[o, , drop = FALSE])
}

# Beam step: every kept panel extended by each marker it lacks, duplicates removed
.nest_extend <- function(parents, p) {
  km1 <- nrow(parents)
  pieces <- lapply(seq_len(ncol(parents)), function(i) {
    par <- parents[, i]; js <- setdiff(seq_len(p), par)
    if (!length(js)) return(NULL)
    M <- apply(rbind(matrix(par, nrow = km1, ncol = length(js)), js), 2L, sort.int)
    if (!is.matrix(M)) M <- matrix(M, nrow = km1 + 1L)
    M
  })
  cand <- do.call(cbind, pieces)
  if (is.null(cand)) return(NULL)
  dimnames(cand) <- NULL; storage.mode(cand) <- "integer"
  cand[, !duplicated(t(cand)), drop = FALSE]
}

# Script 02's search on the design loaded on the workers: all single markers
# and pairs, then beam search; returns the kept panels and scores at each size
.nest_search <- function(cl, p) {
  kept <- list(); parents <- NULL
  for (k in seq_len(.max_k)) {
    cmb <- if (k == 1L) matrix(seq_len(p), nrow = 1L) else if (k == 2L) utils::combn(p, 2L) else .nest_extend(parents, p)
    if (is.null(cmb) || !ncol(cmb)) break
    storage.mode(cmb) <- "integer"
    kept[[k]] <- .nest_keep(cmb, .nest_score_par(cl, cmb), NESTED_BEAM_WIDTH)
    parents <- kept[[k]]$cmb
    if (!ncol(parents)) break
  }
  kept
}

# Choice as in script 02: highest balanced accuracy, then AUC, then fewer
# markers; also the best panel of each size. Threshold: that panel's Youden
# threshold from the same predictions.
.nest_choose <- function(kept) {
  cand <- do.call(rbind, lapply(seq_along(kept), function(k) {
    kk <- kept[[k]]
    if (is.null(kk) || !ncol(kk$cmb)) return(NULL)
    data.frame(k = k, j = seq_len(ncol(kk$cmb)), BA = kk$sc[, 1], AUC = kk$sc[, 2], Threshold = kk$sc[, 3])
  }))
  if (is.null(cand) || !nrow(cand)) stop("No valid panel fits in an outer training set.")
  cand <- cand[order(-cand$BA, -cand$AUC, cand$k), , drop = FALSE]
  pick <- function(r) {
    cols <- kept[[r$k]]$cmb[, r$j]
    list(cols = cols, markers = ND_MARKERS[cols], k = r$k, inner_BA = r$BA, inner_AUC = r$AUC, threshold = r$Threshold)
  }
  list(best = pick(cand[1, ]), by_size = lapply(split(cand, cand$k), function(d) pick(d[1, ])))
}

# Fit on the training rows, predict the test rows
.nest_fit_predict <- function(Btr, Xtr, ytr, Bte, Xte) {
  fit <- .nest_fit(cbind(Btr, Xtr), ytr)
  if (is.null(fit)) return(rep(NA_real_, nrow(Bte)))
  stats::binomial()$linkinv(as.vector(cbind(Bte, Xte) %*% fit$coef))
}

# Stratified folds: LD and HD each spread evenly; list(tr, te) per fold
.nest_folds <- function(y, k, repeats, seed) {
  set.seed(seed); out <- list()
  for (r in seq_len(repeats)) {
    ord <- unlist(lapply(split(seq_along(y), y), function(ii) ii[sample.int(length(ii))]), use.names = FALSE)
    grp <- integer(length(y)); grp[ord] <- rep_len(seq_len(k), length(ord))
    for (g in seq_len(k)) out[[sprintf("Fold%d.Rep%02d", g, r)]] <- list(tr = which(grp != g), te = which(grp == g))
  }
  out
}
.foldid <- function(folds, n) {
  first <- folds[endsWith(names(folds), ".Rep01")]
  id <- integer(n); for (g in seq_along(first)) id[first[[g]]$te] <- g
  id
}


# ==========================================================================
# 5) BENCHMARKS (each tuned on the outer training rows only)
# ==========================================================================

# glmnet: lambda by inner CV (deviance) on the first inner repeat's folds;
# covariates unpenalized
.bench_glmnet <- function(Ztr, Xtr, ytr, Zte, Xte, alpha, foldid) {
  xtr <- cbind(Xtr, Ztr); xte <- cbind(Xte, Zte)
  colnames(xtr) <- colnames(xte) <- c(ND_MARKERS, colnames(Ztr))
  pf <- c(rep(1, ncol(Xtr)), rep(0, ncol(Ztr)))
  cv <- suppressWarnings(glmnet::cv.glmnet(xtr, ytr, family = "binomial", alpha = alpha, penalty.factor = pf,
                                           foldid = foldid, keep = TRUE, type.measure = "deviance"))
  beta_m <- as.matrix(cv$glmnet.fit$beta[seq_len(ncol(Xtr)), , drop = FALSE])
  if (ncol(beta_m) != length(cv$lambda) || length(cv$cvm) != length(cv$lambda))
    stop("Unexpected glmnet output: lambda path lengths differ.")
  nz <- colSums(beta_m != 0)
  i <- match(cv$lambda.min, cv$lambda)
  pin <- stats::plogis(cv$fit.preval[, i]); ok <- is.finite(pin)   # prevalidated fits are on the logit scale
  list(p = as.numeric(predict(cv$glmnet.fit, newx = xte, s = cv$lambda[i], type = "response")),
       threshold = .nest_roc(pin[ok], ytr[ok])[3], n_markers = unname(nz[i]),
       markers = paste(ND_MARKERS[beta_m[, i] != 0], collapse = "; "))
}

# Random forest: threshold from the out-of-bag predictions
.bench_rf <- function(tr, te, covs, seed) {
  vars <- c(ND_MARKERS, covs)
  fit <- ranger::ranger(x = ND[tr, vars, drop = FALSE], y = ND$condL[tr], probability = TRUE,
                        num.trees = RF_TREES, seed = seed, num.threads = NESTED_CORES)
  pin <- fit$predictions[, "HD"]; ok <- is.finite(pin)
  list(p = as.numeric(stats::predict(fit, data = ND[te, vars, drop = FALSE], num.threads = NESTED_CORES)$predictions[, "HD"]),
       threshold = .nest_roc(pin[ok], ND_Y[tr][ok])[3], n_markers = length(ND_MARKERS), markers = "all")
}

# Metrics of pooled held-out predictions; `thr` holds each prediction's own
# threshold (chosen in its training rows)
.nest_metrics <- function(p, y, thr) {
  st <- .nest_roc(p, y)
  sens <- mean(p[y == 1] >= thr[y == 1]); spec <- mean(p[y == 0] < thr[y == 0])
  s5 <- mean(p[y == 1] >= 0.5); c5 <- mean(p[y == 0] < 0.5)
  brier <- mean((p - y)^2); ybar <- mean(y)
  lp <- stats::qlogis(pmin(pmax(p, NESTED_CAL_EPS), 1 - NESTED_CAL_EPS))
  a <- tryCatch(unname(stats::coef(suppressWarnings(stats::glm(y ~ 1, offset = lp, family = stats::binomial())))[1]),
                error = function(e) NA_real_)
  b <- tryCatch(unname(stats::coef(suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())))[2]),
                error = function(e) NA_real_)
  c(AUC = st[2], Balanced_Accuracy = (sens + spec) / 2, Sensitivity = sens, Specificity = spec,
    Balanced_Accuracy_at_0.5 = (s5 + c5) / 2, Sensitivity_at_0.5 = s5, Specificity_at_0.5 = c5,
    Brier = brier, Scaled_Brier = 1 - brier / (ybar * (1 - ybar)),
    Calibration_intercept = a, Calibration_slope = b)
}

.fmt_secs <- function(s) if (s >= 3600) sprintf("%dh %dm", s %/% 3600, round((s %% 3600) / 60)) else
  if (s >= 60) sprintf("%dm %ds", s %/% 60, round(s %% 60)) else sprintf("%ds", round(s))
.fmt_time <- function(t0) .fmt_secs(as.numeric(difftime(Sys.time(), t0, units = "secs")))


# ==========================================================================
# 6) THE NESTED LOOP
# ==========================================================================

.methods <- c(Panel_search = "Panel search", Covariates_only = "Covariates only", Lasso = "Lasso",
              Elastic_net = "Elastic net (alpha 0.5)", Random_forest = "Random forest")
.size_method <- function(k) sprintf("Best %d-marker panel", k)
NESTED_PRED <- list(); NESTED_PANELS <- list(); NESTED_BENCH <- list()
.add_pred <- function(set, method, r, f, te, res) {
  if (is.null(res)) return(invisible(NULL))
  NESTED_PRED[[length(NESTED_PRED) + 1L]] <<- data.frame(
    Set = set, Method = method, Outer_repeat = r, Outer_fold = f, Row = te,
    Source_row = ND_SOURCE_ROWS[te], Outcome = as.character(ND$condL[te]), y = ND_Y[te],
    p = res$p, Threshold = res$threshold, stringsAsFactors = FALSE)
}
.bases <- lapply(NESTED_SET_LIST, function(s) .base_matrix(s$covariates))
.projected <- FALSE

.t_all <- Sys.time()
.cl <- parallel::makeCluster(NESTED_CORES)
tryCatch({
  parallel::clusterExport(.cl, c(".nest_glmfit", ".nest_fit", ".nest_roc", ".nest_cv_pred", ".nest_score_one",
                                 ".nest_score_cmb"), envir = globalenv())
  for (.r in seq_len(NESTED_OUTER_REPEATS)) {
    .outer <- .nest_folds(ND$condL, NESTED_OUTER_FOLDS, 1L, .SEED + 1000L * .r)
    for (.f in seq_along(.outer)) {
      .te <- .outer[[.f]]$te; .tr <- .outer[[.f]]$tr; .ytr <- ND_Y[.tr]
      .seed  <- .SEED + 100L * .r + .f
      .inner <- .nest_folds(ND$condL[.tr], NESTED_INNER_FOLDS, NESTED_INNER_REPEATS, .seed)
      .fid   <- .foldid(.inner, length(.tr))
      cat(sprintf("\n--- Outer repeat %d, fold %d: train %d (HD %d), test %d (HD %d) ---\n", .r, .f,
                  length(.tr), sum(.ytr), length(.te), sum(ND_Y[.te])))
      for (.s in NESTED_SET_LIST) {
        .t0 <- Sys.time(); .B <- .bases[[.s$name]]
        .Btr <- .B[.tr, , drop = FALSE]; .Bte <- .B[.te, , drop = FALSE]
        if (.s$biomarkers) {
          parallel::clusterCall(.cl, function(D) { assign(".NW", D, envir = globalenv()); invisible(NULL) },
                                .nest_prepare(.Btr, XM[.tr, , drop = FALSE], .ytr, .inner))
          .ch <- .nest_choose(.nest_search(.cl, ncol(XM)))
          .choices <- c(list(.ch$best), .ch$by_size)
          names(.choices) <- c(.methods[["Panel_search"]], .size_method(as.integer(names(.ch$by_size))))
          for (.lab in names(.choices)) {
            .sel <- .choices[[.lab]]
            .p <- .nest_fit_predict(.Btr, XM[.tr, .sel$cols, drop = FALSE], .ytr, .Bte, XM[.te, .sel$cols, drop = FALSE])
            .add_pred(.s$name, .lab, .r, .f, .te, list(p = .p, threshold = .sel$threshold))
            NESTED_PANELS[[length(NESTED_PANELS) + 1L]] <- data.frame(
              Set = .s$name, Choice = .lab, Outer_repeat = .r, Outer_fold = .f, N_train = length(.tr),
              N_test = length(.te), Panel_size = .sel$k, Panel = paste(.sel$markers, collapse = "; "),
              Inner_Balanced_Accuracy = .sel$inner_BA, Inner_AUC = .sel$inner_AUC, Inner_Threshold = .sel$threshold,
              Heldout_AUC = .nest_roc(.p, ND_Y[.te])[2], stringsAsFactors = FALSE)
          }
          if (NESTED_BENCHMARKS) {
            .Ztr <- .Btr[, -1, drop = FALSE]; .Zte <- .Bte[, -1, drop = FALSE]
            .bench <- list(
              Lasso       = .bench_glmnet(.Ztr, XM[.tr, , drop = FALSE], .ytr, .Zte, XM[.te, , drop = FALSE], 1, .fid),
              Elastic_net = .bench_glmnet(.Ztr, XM[.tr, , drop = FALSE], .ytr, .Zte, XM[.te, , drop = FALSE], 0.5, .fid),
              Random_forest = .bench_rf(.tr, .te, .s$covariates, .seed))
            for (.m in names(.bench)) {
              .add_pred(.s$name, .methods[[.m]], .r, .f, .te, .bench[[.m]])
              if (!is.null(.bench[[.m]]))
                NESTED_BENCH[[length(NESTED_BENCH) + 1L]] <- data.frame(
                  Set = .s$name, Method = .methods[[.m]], Outer_repeat = .r, Outer_fold = .f,
                  N_markers = .bench[[.m]]$n_markers, Markers = .bench[[.m]]$markers, stringsAsFactors = FALSE)
            }
          }
          .last <- NESTED_PANELS[[length(NESTED_PANELS) - length(.choices) + 1L]]
          cat(sprintf("  %-28s %d marker(s) [%s] | inner BA %.3f, AUC %.3f | held-out AUC %.3f | %s\n", .s$name,
                      .ch$best$k, paste(.ch$best$markers, collapse = ", "), .ch$best$inner_BA, .ch$best$inner_AUC,
                      .last$Heldout_AUC, .fmt_time(.t0)))
          if (!.projected) {
            .projected <- TRUE
            .proj <- as.numeric(difftime(Sys.time(), .t0, units = "secs")) * length(.bm_sets) *
              NESTED_OUTER_FOLDS * NESTED_OUTER_REPEATS
            cat(sprintf("  Projected run time: about %s in all, from this set\n", .fmt_secs(.proj)))
          }
        } else {
          .thr <- .nest_score_one(integer(0), .nest_prepare(.Btr, XM[.tr, 0, drop = FALSE], .ytr, .inner))[3]
          .p   <- .nest_fit_predict(.Btr, XM[.tr, 0, drop = FALSE], .ytr, .Bte, XM[.te, 0, drop = FALSE])
          .add_pred(.s$name, .methods[["Covariates_only"]], .r, .f, .te, list(p = .p, threshold = .thr))
          cat(sprintf("  %-28s covariates only | held-out AUC %.3f\n", .s$name, .nest_roc(.p, ND_Y[.te])[2]))
        }
      }
    }
  }
}, finally = { parallel::stopCluster(.cl) })
NESTED_PRED   <- do.call(rbind, NESTED_PRED)
NESTED_PANELS <- if (length(NESTED_PANELS)) do.call(rbind, NESTED_PANELS) else data.frame()
NESTED_BENCH  <- if (length(NESTED_BENCH)) do.call(rbind, NESTED_BENCH) else data.frame()


# ==========================================================================
# 7) RESULTS
# ==========================================================================

.metric_names <- c("AUC", "Balanced_Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy_at_0.5",
                   "Sensitivity_at_0.5", "Specificity_at_0.5", "Brier", "Scaled_Brier",
                   "Calibration_intercept", "Calibration_slope")
.combos <- unique(NESTED_PRED[, c("Set", "Method")])

# per outer repeat (each participant predicted once per repeat); a method with
# a failed fit in any outer fold is left out of that repeat, with a warning
.incomplete <- character(0)
NESTED_BY_REPEAT <- do.call(rbind, lapply(seq_len(nrow(.combos)), function(i) {
  do.call(rbind, lapply(seq_len(NESTED_OUTER_REPEATS), function(r) {
    d <- NESTED_PRED[NESTED_PRED$Set == .combos$Set[i] & NESTED_PRED$Method == .combos$Method[i] &
                       NESTED_PRED$Outer_repeat == r, ]
    if (nrow(d) != .n || anyDuplicated(d$Row) || anyNA(d$p) || anyNA(d$Threshold)) {
      .incomplete <<- c(.incomplete, sprintf("%s / %s (outer repeat %d)", .combos$Set[i], .combos$Method[i], r))
      return(NULL)
    }
    data.frame(Set = .combos$Set[i], Method = .combos$Method[i], Outer_repeat = r,
               t(.nest_metrics(d$p, d$y, d$Threshold)), check.names = FALSE)
  }))
}))
if (is.null(NESTED_BY_REPEAT)) stop("No method produced predictions for every participant.", call. = FALSE)
if (length(.incomplete))
  warning("Left out (a fit failed or no threshold in at least one outer fold): ", paste(.incomplete, collapse = "; "),
          call. = FALSE)

.size_of <- function(method) suppressWarnings(as.integer(sub("^Best ([0-9]+)-marker panel$", "\\1", method)))
.markers_used <- function(set, method) {
  if (method == .methods[["Panel_search"]])
    return(NESTED_PANELS$Panel_size[NESTED_PANELS$Set == set & NESTED_PANELS$Choice == method])
  if (method == .methods[["Covariates_only"]]) return(0L)
  if (!is.na(.size_of(method))) return(.size_of(method))
  if (nrow(NESTED_BENCH)) NESTED_BENCH$N_markers[NESTED_BENCH$Set == set & NESTED_BENCH$Method == method] else NA_integer_
}
.cov_text <- function(set) { s <- .all_sets[[set]]; if (length(s$covariates)) paste(s$covariates, collapse = ", ") else "none" }

# nested results (mean over outer repeats), long format
NESTED_PERF <- do.call(rbind, lapply(seq_len(nrow(.combos)), function(i) {
  set <- .combos$Set[i]; method <- .combos$Method[i]
  d <- NESTED_BY_REPEAT[NESTED_BY_REPEAT$Set == set & NESTED_BY_REPEAT$Method == method, .metric_names, drop = FALSE]
  if (!nrow(d)) return(NULL)
  mk <- .markers_used(set, method)
  out <- data.frame(Set = set, Covariates = .cov_text(set), Method = method, Estimate = "Nested CV",
                    Markers_median = stats::median(mk), Markers_min = min(mk), Markers_max = max(mk),
                    t(colMeans(d)), check.names = FALSE, stringsAsFactors = FALSE)
  if (NESTED_OUTER_REPEATS > 1) for (m in c("AUC", "Balanced_Accuracy")) {
    out[[paste0(m, "_repeat_min")]] <- min(d[[m]]); out[[paste0(m, "_repeat_max")]] <- max(d[[m]])
  }
  out
}))

# apparent results copied from script 02's Excel file
.app_rows <- function(tab, set_col, method, size_col) {
  if (!nrow(tab)) return(NULL)
  tab <- tab[tab[[set_col]] %in% names(NESTED_SET_LIST), , drop = FALSE]
  if (!nrow(tab)) return(NULL)
  keep <- intersect(c("AUC", "Balanced_Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy_at_0.5",
                      "Sensitivity_at_0.5", "Specificity_at_0.5"), names(tab))
  data.frame(Set = tab[[set_col]], Covariates = vapply(tab[[set_col]], .cov_text, character(1)), Method = method,
             Estimate = "Apparent (script 02)", Markers_median = if (is.null(size_col)) 0 else tab[[size_col]],
             tab[, keep, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
}
# (the Overview and Covariates_Only sheets hold different metric columns)
.app_list <- Filter(Negate(is.null), list(.app_rows(.apparent, "Set", .methods[["Panel_search"]], "Best_panel_size"),
                                          .app_rows(.app_cov, "Model", .methods[["Covariates_only"]], NULL)))
.fill <- function(x, cols) { for (cc in setdiff(cols, names(x))) x[[cc]] <- NA; x[, cols, drop = FALSE] }
.cols <- unique(c(names(NESTED_PERF), unlist(lapply(.app_list, names))))
NESTED_ALL <- do.call(rbind, lapply(c(.app_list, list(NESTED_PERF)), .fill, cols = .cols))
.method_order <- c(.methods, .size_method(seq_len(.max_k)))
NESTED_ALL <- NESTED_ALL[order(match(NESTED_ALL$Set, names(NESTED_SET_LIST)), match(NESTED_ALL$Method, .method_order),
                               NESTED_ALL$Estimate != "Apparent (script 02)"), ]
rownames(NESTED_ALL) <- NULL
NESTED_PERFORMANCE <- NESTED_ALL[NESTED_ALL$Method %in% .methods, , drop = FALSE]
rownames(NESTED_PERFORMANCE) <- NULL
.get <- function(set, method, est, metric) {
  v <- NESTED_ALL[[metric]][NESTED_ALL$Set == set & NESTED_ALL$Method == method & NESTED_ALL$Estimate == est]
  if (length(v)) v[1] else NA_real_
}

# script 02's panels: best overall (Overview) and best of each size (Best_By_Size)
.app_size_row <- function(set, k)
  if (nrow(.app_size)) .app_size[.app_size$Model == set & .app_size$N_Biomarkers == k, , drop = FALSE] else data.frame()
.script02_panel <- function(set, choice) {
  txt <- if (choice == .methods[["Panel_search"]]) {
    r <- .apparent[.apparent$Set == set, , drop = FALSE]; if (nrow(r)) r$Best_panel[1] else NA
  } else {
    a <- .app_size_row(set, .size_of(choice)); if (nrow(a)) a$Panel[1] else NA
  }
  if (is.na(txt)) character(0) else trimws(strsplit(as.character(txt), ";", fixed = TRUE)[[1]])
}
.share_chosen <- function(set, choice) {
  pp <- NESTED_PANELS[NESTED_PANELS$Set == set & NESTED_PANELS$Choice == choice, , drop = FALSE]
  fp <- .script02_panel(set, choice)
  if (!length(fp) || !nrow(pp)) return(NA_real_)
  chosen <- strsplit(pp$Panel, "; ", fixed = TRUE)
  mean(vapply(fp, function(m) mean(vapply(chosen, function(x) m %in% x, logical(1))), numeric(1)))
}
if (nrow(NESTED_PANELS))
  NESTED_PANELS$Markers_shared_with_script02_panel <- vapply(seq_len(nrow(NESTED_PANELS)), function(i)
    sum(strsplit(NESTED_PANELS$Panel[i], "; ", fixed = TRUE)[[1]] %in%
          .script02_panel(NESTED_PANELS$Set[i], NESTED_PANELS$Choice[i])), integer(1))

# how often each marker was in the best panel, per set
.best_panels <- if (nrow(NESTED_PANELS)) NESTED_PANELS[NESTED_PANELS$Choice == .methods[["Panel_search"]], , drop = FALSE] else data.frame()
NESTED_FREQ <- if (nrow(.best_panels)) do.call(rbind, lapply(unique(.best_panels$Set), function(set) {
  pp <- .best_panels[.best_panels$Set == set, ]
  mk <- table(unlist(strsplit(pp$Panel, "; ", fixed = TRUE)))
  out <- data.frame(Set = set, Marker = names(mk), Times_chosen = as.integer(mk),
                    Share_of_outer_folds = as.integer(mk) / nrow(pp),
                    In_script02_panel = names(mk) %in% .script02_panel(set, .methods[["Panel_search"]]),
                    stringsAsFactors = FALSE)
  out[order(-out$Times_chosen, out$Marker), ]
})) else data.frame()

# one row per covariate set with biomarkers (best panel)
.ps <- .methods[["Panel_search"]]; .app_lab <- "Apparent (script 02)"
NESTED_OVERVIEW <- if (length(.bm_sets)) do.call(rbind, lapply(.bm_sets, function(set) {
  ar <- .apparent[.apparent$Set == set, , drop = FALSE]
  cov_model <- if (nrow(ar) && "Covariates_only_model" %in% names(ar)) as.character(ar$Covariates_only_model[1]) else NA_character_
  pp <- .best_panels[.best_panels$Set == set, , drop = FALSE]
  fr <- NESTED_FREQ[NESTED_FREQ$Set == set, , drop = FALSE]
  fp <- .script02_panel(set, .ps)
  data.frame(
    Set = set, Covariates = .cov_text(set), N = .n,
    Script02_panel = if (length(fp)) paste(fp, collapse = "; ") else NA,
    Apparent_AUC = .get(set, .ps, .app_lab, "AUC"),
    Nested_AUC = .get(set, .ps, "Nested CV", "AUC"),
    AUC_optimism = .get(set, .ps, .app_lab, "AUC") - .get(set, .ps, "Nested CV", "AUC"),
    Apparent_Balanced_Accuracy = .get(set, .ps, .app_lab, "Balanced_Accuracy"),
    Nested_Balanced_Accuracy = .get(set, .ps, "Nested CV", "Balanced_Accuracy"),
    Balanced_Accuracy_optimism = .get(set, .ps, .app_lab, "Balanced_Accuracy") - .get(set, .ps, "Nested CV", "Balanced_Accuracy"),
    Nested_Sensitivity = .get(set, .ps, "Nested CV", "Sensitivity"),
    Nested_Specificity = .get(set, .ps, "Nested CV", "Specificity"),
    Nested_Balanced_Accuracy_at_0.5 = .get(set, .ps, "Nested CV", "Balanced_Accuracy_at_0.5"),
    Nested_Brier = .get(set, .ps, "Nested CV", "Brier"),
    Nested_Scaled_Brier = .get(set, .ps, "Nested CV", "Scaled_Brier"),
    Nested_Calibration_slope = .get(set, .ps, "Nested CV", "Calibration_slope"),
    Covariates_only_model = cov_model,
    Nested_Covariates_only_AUC = if (is.na(cov_model)) NA_real_ else .get(cov_model, .methods[["Covariates_only"]], "Nested CV", "AUC"),
    Nested_AUC_added_by_biomarkers = if (is.na(cov_model)) NA_real_ else
      .get(set, .ps, "Nested CV", "AUC") - .get(cov_model, .methods[["Covariates_only"]], "Nested CV", "AUC"),
    Nested_AUC_Lasso = .get(set, .methods[["Lasso"]], "Nested CV", "AUC"),
    Nested_AUC_Elastic_net = .get(set, .methods[["Elastic_net"]], "Nested CV", "AUC"),
    Nested_AUC_Random_forest = .get(set, .methods[["Random_forest"]], "Nested CV", "AUC"),
    Panel_sizes_chosen = if (!nrow(pp)) NA else if (min(pp$Panel_size) == max(pp$Panel_size)) as.character(pp$Panel_size[1])
                         else paste(range(pp$Panel_size), collapse = " to "),
    Most_chosen_markers = if (nrow(fr)) paste(sprintf("%s (%d/%d)", head(fr$Marker, 3), head(fr$Times_chosen, 3), nrow(pp)),
                                              collapse = "; ") else NA,
    Script02_markers_chosen_share = .share_chosen(set, .ps),
    check.names = FALSE, stringsAsFactors = FALSE)
})) else data.frame()

# one row per covariate set with biomarkers and panel size (best panel of that size)
NESTED_BY_SIZE <- if (length(.bm_sets)) do.call(rbind, lapply(.bm_sets, function(set) do.call(rbind, lapply(seq_len(.max_k), function(k) {
  lab <- .size_method(k); a <- .app_size_row(set, k)
  ga <- function(m) if (nrow(a) && m %in% names(a)) as.numeric(a[[m]][1]) else NA_real_
  data.frame(Set = set, Covariates = .cov_text(set), Panel_size = k,
             Script02_panel = if (nrow(a)) as.character(a$Panel[1]) else NA_character_,
             Apparent_AUC = ga("AUC"), Nested_AUC = .get(set, lab, "Nested CV", "AUC"),
             AUC_optimism = ga("AUC") - .get(set, lab, "Nested CV", "AUC"),
             Apparent_Balanced_Accuracy = ga("Balanced_Accuracy"),
             Nested_Balanced_Accuracy = .get(set, lab, "Nested CV", "Balanced_Accuracy"),
             Balanced_Accuracy_optimism = ga("Balanced_Accuracy") - .get(set, lab, "Nested CV", "Balanced_Accuracy"),
             Nested_Sensitivity = .get(set, lab, "Nested CV", "Sensitivity"),
             Nested_Specificity = .get(set, lab, "Nested CV", "Specificity"),
             Nested_Balanced_Accuracy_at_0.5 = .get(set, lab, "Nested CV", "Balanced_Accuracy_at_0.5"),
             Nested_Scaled_Brier = .get(set, lab, "Nested CV", "Scaled_Brier"),
             Nested_Calibration_slope = .get(set, lab, "Nested CV", "Calibration_slope"),
             Script02_markers_chosen_share = .share_chosen(set, lab),
             check.names = FALSE, stringsAsFactors = FALSE)
})))) else data.frame()


# ==========================================================================
# 8) OUTPUT
# ==========================================================================

# Script 1's LD / HD cut-off (DMFT_CUTOFF_USED, set again by every run of script 1)
.cut <- get0("DMFT_CUTOFF_USED")
.group_text <- if (length(.cut) == 1 && is.numeric(.cut) && is.finite(.cut))
  sprintf("DMFT <= %g -> LD (low DMFT) | > %g -> HD (high DMFT) (%s)", .cut, .cut,
          if (is.null(get0("DMFT_CUTOFF"))) "median of the analysed subjects" else "fixed") else "not recorded"
.settings <- data.frame(Parameter = c(
  "Registry key", "LD / HD", "Script 02 file", "Participants (LD / HD)", "Candidate markers", "Covariate sets",
  "In script 02's file, not in script 01 (left out)",
  "Outer folds x repeats", "Search in each outer training set", "Choice", "Threshold", "Benchmarks",
  "Random forest trees", "Seed (script 01's SEED)", "Cores", "Run time"),
  Value = c(registry_key, .group_text, basename(.res_file), sprintf("%d (%d / %d)", .n, sum(ND_Y == 0), sum(ND_Y == 1)),
            length(ND_MARKERS), paste(names(NESTED_SET_LIST), collapse = ", "),
            if (length(.left_out)) paste(.left_out, collapse = ", ") else "none",
            sprintf("%d x %d", NESTED_OUTER_FOLDS, NESTED_OUTER_REPEATS),
            sprintf("every 1- and 2-marker panel, then beam search to %d markers keeping %d panels per size; each panel scored by %d-fold CV x %d repeats",
                    .max_k, NESTED_BEAM_WIDTH, NESTED_INNER_FOLDS, NESTED_INNER_REPEATS),
            "highest balanced accuracy, then AUC, then fewer markers (as script 02); also the best panel of each size",
            "Youden threshold of each model's inner (training-only) predictions; also a fixed 0.5",
            if (NESTED_BENCHMARKS) paste(.methods[c("Lasso", "Elastic_net", "Random_forest")], collapse = "; ") else "none",
            if (NESTED_BENCHMARKS) RF_TREES else NA, .SEED, NESTED_CORES, .fmt_time(.t_all)),
  stringsAsFactors = FALSE)
.notes <- data.frame(Note = c(
  "Nested CV: in each outer fold the panel search, panel choice, threshold and all tuning used the training participants only; the held-out participants were then predicted.",
  "Every covariate set was analysed the same way, on the same outer and inner folds.",
  "The search inside each outer fold is script 02's search, with the settings on the Settings sheet (the defaults are script 02's settings).",
  "Apparent results are copied from script 02's Excel file: panels and thresholds there were chosen on the reported predictions.",
  "Optimism = apparent minus nested estimate.",
  "Panels chosen in the outer folds can differ from each other and from script 02's panel; Marker_frequency and Script02_markers_chosen_share show how stable the choice is.",
  "Script 02's panel, built from all participants, is the panel to report; the nested estimates describe how a panel found this way performs in participants not used to find it."))
.out_file <- sprintf("DMFT_Nested_CV_%s.xlsx", registry_key)
writexl::write_xlsx(list(Overview = NESTED_OVERVIEW, Performance = NESTED_PERFORMANCE, By_size = NESTED_BY_SIZE,
                         Selected_panels = NESTED_PANELS, Marker_frequency = NESTED_FREQ,
                         Benchmark_markers = NESTED_BENCH, Performance_by_repeat = NESTED_BY_REPEAT,
                         Predictions = NESTED_PRED, Settings = .settings, Notes = .notes), .out_file)

.p1 <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))
.f3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
.text <- c(
  "METHODS:",
  paste0("Because panels and thresholds were chosen using the same out-of-fold predictions that were reported, ",
         "performance was also estimated with nested cross-validation [ref]. Participants were split into ",
         NESTED_OUTER_FOLDS, " outer folds", if (NESTED_OUTER_REPEATS > 1) sprintf(", repeated %d times", NESTED_OUTER_REPEATS),
         ". Within each outer training set, the whole selection procedure was repeated: all one- and two-marker panels ",
         "and beam search to ", .max_k, " markers (", NESTED_BEAM_WIDTH, " panels kept per size), each panel scored by ",
         NESTED_INNER_FOLDS, "-fold cross-validation repeated ", NESTED_INNER_REPEATS, " times; the panel with the ",
         "highest balanced accuracy was chosen, with its Youden threshold from the same inner predictions. ",
         "The chosen panel was fitted to the outer training set and applied to the held-out participants. ",
         "This was done for every covariate set on the same folds",
         if (NESTED_BENCHMARKS) paste0(", together with the lasso and elastic net (alpha 0.5) with the covariates ",
                                       "unpenalized and the penalty chosen by inner cross-validation [refs], and a ",
                                       "random forest with ", RF_TREES,
                                       " trees [ref], each tuned within the outer training sets") else "", "."),
  "",
  "RESULTS (one line per covariate set):",
  if (nrow(NESTED_OVERVIEW)) vapply(seq_len(nrow(NESTED_OVERVIEW)), function(i) {
    o <- NESTED_OVERVIEW[i, ]
    paste0(o$Set, " (covariates: ", o$Covariates, "): nested AUC ", .f3(o$Nested_AUC), " and balanced accuracy ",
           .p1(o$Nested_Balanced_Accuracy), " (apparent ", .f3(o$Apparent_AUC), " and ", .p1(o$Apparent_Balanced_Accuracy),
           ")", if (!is.na(o$Nested_Covariates_only_AUC)) paste0("; covariates alone AUC ", .f3(o$Nested_Covariates_only_AUC)) else "",
           if (NESTED_BENCHMARKS) paste0("; lasso ", .f3(o$Nested_AUC_Lasso), ", elastic net ", .f3(o$Nested_AUC_Elastic_net),
                                         ", random forest ", .f3(o$Nested_AUC_Random_forest)) else "",
           "; chosen panels had ", o$Panel_sizes_chosen, " markers; most often chosen: ", o$Most_chosen_markers, ".")
  }, character(1)) else "(no covariate sets with biomarkers)")
writeLines(.text, sprintf("DMFT_Nested_CV_text_%s.txt", registry_key), useBytes = TRUE)

cat("\n===== Nested cross-validation: best panel per covariate set =====\n")
if (nrow(NESTED_OVERVIEW))
  print(NESTED_OVERVIEW[, c("Set", "Apparent_AUC", "Nested_AUC", "Apparent_Balanced_Accuracy", "Nested_Balanced_Accuracy",
                            "Nested_Covariates_only_AUC", "Nested_AUC_Lasso", "Nested_AUC_Random_forest")],
        row.names = FALSE, digits = 3)
cat(sprintf("\n  Total time %s\n  Saved in %s:\n    %s\n    DMFT_Nested_CV_text_%s.txt\n", .fmt_time(.t_all),
            normalizePath(getwd()), .out_file, registry_key))
