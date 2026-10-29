# ==========================================================================
# DMFT SALIVARY BIOMARKER CALCULATOR (R Shiny, web + desktop)
# ==========================================================================
# Two tasks, chosen on the start page:
#
#   1. Data analysis - upload participant information and raw biomarker
#      measurements, match them by participant ID, define the two comparison
#      groups, search logistic regression panels with repeated cross-validation,
#      train a model from any panel and download the results.
#   2. Patient prediction - enter absolute raw concentration values to estimate
#      the risk of high DMFT with the selected age/BMI model, or
#      with a model trained in task 1 (or loaded from a saved model file).
#      The page shows the equation with its coefficients, the calculation for
#      the values entered, a speedometer gauge of the estimate and, for models
#      with age and BMI, risk curves across age at two BMI values.
#
# Run:   shiny::runApp("app.R")      (or open in RStudio and click Run App)
# Needs: shiny. Optional: readxl (Excel uploads) and writexl (Excel downloads);
#        without them, upload and download CSV files.
#
# Two versions are built from this file:
#   Desktop - this file, run in a local R session. A dmft_age_bmi_model.csv
#             (exported by 03_Figures_DMFT_Age_BMI.R) placed next to it replaces the
#             bundled coefficients after it passes the checks below.
#   Browser - app/app.R in the dmft-calculator repository: this file with the
#             bundled coefficients only (see the PUB line in section 11),
#             compiled by build.R with Shinylive to run R in the browser (webR).
#
# Study specification: age + BMI; exhaustive sizes 1–2, beam from 3.
# The previous age/sex model is intentionally not bundled. Load the new export.
# ==========================================================================

suppressPackageStartupMessages(library(shiny))
HAS_READXL  <- requireNamespace("readxl", quietly = TRUE)
HAS_WRITEXL <- requireNamespace("writexl", quietly = TRUE)
options(shiny.maxRequestSize = 100 * 1024^2)

MINUS <- "−"
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


# ==========================================================================
# 1) THE PUBLISHED MODEL
# ==========================================================================
# Raw-input form: logit p(HD) = intercept + sum_j coef_j * log2(x_j + delta)
#                               + coef_age * Age + coef_bmi * BMI
# (the log2 transform is applied to the raw value; the LD-group centring and
# scaling used in fitting are absorbed into the intercept and slopes).

PUBLISHED_MODEL <- list(
  source = "published", name = "DMFT model (age + BMI; rerun required)",
  description = "Load the fitted age/BMI model exported after the revised panel search.",
  case_label = "high DMFT (DMFT > 13)", ref_label = "low DMFT (DMFT ≤ 13)",
  case_short = "HD", ref_short = "LD", transform = "log2", delta = 5e-7,
  intercept = NA_real_,
  markers = data.frame(key=character(), abbr=character(), label=character(), coef=numeric()),
  covariates = data.frame(var=c("Age_DOS","BMI"), term=c("Age_DOS","BMI"),
    abbr=c("Age","BMI"), label=c("Age","Body mass index"), unit=c("years","kg/m²"),
    type=c("numeric","numeric"), level=c(NA_character_,NA_character_), coef=c(NA_real_,NA_real_)),
  cov_levels=list(), level_labels=list(), threshold=NA_real_,
  cohort="Age/BMI model results are pending the revised analysis.",
  performance=NULL, performance_note=NULL,
  caution="Load dmft_age_bmi_model.csv from the rerun. Previous age/sex coefficients and performance do not apply."
)

# Display label of a factor level (the level itself when the model has no labels)
level_label <- function(m, v, lev) {
  labs <- m$level_labels[[v]]
  out <- if (!is.null(labs)) unname(labs[as.character(lev)]) else rep(NA_character_, length(lev))
  ifelse(is.na(out), as.character(lev), out)
}
# Choices for a factor input: values are the model's levels, names their labels
level_choices <- function(m, v, lv) stats::setNames(lv, level_label(m, v, lv))

EXAMPLE <- list(markers = c(1, 1, 1, 0, 0, 0), covariates = c(40, NA))

model_ready <- function(m) {
  !is.null(m) && nrow(m$markers) > 0L && is.finite(m$threshold) && is.finite(m$intercept) && all(is.finite(m$markers$coef)) &&
    (nrow(m$covariates) == 0 || all(is.finite(m$covariates$coef)))
}

# Load only an explicitly named, complete model export; never reuse a legacy file.
selected_age_bmi_model <- function(tab, filename = "dmft_age_bmi_model.csv") {
  m <- model_from_table(tab, filename)
  meta <- function(key) tab$Label[tab$Term == key]
  if (sum(tab$Term == "(Threshold)") != 1L ||
      !is.finite(as.numeric(tab$Coefficient_raw_input[tab$Term == "(Threshold)"])) ||
      m$threshold <= 0 || m$threshold >= 1 ||
      !model_ready(m) || !setequal(unique(m$covariates$var), c("Age_DOS","BMI")) ||
      !all(m$covariates$type == "numeric") || m$transform != "log2" ||
      abs(m$delta - 5e-7) > 1e-15 ||
      !identical(meta("(Model configuration)"), "Biomarkers_Age_BMI") ||
      !identical(meta("(Search design)"), "exhaustive12_beam3"))
    stop("Use dmft_age_bmi_model.csv exported after the new age/BMI search (exhaustive sizes 1–2, beam from 3).")
  m$analysis_configuration <- "Biomarkers_Age_BMI"; m$search_design <- "exhaustive12_beam3"
  m$source <- "published"; m$name <- "Selected DMFT model (age + BMI)"
  m$description <- sprintf("%d salivary biomarkers with age and BMI.", nrow(m$markers))
  m$covariates$label <- ifelse(m$covariates$var=="BMI","Body mass index","Age")
  m$covariates$unit <- ifelse(m$covariates$var=="BMI","kg/m²","years")
  m$cohort <- "See the rerun analysis workbook for cohort size and cross-validated performance."
  m$caution <- "Development model; independent validation is required."
  m
}
load_published_coefficients <- function(model, dir = getwd()) {
  f <- file.path(dir, "dmft_age_bmi_model.csv")
  if (!file.exists(f)) return(list(model=model, msg=NULL))
  tryCatch(list(model=selected_age_bmi_model(utils::read.csv(f, stringsAsFactors=FALSE),basename(f)),
                msg="Age/BMI model loaded."),
           error=function(e) list(model=model,msg=conditionMessage(e)))
}

# ==========================================================================
# 2) SMALL HELPERS
# ==========================================================================

# Number for display: 4 significant digits, Unicode minus, x 10^n for tiny or huge values
fmt_num <- function(x, sig = 4, html = TRUE) {
  if (length(x) != 1 || is.na(x)) return("NA")
  if (!is.finite(x)) return(if (x > 0) "∞" else paste0(MINUS, "∞"))
  if (x == 0) return("0")
  ax <- abs(x)
  s <- if (ax >= 1e-3 && ax < 1e7) {
    formatC(signif(x, sig), format = "fg", digits = sig, big.mark = ",")
  } else {
    e <- floor(log10(ax)); m <- signif(x / 10^e, sig)
    ex <- sub("-", MINUS, as.character(e))
    if (html) sprintf("%s × 10<sup>%s</sup>", formatC(m, format = "fg", digits = sig), ex)
    else sprintf("%s × 10^%s", formatC(m, format = "fg", digits = sig), ex)
  }
  sub("^-", MINUS, trimws(s))
}
fmt_pct <- function(p, d = 1) ifelse(is.na(p), "NA", sprintf("%.*f%%", d, 100 * p))
# A probability for display: never shows a misleading 0.0% or 100.0%
fmt_prob <- function(p) if (is.na(p)) "NA" else if (p > 0.999) "> 99.9%" else if (p < 0.001) "< 0.1%" else fmt_pct(p)
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_time <- function(sec) {
  sec <- max(0, round(sec))
  if (sec >= 3600) sprintf("%dh %02dm", sec %/% 3600, (sec %% 3600) %/% 60)
  else if (sec >= 60) sprintf("%dm %02ds", sec %/% 60, sec %% 60) else sprintf("%ds", sec)
}
# *text* -> <i>text</i>, everything else escaped
md_html <- function(x) gsub("\\*([^*]+)\\*", "<i>\\1</i>", htmltools::htmlEscape(x))
# Readable label from a column name
pretty_name <- function(x) gsub("_", " ", x, fixed = TRUE)
# "a", "a and b", "a, b and c"
and_list <- function(x) if (length(x) <= 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])


# ==========================================================================
# 3) READING AND MATCHING FILES
# ==========================================================================

read_table_file <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  na  <- c("", "NA", "N/A", "NaN", "NULL", ".", "-")
  df <- if (ext %in% c("xlsx", "xls")) {
    if (!HAS_READXL) stop("Reading Excel files needs the readxl package; install it or upload a CSV file.")
    as.data.frame(readxl::read_excel(path, guess_max = 100000, na = na), check.names = FALSE,
                  stringsAsFactors = FALSE)
  } else if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = na)
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = na)
  } else stop("Use a .csv, .tsv, .txt, .xlsx or .xls file.")
  if (!nrow(df) || !ncol(df)) stop("The file has no rows or columns.")
  names(df) <- make.unique(trimws(as.character(names(df))))
  df
}

guess_id_col <- function(df) {
  nm  <- names(df)
  hit <- grep("^(study|participant|subject|patient|sample)?[ _.]?id$", nm, ignore.case = TRUE)
  if (length(hit)) nm[hit[1]] else nm[1]
}

# Columns that are numeric, or character columns that are entirely numeric text
numeric_cols <- function(df, exclude = character(0)) {
  keep <- vapply(names(df), function(v) {
    x <- df[[v]]
    if (is.numeric(x)) return(any(is.finite(x)))
    if (is.character(x)) {
      y <- suppressWarnings(as.numeric(x[!is.na(x) & nzchar(trimws(x))]))
      return(length(y) > 0 && !anyNA(y))
    }
    FALSE
  }, logical(1))
  setdiff(names(df)[keep], exclude)
}
as_num <- function(x) if (is.numeric(x)) as.numeric(x) else suppressWarnings(as.numeric(as.character(x)))

match_files <- function(pinfo, pid, bio, bid) {
  a <- trimws(as.character(pinfo[[pid]])); b <- trimws(as.character(bio[[bid]]))
  dup_a <- unique(a[duplicated(a) & !is.na(a) & nzchar(a)])
  dup_b <- unique(b[duplicated(b) & !is.na(b) & nzchar(b)])
  if (length(dup_a) || length(dup_b))
    return(list(ok = FALSE, msg = paste0("Participant IDs must be unique. Duplicated: ",
                                         paste(head(c(dup_a, dup_b), 10), collapse = ", "))))
  common <- intersect(a[!is.na(a) & nzchar(a)], b)
  if (!length(common)) return(list(ok = FALSE, msg = "No participant ID appears in both files."))
  left  <- pinfo[match(common, a), , drop = FALSE]
  right <- bio[match(common, b), setdiff(names(bio), bid), drop = FALSE]
  clash <- intersect(names(right), names(left))
  if (length(clash)) names(right)[match(clash, names(right))] <- paste0(clash, " (biomarker file)")
  merged <- cbind(left, right)
  rownames(merged) <- NULL
  list(ok = TRUE, data = merged, id = pid, bio_cols = names(right),
       n_info = sum(!is.na(a) & nzchar(a)), n_bio = sum(!is.na(b) & nzchar(b)), n_matched = length(common),
       only_info = setdiff(a[!is.na(a) & nzchar(a)], b), only_bio = setdiff(b[!is.na(b) & nzchar(b)], a))
}


# ==========================================================================
# 4) COMPARISON GROUPS AND PREPROCESSING
# ==========================================================================

# Outcome value: one column, or the sum of several numeric columns (e.g. D + M + F = DMFT)
outcome_value <- function(df, cols) {
  if (length(cols) == 1) return(df[[cols]])
  m <- vapply(cols, function(v) as_num(df[[v]]), numeric(nrow(df)))
  if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
  rowSums(m)                                   # NA if any component is missing
}

# y = 1 (case), 0 (reference) or NA (in neither group)
define_groups <- function(value, s) {
  if (s$mode == "levels") {
    v <- trimws(as.character(value))
    y <- ifelse(v == s$case_level, 1L, ifelse(v == s$ref_level, 0L, NA_integer_))
    desc <- sprintf("%s = \"%s\" vs %s = \"%s\"", s$ref_name, s$ref_level, s$case_name, s$case_level)
  } else {
    v <- as_num(value)
    if (s$mode == "cutoff") {
      y <- ifelse(is.na(v), NA_integer_, as.integer(v > s$cutoff))
      desc <- sprintf("%s ≤ %s vs %s > %s", s$ref_name, format(s$cutoff), s$case_name, format(s$cutoff))
    } else {
      y <- ifelse(is.na(v), NA_integer_, ifelse(v <= s$ref_max, 0L, ifelse(v >= s$case_min, 1L, NA_integer_)))
      desc <- sprintf("%s ≤ %s vs %s ≥ %s (values in between left out)", s$ref_name, format(s$ref_max),
                      s$case_name, format(s$case_min))
    }
  }
  list(y = y, desc = desc)
}

# Missing biomarker values: k-nearest neighbours (Gower distance on the
# available markers, median of the k donors), as in VIM::kNN
impute_knn <- function(X, k = 5) {
  miss <- is.na(X)
  if (!any(miss)) return(X)
  rng <- apply(X, 2, function(v) { r <- suppressWarnings(diff(range(v, na.rm = TRUE))); if (is.finite(r) && r > 0) r else 1 })
  Xs  <- sweep(X, 2, rng, "/")
  for (i in which(rowSums(miss) > 0)) {
    obs <- !miss[i, ]
    if (any(obs)) {
      D <- abs(sweep(Xs[, obs, drop = FALSE], 2, Xs[i, obs], "-"))
      d <- rowMeans(D, na.rm = TRUE)
    } else d <- rep(0, nrow(X))
    d[i] <- Inf; d[!is.finite(d)] <- Inf
    for (j in which(miss[i, ])) {
      donors <- which(!miss[, j])
      if (!length(donors)) next
      nn <- donors[order(d[donors])][seq_len(min(k, length(donors)))]
      X[i, j] <- stats::median(X[nn, j])
    }
  }
  X
}

# Builds the analysis data: outcome, covariate design, transformed and
# reference-standardised biomarkers
prepare_analysis <- function(data, id_col, outcome_cols, groups, markers, covariates, s) {
  notes <- character(0)
  g <- define_groups(outcome_value(data, outcome_cols), groups)
  keep <- !is.na(g$y)
  if (sum(keep) < 10) stop("Fewer than 10 participants fall into the two comparison groups.")
  data <- data[keep, , drop = FALSE]; y <- g$y[keep]
  if (sum(!keep)) notes <- c(notes, sprintf("%d participant(s) outside both groups or with a missing outcome were left out.",
                                            sum(!keep)))

  # covariates: rows with a missing value are left out once, before the folds
  covariates <- setdiff(covariates, c(id_col, outcome_cols, markers))
  Z <- matrix(numeric(0), nrow(data), 0); cov_info <- data.frame()
  cov_levels <- list(); cov_numeric <- list()
  if (length(covariates)) {
    cdf <- data[, covariates, drop = FALSE]
    for (v in covariates) {
      x <- cdf[[v]]
      if (tolower(v) == "sex") {
        sx <- trimws(as.character(x)); sx[is.na(x) | sx == ""] <- NA_character_
        cdf[[v]] <- factor(sx)
      } else if (is.numeric(x) || all(is.na(x) | !is.na(suppressWarnings(as.numeric(as.character(x)))))) {
        cdf[[v]] <- as_num(x)
      } else cdf[[v]] <- factor(trimws(as.character(x)))
    }
    ok <- stats::complete.cases(cdf)
    if (sum(!ok)) notes <- c(notes, sprintf("%d participant(s) with a missing covariate were left out.", sum(!ok)))
    cdf <- cdf[ok, , drop = FALSE]; data <- data[ok, , drop = FALSE]; y <- y[ok]
    for (v in covariates) if (is.factor(cdf[[v]])) cdf[[v]] <- droplevels(cdf[[v]])
    const <- vapply(covariates, function(v) length(unique(cdf[[v]])) < 2, logical(1))
    if (any(const)) {
      notes <- c(notes, sprintf("Covariate(s) with a single value were dropped: %s.", paste(covariates[const], collapse = ", ")))
      covariates <- covariates[!const]; cdf <- cdf[, covariates, drop = FALSE]
    }
    for (v in covariates) {
      if (is.factor(cdf[[v]])) {
        lv <- levels(cdf[[v]]); cov_levels[[v]] <- lv
        for (l in lv[-1]) {
          Z <- cbind(Z, as.numeric(cdf[[v]] == l))
          cov_info <- rbind(cov_info, data.frame(var = v, term = paste0(v, "=", l), type = "factor", level = l,
                                                 stringsAsFactors = FALSE))
        }
      } else {
        Z <- cbind(Z, cdf[[v]])
        cov_info <- rbind(cov_info, data.frame(var = v, term = v, type = "numeric", level = NA_character_,
                                               stringsAsFactors = FALSE))
        cov_numeric[[v]] <- stats::median(cdf[[v]])
      }
    }
  }

  # biomarkers
  X <- vapply(markers, function(v) as_num(data[[v]]), numeric(nrow(data)))
  if (is.null(dim(X))) X <- matrix(X, nrow = nrow(data), dimnames = list(NULL, markers))
  empty <- colSums(!is.na(X)) == 0
  if (any(empty)) { notes <- c(notes, sprintf("%d biomarker column(s) with no values were dropped.", sum(empty)))
                    X <- X[, !empty, drop = FALSE] }
  n_miss <- sum(is.na(X))
  if (n_miss) {
    if (s$missing == "knn") {
      X <- impute_knn(X, 5); notes <- c(notes, sprintf("%s missing biomarker value(s) imputed by k-nearest neighbours (k = 5).", fmt_int(n_miss)))
    } else if (s$missing == "median") {
      for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- stats::median(X[, j], na.rm = TRUE)
      notes <- c(notes, sprintf("%s missing biomarker value(s) replaced by the column median.", fmt_int(n_miss)))
    } else {
      ok <- stats::complete.cases(X)
      notes <- c(notes, sprintf("%d participant(s) with a missing biomarker value were left out.", sum(!ok)))
      X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]; y <- y[ok]; data <- data[ok, , drop = FALSE]
    }
  }
  raw_median <- apply(X, 2, stats::median)
  if (s$transform == "log2") {
    if (any(X + s$delta <= 0, na.rm = TRUE))
      stop("The log2 transform needs raw values above −offset; some values are negative. Turn the transform off or check the data.")
    Tx <- log2(X + s$delta)
  } else Tx <- X
  mu <- rep(0, ncol(Tx)); sdv <- rep(1, ncol(Tx))
  if (isTRUE(s$standardize)) {
    ref <- y == 0
    mu  <- colMeans(Tx[ref, , drop = FALSE])
    sdv <- apply(Tx[ref, , drop = FALSE], 2, stats::sd)
    bad <- !is.finite(sdv) | sdv == 0
    if (any(bad)) {   # no spread in the reference group: use the SD of all participants
      sd_all <- apply(Tx, 2, stats::sd); sdv[bad] <- sd_all[bad]
      notes <- c(notes, sprintf("%d biomarker(s) had no spread in the reference group and were scaled by the SD of all participants.", sum(bad)))
    }
  }
  const <- !is.finite(sdv) | sdv == 0 | apply(Tx, 2, function(v) length(unique(v)) < 2)
  if (any(const)) {
    notes <- c(notes, sprintf("%d biomarker(s) with a single value were dropped.", sum(const)))
    X <- X[, !const, drop = FALSE]; Tx <- Tx[, !const, drop = FALSE]; mu <- mu[!const]; sdv <- sdv[!const]
    raw_median <- raw_median[!const]
  }
  if (!ncol(Tx)) stop("No usable biomarkers remain.")
  M <- sweep(sweep(Tx, 2, mu, "-"), 2, sdv, "/")
  n_ref <- sum(y == 0); n_case <- sum(y == 1)
  if (min(n_ref, n_case) < s$folds)
    stop(sprintf("Each group needs at least %d participants for %d-fold cross-validation (now %d and %d).",
                 s$folds, s$folds, n_ref, n_case))
  list(y = y, ids = as.character(data[[id_col]]), M = M, markers = colnames(Tx), mu = unname(mu), sd = unname(sdv),
       raw_median = unname(raw_median), Z = Z, cov_info = cov_info, cov_levels = cov_levels,
       cov_numeric = cov_numeric, covariates = covariates, n_ref = n_ref, n_case = n_case,
       group_desc = g$desc, groups = groups, outcome_cols = outcome_cols, notes = notes, s = s)
}


# ==========================================================================
# 5) LOGISTIC REGRESSION AND CROSS-VALIDATION
# ==========================================================================

# Logistic regression by iteratively reweighted least squares; X includes the
# intercept column. Same convergence rule as glm (relative deviance change < tol).
irls_fit <- function(X, y, start = NULL, maxit = 25L, tol = 1e-8) {
  p <- ncol(X)
  b <- if (is.null(start) || length(start) != p) numeric(p) else start
  dev_old <- Inf; conv <- FALSE
  for (it in seq_len(maxit)) {
    eta <- drop(X %*% b)
    mu  <- 1 / (1 + exp(-pmin(pmax(eta, -30), 30)))
    dev <- -2 * sum(y * log(mu) + (1 - y) * log1p(-mu))
    if (!is.finite(dev)) return(NULL)
    if (abs(dev - dev_old) < tol * (abs(dev) + 0.1)) { conv <- TRUE; break }
    dev_old <- dev
    w  <- pmax(mu * (1 - mu), 1e-10)
    R  <- tryCatch(chol(crossprod(X, X * w)), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    b  <- drop(backsolve(R, backsolve(R, crossprod(X, w * eta + (y - mu)), transpose = TRUE)))
  }
  list(coef = b, conv = conv)
}

# Stratified folds, repeated R times. Repeat r uses its own random split.
make_folds <- function(y, k, R, seed = 123) {
  set.seed(seed)
  n <- length(y); tests <- list(); rep_id <- integer(0)
  for (r in seq_len(R)) {
    fold <- integer(n)
    for (cl in c(0L, 1L)) {
      idx <- which(y == cl); idx <- idx[sample.int(length(idx))]
      fold[idx] <- (seq_along(idx) + sample.int(k, 1)) %% k + 1L
    }
    for (f in seq_len(k)) { tests[[length(tests) + 1]] <- which(fold == f); rep_id <- c(rep_id, r) }
  }
  list(tests = tests, trains = lapply(tests, function(t) setdiff(seq_len(n), t)), rep = rep_id, k = k, R = R,
       pooled_y = y[unlist(tests)], rep_pos = rep(rep_id, lengths(tests)))
}

# Exact cross-validation: one fit per training set, each warm-started from the previous one
cv_exact <- function(X, y, folds) {
  out <- numeric(length(folds$pooled_y)); pos <- 0L; b <- NULL
  for (f in seq_along(folds$tests)) {
    tr  <- folds$trains[[f]]; te <- folds$tests[[f]]
    fit <- irls_fit(X[tr, , drop = FALSE], y[tr], start = b)
    if (is.null(fit) && !is.null(b)) fit <- irls_fit(X[tr, , drop = FALSE], y[tr])
    if (is.null(fit)) return(NULL)
    b <- if (fit$conv) fit$coef else NULL
    out[pos + seq_along(te)] <- 1 / (1 + exp(-pmin(pmax(drop(X[te, , drop = FALSE] %*% fit$coef), -30), 30)))
    pos <- pos + length(te)
  }
  out
}

# One-step approximation: one fit on all participants, then each training
# set's coefficients from a single Newton step that removes the held-out fold S,
#   b(-S) = b - (X'WX - X_S'W_S X_S)^-1 X_S'(y_S - mu_S)
# Returns NULL near separation (no convergence, |linear predictor| > eta_max,
# or a training information matrix that cannot be factorised); the caller then
# refits every fold exactly.
cv_onestep <- function(X, y, folds, eta_max = 10) {
  full <- irls_fit(X, y)
  if (is.null(full) || !full$conv) return(NULL)
  b <- full$coef; eta <- drop(X %*% b)
  if (max(abs(eta)) > eta_max) return(NULL)
  mu <- 1 / (1 + exp(-eta)); w <- mu * (1 - mu); r <- y - mu
  H  <- crossprod(X, X * w)
  out <- numeric(length(folds$pooled_y)); pos <- 0L
  for (S in folds$tests) {
    XS <- X[S, , drop = FALSE]
    R  <- tryCatch(chol(H - crossprod(XS, XS * w[S])), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    bf <- b - drop(backsolve(R, backsolve(R, crossprod(XS, r[S]), transpose = TRUE)))
    out[pos + seq_along(S)] <- 1 / (1 + exp(-drop(XS %*% bf)))
    pos <- pos + length(S)
  }
  out
}

# Pooled ROC metrics. AUC is the trapezoid under the tie-collapsed ROC curve
# (pROC's value). The threshold maximises Youden's J (sensitivity + specificity
# - 1); ties go to the threshold closest to 0.5. A participant is called a case
# when the predicted probability is >= the threshold.
roc_stats <- function(p, y) {
  np <- sum(y == 1); nn <- sum(y == 0)
  o  <- order(p, decreasing = TRUE, method = "radix")
  ps <- p[o]; ys <- y[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  m  <- length(ps)
  keep <- c(ps[-1L] != ps[-m], TRUE)
  tp <- tp[keep]; fp <- fp[keep]; s <- ps[keep]; nk <- length(tp)
  auc <- sum((c(0, tp[-nk]) + tp) * (fp - c(0, fp[-nk]))) / (2 * np * nn)
  tp <- c(0, tp); fp <- c(0, fp)
  thr  <- c(Inf, if (nk > 1) (s[-nk] + s[-1L]) / 2, -Inf)
  sens <- tp / np; spec <- 1 - fp / nn
  J    <- sens + spec
  tied <- which(J >= max(J) - 0.5 / (np * nn))
  j    <- tied[order(abs(thr[tied] - 0.5), -thr[tied])[1]]
  TP <- tp[j]; FP <- fp[j]; FN <- np - TP; TN <- nn - FP
  c(acc = (TP + TN) / (np + nn), sens = sens[j], spec = spec[j],
    ppv = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
    npv = if (TN + FN > 0) TN / (TN + FN) else NA_real_,
    thr = thr[j], auc = auc)
}

roc_curve_points <- function(p, y) {
  o <- order(p, decreasing = TRUE); ys <- y[o]; ps <- p[o]
  keep <- c(ps[-1L] != ps[-length(ps)], TRUE)
  tp <- cumsum(ys)[keep]; fp <- cumsum(1 - ys)[keep]
  data.frame(fpr = c(0, fp / sum(y == 0)), tpr = c(0, tp / sum(y == 1)))
}

# 2.5th and 97.5th percentiles of the repeat-specific balanced accuracies at a fixed threshold
ba_percentiles <- function(p, folds, thr) {
  v <- vapply(split(seq_along(p), folds$rep_pos), function(ix) {
    pp <- p[ix]; yy <- folds$pooled_y[ix]
    (mean(pp[yy == 1] >= thr) + mean(pp[yy == 0] < thr)) / 2
  }, numeric(1))
  if (length(v) >= 2) stats::quantile(v, c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
}


# ==========================================================================
# 6) PANEL SEARCH (runs in short steps so the page stays responsive)
# ==========================================================================
# Level 1 scores every single marker. Each later level is exhaustive while the
# number of combinations is within the limit, otherwise a beam search that
# extends each of the best `beam` panels of the previous level by one marker.
# Panels are ranked by balanced accuracy (ties: AUC) from rapid scoring on the
# first `search_repeats` repetitions; the best `beam` panels per size are then
# refitted exactly on all `final_repeats` repetitions of the same folds, and all
# reported metrics come from those exact refits.

HARD_EXHAUSTIVE_MAX <- 250000

new_search <- function(prep, s) {
  job <- new.env()
  job$prep <- prep; job$s <- s
  job$FULL <- cbind(1, prep$Z, prep$M)
  job$nb   <- 1L + ncol(prep$Z)
  job$p    <- ncol(prep$M)
  job$folds_final  <- make_folds(prep$y, s$folds, s$final_repeats, s$seed)
  n_s <- s$folds * s$search_repeats
  ff  <- job$folds_final
  job$folds_search <- list(tests = ff$tests[seq_len(n_s)], trains = ff$trains[seq_len(n_s)],
                           rep = ff$rep[seq_len(n_s)], pooled_y = prep$y[unlist(ff$tests[seq_len(n_s)])])
  job$k_max <- min(s$max_size, job$p)
  job$retained <- list(); job$level_info <- list(); job$log <- character(0)
  job$phase <- "search"; job$t0 <- Sys.time(); job$n_scored <- 0
  start_level(job, 1L)
  job
}

level_candidates <- function(job, k) {
  p <- job$p
  if (k == 1L) return(list(cand = matrix(seq_len(p), nrow = 1), mode = "exhaustive"))
  n_all <- choose(p, k)
  if (k == 2L) {
    if (n_all > HARD_EXHAUSTIVE_MAX) stop("Exhaustive two-marker search exceeds the browser limit; use the desktop analysis script or fewer candidate markers.")
    return(list(cand=utils::combn(p, k), mode="exhaustive"))
  }
  want_exh <- job$s$search_mode == "exhaustive" ||
              (job$s$search_mode == "auto" && n_all <= job$s$exhaustive_limit)
  if (want_exh && n_all <= HARD_EXHAUSTIVE_MAX) return(list(cand = utils::combn(p, k), mode = "exhaustive"))
  parents <- job$parents
  pieces <- lapply(seq_len(ncol(parents)), function(i) {
    par <- parents[, i]; js <- setdiff(seq_len(p), par)
    if (!length(js)) return(NULL)
    M <- rbind(matrix(par, nrow = k - 1L, ncol = length(js)), js)
    matrix(apply(M, 2, sort.int), nrow = k)
  })
  cand <- do.call(cbind, pieces)
  cand <- cand[, !duplicated(apply(cand, 2, paste, collapse = ",")), drop = FALSE]
  list(cand = cand, mode = if (want_exh) "beam (exhaustive too large)" else "beam")
}

start_level <- function(job, k) {
  lc <- level_candidates(job, k)
  job$k <- k; job$cand <- lc$cand; job$mode <- lc$mode
  job$met <- matrix(NA_real_, ncol(lc$cand), 9,
                    dimnames = list(NULL, c("acc", "sens", "spec", "ppv", "npv", "thr", "auc", "ba", "exact")))
  job$pos <- 1L; job$t_level <- Sys.time()
  job$log <- c(job$log, sprintf("[%s] Panel size %d: %s, %s candidate panels.", format(Sys.time(), "%H:%M:%S"),
                                k, lc$mode, fmt_int(ncol(lc$cand))))
}

score_rapid <- function(job, combo) {
  X <- job$FULL[, c(seq_len(job$nb), job$nb + combo), drop = FALSE]
  y <- job$prep$y
  pred <- if (isTRUE(job$s$onestep)) cv_onestep(X, y, job$folds_search, job$s$eta_max) else NULL
  exact <- is.null(pred)
  if (exact) pred <- cv_exact(X, y, job$folds_search)
  if (is.null(pred) || anyNA(pred)) return(c(rep(NA_real_, 8), as.numeric(exact)))
  st <- roc_stats(pred, job$folds_search$pooled_y)
  c(st, (st[["sens"]] + st[["spec"]]) / 2, as.numeric(exact))
}

finish_level <- function(job) {
  k  <- job$k
  ok <- which(!is.na(job$met[, "ba"]))
  ord  <- ok[order(-job$met[ok, "ba"], -job$met[ok, "auc"])]
  keep <- head(ord, job$s$beam)
  job$retained[[k]] <- list(cand = job$cand[, keep, drop = FALSE], rapid = job$met[keep, , drop = FALSE])
  job$parents <- job$cand[, keep, drop = FALSE]
  n <- ncol(job$cand)
  job$level_info[[k]] <- list(mode = job$mode, n = n, exact_share = mean(job$met[, "exact"] == 1, na.rm = TRUE),
                              best = if (length(ord)) job$met[ord[1], "ba"] else NA_real_)
  job$log <- c(job$log, sprintf("[%s] Panel size %d done in %s: best rapid balanced accuracy %s%s.",
                                format(Sys.time(), "%H:%M:%S"), k,
                                fmt_time(as.numeric(difftime(Sys.time(), job$t_level, units = "secs"))),
                                fmt_pct(job$level_info[[k]]$best),
                                if (isTRUE(job$s$onestep)) sprintf("; %s of panels refitted fold by fold (near separation)",
                                                                   fmt_pct(job$level_info[[k]]$exact_share)) else ""))
  if (k < job$k_max) {
    start_level(job, k + 1L)
  } else {
    job$phase <- "final"
    job$fq <- do.call(rbind, lapply(seq_along(job$retained), function(kk)
      cbind(kk, seq_len(ncol(job$retained[[kk]]$cand)))))
    job$fpos  <- 1L
    job$final <- lapply(job$retained, function(r) matrix(NA_real_, ncol(r$cand), 10,
      dimnames = list(NULL, c("acc", "sens", "spec", "ppv", "npv", "thr", "auc", "ba", "ba_lo", "ba_hi"))))
    job$final_pred <- lapply(job$retained, function(r) vector("list", ncol(r$cand)))
    job$t_level <- Sys.time()
    job$log <- c(job$log, sprintf("[%s] Exact refits of the %s retained panels on %d × %d cross-validation...",
                                  format(Sys.time(), "%H:%M:%S"), fmt_int(nrow(job$fq)), job$s$folds,
                                  job$s$final_repeats))
  }
}

search_tick <- function(job, budget = 0.3) {
  t_end <- Sys.time() + budget
  if (job$phase == "search") {
    while (job$pos <= ncol(job$cand) && Sys.time() < t_end) {
      job$met[job$pos, ] <- score_rapid(job, job$cand[, job$pos])
      job$pos <- job$pos + 1L; job$n_scored <- job$n_scored + 1
    }
    done <- job$pos - 1L; tot <- ncol(job$cand)
    rate <- done / max(1e-6, as.numeric(difftime(Sys.time(), job$t_level, units = "secs")))
    job$msg <- sprintf("Panel size %d of %d (%s): %s of %s panels (%.0f%%), about %s left at this size.",
                       job$k, job$k_max, job$mode, fmt_int(done), fmt_int(tot), 100 * done / tot,
                       fmt_time((tot - done) / max(rate, 1e-6)))
    job$frac <- ((job$k - 1) + done / tot) / (job$k_max + 1)
    if (job$pos > ncol(job$cand)) finish_level(job)
  } else if (job$phase == "final") {
    y <- job$prep$y; ff <- job$folds_final
    while (job$fpos <= nrow(job$fq) && Sys.time() < t_end) {
      kk <- job$fq[job$fpos, 1]; j <- job$fq[job$fpos, 2]
      combo <- job$retained[[kk]]$cand[, j]
      X <- job$FULL[, c(seq_len(job$nb), job$nb + combo), drop = FALSE]
      pred <- cv_exact(X, y, ff)
      if (!is.null(pred) && !anyNA(pred)) {
        st <- roc_stats(pred, ff$pooled_y)
        job$final[[kk]][j, ] <- c(st, (st[["sens"]] + st[["spec"]]) / 2, ba_percentiles(pred, ff, st[["thr"]]))
        job$final_pred[[kk]][[j]] <- pred
      }
      job$fpos <- job$fpos + 1L
    }
    job$msg <- sprintf("Exact refits: %s of %s retained panels.", fmt_int(job$fpos - 1L), fmt_int(nrow(job$fq)))
    job$frac <- (job$k_max + (job$fpos - 1) / nrow(job$fq)) / (job$k_max + 1)
    if (job$fpos > nrow(job$fq)) assemble_results(job)
  }
  invisible(job$phase)
}

assemble_results <- function(job) {
  mk <- job$prep$markers
  tabs <- list(); best <- list()
  for (k in seq_along(job$retained)) {
    r <- job$retained[[k]]; f <- job$final[[k]]
    ok <- which(!is.na(f[, "ba"]))
    if (!length(ok)) next
    ord <- ok[order(-f[ok, "ba"], -f[ok, "auc"])]
    cand <- r$cand[, ord, drop = FALSE]
    bm <- t(apply(cand, 2, function(ix) mk[ix])); if (k == 1) bm <- matrix(mk[cand], ncol = 1)
    colnames(bm) <- paste0("Biomarker", seq_len(k))
    tabs[[k]] <- data.frame(
      Rank = seq_along(ord),
      Panel = apply(bm, 1, paste, collapse = " + "),
      bm,
      Balanced_Accuracy = f[ord, "ba"], BA_2.5th_pct = f[ord, "ba_lo"], BA_97.5th_pct = f[ord, "ba_hi"],
      Sensitivity = f[ord, "sens"], Specificity = f[ord, "spec"], Accuracy = f[ord, "acc"],
      PPV = f[ord, "ppv"], NPV = f[ord, "npv"], AUC = f[ord, "auc"], Threshold = f[ord, "thr"],
      Rapid_Balanced_Accuracy = r$rapid[ord, "ba"], Rapid_Refitted_Per_Fold = r$rapid[ord, "exact"],
      check.names = FALSE, stringsAsFactors = FALSE)
    attr(tabs[[k]], "idx") <- cand
    best[[k]] <- list(pred = job$final_pred[[k]][[ord[1]]], panel = tabs[[k]]$Panel[1])
  }
  job$results <- list(tables = tabs, best = best, folds = job$folds_final, prep = job$prep, s = job$s,
                      level_info = job$level_info, n_scored = job$n_scored,
                      elapsed = as.numeric(difftime(Sys.time(), job$t0, units = "secs")))
  job$phase <- "done"; job$frac <- 1
  job$log <- c(job$log, sprintf("[%s] Finished: %s panels scored in %s.", format(Sys.time(), "%H:%M:%S"),
                                fmt_int(job$n_scored), fmt_time(job$results$elapsed)))
  job$msg <- "Finished."
}

# Run-time projection from a timed sample at the largest panel size
project_runtime <- function(job) {
  p <- job$p; K <- job$k_max; nb <- job$s$beam
  counts <- vapply(seq_len(K), function(k) {
    if (k <= 2) return(choose(p, k))
    n_all <- choose(p, k)
    exh <- job$s$search_mode == "exhaustive" || (job$s$search_mode == "auto" && n_all <= job$s$exhaustive_limit)
    if (exh && n_all <= HARD_EXHAUSTIVE_MAX) n_all else min(n_all, min(nb, choose(p, k - 1)) * (p - k + 1))
  }, numeric(1))
  set.seed(1)
  smp <- replicate(12, sort(sample.int(p, K)))
  if (is.null(dim(smp))) smp <- matrix(smp, nrow = K)
  timed <- function(onestep) {
    old <- job$s$onestep; job$s$onestep <- onestep; on.exit(job$s$onestep <- old)
    system.time(for (j in seq_len(ncol(smp))) score_rapid(job, smp[, j]))[["elapsed"]] / ncol(smp)
  }
  per_exact <- timed(FALSE)
  per_fast  <- if (isTRUE(job$s$onestep)) timed(TRUE) else per_exact
  # one-step scoring is skipped for panels near separation, so the run time
  # lies between the two projections
  list(counts = counts, seconds = sum(counts) * min(per_fast, per_exact), seconds_max = sum(counts) * per_exact)
}

# Panel names one per line (for tables rendered with HTML allowed)
panel_html <- function(t, rows) {
  k <- sum(grepl("^Biomarker[0-9]+$", names(t)))
  vapply(rows, function(i) paste(htmltools::htmlEscape(unlist(t[i, paste0("Biomarker", seq_len(k))])), collapse = "<br>"),
         character(1))
}

summary_table <- function(res) {
  rows <- lapply(seq_along(res$tables), function(k) {
    t <- res$tables[[k]]; if (is.null(t)) return(NULL)
    data.frame(`Panel size` = k, `Best panel` = panel_html(t, 1), `Balanced accuracy` = fmt_pct(t$Balanced_Accuracy[1]),
               Sensitivity = fmt_pct(t$Sensitivity[1]), Specificity = fmt_pct(t$Specificity[1]),
               AUC = sprintf("%.3f", t$AUC[1]), check.names = FALSE)
  })
  do.call(rbind, rows)
}


# ==========================================================================
# 7) TRAINING A MODEL FROM A PANEL, AND PREDICTION
# ==========================================================================

fit_panel_model <- function(res, k, rank) {
  prep <- res$prep; tab <- res$tables[[k]]
  idx  <- attr(tab, "idx")[, rank]
  Mz   <- prep$M[, idx, drop = FALSE]
  X    <- cbind(prep$Z, Mz)
  nm   <- c(prep$cov_info$term, prep$markers[idx]) %||% character(0)
  safe <- make.names(c(paste0("cov", seq_len(ncol(prep$Z))), paste0("bm", seq_along(idx))))
  df   <- as.data.frame(X); names(df) <- safe; df$.y <- prep$y
  warn <- character(0)
  fit  <- withCallingHandlers(stats::glm(.y ~ ., family = stats::binomial(), data = df),
                              warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  co <- summary(fit)$coefficients
  b  <- stats::coef(fit)
  se <- co[, "Std. Error"][names(b)]; pv <- co[, "Pr(>|z|)"][names(b)]
  b_bm  <- b[paste0("bm", seq_along(idx))]; mu <- prep$mu[idx]; sdv <- prep$sd[idx]
  b_raw <- b_bm / sdv
  b0    <- b[["(Intercept)"]] - sum(b_bm * mu / sdv)
  cov_df <- if (nrow(prep$cov_info)) {
    data.frame(var = prep$cov_info$var, term = prep$cov_info$term,
               abbr = ifelse(prep$cov_info$type == "factor", paste0(prep$cov_info$var, " = ", prep$cov_info$level),
                             prep$cov_info$var),
               label = pretty_name(prep$cov_info$var), unit = "", type = prep$cov_info$type,
               level = prep$cov_info$level, coef = unname(b[paste0("cov", seq_len(ncol(prep$Z)))]),
               stringsAsFactors = FALSE)
  } else data.frame(var = character(0), term = character(0), abbr = character(0), label = character(0),
                    unit = character(0), type = character(0), level = character(0), coef = numeric(0))
  sep <- any(grepl("fitted probabilities numerically 0 or 1|did not converge", warn)) || max(abs(b[-1])) > 15
  s <- prep$s; g <- prep$groups
  list(source = "trained",
       name = sprintf("Model trained on your data (%d-marker panel, rank %d)", k, rank),
       description = paste("Panel:", tab$Panel[rank]),
       case_label = g$case_name, ref_label = g$ref_name, case_short = g$case_name, ref_short = g$ref_name,
       transform = s$transform, delta = s$delta, intercept = unname(b0),
       markers = data.frame(key = prep$markers[idx], abbr = paste0("B", seq_along(idx)),
                            label = pretty_name(prep$markers[idx]), coef = unname(b_raw),
                            ref_mean = mu, ref_sd = sdv, raw_median = prep$raw_median[idx],
                            stringsAsFactors = FALSE),
       covariates = cov_df, cov_levels = prep$cov_levels, cov_numeric = prep$cov_numeric,
       threshold = tab$Threshold[rank],
       cohort = sprintf("Trained on %d participants (%d %s, %d %s): %s.", length(prep$y), prep$n_ref, g$ref_name,
                        prep$n_case, g$case_name, prep$group_desc),
       performance = c("Balanced accuracy" = tab$Balanced_Accuracy[rank], "Sensitivity" = tab$Sensitivity[rank],
                       "Specificity" = tab$Specificity[rank], "AUC" = tab$AUC[rank]),
       performance_note = sprintf("%d-fold cross-validation repeated %d times in the uploaded data; selection and evaluation used the same data.",
                                  s$folds, s$final_repeats),
       caution = if (sep) "Possible separation: at least one coefficient is very large or the fit did not converge, so estimates for this model are unstable." else NULL,
       fit_table = data.frame(Term = c("(Intercept)", nm),
                              Estimate_model_scale = unname(b), SE = unname(se), p_value = unname(pv),
                              Coefficient_raw_input = c(b0, cov_df$coef, b_raw),
                              stringsAsFactors = FALSE),
       warnings = unique(warn))
}

# Coefficient table that can be loaded back into Patient prediction
model_to_table <- function(m) {
  tr <- if (m$transform == "log2") sprintf("log2(raw value + %s)", format(m$delta)) else "raw value"
  rbind(
    if (!is.null(m$analysis_configuration)) data.frame(
      Term = c("(Model configuration)", "(Search design)"),
      Label = c(m$analysis_configuration, m$search_design), Coefficient_raw_input = NA_real_, Input = ""),
    data.frame(Term = "(Intercept)", Label = "Intercept", Coefficient_raw_input = m$intercept, Input = "", stringsAsFactors = FALSE),
    if (nrow(m$covariates)) data.frame(Term = m$covariates$term, Label = m$covariates$label, Coefficient_raw_input = m$covariates$coef,
                                       Input = ifelse(m$covariates$type == "factor",
                                                      sprintf("indicator: %s = %s", m$covariates$var, m$covariates$level),
                                                      "value in original units"), stringsAsFactors = FALSE),
    data.frame(Term = m$markers$key, Label = m$markers$label, Coefficient_raw_input = m$markers$coef, Input = tr,
               stringsAsFactors = FALSE),
    if (length(m$cov_levels)) data.frame(Term = paste0("(Reference level: ", names(m$cov_levels), ")"),
      Label = vapply(m$cov_levels, function(x) x[1], character(1)), Coefficient_raw_input = NA_real_, Input = ""),
    if (length(m$level_labels)) do.call(rbind, lapply(names(m$level_labels), function(v)
      data.frame(Term = sprintf("(Level label: %s = %s)", v, names(m$level_labels[[v]])),
                 Label = unname(m$level_labels[[v]]), Coefficient_raw_input = NA_real_, Input = "",
                 stringsAsFactors = FALSE))),
    data.frame(Term = c("(Threshold)", "(Delta)", "(Case group)", "(Reference group)"),
               Label = c("Classification threshold", "log2 offset", m$case_label, m$ref_label),
               Coefficient_raw_input = c(m$threshold, if (m$transform == "log2") m$delta else NA, NA, NA),
               Input = "", stringsAsFactors = FALSE))
}

model_from_table <- function(tab, name) {
  need <- c("Term", "Coefficient_raw_input", "Input")
  if (!all(need %in% names(tab))) stop("A model file needs the columns Term, Coefficient_raw_input and Input.")
  if (anyDuplicated(trimws(as.character(tab$Term)))) stop("Duplicate model terms are not allowed.")
  tab$Term <- trimws(as.character(tab$Term)); tab$Input <- trimws(as.character(tab$Input %||% ""))
  val <- function(t) { v <- suppressWarnings(as.numeric(tab$Coefficient_raw_input[tab$Term == t])); if (length(v)) v[1] else NA_real_ }
  lab <- function(t, d) { v <- tab$Label[tab$Term == t]; if (length(v) && !is.na(v[1]) && nzchar(v[1])) v[1] else d }
  tab$Input[is.na(tab$Input)] <- ""
  is_bm  <- grepl("^log2\\(|^raw value", tab$Input)
  is_num <- tab$Input == "value in original units"
  is_ind <- grepl("^indicator:", tab$Input)
  log2_in <- any(grepl("^log2\\(", tab$Input[is_bm]))
  delta <- val("(Delta)")
  if (is.na(delta) && log2_in) {
    d <- suppressWarnings(as.numeric(sub("^log2\\(raw value \\+ ([^)]+)\\)$", "\\1", tab$Input[is_bm][1])))
    delta <- if (is.finite(d)) d else 5e-7
  }
  ind <- tab[is_ind, , drop = FALSE]
  ind_var <- sub("^indicator: (.*) = (.*)$", "\\1", ind$Input); ind_lev <- sub("^indicator: (.*) = (.*)$", "\\2", ind$Input)
  covs <- rbind(
    if (any(is_num)) data.frame(var = tab$Term[is_num], term = tab$Term[is_num], abbr = tab$Term[is_num],
                                label = if ("Label" %in% names(tab)) tab$Label[is_num] else tab$Term[is_num],
                                unit = "", type = "numeric", level = NA_character_,
                                coef = as.numeric(tab$Coefficient_raw_input[is_num]), stringsAsFactors = FALSE),
    if (nrow(ind)) data.frame(var = ind_var, term = ind$Term, abbr = paste0(ind_var, " = ", ind_lev), label = ind_var,
                              unit = "", type = "factor", level = ind_lev,
                              coef = as.numeric(ind$Coefficient_raw_input), stringsAsFactors = FALSE))
  if (is.null(covs)) covs <- data.frame(var = character(0), term = character(0), abbr = character(0), label = character(0),
                                        unit = character(0), type = character(0), level = character(0), coef = numeric(0))
  cov_levels <- lapply(split(covs$level[covs$type == "factor"], covs$var[covs$type == "factor"]),
                       function(l) c("(reference level)", l))
  for (v in names(cov_levels)) {
    ref <- lab(paste0("(Reference level: ", v, ")"), "(reference level)")
    cov_levels[[v]][1] <- ref
  }
  # Optional "(Level label: <variable> = <level>)" rows give the meaning of each code
  ll <- tab[grepl("^\\(Level label: .* = .*\\)$", tab$Term), , drop = FALSE]
  level_labels <- list()
  if (nrow(ll) && "Label" %in% names(ll)) {
    lv_var <- sub("^\\(Level label: (.*) = (.*)\\)$", "\\1", ll$Term)
    lv_lev <- sub("^\\(Level label: (.*) = (.*)\\)$", "\\2", ll$Term)
    for (v in unique(lv_var)) level_labels[[v]] <- stats::setNames(as.character(ll$Label[lv_var == v]), lv_lev[lv_var == v])
    for (i in which(covs$type == "factor" & covs$var %in% names(level_labels))) {
      labs <- level_labels[[covs$var[i]]]
      if (!is.na(labs[covs$level[i]])) {
        others <- setdiff(names(labs), covs$level[i])
        covs$abbr[i]  <- unname(labs[covs$level[i]])
        covs$label[i] <- sprintf("%s %s (vs %s)", unname(labs[covs$level[i]]), tolower(pretty_name(covs$var[i])),
                                 and_list(tolower(unname(labs[others]))))
      }
    }
  }
  bm <- tab[is_bm, , drop = FALSE]
  if (!nrow(bm)) stop("No biomarker rows found (Input should start with \"log2(\" or \"raw value\").")
  thr <- val("(Threshold)")
  list(source = "loaded", name = paste("Model loaded from", name), description = "",
       case_label = lab("(Case group)", "case group"), ref_label = lab("(Reference group)", "reference group"),
       case_short = lab("(Case group)", "case"), ref_short = lab("(Reference group)", "reference"),
       transform = if (log2_in) "log2" else "none", delta = if (log2_in) delta else 0,
       intercept = val("(Intercept)"),
       markers = data.frame(key = bm$Term, abbr = paste0("B", seq_len(nrow(bm))),
                            label = if ("Label" %in% names(bm)) as.character(bm$Label) else bm$Term,
                            coef = as.numeric(bm$Coefficient_raw_input), stringsAsFactors = FALSE),
       covariates = covs, cov_levels = cov_levels, level_labels = if (length(level_labels)) level_labels else NULL,
       threshold = if (is.finite(thr)) thr else 0.5,
       cohort = "Loaded from a saved model file.", performance = NULL,
       performance_note = if (is.finite(thr)) "" else "No threshold in the file; 0.5 is used.", caution = NULL)
}

# Linear predictor and probability, with each term of the sum
predict_model <- function(m, marker_vals, cov_vals) {
  if (!model_ready(m)) stop("Incomplete or non-finite model coefficients.")
  if (length(marker_vals) != nrow(m$markers) || any(!is.finite(marker_vals)) || any(marker_vals < 0))
    stop("Biomarker inputs must be finite, nonnegative values for every marker.")
  for (v in unique(m$covariates$var)) {
    vv <- cov_vals[[v]]
    if (length(vv) != 1L || is.na(vv)) stop("Missing participant characteristic: ", v)
    if (m$covariates$type[m$covariates$var == v][1] == "factor") {
      if (!(as.character(vv) %in% m$cov_levels[[v]])) stop("Unknown category for ", v)
    } else if (!is.finite(as.numeric(vv))) stop("Non-finite participant characteristic: ", v)
  }
  tx <- if (m$transform == "log2") log2(marker_vals + m$delta) else marker_vals
  terms <- data.frame(Term = m$markers$label, Abbr = m$markers$abbr, Input = marker_vals, Transformed = tx,
                      Coefficient = m$markers$coef, Contribution = m$markers$coef * tx, kind = "marker",
                      Shown = NA_character_, stringsAsFactors = FALSE)
  if (nrow(m$covariates)) {
    cv <- m$covariates
    xin <- vapply(seq_len(nrow(cv)), function(i) {
      v <- cov_vals[[cv$var[i]]]
      if (cv$type[i] == "factor") as.numeric(identical(as.character(v), cv$level[i])) else as.numeric(v)
    }, numeric(1))
    # indicator values are shown with the category entered, e.g. "0 (female)"
    shown <- vapply(seq_len(nrow(cv)), function(i) {
      if (cv$type[i] != "factor") return(NA_character_)
      sprintf("%d (%s)", as.integer(xin[i]), tolower(level_label(m, cv$var[i], cov_vals[[cv$var[i]]])))
    }, character(1))
    terms <- rbind(terms, data.frame(Term = cv$label, Abbr = cv$abbr, Input = xin, Transformed = xin,
                                     Coefficient = cv$coef, Contribution = cv$coef * xin, kind = cv$type,
                                     Shown = shown, stringsAsFactors = FALSE))
  }
  eta <- m$intercept + sum(terms$Contribution)
  list(terms = terms, eta = eta, p = 1 / (1 + exp(-eta)))
}


# ==========================================================================
# 8) DISPLAYS: EQUATION, GAUGE, ROC
# ==========================================================================


is_age_variable <- function(v) tolower(v) %in% c("age", "age_dos")

age_slider_limits <- function(m, v) {
  middle <- m$cov_numeric[[v]] %||% 40
  c(if (identical(m$source, "published")) 18 else 0, max(100, ceiling(middle)))
}

age_bmi_variables <- function(m) {
  if (is.null(m) || !nrow(m$covariates)) return(NULL)
  age <- unique(m$covariates$var[is_age_variable(m$covariates$var) & m$covariates$type == "numeric"])
  bmi <- unique(m$covariates$var[tolower(m$covariates$var) == "bmi" & m$covariates$type == "numeric"])
  if (length(age) != 1 || length(bmi) != 1) return(NULL)
  list(age = age, bmi = bmi)
}

age_bmi_curve_data <- function(m, marker_vals, cov_vals, ages) {
  vars <- age_bmi_variables(m)
  if (is.null(vars)) stop("Comparison requires numeric age and numeric BMI in the model.")
  # Validate all current inputs before generating any scenarios.
  predict_model(m, marker_vals, cov_vals)
  do.call(rbind, lapply(c(as.numeric(cov_vals[[vars$bmi]]), as.numeric(cov_vals[[vars$bmi]]) + 5), function(level) {
    risk <- vapply(ages, function(age) {
      scenario <- cov_vals
      scenario[[vars$age]] <- age
      scenario[[vars$bmi]] <- level
      predict_model(m, marker_vals, scenario)$p
    }, numeric(1))
    data.frame(age = ages, bmi = level, risk = risk, stringsAsFactors = FALSE)
  }))
}

draw_age_bmi_curves <- function(curves, selected_age, selected_bmi, threshold, labels = NULL) {
  levels <- unique(curves$bmi)
  if (is.null(labels)) labels <- paste("BMI", levels)
  colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), length.out = length(levels))
  old <- par(mar = c(4.5, 4.5, 2, 1)); on.exit(par(old))
  plot(range(curves$age), c(0, 100), type = "n", xlab = "Age (years)",
       ylab = "Estimated risk (%)", xaxs = "i", yaxs = "i")
  abline(h = seq(0, 100, 20), col = "#E6E6E6", lwd = 0.8)
  abline(h = 100 * threshold, col = "#666666", lty = 3)
  abline(v = selected_age, col = "#555555", lty = 3)
  for (i in seq_along(levels)) {
    d <- curves[curves$bmi == levels[i], ]
    lines(d$age, 100 * d$risk, col = colors[i], lwd = 2.5, lty = i)
    p <- d$risk[which.min(abs(d$age - selected_age))]
    points(selected_age, 100 * p, col = colors[i], pch = if (as.character(levels[i]) == as.character(selected_bmi)) 19 else 1, cex = 1.4, lwd = 2)
  }
  legend("topleft", legend = labels, col = colors, lty = seq_along(levels),
         lwd = 2.5, bg = "white", cex = 0.9)
  mtext(sprintf("Selected age: %s years | Classification threshold: %.1f%%", selected_age, 100 * threshold),
        side = 3, line = 0.3, cex = 0.85)
}


equation_html <- function(m) {
  co <- function(x, sym) if (is.finite(x)) fmt_num(abs(x)) else sym
  sg <- function(x) if (is.finite(x) && x < 0) paste0(MINUS, " ") else "+ "
  inp <- function(a) if (m$transform == "log2") sprintf("log<sub>2</sub>(%s + &delta;)", a) else a
  b0 <- if (is.finite(m$intercept)) fmt_num(m$intercept) else "&beta;<sub>0</sub>"
  terms <- character(0)
  for (j in seq_len(nrow(m$markers))) {
    terms <- c(terms, sprintf("%s%s&thinsp;%s", sg(m$markers$coef[j]), co(m$markers$coef[j], sprintf("b<sub>%d</sub>", j)),
                              inp(m$markers$abbr[j])))
  }
  if (nrow(m$covariates)) for (i in seq_len(nrow(m$covariates))) {
    terms <- c(terms, sprintf("%s%s&thinsp;%s", sg(m$covariates$coef[i]),
                              co(m$covariates$coef[i], sprintf("&gamma;<sub>%d</sub>", i)),
                              htmltools::htmlEscape(m$covariates$abbr[i])))
  }
  # each term stays on one line; lines break only between terms, and
  # continuation lines are indented (hanging indent), whatever the window width
  lead <- sprintf("logit&thinsp;<i>p</i>(%s) = %s", htmltools::htmlEscape(m$case_short), b0)
  body <- paste(sprintf("<span class='eq-term'>%s</span>", c(lead, terms)), collapse = " ")
  prob <- sprintf("<i>p</i>(%s) = 1 / [1 + exp(%s&eta;)], where &eta; is the right-hand side above",
                  htmltools::htmlEscape(m$case_short), MINUS)
  legend <- paste(sprintf("<b>%s</b> = %s", htmltools::htmlEscape(m$markers$abbr), md_html(m$markers$label)), collapse = "; ")
  delta <- if (m$transform == "log2") sprintf("; &delta; = %s, in the same units as the raw values", fmt_num(m$delta)) else ""
  covtxt <- ""
  if (nrow(m$covariates)) {
    cv <- m$covariates; esc <- htmltools::htmlEscape
    parts <- character(0); plain <- character(0)
    for (i in seq_len(nrow(cv))) {
      if (cv$type[i] == "factor") {
        labs <- m$level_labels[[cv$var[i]]]
        if (!is.null(labs) && !is.na(labs[cv$level[i]])) {
          ref <- level_label(m, cv$var[i], m$cov_levels[[cv$var[i]]][1])
          parts <- c(parts, sprintf("<b>%s</b> is 1 for %s and 0 for %s (the reference category)", esc(cv$abbr[i]),
                                    esc(tolower(labs[[cv$level[i]]])), esc(tolower(ref))))
        } else {
          parts <- c(parts, sprintf("<b>%s</b> is 1 if %s is %s and 0 otherwise", esc(cv$abbr[i]), esc(pretty_name(cv$var[i])),
                                    esc(cv$level[i])))
        }
      } else if (nzchar(cv$unit[i])) {
        parts <- c(parts, sprintf("<b>%s</b> in %s", esc(cv$abbr[i]), esc(cv$unit[i])))
      } else if (identical(m$source, "published")) {
        parts <- c(parts, sprintf("<b>%s</b> is the %s", esc(cv$abbr[i]), esc(tolower(cv$label[i]))))
      } else plain <- c(plain, sprintf("<b>%s</b>", esc(cv$abbr[i])))
    }
    if (length(plain)) parts <- c(parts, sprintf("%s in the units of the data the model was trained on", and_list(plain)))
    covtxt <- paste0("; ", paste(parts, collapse = "; "))
  }
  HTML(sprintf("<div class='equation eq-wrap'>%s</div><div class='equation eq-prob'>%s</div><p class='eq-legend'>%s%s%s.</p>",
               body, prob, legend, delta, covtxt))
}

# Speedometer: diverging blue (low) to red (high), neutral at the threshold
draw_gauge <- function(p, thr, lo_lab = "Lower", hi_lab = "Higher") {
  op <- graphics::par(mar = c(0, 0, 0, 0), bg = "white"); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(xlim = c(-1.3, 1.3), ylim = c(-0.3, 1.36), asp = 1)
  ang <- function(f) pi * (1 - f)
  arm_lo <- grDevices::colorRampPalette(c("#184f95", "#3987e5", "#9ec5f4", "#e4e2dc"))
  arm_hi <- grDevices::colorRampPalette(c("#e4e2dc", "#f2b3a6", "#e34948", "#a32727"))
  n <- 240; fr <- seq(0, 1, length.out = n + 1); mid <- (fr[-1] + fr[-(n + 1)]) / 2
  t0 <- if (is.finite(thr)) thr else 0.5
  n_lo <- sum(mid < t0); cols <- c(arm_lo(max(n_lo, 1))[seq_len(n_lo)], arm_hi(max(n - n_lo, 1))[seq_len(n - n_lo)])
  r_out <- 1; r_in <- 0.70
  for (i in seq_len(n)) {
    th <- seq(ang(fr[i]), ang(fr[i + 1]), length.out = 4)
    graphics::polygon(c(r_out * cos(th), rev(r_in * cos(th))), c(r_out * sin(th), rev(r_in * sin(th))),
                      col = cols[i], border = cols[i], lwd = 0.6)
  }
  # ticks every 5%, labelled every 25% (a label next to the threshold is left out)
  for (f in seq(0, 1, by = 0.05)) {
    a <- ang(f); major <- abs(f * 4 - round(f * 4)) < 1e-9
    graphics::segments(1.02 * cos(a), 1.02 * sin(a), (if (major) 1.09 else 1.055) * cos(a),
                       (if (major) 1.09 else 1.055) * sin(a), col = "#52514e", lwd = if (major) 1.4 else 0.8)
  }
  for (f in c(0.25, 0.5, 0.75)) {
    if (is.finite(thr) && abs(f - thr) < 0.12) next
    a <- ang(f)
    graphics::text(1.19 * cos(a), 1.19 * sin(a), sprintf("%d%%", round(100 * f)), cex = 0.82, col = "#52514e")
  }
  graphics::text(-1.12, -0.06, "0%", cex = 0.82, col = "#52514e", adj = c(1, 0.5))
  graphics::text(1.12, -0.06, "100%", cex = 0.82, col = "#52514e", adj = c(0, 0.5))
  # threshold marker across the band, labelled outside it on its own side
  if (is.finite(thr)) {
    a <- ang(thr)
    graphics::segments(0.66 * cos(a), 0.66 * sin(a), 1.09 * cos(a), 1.09 * sin(a), col = "#0b0b0b", lwd = 2.4)
    hx <- if (cos(a) > 0.2) 0 else if (cos(a) < -0.2) 1 else 0.5
    graphics::text(1.13 * cos(a), 1.13 * sin(a) + 0.04, sprintf("threshold %s", fmt_pct(thr)), cex = 0.78,
                   col = "#0b0b0b", font = 2, adj = c(hx, 0))
  }
  graphics::text(-0.85, -0.17, lo_lab, cex = 0.85, col = "#52514e", font = 2)
  graphics::text(0.85, -0.17, hi_lab, cex = 0.85, col = "#52514e", font = 2)
  # needle
  if (is.finite(p)) {
    a <- ang(min(max(p, 0), 1)); L <- 0.9; w <- 0.045
    tip <- c(L * cos(a), L * sin(a)); bl <- c(w * cos(a + pi / 2), w * sin(a + pi / 2)); br <- -bl
    graphics::polygon(c(tip[1], bl[1], -0.12 * cos(a), br[1]), c(tip[2], bl[2], -0.12 * sin(a), br[2]),
                      col = "#1f2937", border = "#1f2937")
  }
  th <- seq(0, 2 * pi, length.out = 60)
  graphics::polygon(0.075 * cos(th), 0.075 * sin(th), col = "#1f2937", border = "white", lwd = 2)
}

ROC_COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
draw_roc <- function(res) {
  op <- graphics::par(mar = c(4.2, 4.4, 1, 1), las = 1); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "1 − specificity",
                 ylab = "Sensitivity", bty = "l", col.axis = "#52514e", col.lab = "#0b0b0b")
  graphics::abline(h = seq(0.2, 0.8, 0.2), v = seq(0.2, 0.8, 0.2), col = "#ebeae6", lwd = 0.8)
  graphics::abline(0, 1, col = "#b5b3ac", lty = 2)
  lab <- character(0); cols <- character(0); pchs <- integer(0)
  for (k in seq_along(res$best)) {
    b <- res$best[[k]]; t <- res$tables[[k]]
    if (is.null(b) || is.null(b$pred)) next
    pts <- roc_curve_points(b$pred, res$folds$pooled_y)
    col <- ROC_COLS[(k - 1) %% length(ROC_COLS) + 1]
    graphics::lines(pts$fpr, pts$tpr, col = col, lwd = 2)
    graphics::points(1 - t$Specificity[1], t$Sensitivity[1], pch = c(16, 17, 15, 18, 4, 8)[(k - 1) %% 6 + 1],
                     col = col, cex = 1.5)
    lab <- c(lab, sprintf("%d marker%s (AUC %.3f)", k, if (k > 1) "s" else "", t$AUC[1]))
    cols <- c(cols, col); pchs <- c(pchs, c(16, 17, 15, 18, 4, 8)[(k - 1) %% 6 + 1])
  }
  if (length(lab)) graphics::legend("bottomright", legend = lab, col = cols, pch = pchs, lwd = 2, bty = "n",
                                    cex = 0.9, inset = 0.02)
}


# ==========================================================================
# 9) SIMULATED EXAMPLE DATA (for trying the analysis; not study data)
# ==========================================================================

make_example_data <- function(seed = 2026) {
  set.seed(seed)
  n <- 120
  id <- sprintf("P%03d", seq_len(n))
  age <- round(stats::rnorm(n, 40, 13)); age <- pmin(pmax(age, 18), 85)
  high <- stats::rbinom(n, 1, stats::plogis(-2.4 + 0.06 * age))
  missing <- pmin(pmax(round(stats::rnorm(n, 3 + 2.5 * high, 2)), 0), 20)
  decayed <- stats::rpois(n, 0.8 + 2 * high)
  filled  <- stats::rpois(n, 3.5 + 6.5 * high)
  pinfo <- data.frame(StudyID = id, Age = age, Sex = sample(c("Female", "Male"), n, TRUE),
                      BMI = round(stats::rnorm(n, 27 + 2.5 * high, 6), 1), Teeth = 32 - missing,
                      Decayed = decayed, Missing = missing, Filled = filled, stringsAsFactors = FALSE)
  p <- 24; nm <- sprintf("Taxon_%02d", seq_len(p))
  eff <- c(1.2, -1.0, 0.9, 0.8, -0.7, rep(0, p - 5))
  lin <- 0.4 * (pinfo$Decayed + pinfo$Filled + pinfo$Missing - 12) / 4
  bio <- vapply(seq_len(p), function(j) {
    det <- stats::rbinom(n, 1, stats::plogis(1 + 0.8 * eff[j] * (lin > 0)))
    round(det * exp(stats::rnorm(n, 5 + eff[j] * lin, 1.3)), 2)
  }, numeric(n))
  colnames(bio) <- nm
  bio <- data.frame(StudyID = id, bio, check.names = FALSE)
  bio[sample.int(n, 6), sample(nm, 1)] <- NA            # a few missing values
  bio <- bio[sample.int(n), ]                           # different row order
  bio <- bio[-c(3, 7), ]                                # two participants without biomarker data
  pinfo <- pinfo[-11, ]                                 # one participant without clinical data
  list(pinfo = pinfo, bio = bio)
}


# ==========================================================================
# 10) USER INTERFACE
# ==========================================================================

CSS <- "
:root { --accent:#0b7a75; --accent-dark:#075e5a; --ink:#0b0b0b; --ink2:#52514e; --card:#f6f7f8; --line:#dde1e5; }
body { color: var(--ink); font-family: system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }
.navbar-default { background:#ffffff; border-bottom:1px solid var(--line); }
.navbar-default .navbar-brand { color: var(--accent-dark); font-weight:600; }
.navbar-default .navbar-nav > .active > a, .navbar-default .navbar-nav > .active > a:focus,
.navbar-default .navbar-nav > .active > a:hover { background:#eef6f5; color: var(--accent-dark); }
.card { background: var(--card); border:1px solid var(--line); border-radius:10px; padding:18px 20px; margin-bottom:16px; }
.card h3, .card h4 { margin-top:0; }
.btn-accent { background: var(--accent); color:#fff; border:none; font-size:16px; padding:9px 18px; border-radius:6px; }
.btn-accent:hover, .btn-accent:focus { background: var(--accent-dark); color:#fff; }
.btn-soft { background:#fff; border:1px solid var(--line); border-radius:6px; }
.lead-muted { color: var(--ink2); font-size:16px; }
.step-title { font-weight:600; margin:14px 0 6px; }
.help-note { color: var(--ink2); font-size:13px; margin-top:-4px; }
.raw-note { background:#fff8e6; border:1px solid #f1d58a; border-radius:8px; padding:10px 12px; margin-bottom:12px; }
.equation { font-family: 'Cambria Math', 'STIX Two Math', 'Times New Roman', serif; font-size:17px; line-height:1.9;
            background:#fff; border:1px solid var(--line); border-radius:8px; padding:10px 14px; overflow-x:auto; }
.eq-wrap { padding-left:calc(14px + 2.2em); text-indent:-2.2em; }
.eq-term { white-space:nowrap; }
.eq-prob { margin-top:8px; }
.eq-legend { color: var(--ink2); font-size:13px; margin-top:8px; }
.hero { font-size:56px; font-weight:600; line-height:1.05; color: var(--ink); }
.hero-sub { color: var(--ink2); font-size:15px; }
.status { font-size:16px; font-weight:600; margin-top:8px; }
.status.hi { color:#a32727; } .status.lo { color:#184f95; }
.calc-table { width:100%; border-collapse:collapse; background:#fff; font-size:14px; }
.calc-table th, .calc-table td { border-bottom:1px solid var(--line); padding:6px 8px; text-align:right; }
.calc-table th:first-child, .calc-table td:first-child { text-align:left; }
.calc-table tr.total td { font-weight:600; border-top:2px solid #b5b3ac; }
.progress-wrap { background:#e8ebee; border-radius:6px; height:12px; overflow:hidden; margin:8px 0; }
.progress-fill { background: var(--accent); height:12px; }
pre.config { background:#fff; border:1px solid var(--line); font-size:13px; white-space:pre-wrap; }
.notice { background:#fdecec; border:1px solid #f3b4b4; border-radius:8px; padding:10px 12px; margin-bottom:12px; }
.okmsg { background:#ecf7ef; border:1px solid #b9dfc3; border-radius:8px; padding:8px 12px; margin-bottom:12px; }
table.table { background:#fff; }
.card .shiny-input-container:not(.shiny-input-container-inline) { width:100% !important; max-width:100%; }
.res-wrap { overflow-x:auto; }
.res-wrap td { vertical-align:top; }
.footer-note { color: var(--ink2); font-size:12px; margin:24px 0 12px; }
#config_text, #run_log { white-space:pre-wrap; overflow-wrap:break-word; }
"

start_tab <- tabPanel(
  "Start", value = "start",
  fluidRow(column(12,
    h2("Salivary biomarker DMFT calculator"),
    p(class = "lead-muted", "Search for salivary biomarker panels that separate low and high DMFT in your own data,",
      "or estimate an individual's risk of high DMFT from absolute raw concentration values.")
  )),
  fluidRow(
    column(6, div(class = "card",
      h3("Analyze uploaded data"),
      p(class = "lead-muted", "Upload participant information and raw biomarker measurements, join them by participant ID,",
        "define the comparison groups, search for biomarker panels and train your own model."),
      actionButton("go_analysis", "Go to Train Your Own Data", class = "btn-accent"))),
    column(6, div(class = "card",
      h3("Calculate a patient's risk"),
      p(class = "lead-muted", "Load the revised model, then enter raw biomarker values with age and BMI to estimate the",
        "risk of high DMFT (DMFT > 13) with the selected model, or use a model you trained."),
      actionButton("go_predict", "Go to Risk Calculator", class = "btn-accent")))
  ),
  fluidRow(column(12, p(class = "footer-note",
    "Research use only. Estimates come from a cross-sectional development cohort and have not been externally validated;",
    "they are not a diagnosis and do not replace a dental examination.", sprintf("Calculator version %s.", APP_VERSION))))
)

analysis_tab <- tabPanel(
  "Train Your Own Data", value = "analysis",
  fluidRow(column(12, h3("Biomarker panel analysis"),
                  p(class = "lead-muted", "Upload participant information and raw biomarker measurements, then match records by participant ID."))),
  fluidRow(
    column(4, div(class = "card",
      div(style = "display:flex; gap:8px; flex-wrap:wrap;",
          actionButton("load_example", "Load example data", class = "btn-soft"),
          downloadButton("dl_example_info", "Example participant file", class = "btn-soft"),
          downloadButton("dl_example_bio", "Example biomarker file", class = "btn-soft")),
      p(class = "help-note", style = "margin-top:6px;", "The example data are simulated."),
      p(class = "help-note", "Browser version: CSV/TSV inputs and CSV exports work without extra packages. Excel files require readxl/writexl, which may not be available in WebR. Large exhaustive searches may be slow or run out of browser memory."),
      div(class = "step-title", "1. Participant information file"),
      fileInput("pinfo_file", NULL, accept = if (HAS_READXL) c(".csv", ".tsv", ".txt", ".xlsx", ".xls") else c(".csv", ".tsv", ".txt")),
      selectInput("pid", "Participant ID column", choices = NULL),
      div(class = "step-title", "2. Raw biomarker file"),
      checkboxInput("one_file", "Biomarker values are in the participant information file", FALSE),
      conditionalPanel("!input.one_file",
        fileInput("bio_file", NULL, accept = if (HAS_READXL) c(".csv", ".tsv", ".txt", ".xlsx", ".xls") else c(".csv", ".tsv", ".txt")),
        selectInput("bid", "Participant ID column", choices = NULL)),
      div(class = "step-title", "3. Comparison groups"),
      selectizeInput("outcome_cols", "Outcome column(s)", choices = NULL, multiple = TRUE),
      p(class = "help-note", "Several numeric columns are added together, e.g. decayed + missing + filled teeth = DMFT."),
      radioButtons("grp_mode", NULL, inline = FALSE,
                   choices = c("One cutoff: reference ≤ cutoff, case > cutoff" = "cutoff",
                               "Two cutoffs: leave out the values in between" = "two",
                               "Two categories" = "levels")),
      conditionalPanel("input.grp_mode == 'cutoff'",
        numericInput("cutoff", "Cutoff (default: the median)", value = NA),
        uiOutput("cutoff_hint")),
      conditionalPanel("input.grp_mode == 'two'",
        fluidRow(column(6, numericInput("ref_max", "Reference ≤", value = 0)),
                 column(6, numericInput("case_min", "Case ≥", value = 3)))),
      conditionalPanel("input.grp_mode == 'levels'",
        fluidRow(column(6, selectInput("ref_level", "Reference category", choices = NULL)),
                 column(6, selectInput("case_level", "Case category", choices = NULL)))),
      fluidRow(column(6, textInput("ref_name", "Reference label", "LD")),
               column(6, textInput("case_name", "Case label", "HD"))),
      div(class = "step-title", "4. Biomarkers and adjustment variables"),
      checkboxInput("all_markers", "Use every numeric column of the biomarker file", TRUE),
      conditionalPanel("!input.all_markers", selectizeInput("markers", "Biomarkers", choices = NULL, multiple = TRUE)),
      selectizeInput("covariates", "Adjustment covariates (optional)", choices = NULL, multiple = TRUE),
      div(class = "step-title", "5. Options"),
      fluidRow(column(6, numericInput("max_size", "Largest panel size", value = 3, min = 1, max = 6, step = 1)),
               column(6, selectInput("search_mode", "Search", choices = c("Beam from 3 (exhaustive 1–2)" = "beam", "Auto" = "auto", "Exhaustive" = "exhaustive")))),
      fluidRow(column(6, numericInput("beam", "Panels kept per size", value = 50, min = 5, max = 500, step = 5)),
               column(6, numericInput("exh_limit", "Exhaustive up to (panels)", value = 20000, min = 100, step = 1000))),
      fluidRow(column(4, numericInput("folds", "CV folds", value = 5, min = 2, max = 10, step = 1)),
               column(4, numericInput("search_repeats", "Search repeats", value = 5, min = 1, max = 50, step = 1)),
               column(4, numericInput("final_repeats", "Final repeats", value = 20, min = 1, max = 200, step = 1))),
      checkboxInput("onestep", "Rapid scoring: one Newton step per fold (exact refits near separation)", TRUE),
      fluidRow(column(6, checkboxInput("log2", "log2 transform", TRUE)),
               column(6, numericInput("delta", "Offset added first", value = 5e-7, min = 0))),
      checkboxInput("standardize", "Standardize to the reference group's mean and SD", TRUE),
      fluidRow(column(8, selectInput("missing", "Missing biomarker values",
                                     choices = c("k-nearest neighbours (k = 5)" = "knn", "Column median" = "median",
                                                 "Leave the participant out" = "complete"))),
               column(4, numericInput("seed", "Random seed", value = 123, step = 1))),
      div(style = "display:flex; gap:8px; margin-top:8px;",
          actionButton("run", "Run search", class = "btn-accent", icon = icon("play")),
          actionButton("stop", "Stop", class = "btn-soft", icon = icon("stop")))
    )),
    column(8,
      tabsetPanel(id = "an_tabs",
        tabPanel("Matching", br(), uiOutput("match_ui")),
        tabPanel("Search settings", br(), verbatimTextOutput("config_text")),
        tabPanel("Run log", br(), uiOutput("progress_ui"), verbatimTextOutput("run_log")),
        tabPanel("Results", br(), uiOutput("results_ui")),
        tabPanel("Train a model", br(), uiOutput("train_ui"))
      ))
  )
)

predict_tab <- tabPanel(
  "Risk Calculator", value = "predict",
  fluidRow(column(12,
    h3("Calculate individual risk"),
    radioButtons("model_choice", NULL, inline = TRUE,
                 choices = c("Selected DMFT model" = "published", "Model trained on your data" = "trained",
                             "Load a saved model" = "loaded")),
    conditionalPanel("input.model_choice == 'loaded'",
      fileInput("model_file", "Saved model file (.csv from \"Download model\")", accept = ".csv")),
    conditionalPanel("input.model_choice == 'published'",
      fileInput("selected_model_file", "Age/BMI model export from the rerun (.csv)", accept=".csv")),
    uiOutput("model_header"))),
  fluidRow(
    column(5, div(class = "card",
      h4("Measurements"),
      div(class = "raw-note", HTML(paste(
        "<b>Enter absolute raw concentration values</b> for each biomarker, exactly as reported by the laboratory and",
        "in the same units used to develop the model. Do not log-transform or standardize them, and do not convert",
        "them to relative abundance or percentages. Enter 0 if the biomarker was not detected."))),
      uiOutput("pred_inputs"),
      div(style = "display:flex; gap:8px; margin-top:6px; flex-wrap:wrap;",
          actionButton("calc", "Calculate risk", class = "btn-accent"),
          actionButton("clear", "Clear", class = "btn-soft"),
          uiOutput("fill_btn", inline = TRUE)))),
    column(7,
      div(class = "card",
        h4(textOutput("gauge_title", inline = TRUE)),
        plotOutput("gauge", height = "250px"),
        uiOutput("hero_ui")),
      div(class = "card", h4("What the risk estimate represents"), uiOutput("represents_ui")))
  ),
  fluidRow(column(12, uiOutput("age_bmi_card"))),
  fluidRow(column(12, div(class = "card",
    h4("How risk is calculated"),
    uiOutput("equation_ui"),
    uiOutput("calc_ui"))))
)

ui <- tagList(
  tags$head(tags$style(HTML(CSS)), tags$title("DMFT biomarker calculator")),
  navbarPage(title = "Biomarker analysis & prediction", id = "nav", collapsible = TRUE,
             start_tab, analysis_tab, predict_tab)
)


# ==========================================================================
# 11) SERVER
# ==========================================================================

PUB <- list(model = PUBLISHED_MODEL, msg = NULL)  # Browser version: load the new model CSV in Patient prediction

server <- function(input, output, session) {

  pub <- PUB
  rv  <- reactiveValues(pinfo = NULL, pinfo_name = NULL, bio = NULL, bio_name = NULL,
                        trained_candidate = NULL, trained = NULL, loaded = NULL, published = pub$model,
                        last_pred = NULL)
  job <- reactiveValues(state = "idle", env = NULL, log = character(0), msg = "", frac = 0, results = NULL)

  observeEvent(input$go_analysis, updateNavbarPage(session, "nav", selected = "analysis"))
  observeEvent(input$go_predict,  updateNavbarPage(session, "nav", selected = "predict"))

  # ---- data input --------------------------------------------------------
  read_upload <- function(f) tryCatch(read_table_file(f$datapath, f$name), error = function(e) {
    showNotification(conditionMessage(e), type = "error", duration = 8); NULL })
  observeEvent(input$pinfo_file, { d <- read_upload(input$pinfo_file); if (!is.null(d)) { rv$pinfo <- d; rv$pinfo_name <- input$pinfo_file$name } })
  observeEvent(input$bio_file,   { d <- read_upload(input$bio_file);   if (!is.null(d)) { rv$bio <- d; rv$bio_name <- input$bio_file$name } })
  observeEvent(input$load_example, {
    ex <- make_example_data()
    rv$pinfo <- ex$pinfo; rv$pinfo_name <- "example_participants.csv"
    rv$bio <- ex$bio; rv$bio_name <- "example_biomarkers.csv"
    updateCheckboxInput(session, "one_file", value = FALSE)
    showNotification("Simulated example data loaded.", type = "message")
  })
  output$dl_example_info <- downloadHandler(
    filename = function() "example_participants.csv",
    content = function(file) utils::write.csv(make_example_data()$pinfo, file, row.names = FALSE, na = ""))
  output$dl_example_bio <- downloadHandler(
    filename = function() "example_biomarkers.csv",
    content = function(file) utils::write.csv(make_example_data()$bio, file, row.names = FALSE, na = ""))

  observeEvent(rv$pinfo, {
    d <- rv$pinfo; nm <- names(d)
    updateSelectInput(session, "pid", choices = nm, selected = guess_id_col(d))
    dm <- grep("dmft", nm, ignore.case = TRUE, value = TRUE)
    comp <- grep("^(decayed|missing|filled|missing_teeth|restored_teeth|teeth_wcarries)$", nm, ignore.case = TRUE, value = TRUE)
    sel <- if (length(dm)) dm[1] else if (length(comp) == 3) comp else character(0)
    updateSelectizeInput(session, "outcome_cols", choices = nm, selected = sel)
    cov_guess <- c(intersect(c("Age_DOS", "Age"), nm)[1], intersect("BMI", nm))
    cov_guess <- cov_guess[!is.na(cov_guess)]
    updateSelectizeInput(session, "covariates", choices = nm, selected = cov_guess)
  })
  observeEvent(rv$bio, updateSelectInput(session, "bid", choices = names(rv$bio), selected = guess_id_col(rv$bio)))

  matched <- reactive({
    req(rv$pinfo, input$pid)
    if (isTRUE(input$one_file)) {
      d <- rv$pinfo
      return(list(ok = TRUE, data = d, id = input$pid, bio_cols = numeric_cols(d, c(input$pid, input$outcome_cols, input$covariates)),
                  n_info = nrow(d), n_bio = nrow(d), n_matched = nrow(d), only_info = character(0), only_bio = character(0)))
    }
    req(rv$bio, input$bid, input$pid %in% names(rv$pinfo), input$bid %in% names(rv$bio))
    match_files(rv$pinfo, input$pid, rv$bio, input$bid)
  })

  marker_choices <- reactive({
    m <- matched(); req(isTRUE(m$ok))
    numeric_cols(m$data[, m$bio_cols, drop = FALSE])
  })
  observe(updateSelectizeInput(session, "markers", choices = marker_choices(), selected = marker_choices(), server = TRUE))

  outcome_vals <- reactive({
    m <- matched(); req(isTRUE(m$ok), length(input$outcome_cols) > 0, all(input$outcome_cols %in% names(m$data)))
    outcome_value(m$data, input$outcome_cols)
  })
  observe({
    v <- outcome_vals()
    if (length(input$outcome_cols) == 1 && !is.numeric(v)) {
      lv <- sort(unique(trimws(as.character(v[!is.na(v)]))))
      updateSelectInput(session, "ref_level", choices = lv, selected = lv[1])
      updateSelectInput(session, "case_level", choices = lv, selected = lv[min(2, length(lv))])
      if (input$grp_mode != "levels") updateRadioButtons(session, "grp_mode", selected = "levels")
    }
  })
  observeEvent(outcome_vals(), {
    v <- as_num(outcome_vals())
    if (any(is.finite(v))) updateNumericInput(session, "cutoff", value = stats::median(v, na.rm = TRUE))
  })
  output$cutoff_hint <- renderUI({
    v <- tryCatch(as_num(outcome_vals()), error = function(e) NULL); req(v, any(is.finite(v)), is.finite(input$cutoff))
    p(class = "help-note", sprintf("Median %s; %d participants ≤ %s and %d > %s.", format(stats::median(v, na.rm = TRUE)),
                                   sum(v <= input$cutoff, na.rm = TRUE), format(input$cutoff), sum(v > input$cutoff, na.rm = TRUE),
                                   format(input$cutoff)))
  })

  output$match_ui <- renderUI({
    if (is.null(rv$pinfo)) return(p(class = "lead-muted", "Upload the participant information file (or load the example data) to begin."))
    if (!isTRUE(input$one_file) && is.null(rv$bio)) return(p(class = "lead-muted", "Now upload the raw biomarker file."))
    m <- tryCatch(matched(), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (!isTRUE(m$ok)) return(div(class = "notice", m$msg))
    tagList(
      div(class = "okmsg", sprintf("%d participants matched by ID.", m$n_matched)),
      tags$table(class = "table table-condensed", style = "max-width:560px;",
        tags$tr(tags$td("Participant information file"), tags$td(rv$pinfo_name %||% ""), tags$td(sprintf("%d participants", m$n_info))),
        if (!isTRUE(input$one_file)) tags$tr(tags$td("Biomarker file"), tags$td(rv$bio_name %||% ""), tags$td(sprintf("%d participants", m$n_bio))),
        tags$tr(tags$td("Matched"), tags$td(""), tags$td(sprintf("%d participants", m$n_matched))),
        tags$tr(tags$td("Biomarker columns"), tags$td(""), tags$td(length(marker_choices())))),
      if (length(m$only_info)) p(class = "help-note", sprintf("Only in the participant file (%d): %s", length(m$only_info),
                                                              paste(head(m$only_info, 12), collapse = ", "))),
      if (length(m$only_bio)) p(class = "help-note", sprintf("Only in the biomarker file (%d): %s", length(m$only_bio),
                                                             paste(head(m$only_bio, 12), collapse = ", ")))
    )
  })

  # ---- settings ----------------------------------------------------------
  settings <- reactive({
    list(folds = as.integer(input$folds %||% 5), search_repeats = as.integer(input$search_repeats %||% 5),
         final_repeats = max(as.integer(input$final_repeats %||% 20), as.integer(input$search_repeats %||% 5)),
         max_size = as.integer(input$max_size %||% 3), search_mode = input$search_mode, beam = as.integer(input$beam %||% 50),
         exhaustive_limit = as.numeric(input$exh_limit %||% 20000), onestep = isTRUE(input$onestep), eta_max = 10,
         transform = if (isTRUE(input$log2)) "log2" else "none", delta = as.numeric(input$delta %||% 5e-7),
         standardize = isTRUE(input$standardize), missing = input$missing, seed = as.integer(input$seed %||% 123))
  })
  groups <- reactive({
    list(mode = input$grp_mode, cutoff = input$cutoff, ref_max = input$ref_max, case_min = input$case_min,
         ref_level = input$ref_level, case_level = input$case_level,
         ref_name = trimws(input$ref_name %||% "Reference"), case_name = trimws(input$case_name %||% "Case"))
  })
  selected_markers <- reactive(if (isTRUE(input$all_markers)) marker_choices() else intersect(input$markers, marker_choices()))

  output$config_text <- renderText({
    m <- tryCatch(matched(), error = function(e) NULL)
    if (is.null(m) || !isTRUE(m$ok) || !length(input$outcome_cols)) return("Upload and match the files, and choose the outcome, to see the configuration.")
    s <- settings(); g <- groups()
    gg <- tryCatch(define_groups(outcome_vals(), g), error = function(e) NULL)
    nref <- if (is.null(gg)) NA else sum(gg$y == 0, na.rm = TRUE); ncase <- if (is.null(gg)) NA else sum(gg$y == 1, na.rm = TRUE)
    p <- length(selected_markers()); K <- min(s$max_size, p)
    lv <- vapply(seq_len(max(K, 1)), function(k) {
      if (k <= 2) return(sprintf("k = %d: exhaustive | %s candidate models", k, fmt_int(choose(p,k))))
      n_all <- choose(p, k)
      exh <- s$search_mode == "exhaustive" || (s$search_mode == "auto" && n_all <= s$exhaustive_limit)
      if (exh && n_all <= HARD_EXHAUSTIVE_MAX) sprintf("k = %d: exhaustive | %s candidate models", k, fmt_int(n_all))
      else sprintf("k = %d: expand the %d retained panels | up to %s candidate models", k, s$beam,
                   fmt_int(min(n_all, s$beam * (p - k + 1))))
    }, character(1))
    paste(c(
      sprintf("Participants: %d | %s: %s | %s: %s", m$n_matched, g$ref_name, nref, g$case_name, ncase),
      sprintf("Outcome: %s | %s", paste(input$outcome_cols, collapse = " + "), if (is.null(gg)) "" else gg$desc),
      sprintf("Candidate biomarkers: %d | Panel sizes: 1–%d", p, K),
      sprintf("Adjustment covariates: %s", if (length(input$covariates)) paste(input$covariates, collapse = ", ") else "none"),
      sprintf("Cross-validation: %d folds × %d repetitions (search), %d repetitions (final metrics) | seed %d",
              s$folds, s$search_repeats, s$final_repeats, s$seed),
      sprintf("Positive class: %s | Selection metric: balanced accuracy; AUC breaks ties", g$case_name),
      lv,
      sprintf("Rapid scoring: %s", if (s$onestep) "one Newton step per fold from the full-data fit; exact refits near separation" else "exact refit in every fold"),
      sprintf("Raw values → %s%s | missing values: %s",
              if (s$transform == "log2") sprintf("log2(x + %s)", format(s$delta)) else "no transform",
              if (s$standardize) sprintf(" → standardized to the %s mean and SD", g$ref_name) else "",
              c(knn = "k-nearest neighbours (k = 5)", median = "column median", complete = "participant left out")[[s$missing]]),
      "Outputs: ranked panels, performance metrics, ROC curves and trainable models"), collapse = "\n")
  })

  # ---- running the search ------------------------------------------------
  observeEvent(input$run, {
    if (identical(job$state, "running")) return()
    m <- tryCatch(matched(), error = function(e) NULL)
    if (is.null(m) || !isTRUE(m$ok)) { showNotification("Upload and match the files first.", type = "error"); return() }
    if (!length(input$outcome_cols)) { showNotification("Choose the outcome column(s).", type = "error"); return() }
    mk <- selected_markers()
    if (!length(mk)) { showNotification("Choose at least one biomarker.", type = "error"); return() }
    s <- settings()
    prep <- tryCatch(withProgress(message = "Preparing the data", value = 0.5,
                                  prepare_analysis(m$data, input$pid, input$outcome_cols, groups(), mk, input$covariates, s)),
                     error = function(e) e)
    if (inherits(prep, "error")) { showNotification(conditionMessage(prep), type = "error", duration = 10); return() }
    env <- new_search(prep, s)
    proj <- project_runtime(env)
    env$log <- c(sprintf("Data: %d participants (%d %s, %d %s), %d biomarkers, %d covariate term(s).", length(prep$y),
                         prep$n_ref, prep$groups$ref_name, prep$n_case, prep$groups$case_name, ncol(prep$M), ncol(prep$Z)),
                 prep$notes,
                 sprintf(paste("Projected search: %s candidate panels at most, taking about %s with rapid scoring and up to %s",
                               "if most panels need exact refits (common with sparse taxa)."),
                         fmt_int(sum(proj$counts)), fmt_time(proj$seconds), fmt_time(proj$seconds_max)),
                 env$log)
    job$env <- env; job$log <- env$log; job$msg <- "Starting..."; job$frac <- 0; job$results <- NULL
    job$state <- "running"
    updateTabsetPanel(session, "an_tabs", selected = "Run log")
  })
  observeEvent(input$stop, {
    if (identical(job$state, "running")) { job$state <- "stopped"; job$log <- c(isolate(job$log), "Stopped by the user."); job$msg <- "Stopped." }
  })
  observe({
    req(identical(job$state, "running"))
    env <- isolate(job$env)
    st <- tryCatch(search_tick(env, 0.3), error = function(e) e)
    if (inherits(st, "error")) {
      job$state <- "error"; job$log <- c(env$log, paste("Error:", conditionMessage(st))); job$msg <- "Stopped by an error."
      return()
    }
    job$log <- env$log; job$msg <- env$msg; job$frac <- env$frac %||% 0
    if (identical(env$phase, "done")) {
      job$results <- env$results; job$state <- "done"
      updateTabsetPanel(session, "an_tabs", selected = "Results")
    } else invalidateLater(10)
  })
  output$progress_ui <- renderUI({
    if (identical(job$state, "idle")) return(p(class = "lead-muted", "Set the options on the left and click Run search."))
    tagList(div(class = "progress-wrap", div(class = "progress-fill", style = sprintf("width:%.1f%%;", 100 * (job$frac %||% 0)))),
            p(job$msg))
  })
  output$run_log <- renderText(paste(job$log, collapse = "\n"))

  # ---- results -------------------------------------------------------------
  output$results_ui <- renderUI({
    res <- job$results
    if (is.null(res)) return(p(class = "lead-muted", "Results appear here when a search finishes."))
    ks <- which(!vapply(res$tables, is.null, logical(1)))
    tagList(
      h4("Best panel at each size"),
      div(class = "res-wrap", tableOutput("summary_tab")),
      fluidRow(
        column(7, plotOutput("roc_plot", height = "380px")),
        column(5, p(class = "help-note", style = "margin-top:10px;",
                    sprintf(paste("ROC curves of the best panel at each size from pooled out-of-fold predictions (%d folds × %d",
                                  "repetitions). Symbols mark each panel's Youden threshold. Panels were selected and evaluated on the",
                                  "same data, so these estimates are optimistic."), res$s$folds, res$s$final_repeats)),
               if (HAS_WRITEXL) downloadButton("dl_results", "Download Excel results", class = "btn-accent")
               else downloadButton("dl_results", "Download results (CSV)", class = "btn-accent"))),
      hr(),
      fluidRow(column(4, selectInput("res_k", "Panel size", choices = ks, selected = max(ks))),
               column(4, numericInput("res_n", "Rows shown", value = 15, min = 1, max = 500))),
      div(class = "res-wrap", tableOutput("res_table")),
      p(class = "help-note", "Balanced accuracy range: 2.5th–97.5th percentiles across repetitions (variation from the fold split, not a confidence interval).")
    )
  })
  output$summary_tab <- renderTable(summary_table(job$results), striped = TRUE, bordered = FALSE,
                                    sanitize.text.function = function(x) x)
  output$roc_plot <- renderPlot(draw_roc(job$results), res = 96)
  output$res_table <- renderTable({
    res <- job$results; req(res, input$res_k)
    t <- res$tables[[as.integer(input$res_k)]]; req(t)
    t <- head(t, max(1, input$res_n %||% 15))
    data.frame(Rank = t$Rank, Panel = panel_html(t, seq_len(nrow(t))), `Balanced accuracy` = fmt_pct(t$Balanced_Accuracy),
               `BA range` = sprintf("%s–%s", fmt_pct(t$BA_2.5th_pct), fmt_pct(t$BA_97.5th_pct)),
               Sensitivity = fmt_pct(t$Sensitivity), Specificity = fmt_pct(t$Specificity), Accuracy = fmt_pct(t$Accuracy),
               PPV = fmt_pct(t$PPV), NPV = fmt_pct(t$NPV), AUC = sprintf("%.3f", t$AUC), Threshold = sprintf("%.3f", t$Threshold),
               check.names = FALSE)
  }, striped = TRUE, sanitize.text.function = function(x) x)

  results_params <- function(res) {
    s <- res$s; pr <- res$prep; g <- pr$groups
    data.frame(Parameter = c("Participants", paste0("  n (", g$ref_name, ")"), paste0("  n (", g$case_name, ")"),
                             "Groups", "Outcome column(s)", "Candidate biomarkers", "Adjustment covariates",
                             "Cross-validation", "Rapid scoring", "Search", "Panels kept per size", "Transform",
                             "Standardization", "Missing biomarker values", "Data notes", "Panels scored", "Run time",
                             "Positive class", "Selection"),
               Value = c(length(pr$y), pr$n_ref, pr$n_case, pr$group_desc, paste(pr$outcome_cols, collapse = " + "),
                         ncol(pr$M), if (length(pr$covariates)) paste(pr$covariates, collapse = ", ") else "none",
                         sprintf("%d folds; %d repetitions for the search, %d for the final metrics; seed %d", s$folds,
                                 s$search_repeats, s$final_repeats, s$seed),
                         if (s$onestep) "one Newton step per fold from the full-data fit; exact refits when the full fit did not converge or |linear predictor| > 10" else "exact refit in every fold",
                         paste(vapply(seq_along(res$level_info), function(k) sprintf("k = %d: %s (%s panels)", k, res$level_info[[k]]$mode,
                                                                                       fmt_int(res$level_info[[k]]$n)), character(1)), collapse = "; "),
                         s$beam, if (s$transform == "log2") sprintf("log2(raw + %s)", format(s$delta)) else "none",
                         if (s$standardize) paste("reference group (", g$ref_name, ") mean and SD") else "none",
                         c(knn = "k-nearest neighbours (k = 5)", median = "column median", complete = "participant left out")[[s$missing]],
                         if (length(pr$notes)) paste(pr$notes, collapse = " ") else "none",
                         res$n_scored, fmt_time(res$elapsed), g$case_name,
                         "Balanced accuracy at the Youden threshold of pooled out-of-fold predictions; ties broken by AUC"),
               stringsAsFactors = FALSE)
  }
  output$dl_results <- downloadHandler(
    filename = function() if (HAS_WRITEXL) "biomarker_panel_results.xlsx" else "biomarker_panel_results.csv",
    content = function(file) {
      res <- job$results
      tabs <- lapply(res$tables, function(t) { if (is.null(t)) return(NULL); attr(t, "idx") <- NULL; t })
      if (HAS_WRITEXL) {
        sheets <- c(list(Analysis_Parameters = results_params(res)),
                    stats::setNames(tabs[!vapply(tabs, is.null, logical(1))],
                                    paste0("Results_", which(!vapply(tabs, is.null, logical(1))), "_Biomarkers")))
        writexl::write_xlsx(sheets, file)
      } else {
        out <- do.call(rbind, lapply(seq_along(tabs), function(k) {
          t <- tabs[[k]]; if (is.null(t)) return(NULL)
          data.frame(Panel_size = k, t[, setdiff(names(t), paste0("Biomarker", seq_len(k)))], check.names = FALSE)
        }))
        utils::write.csv(out, file, row.names = FALSE)
      }
    })

  # ---- training a model --------------------------------------------------
  output$train_ui <- renderUI({
    res <- job$results
    if (is.null(res)) return(p(class = "lead-muted", "Run a search first; then choose a panel here to train a model on all participants."))
    ks <- which(!vapply(res$tables, is.null, logical(1)))
    tagList(
      p("Choose a panel. The model is refitted to all participants, and its coefficients are converted so that it can",
        "be applied directly to absolute raw concentration values."),
      fluidRow(column(4, selectInput("tr_k", "Panel size", choices = ks, selected = max(ks))),
               column(4, numericInput("tr_rank", "Rank within that size", value = 1, min = 1, step = 1)),
               column(4, br(), actionButton("tr_fit", "Train model", class = "btn-accent"))),
      fluidRow(column(12, uiOutput("tr_panel_name"))),
      fluidRow(column(12, uiOutput("tr_result"))))
  })
  output$tr_panel_name <- renderUI({
    res <- job$results; req(res, input$tr_k, input$tr_rank)
    t <- res$tables[[as.integer(input$tr_k)]]; req(t, input$tr_rank <= nrow(t))
    p(class = "help-note", HTML(sprintf("Panel: <b>%s</b> (cross-validated balanced accuracy %s).",
                                        htmltools::htmlEscape(t$Panel[input$tr_rank]), fmt_pct(t$Balanced_Accuracy[input$tr_rank]))))
  })
  observeEvent(job$results, rv$trained_candidate <- NULL, ignoreNULL = FALSE)
  observeEvent(input$tr_fit, {
    res <- job$results; req(res)
    k <- as.integer(input$tr_k); r <- as.integer(input$tr_rank)
    if (is.na(r) || r < 1 || r > nrow(res$tables[[k]])) { showNotification("That rank does not exist for this panel size.", type = "error"); return() }
    m <- tryCatch(fit_panel_model(res, k, r), error = function(e) e)
    if (inherits(m, "error")) { showNotification(conditionMessage(m), type = "error"); return() }
    rv$trained_candidate <- m
  })
  output$tr_result <- renderUI({
    m <- rv$trained_candidate; req(m)
    ft <- m$fit_table
    tagList(
      h4(m$name),
      if (!is.null(m$caution)) div(class = "notice", m$caution),
      equation_html(m),
      tags$table(class = "calc-table", style = "margin-top:10px;",
        tags$tr(tags$th("Term"), tags$th("Coefficient (model scale)"), tags$th("SE"), tags$th("p"), tags$th("Coefficient for raw input")),
        lapply(seq_len(nrow(ft)), function(i) tags$tr(tags$td(ft$Term[i]), tags$td(HTML(fmt_num(ft$Estimate_model_scale[i]))),
                                                      tags$td(HTML(fmt_num(ft$SE[i]))), tags$td(HTML(fmt_num(ft$p_value[i], 3))),
                                                      tags$td(HTML(fmt_num(ft$Coefficient_raw_input[i])))))),
      p(class = "help-note", style = "margin-top:8px;",
        sprintf("Classification threshold %s (Youden, cross-validated). Model scale: biomarkers as transformed and standardized for fitting.",
                format(round(m$threshold, 4)))),
      div(style = "display:flex; gap:8px;",
          actionButton("tr_use", "Use this model in patient prediction", class = "btn-accent"),
          downloadButton("tr_download", "Download model (CSV)", class = "btn-soft")))
  })
  observeEvent(input$tr_use, {
    req(rv$trained_candidate)
    rv$trained <- rv$trained_candidate
    updateRadioButtons(session, "model_choice", selected = "trained")
    updateNavbarPage(session, "nav", selected = "predict")
  })
  output$tr_download <- downloadHandler(
    filename = function() "trained_model.csv",
    content = function(file) utils::write.csv(model_to_table(rv$trained_candidate), file, row.names = FALSE, na = ""))

  # ---- patient prediction ------------------------------------------------
  observeEvent(input$model_file, {
    tab <- tryCatch(utils::read.csv(input$model_file$datapath, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) e)
    m <- if (inherits(tab, "error")) tab else tryCatch(model_from_table(tab, input$model_file$name), error = function(e) e)
    if (inherits(m, "error")) { showNotification(conditionMessage(m), type = "error", duration = 8); return() }
    rv$loaded <- m
  })
  observeEvent(input$selected_model_file, {
    req(input$selected_model_file)
    m <- tryCatch(selected_age_bmi_model(utils::read.csv(input$selected_model_file$datapath,
                        stringsAsFactors=FALSE),input$selected_model_file$name), error=function(e)e)
    if (inherits(m,"error")) {showNotification(conditionMessage(m),type="error",duration=10);return()}
    rv$published <- m
  })
  current_model <- reactive(switch(input$model_choice %||% "published", published = rv$published,
                                   trained = rv$trained, loaded = rv$loaded))
  observeEvent(current_model(), rv$last_pred <- NULL, ignoreNULL = FALSE)

  output$model_header <- renderUI({
    m <- current_model()
    if (is.null(m)) return(div(class = "notice", switch(input$model_choice,
      trained = "No trained model yet: run a search in Train Your Own Data, train a model from a panel and click “Use this model”.",
      loaded = "Choose a saved model file.")))
    tagList(
      p(class = "lead-muted", HTML(paste0("<b>", htmltools::htmlEscape(m$name), "</b>. ", md_html(m$description %||% "")))),
      if (identical(m$source, "published") && !is.null(pub$msg)) div(class = if (model_ready(m)) "okmsg" else "notice", pub$msg),
      if (!model_ready(m)) div(class = "notice", HTML(paste(
        "The age/BMI model must be rerun and its complete coefficient export loaded before calculating risk.",
        "Use the <b>Age/BMI model export</b> input above to upload <code>dmft_age_bmi_model.csv</code>, generated by 03_Figures_DMFT_Age_BMI.R.",
        "For desktop use, you can also place that file next to app.R and restart."))))
  })

  output$pred_inputs <- renderUI({
    m <- current_model(); req(m)
    mk <- lapply(seq_len(nrow(m$markers)), function(j) {
      numericInput(paste0("mk_", j), HTML(sprintf("%s <span class='help-note'>(%s) absolute raw concentration</span>",
                                                  md_html(m$markers$label[j]), htmltools::htmlEscape(m$markers$abbr[j]))),
                   value = NA, min = 0)
    })
    cv <- if (nrow(m$covariates)) {
      vars <- unique(m$covariates$var)
      lapply(seq_along(vars), function(i) {
        v <- vars[i]; row <- m$covariates[m$covariates$var == v, ][1, ]
        if (row$type == "factor") {
          lv <- m$cov_levels[[v]] %||% c("(reference level)", m$covariates$level[m$covariates$var == v])
          selectInput(paste0("cv_", i), pretty_name(v), choices = level_choices(m, v, lv))
        } else if (is_age_variable(v)) {
          limits <- age_slider_limits(m, v)
          sliderInput(paste0("cv_", i), "Age (years)", min = limits[1], max = limits[2],
                      value = max(limits[1], min(limits[2], m$cov_numeric[[v]] %||% 40)), step = 1, sep = "",
                      ticks = FALSE)
        } else numericInput(paste0("cv_", i), sprintf("%s%s", row$label, if (nzchar(row$unit)) sprintf(" (%s)", row$unit) else ""),
                            value = NA)
      })
    }
    tagList(div(class = "step-title", "Biomarkers"), mk,
            if (length(cv)) tagList(div(class = "step-title", "Participant characteristics"),
                                    if (!identical(m$source, "published"))
                                      p(class = "help-note", "Use the same units as the data the model was trained on."),
                                    cv))
  })
  output$fill_btn <- renderUI({
    m <- current_model()
    if (is.null(m) || !model_ready(m)) return(NULL)
    
    if (!is.null(m$markers$raw_median))
      actionButton("fill", "Fill with training-data medians", class = "btn-soft")
  })
  observeEvent(input$example, {
    m <- current_model(); req(m, identical(m$source, "published"))
    for (j in seq_len(nrow(m$markers))) updateNumericInput(session, paste0("mk_", j), value = EXAMPLE$markers[j])
    for (i in seq_along(unique(m$covariates$var))) {
      v <- unique(m$covariates$var)[i]
      if (is_age_variable(v)) updateSliderInput(session, paste0("cv_", i), value = 40)
      if (v == "BMI") updateNumericInput(session, paste0("cv_", i), value = 25)
    }
    rv$last_pred <- NULL
    showNotification("Illustrative values filled in (not a real participant). Click Calculate risk.",
                     type = "message", duration = 6)
  })
  observeEvent(input$fill, {
    m <- current_model(); req(m, m$markers$raw_median)
    for (j in seq_len(nrow(m$markers))) updateNumericInput(session, paste0("mk_", j), value = signif(m$markers$raw_median[j], 4))
    vars <- unique(m$covariates$var)
    for (i in seq_along(vars)) if (!is.null(m$cov_numeric[[vars[i]]])) {
      if (is_age_variable(vars[i])) updateSliderInput(session, paste0("cv_", i), value = m$cov_numeric[[vars[i]]])
      else updateNumericInput(session, paste0("cv_", i), value = m$cov_numeric[[vars[i]]])
    }
  })
  observeEvent(input$clear, {
    m <- current_model(); req(m)
    for (j in seq_len(nrow(m$markers))) updateNumericInput(session, paste0("mk_", j), value = NA)
    vars <- unique(m$covariates$var)
    for (i in seq_along(vars)) if (m$covariates$type[m$covariates$var == vars[i]][1] == "numeric") {
      if (is_age_variable(vars[i])) updateSliderInput(session, paste0("cv_", i), value = 40)
      else updateNumericInput(session, paste0("cv_", i), value = NA)
    }
    rv$last_pred <- NULL
  })
  observeEvent(input$calc, {
    m <- current_model(); req(m)
    if (!model_ready(m)) { showNotification("This model's coefficients are incomplete, so no risk estimate can be calculated.", type = "error"); return() }
    x <- vapply(seq_len(nrow(m$markers)), function(j) as.numeric(input[[paste0("mk_", j)]] %||% NA), numeric(1))
    if (any(!is.finite(x))) { showNotification("Enter a finite value for every biomarker (0 if not detected).", type = "error"); return() }
    if (any(x < 0)) { showNotification("Absolute raw concentration values cannot be negative.", type = "error"); return() }
    vars <- unique(m$covariates$var)
    cv <- stats::setNames(lapply(seq_along(vars), function(i) input[[paste0("cv_", i)]]), vars)
    num <- vars[vapply(vars, function(v) m$covariates$type[m$covariates$var == v][1] == "numeric", logical(1))]
    if (length(num) && any(vapply(num, function(v) is.null(cv[[v]]) || is.na(cv[[v]]), logical(1)))) {
      showNotification("Enter every participant characteristic.", type = "error"); return()
    }
    warn <- character(0)
    if ("Teeth" %in% vars && is.finite(cv$Teeth) && (cv$Teeth < 0 || cv$Teeth > 32)) warn <- c(warn, "Number of teeth should be between 0 and 32.")
    if ("Age_DOS" %in% vars && is.finite(cv$Age_DOS) && (cv$Age_DOS < 18 || cv$Age_DOS > 100)) warn <- c(warn, "The model was developed in adults.")
    if ("BMI" %in% vars && is.finite(cv$BMI) && (cv$BMI < 12 || cv$BMI > 70)) warn <- c(warn, "BMI is outside the usual range; check the value.")
    out <- tryCatch(predict_model(m, x, cv), error = function(e) e)
    if (inherits(out, "error")) { showNotification(conditionMessage(out), type = "error"); return() }
    out$warn <- warn
    rv$last_pred <- out
  })


  comparison_inputs <- reactive({
    m <- current_model(); req(m)
    x <- vapply(seq_len(nrow(m$markers)), function(j) as.numeric(input[[paste0("mk_", j)]] %||% NA), numeric(1))
    vars <- unique(m$covariates$var)
    cv <- stats::setNames(lapply(seq_along(vars), function(i) input[[paste0("cv_", i)]]), vars)
    list(model = m, markers = x, covariates = cv)
  })
  live_prediction <- reactive({
    d <- comparison_inputs()
    tryCatch(predict_model(d$model, d$markers, d$covariates), error = function(e) NULL)
  })
  # Once calculated, keep the gauge and equation synchronized with slider changes.
  observeEvent(live_prediction(), {
    if (!is.null(isolate(rv$last_pred))) rv$last_pred <- live_prediction()
  }, ignoreNULL = FALSE)

  output$age_bmi_card <- renderUI({
    m <- current_model()
    if (is.null(age_bmi_variables(m))) return(NULL)
    div(class = "card", h4("Compare risk across age and BMI"),
        p(class = "help-note", "Move the age slider or change BMI above. The two curves compare the entered BMI with BMI five units higher, holding biomarker values fixed. The filled point marks the entered age and BMI."),
        plotOutput("age_bmi_plot", height = "360px"),
        uiOutput("age_bmi_values"),
        p(class = "help-note", "These curves show model scenarios, not changes observed over time or the effect of aging. The displayed age range is for exploration and is not a validated prediction range."))
  })
  comparison_curves <- reactive({
    d <- comparison_inputs(); vars <- age_bmi_variables(d$model); req(vars)
    validate(need(model_ready(d$model), "Load complete model coefficients to show the comparison."),
             need(!is.null(live_prediction()), "Enter all biomarker values and participant characteristics to show the comparison."))
    age <- as.numeric(d$covariates[[vars$age]])
    limits <- age_slider_limits(d$model, vars$age)
    curves <- age_bmi_curve_data(d$model, d$markers, d$covariates,
                                sort(unique(c(seq(limits[1], limits[2], by = 1), age))))
    list(curves = curves, age = age, bmi = as.character(d$covariates[[vars$bmi]]), threshold = d$model$threshold,
         labels = paste("BMI", unique(curves$bmi)), model = d$model, bmi_var = vars$bmi)
  })
  output$age_bmi_plot <- renderPlot({
    d <- comparison_curves()
    draw_age_bmi_curves(d$curves, d$age, d$bmi, d$threshold, d$labels)
  }, res = 110)
  output$age_bmi_values <- renderUI({
    d <- comparison_curves()
    rows <- d$curves[d$curves$age == d$age, ]
    p(sprintf("At age %s: %s.", d$age, paste(sprintf("%s %.2f%%", paste("BMI", rows$bmi),
                                                     100 * rows$risk), collapse = "; ")))
  })

  output$gauge_title <- renderText({
    m <- current_model(); if (is.null(m)) "Estimated risk" else paste("Estimated risk of", m$case_label)
  })
  output$gauge <- renderPlot({
    m <- current_model(); pr <- rv$last_pred
    req(m, model_ready(m))
    draw_gauge(if (is.null(pr)) NA else pr$p, if (is.null(m)) NA else m$threshold,
               lo_lab = "Lower", hi_lab = "Higher")
  }, res = 110)
  output$hero_ui <- renderUI({
    m <- current_model(); pr <- rv$last_pred
    if (is.null(m)) return(NULL)
    if (is.null(pr)) return(p(class = "hero-sub", "Enter the values and click Calculate risk."))
    above <- pr$p >= m$threshold
    tagList(
      div(class = "hero", fmt_prob(pr$p)),
      div(class = "hero-sub", sprintf("estimated risk of %s versus %s", m$case_label, m$ref_label)),
      div(class = paste("status", if (above) "hi" else "lo"),
          HTML(sprintf("%s %s the %s classification threshold: classified as %s", if (above) "&#9650;" else "&#9660;",
                       if (above) "At or above" else "Below", fmt_pct(m$threshold),
                       htmltools::htmlEscape(if (above) m$case_label else m$ref_label)))),
      if (length(pr$warn)) div(class = "notice", style = "margin-top:10px;", HTML(md_html(paste(pr$warn, collapse = " ")))))
  })
  output$represents_ui <- renderUI({
    m <- current_model(); req(m)
    perf <- if (!is.null(m$performance)) paste(sprintf("%s %s", names(m$performance),
                                                      ifelse(names(m$performance) == "AUC", sprintf("%.3f", m$performance), fmt_pct(m$performance))),
                                              collapse = ", ")
    tagList(
      p(m$cohort),
      if (identical(m$source, "published")) p("Risk means the estimated likelihood of having high DMFT at assessment (DMFT > 13); it does not predict future disease."),
      if (!is.null(perf)) p(sprintf("Cross-validated performance at the %s threshold: %s.", fmt_pct(m$threshold), perf),
                            span(class = "help-note", m$performance_note)),
      if (!is.null(m$caution)) p(HTML(md_html(m$caution))),
      p(class = "help-note", "Research estimate only: calibration and performance in new participants have not been established,",
        "and the risk estimate is not a diagnosis."))
  })
  output$equation_ui <- renderUI({
    m <- current_model(); req(m, model_ready(m))
    tagList(p("The model is a logistic regression. Its linear predictor, with the fitted coefficients, is"),
            equation_html(m),
            p(class = "help-note", "Coefficients are shown to four significant digits; the calculation uses the full values."))
  })
  output$calc_ui <- renderUI({
    m <- current_model(); pr <- rv$last_pred
    if (is.null(m) || is.null(pr)) return(NULL)
    t <- pr$terms
    rows <- lapply(seq_len(nrow(t)), function(i) {
      tx <- if (t$kind[i] == "marker" && m$transform == "log2")
        sprintf("log<sub>2</sub>(%s + &delta;) = %s", fmt_num(t$Input[i], 7), fmt_num(t$Transformed[i]))
      else if (!is.null(t$Shown) && !is.na(t$Shown[i])) htmltools::htmlEscape(t$Shown[i]) else fmt_num(t$Input[i], 7)
      nm <- if (grepl(t$Abbr[i], t$Term[i], fixed = TRUE)) md_html(t$Term[i])
            else sprintf("%s (%s)", md_html(t$Term[i]), htmltools::htmlEscape(t$Abbr[i]))
      tags$tr(tags$td(HTML(nm)), tags$td(HTML(tx)),
              tags$td(HTML(fmt_num(t$Coefficient[i]))), tags$td(HTML(fmt_num(t$Contribution[i]))))
    })
    tagList(
      h4("Calculation for the values entered", style = "margin-top:18px;"),
      tags$table(class = "calc-table",
        tags$tr(tags$th("Term"), tags$th("Value used"), tags$th("Coefficient"), tags$th("Coefficient × value")),
        tags$tr(tags$td("Intercept"), tags$td(""), tags$td(""), tags$td(HTML(fmt_num(m$intercept)))),
        rows,
        tags$tr(class = "total", tags$td(HTML("Linear predictor &eta; (sum)")), tags$td(""), tags$td(""), tags$td(HTML(fmt_num(pr$eta))))),
      div(class = "equation", style = "margin-top:12px;",
          HTML(sprintf("<i>p</i> = 1 / [1 + exp(%s%s)] = %s, i.e. <b>%s</b>", MINUS,
                       if (pr$eta < 0) sprintf("(%s)", fmt_num(pr$eta)) else fmt_num(pr$eta),
                       if (pr$p > 0.999) sprintf("1 %s %s", MINUS, fmt_num(1 - pr$p, 3)) else fmt_num(pr$p, 4),
                       fmt_prob(pr$p)))))
  })
}

if (!isTRUE(getOption("dmft.app.test"))) shinyApp(ui, server)


