# ==========================================================================
# CALIBRATION OF THE DMFT PANELS (TRIPOD+AI items 12e and 23a)
# ==========================================================================
# Run after Figures_DMFT.V3.R in the same R session:
#
#   source("Script1_Preprocessing_DMFT_Four_Models.R")   # PROJECT = "DMFT"
#   source("Figures_DMFT.V3.R")
#   source("Calibration_DMFT.R")                          # this file
#
# It uses the refits Figures_DMFT.V3.R leaves in the session (FIG_REFITS,
# FIG_BEST, FILE_TAG, draw_roc). Run on its own, it reads the
# OOF_predictions_k<k>_<tag>.csv files that Figures_DMFT.V3.R writes, from the
# working directory. Only base R is needed (readxl for the baseline workbook).
#
# For the best panel at each size, and for the age + sex baseline when its
# workbook (Results_LogReg_HD_vs_LD_Age_Sex_Only_<key>.xlsx, sheet
# CV_Predictions) is in the working directory, from the pooled out-of-fold
# predictions of all repeats:
#   Calibration intercept  logistic regression of HD on an intercept with
#                          logit(p) as offset; 0 is ideal, > 0 means predicted
#                          risks are too low on average, < 0 too high
#   Calibration slope      logistic regression of HD on logit(p); 1 is ideal,
#                          < 1 means predictions are too extreme
#   Brier score            mean of (p - y)^2; 0 is perfect
#   Scaled Brier score     1 - Brier / [ybar (1 - ybar)]; 1 is perfect, 0 is no
#                          better than giving everyone the observed HD proportion
#   AUC, and balanced accuracy, sensitivity and specificity at the panel's
#   pooled Youden threshold (same rule as LogReg_Fast; checked against the
#   workbook), so that they get the same intervals.
# Intervals
#   95% CI   participant-level bootstrap: CAL_BOOT resamples of the participants
#            (within the HD and LD groups), each keeping all of its out-of-fold
#            predictions. Predictions are held fixed, so the interval omits the
#            uncertainty from panel, configuration and threshold selection.
#   Repeats  2.5th and 97.5th percentiles of the repeat-specific values (each
#            participant predicted once per repeat): fold-assignment variability,
#            as reported for balanced accuracy in the manuscript.
# Outputs (OUT_DIR)
#   Calibration_metrics_<tag>.csv          all models, all metrics, both intervals
#   Figure_Calibration_<tag>.(png|pdf|tiff) calibration plot of the selected model
#   Figure_1_ROC_Calibration_<tag>.(...)   ROC curves (A) and calibration (B),
#                                          when draw_roc() is in the session
#   Calibration_text_<tag>.txt             ready-to-paste manuscript sentences
#                                          and the package versions of this session
# ==========================================================================


# ==========================================================================
# 1) SETTINGS
# ==========================================================================

CAL_SEED  <- 2026
CAL_BOOT  <- 2000      # bootstrap resamples of participants (about 30 s per model)
CAL_SPAN  <- 0.75      # loess span of the calibration curves
CAL_BINS  <- 10        # groups (by predicted risk) for the observed-proportion points
CAL_EPS   <- 1e-12     # predictions clipped to [CAL_EPS, 1 - CAL_EPS] before the logit
CAL_K     <- NULL      # panel size shown in the figures; NULL = best overall panel
CAL_HIST_BINS <- 40    # bins of the prediction histogram under the calibration plot

if (!exists("OUT_DIR"))     OUT_DIR     <- getwd()
if (!exists("FIG_FORMATS")) FIG_FORMATS <- c("png", "pdf")
if (!exists("PNG_DPI"))     PNG_DPI     <- 300
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


# ==========================================================================
# 2) THE OUT-OF-FOLD PREDICTIONS
# ==========================================================================
# One entry per model: label, panel size, data frame (id, rep, y, p) and the
# workbook's threshold when known.

.cal_case <- if (exists(".lev")) .lev[2] else "HD"
.cal_rep  <- function(resample) sub("^.*\\.(Rep[0-9]+)$", "\\1", resample)
.cal_frame <- function(rowIndex, resample, obs, p) {
  d <- data.frame(id = as.integer(rowIndex), rep = .cal_rep(as.character(resample)),
                  y = as.integer(as.character(obs) == .cal_case), p = as.numeric(p))
  if (anyNA(d$p) || anyNA(d$id)) stop("Missing predictions or row indices in the out-of-fold predictions.")
  if (length(unique(d$y)) != 2) stop("The out-of-fold predictions need both outcome groups.")
  d
}
.cal_name <- function(k) if (k == 0) "Age + sex only" else sprintf("%d biomarker%s + age + sex", k, if (k == 1) "" else "s")

CAL_MODELS <- list()
if (exists("FIG_REFITS") && length(FIG_REFITS)) {
  for (kk in names(FIG_REFITS)) {
    r <- FIG_REFITS[[kk]]
    CAL_MODELS[[kk]] <- list(k = as.integer(kk), label = .cal_name(as.integer(kk)),
                             data = .cal_frame(r$pred$rowIndex, r$pred$Resample, r$pred$obs, r$pred$p),
                             threshold = if (exists("FIG_BEST")) FIG_BEST[[kk]]$row$Threshold else NA_real_)
  }
  CAL_TAG <- if (exists("FILE_TAG")) FILE_TAG else "DMFT"
} else {
  .files <- list.files(OUT_DIR, pattern = "^OOF_predictions_k([0-9]+)_(.+)\\.csv$")
  if (!length(.files))
    stop("No refits in this session and no OOF_predictions_k<k>_<tag>.csv files in ", OUT_DIR,
         ". Run Figures_DMFT.V3.R first.", call. = FALSE)
  .tags <- unique(sub("^OOF_predictions_k[0-9]+_(.+)\\.csv$", "\\1", .files))
  if (length(.tags) > 1)
    stop("Out-of-fold prediction files from more than one run: ", paste(.tags, collapse = ", "),
         ". Keep one run's files in the folder.", call. = FALSE)
  CAL_TAG <- .tags
  for (f in .files) {
    kk <- sub("^OOF_predictions_k([0-9]+)_.*$", "\\1", f)
    d  <- utils::read.csv(file.path(OUT_DIR, f), stringsAsFactors = FALSE)
    CAL_MODELS[[kk]] <- list(k = as.integer(kk), label = .cal_name(as.integer(kk)),
                             data = .cal_frame(d$rowIndex, d$Resample, d$obs, d$p), threshold = NA_real_)
  }
  CAL_MODELS <- CAL_MODELS[order(as.integer(names(CAL_MODELS)))]
}

# Age + sex baseline from its workbook, if present
.base_file <- list.files(OUT_DIR, pattern = "^Results_LogReg_HD_vs_LD_Age_Sex_Only_.*\\.xlsx$", full.names = TRUE)
if (length(.base_file) == 1 && requireNamespace("readxl", quietly = TRUE)) {
  .bp <- tryCatch(as.data.frame(readxl::read_excel(.base_file, sheet = "CV_Predictions")), error = function(e) NULL)
  if (!is.null(.bp) && all(c("rowIndex", "Resample", "obs", .cal_case) %in% names(.bp))) {
    .bt <- tryCatch(as.data.frame(readxl::read_excel(.base_file, sheet = "Results_0_Biomarkers")), error = function(e) NULL)
    CAL_MODELS <- c(list(`0` = list(k = 0L, label = .cal_name(0),
                                    data = .cal_frame(.bp$rowIndex, .bp$Resample, .bp$obs, .bp[[.cal_case]]),
                                    threshold = if (!is.null(.bt)) as.numeric(.bt$Threshold[1]) else NA_real_)),
                    CAL_MODELS)
  }
} else if (length(.base_file) > 1) {
  message("Several age + sex baseline workbooks found; the baseline is left out.")
}

.n_ids <- vapply(CAL_MODELS, function(m) length(unique(m$data$id)), integer(1))
if (length(unique(.n_ids)) != 1)
  stop("The models do not cover the same participants: ", paste(.n_ids, collapse = ", "), call. = FALSE)


# ==========================================================================
# 3) METRICS
# ==========================================================================

# Youden threshold on pooled predictions as in LogReg_Fast: cut-offs at the
# midpoints between distinct predictions (case when p >= cut); among tied
# cut-offs, the one closest to 0.5 (equally close: the higher)
cal_youden <- function(p, y) {
  np <- sum(y == 1); nn <- sum(y == 0)
  o <- order(p, decreasing = TRUE); ps <- p[o]; ys <- y[o]
  keep <- c(ps[-1] != ps[-length(ps)], TRUE)
  tp <- c(0, cumsum(ys)[keep]); fp <- c(0, cumsum(1 - ys)[keep]); s <- ps[keep]; nk <- length(s)
  thr <- c(Inf, if (nk > 1) (s[-nk] + s[-1]) / 2, -Inf)
  J <- tp / np + 1 - fp / nn
  tied <- which(J >= max(J) - 0.5 / (np * nn))
  thr[tied[order(abs(thr[tied] - 0.5), -thr[tied])[1]]]
}

cal_auc <- function(p, y) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Logistic fits without glm's overhead (used thousands of times in the bootstrap)
.cal_logit_fit <- function(X, y, offset = NULL) {
  f <- suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(),
                                       offset = if (is.null(offset)) rep(0, length(y)) else offset,
                                       control = list(maxit = 50)))
  stats::setNames(f$coefficients, NULL)
}

cal_metrics <- function(p, y, thr) {
  pc <- pmin(pmax(p, CAL_EPS), 1 - CAL_EPS); lp <- stats::qlogis(pc)
  ones <- matrix(1, length(y), 1)
  a <- .cal_logit_fit(ones, y, offset = lp)
  b <- .cal_logit_fit(cbind(1, lp), y)[2]
  brier <- mean((p - y)^2); ybar <- mean(y)
  sens <- mean(p[y == 1] >= thr); spec <- mean(p[y == 0] < thr)
  c(Calibration_intercept = a, Calibration_slope = b, Brier = brier,
    Scaled_Brier = 1 - brier / (ybar * (1 - ybar)), AUC = cal_auc(p, y),
    Balanced_accuracy = (sens + spec) / 2, Sensitivity = sens, Specificity = spec,
    Mean_predicted = mean(p), Observed = ybar)
}

cal_evaluate <- function(m, B = CAL_BOOT, seed = CAL_SEED) {
  d <- m$data
  thr <- cal_youden(d$p, d$y)
  if (is.finite(m$threshold) && abs(thr - m$threshold) > 1e-6)
    warning(sprintf("%s: Youden threshold %.6f differs from the workbook's %.6f; the workbook's is used.",
                    m$label, thr, m$threshold))
  if (is.finite(m$threshold)) thr <- m$threshold
  est <- cal_metrics(d$p, d$y, thr)

  # repeat-specific values
  reps <- split(seq_len(nrow(d)), d$rep)
  rep_vals <- t(vapply(reps, function(ix) cal_metrics(d$p[ix], d$y[ix], thr), numeric(length(est))))
  rep_q <- apply(rep_vals, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE)

  # participant-level bootstrap, predictions held fixed
  by_id <- split(seq_len(nrow(d)), d$id); ids <- names(by_id)
  y_id <- vapply(by_id, function(ix) d$y[ix[1]], integer(1))
  set.seed(seed)
  boot <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {
    # resample within outcome groups so every resample has both groups
    s <- c(sample(ids[y_id == 1], replace = TRUE), sample(ids[y_id == 0], replace = TRUE))
    ix <- unlist(by_id[s], use.names = FALSE)
    boot[b, ] <- cal_metrics(d$p[ix], d$y[ix], thr)
  }
  boot_q <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)

  extreme <- mean(d$p < 0.001 | d$p > 0.999)
  data.frame(Model = m$label, Panel_size = m$k, Participants = length(by_id),
             HD = sum(y_id == 1), LD = sum(y_id == 0), Predictions = nrow(d), Repeats = length(reps),
             Threshold = thr, Metric = names(est), Estimate = unname(est),
             Bootstrap_CI_Lower = boot_q[1, ], Bootstrap_CI_Upper = boot_q[2, ],
             Repeat_P2.5 = rep_q[1, ], Repeat_P97.5 = rep_q[2, ],
             Share_predictions_below_0.001_or_above_0.999 = extreme,
             Clipped_predictions = sum(d$p < CAL_EPS | d$p > 1 - CAL_EPS),
             Bootstrap_resamples = B, stringsAsFactors = FALSE, row.names = NULL)
}

cat(sprintf("\n===== Calibration: %d model(s), %d participants, %d bootstrap resamples =====\n",
            length(CAL_MODELS), .n_ids[1], CAL_BOOT))
CAL_RESULTS <- do.call(rbind, lapply(CAL_MODELS, function(m) {
  t0 <- Sys.time(); r <- cal_evaluate(m)
  g <- function(x) r$Estimate[r$Metric == x]
  cat(sprintf("  %-28s intercept %6.2f  slope %5.2f  Brier %.3f  AUC %.3f  (%.0fs)\n", m$label,
              g("Calibration_intercept"), g("Calibration_slope"), g("Brier"), g("AUC"),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  r
}))
rownames(CAL_RESULTS) <- NULL
utils::write.csv(CAL_RESULTS, file.path(OUT_DIR, sprintf("Calibration_metrics_%s.csv", CAL_TAG)),
                 row.names = FALSE, fileEncoding = "UTF-8")


# ==========================================================================
# 4) CALIBRATION PLOT
# ==========================================================================

.ks_avail <- vapply(CAL_MODELS, `[[`, integer(1), "k")
.k_fig <- if (!is.null(CAL_K)) CAL_K else if (exists("BEST_OVERALL_K")) BEST_OVERALL_K else max(.ks_avail)
if (!.k_fig %in% .ks_avail) stop("CAL_K = ", .k_fig, " is not among the models: ", paste(.ks_avail, collapse = ", "))
CAL_FIG_MODEL <- CAL_MODELS[[as.character(.k_fig)]]

.cal_loess <- function(p, y, grid) {
  fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = CAL_SPAN, degree = 2,
                                                control = stats::loess.control(surface = "direct"))),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, length(grid)))
  pmin(pmax(suppressWarnings(stats::predict(fit, newdata = data.frame(p = grid))), 0), 1)
}

prepare_calibration <- function(m) {
  d <- m$data
  rng  <- range(d$p)
  grid <- seq(rng[1], rng[2], length.out = 200)
  curve <- .cal_loess(d$p, d$y, grid)
  rep_curves <- vapply(split(seq_len(nrow(d)), d$rep), function(ix) {
    g <- grid; inr <- g >= min(d$p[ix]) & g <= max(d$p[ix])
    out <- rep(NA_real_, length(g)); out[inr] <- .cal_loess(d$p[ix], d$y[ix], g[inr]); out
  }, numeric(length(grid)))
  band <- t(apply(rep_curves, 1, function(v) if (sum(is.finite(v)) >= 10)
    stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)))
  br  <- unique(stats::quantile(d$p, seq(0, 1, length.out = CAL_BINS + 1), names = FALSE))
  grp <- cut(d$p, br, include.lowest = TRUE)
  pts <- data.frame(pred = tapply(d$p, grp, mean), obs = tapply(d$y, grp, mean))
  res <- CAL_RESULTS[CAL_RESULTS$Model == m$label, ]
  list(d = d, grid = grid, curve = curve, band = band, pts = pts[stats::complete.cases(pts), ],
       res = res, label = m$label)
}

.fmt_ci <- function(res, metric, digits) {
  r <- res[res$Metric == metric, ]
  f <- function(x) formatC(x, format = "f", digits = digits)
  sprintf("%s (%s to %s)", f(r$Estimate), f(r$Bootstrap_CI_Lower), f(r$Bootstrap_CI_Upper))
}

draw_calibration <- function(cal, main = NULL, hist = TRUE) {
  d <- cal$d
  op <- par(mar = c(if (hist) 0.6 else 4.6, 4.6, 3.2, 1.6), pty = if (hist) "m" else "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = if (hist) "n" else "s",
       xlab = if (hist) "" else "Predicted probability of HD", ylab = "Observed proportion HD",
       main = main %||% "", cex.main = 0.98)
  abline(h = seq(0.2, 0.8, 0.2), v = seq(0.2, 0.8, 0.2), col = "grey92", lwd = 0.8)
  abline(0, 1, col = "grey55", lty = 2, lwd = 1.5)
  ok <- is.finite(cal$band[, 1]) & is.finite(cal$band[, 2])
  if (any(ok)) polygon(c(cal$grid[ok], rev(cal$grid[ok])), c(cal$band[ok, 1], rev(cal$band[ok, 2])),
                       col = grDevices::adjustcolor("#0072B2", 0.18), border = NA)
  lines(cal$grid, cal$curve, col = "#0072B2", lwd = 2.5)
  points(cal$pts$pred, cal$pts$obs, pch = 21, bg = "white", col = "#0b0b0b", cex = 1.2, lwd = 1.4)
  txt <- c(sprintf("Calibration intercept %s", .fmt_ci(cal$res, "Calibration_intercept", 2)),
           sprintf("Calibration slope %s", .fmt_ci(cal$res, "Calibration_slope", 2)),
           sprintf("Brier score %s", .fmt_ci(cal$res, "Brier", 3)))
  legend("bottomright", legend = txt, bty = "n", cex = 0.78, inset = c(0.01, 0.02),
         title = "Estimate (95% bootstrap CI)", title.adj = 0)
  legend("topleft", bty = "n", cex = 0.78, inset = c(0.01, 0.01),
         legend = c("Ideal", "Loess calibration curve", "Middle 95% across repeats", sprintf("Observed in %d risk groups", nrow(cal$pts))),
         lty = c(2, 1, NA, NA), lwd = c(1.5, 2.5, NA, 1.4), pch = c(NA, NA, 15, 21),
         col = c("grey55", "#0072B2", grDevices::adjustcolor("#0072B2", 0.25), "#0b0b0b"),
         pt.bg = c(NA, NA, NA, "white"), pt.cex = c(1, 1, 2, 1.2))
  if (hist) {
    par(mar = c(4.6, 4.6, 0.3, 1.6))
    brk <- seq(0, 1, length.out = CAL_HIST_BINS + 1)
    h1 <- graphics::hist(d$p[d$y == 1], breaks = brk, plot = FALSE)$counts / length(unique(d$rep))
    h0 <- graphics::hist(d$p[d$y == 0], breaks = brk, plot = FALSE)$counts / length(unique(d$rep))
    ym <- max(c(h1, h0), 1)
    plot(NA, xlim = c(0, 1), ylim = c(-ym, ym), xaxs = "i", yaxt = "n", xlab = "Predicted probability of HD",
         ylab = "")
    rect(brk[-length(brk)], 0, brk[-1], h1, col = "#D55E00", border = "white", lwd = 0.4)
    rect(brk[-length(brk)], -h0, brk[-1], 0, col = "#56B4E9", border = "white", lwd = 0.4)
    abline(h = 0, col = "grey40", lwd = 0.8)
    mtext("HD", side = 2, line = 0.6, at = ym / 2, las = 1, cex = 0.7, col = "#D55E00")
    mtext("LD", side = 2, line = 0.6, at = -ym / 2, las = 1, cex = 0.7, col = "#0072B2")
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

.cal_save <- function(stem, width, height, draw) {
  for (fmt in FIG_FORMATS) {
    f <- file.path(OUT_DIR, paste0(stem, ".", fmt))
    switch(fmt,
           png  = grDevices::png(f, width = width, height = height, units = "in", res = PNG_DPI, pointsize = 12),
           pdf  = grDevices::pdf(f, width = width, height = height, pointsize = 12, useDingbats = FALSE),
           tiff = grDevices::tiff(f, width = width, height = height, units = "in", res = PNG_DPI, pointsize = 12,
                                  compression = "lzw"),
           stop("FIG_FORMATS: use png, pdf or tiff."))
    tryCatch(draw(), finally = grDevices::dev.off())
    cat("  Saved:", f, "\n")
  }
}

CAL_FIG <- prepare_calibration(CAL_FIG_MODEL)
.title <- sprintf("Calibration: %s", CAL_FIG_MODEL$label)
.cal_save(sprintf("Figure_Calibration_%s", CAL_TAG), 6.5, 7.6, function() {
  layout(matrix(1:2, 2), heights = c(4.6, 1.4))
  draw_calibration(CAL_FIG, main = .title)
})

# ROC curves (A) beside calibration (B), as one manuscript figure
if (exists("draw_roc") && is.function(draw_roc) && exists("FIG_REFITS")) {
  .cal_save(sprintf("Figure_1_ROC_Calibration_%s", CAL_TAG), 12.5, 6.9, function() {
    layout(matrix(c(1, 1, 2, 3), 2, 2), widths = c(1, 1), heights = c(4.6, 1.4))
    par(oma = c(0, 0, 1.6, 0))
    draw_roc()
    draw_calibration(CAL_FIG, main = .title)
    mtext("A", outer = TRUE, side = 3, at = 0.012, adj = 0, font = 2, cex = 1.5, line = 0)
    mtext("B", outer = TRUE, side = 3, at = 0.512, adj = 0, font = 2, cex = 1.5, line = 0)
  })
}


# ==========================================================================
# 5) READY-TO-PASTE TEXT
# ==========================================================================

.r  <- CAL_RESULTS[CAL_RESULTS$Model == CAL_FIG_MODEL$label, ]
.g  <- function(m, col = "Estimate") .r[[col]][.r$Metric == m]
.f2 <- function(x) formatC(x, format = "f", digits = 2)
.f3 <- function(x) formatC(x, format = "f", digits = 3)
.ci <- function(m, f) sprintf("%s (95%% CI %s to %s)", f(.g(m)), f(.g(m, "Bootstrap_CI_Lower")), f(.g(m, "Bootstrap_CI_Upper")))
.rep <- function(m, f) sprintf("%s to %s", f(.g(m, "Repeat_P2.5")), f(.g(m, "Repeat_P97.5")))

.slope_lo <- .g("Calibration_slope", "Bootstrap_CI_Lower"); .slope_hi <- .g("Calibration_slope", "Bootstrap_CI_Upper")
.int_lo <- .g("Calibration_intercept", "Bootstrap_CI_Lower"); .int_hi <- .g("Calibration_intercept", "Bootstrap_CI_Upper")
.interp_slope <- if (.slope_hi < 1) {
  "The calibration slope below 1 indicates predictions that were too extreme, with high risks overestimated and low risks underestimated, consistent with overfitting."
} else if (.slope_lo > 1) {
  "The calibration slope above 1 indicates predictions that were too moderate, with high risks underestimated and low risks overestimated."
} else "The calibration slope was compatible with 1."
.interp_int <- if (.int_lo > 0) {
  "The positive calibration intercept indicates that predicted risks were, on average, lower than the observed HD proportion."
} else if (.int_hi < 0) {
  "The negative calibration intercept indicates that predicted risks were, on average, higher than the observed HD proportion."
} else "The calibration intercept was compatible with 0, indicating no systematic over- or underestimation."
.extreme <- .r$Share_predictions_below_0.001_or_above_0.999[1]

.base <- CAL_RESULTS[CAL_RESULTS$Panel_size == 0, ]
.base_txt <- if (nrow(.base)) {
  b <- function(m) .base$Estimate[.base$Metric == m]
  sprintf(" For comparison, the age-and-sex baseline had a calibration intercept of %s, slope of %s, and Brier score of %s (AUC %s).",
          .f2(b("Calibration_intercept")), .f2(b("Calibration_slope")), .f3(b("Brier")), .f3(b("AUC")))
} else ""

.versions <- vapply(c("caret", "pROC", "VIM", "doSNOW", "foreach", "fastglm", "readxl", "writexl", "shiny"),
                    function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_,
                    character(1))

.text <- c(
  "RESULTS (Section 3.2, after the ROC paragraph):",
  paste0("Figure 1B shows the calibration of the selected model. The calibration intercept was ",
         .ci("Calibration_intercept", .f2), ", the calibration slope ", .ci("Calibration_slope", .f2),
         ", and the Brier score ", .ci("Brier", .f3), "; the scaled Brier score was ", .f2(.g("Scaled_Brier")),
         ". Across repeats, the slope ranged from ", .rep("Calibration_slope", .f2), " and the Brier score from ",
         .rep("Brier", .f3), " (2.5th to 97.5th percentiles). ", .interp_slope, " ", .interp_int,
         sprintf(" Of the pooled predictions, %s%% were below 0.001 or above 0.999.", formatC(100 * .extreme, format = "f", digits = 1)),
         .base_txt),
  "",
  "ABSTRACT (Results, after the AUC sentence):",
  paste0("Calibration intercept and slope were ", .f2(.g("Calibration_intercept")), " and ", .f2(.g("Calibration_slope")),
         ", and the Brier score was ", .f3(.g("Brier")), "."),
  "",
  "BOOTSTRAP INTERVALS FOR DISCRIMINATION (optional, alongside the repeat percentiles):",
  paste0("AUC ", .ci("AUC", .f3), "; balanced accuracy ", .ci("Balanced_accuracy", .f3),
         " at the ", .f3(.r$Threshold[1]), " threshold."),
  "",
  sprintf("Participants %d (HD %d, LD %d); pooled predictions %d from %d repeats; %d bootstrap resamples; loess span %.2f; predictions clipped to [%g, 1 - %g] before the logit (%d clipped).",
          .r$Participants[1], .r$HD[1], .r$LD[1], .r$Predictions[1], .r$Repeats[1], CAL_BOOT, CAL_SPAN,
          CAL_EPS, CAL_EPS, .r$Clipped_predictions[1]),
  "",
  "PACKAGE VERSIONS IN THIS R SESSION (for Section 2.4; check they match the analysis run):",
  paste0(R.version.string, "; ", paste(sprintf("%s %s", names(.versions)[!is.na(.versions)], .versions[!is.na(.versions)]),
                                       collapse = ", "), ".")
)
writeLines(.text, file.path(OUT_DIR, sprintf("Calibration_text_%s.txt", CAL_TAG)), useBytes = TRUE)
cat("\n", paste(.text, collapse = "\n"), "\n", sep = "")
cat(sprintf("\nSaved: Calibration_metrics_%s.csv and Calibration_text_%s.txt\n", CAL_TAG, CAL_TAG))
