# ==========================================================================
# FIGURES FOR THE DMFT MANUSCRIPT: ROC CURVES, RADAR PLOT, HEATMAP
# ==========================================================================
# Run after Script 1 (PROJECT = "DMFT") in the same session, with the
# LogReg_Fast workbook of that run in the working directory:
#
#   source("script1_preprocessing.R")   # PROJECT = "DMFT"
#   source("Figures_DMFT.R")            # this file; runs at the bottom
#
# 1. ROC curves of the best panel at each size (1 to 6 markers) from the
#    pooled out-of-fold predictions of all repeats. The panels are refitted
#    on LogReg_Fast's folds (caret's own folds after set.seed(123); folds,
#    repeats and seed read from the workbook), so the AUCs, thresholds and
#    balanced accuracies reproduce the workbook. The script checks this.
#    No covariates-only reference curve and no fixed-specificity marker.
# 2. Radar plot of the best panel at each size: balanced accuracy,
#    sensitivity, specificity, PPV, NPV and AUC (workbook values).
# 3. Heatmap of the best overall panel (highest balanced accuracy across
#    sizes, ties by AUC): one row per participant, LD above the line and HD
#    below; each column scaled 0-1; DMFT and the mean out-of-fold P(HD) per
#    row; participants misclassified at the cutoff highlighted.
# Blue/orange heatmap; black outlines identify misclassified cells.
#
# Also written (CSV): the refit check, the
# participants' mean P(HD), and the coefficients of the best overall panel
# with the LD reference means and SDs (for the probability equation).
# ==========================================================================


# ==========================================================================
# 1) SETTINGS
# ==========================================================================

RESULTS_FILE <- NULL           # NULL: automatically compare the three biomarker scenario workbooks for this registry key
OUT_DIR      <- getwd()        # where figures and tables are written
FIG_FORMATS  <- c("png", "pdf")   # any of "png", "pdf", "tiff"
PNG_DPI      <- 300

MAKE_ROC     <- TRUE
MAKE_PERFORMANCE <- TRUE  # manuscript Figure 1: balanced accuracy ranges and AUC
MAKE_RADAR   <- TRUE
MAKE_HEATMAP <- TRUE

# ROC
ROC_COVARIATES_ONLY <- FALSE   # intentionally omitted
ROC_MARK_YOUDEN     <- TRUE    # mark each curve's Youden cut; the symbols also tell the
                               # curves apart without relying on colour

# Radar
RADAR_MIN  <- NULL             # value at the centre; NULL = lowest plotted value rounded down to 0.1
RADAR_AREA <- "full"           # area of each panel's polygon (legend and CSV), as a fraction of the
                               #   polygon with every metric at 1. "full": each metric on a 0-1 radius,
                               #   so RADAR_MIN does not change it; "plotted": the polygon as drawn

# Heatmap of the best overall panel
HEATMAP_CUTOFF     <- "pooled" # "pooled": the panel's Youden cut on the pooled predictions (the
                               #   workbook's Threshold); "participant": Youden cut on the
                               #   participants' mean P(HD); or a number, e.g. 0.5
HEATMAP_ORDER      <- "dmft"   # rows within each group: "dmft" (DMFT, then P(HD)) or "probability"
HEATMAP_COVARIATES <- TRUE     # covariate columns after the markers
HEATMAP_PALETTE    <- "blue_orange"     # two endpoint colours: blue (low) and orange (high)
HEATMAP_ROW_IN     <- 0.06     # row height in inches; the figure height follows (111 rows: ~9.3 in)

# Marker names on the figures
ITALIC_TAXA    <- TRUE         # genus and species names in italics
NAME_OVERRIDES <- c()          # e.g. c(Moraxella_Moraxella_unclassified = "*Moraxella* sp. (unclassified)");
                               # *...* marks italics

# Colour and symbol per panel size (Okabe-Ito, as in the caries ROC figure)
SIZE_COLS <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
SIZE_PCH  <- c(16, 17, 15, 18, 4, 8)


# ==========================================================================
# 2) PACKAGES AND PREREQUISITES
# ==========================================================================

.need_pkgs <- c("caret", "pROC", "readxl")
.miss_pkgs <- .need_pkgs[!vapply(.need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkgs) > 0)
  stop("Missing package(s): ", paste(.miss_pkgs, collapse = ", "),
       "\n  install.packages(c(", paste0('"', .miss_pkgs, '"', collapse = ", "), "))", call. = FALSE)
suppressMessages({ library(caret); library(pROC) })

.msgs <- character(0)
if (!exists("PROJECT") || !exists("registry_key") || !exists("data.dfPH")) {
  .msgs <- "Script 1 has not run in this session: source(\"script1_preprocessing.R\") with PROJECT = \"DMFT\" first."
} else {
  if (!identical(PROJECT, "DMFT"))
    .msgs <- c(.msgs, sprintf("PROJECT = \"%s\"; these figures are for PROJECT = \"DMFT\".", PROJECT))
  .grp <- sort(unique(as.character(data.dfPH$condL)))
  if (!identical(.grp, c("HD", "LD")))
    .msgs <- c(.msgs, sprintf("condL in data.dfPH should be LD and HD; found: %s.", paste(.grp, collapse = ", ")))
  .scenario_names <- c("Biomarkers_Only", "Biomarkers_Age_Sex", "Biomarkers_Age_Sex_BMI")
  .scenario_files <- sprintf("Results_LogReg_HD_vs_LD_%s_%s.xlsx", .scenario_names, registry_key)
  .scenario_best <- list()
  if (is.null(RESULTS_FILE)) {
    if (!all(file.exists(.scenario_files)))
      stop("Put all three biomarker scenario workbooks in the working directory, or set RESULTS_FILE explicitly.")
    for (ff in .scenario_files) {
      for (sh in grep("^Results_[0-9]+_Biomarkers$", readxl::excel_sheets(ff), value = TRUE)) {
        dd <- as.data.frame(readxl::read_excel(ff, sheet = sh))
        dd <- dd[order(-dd$Balanced_Accuracy, -dd$AUC, na.last = TRUE), , drop = FALSE]
        if (nrow(dd)) .scenario_best[[length(.scenario_best) + 1L]] <- data.frame(
          file = ff, k = as.integer(gsub("[^0-9]", "", sh)),
          BA = dd$Balanced_Accuracy[1], AUC = dd$AUC[1])
      }
    }
    .comparison <- do.call(rbind, .scenario_best)
    .comparison <- .comparison[order(-.comparison$BA, -.comparison$AUC), ]
    .res_file <- .comparison$file[1]
    .winners <- .comparison[!duplicated(.comparison$k), ]
    if (any(.winners$file != .res_file))
      stop("Different configurations win at different panel sizes. Set RESULTS_FILE for each configuration and run separately.")
    message("Best configuration across all panel sizes: ", basename(.res_file))
  } else .res_file <- RESULTS_FILE
  if (!file.exists(.res_file))
    .msgs <- c(.msgs, sprintf("Results workbook not found: %s\n  Put it in the working directory (%s) or set RESULTS_FILE.",
                              .res_file, getwd()))
}
if (length(.msgs) > 0)
  stop("\n\n===== PREREQUISITES NOT MET =====\n  ", paste(.msgs, collapse = "\n  "), "\n", call. = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!exists("confounders")) confounders <- NULL


# ==========================================================================
# 3) THE WORKBOOK
# ==========================================================================

.params <- as.data.frame(readxl::read_excel(.res_file, sheet = "Analysis_Parameters", col_types = "text"))
.pget <- function(name) {
  v <- .params$Value[which(trimws(.params$Parameter) == name)]   # which(): blank rows are NA
  if (length(v) && !is.na(v[1])) v[1] else NA_character_
}

.cv <- regmatches(.pget("Cross-validation"),
                  regexec("([0-9]+)-fold CV x ([0-9]+) repeats, seed ([0-9]+)", .pget("Cross-validation")))[[1]]
if (length(.cv) != 4)
  stop("Could not read the cross-validation scheme from the workbook's Analysis_Parameters sheet.", call. = FALSE)
FIG_CV_FOLDS   <- as.integer(.cv[2])
FIG_CV_REPEATS <- as.integer(.cv[3])
FIG_CV_SEED    <- as.integer(.cv[4])

.wb_conf <- .pget("Confounders")
.wb_conf <- if (is.na(.wb_conf) || .wb_conf == "none") character(0) else trimws(strsplit(.wb_conf, ",")[[1]])
.ref  <- .pget("Reference group")
.case <- .pget("Case group (positive)")
if (is.na(.ref) || is.na(.case)) { .ref <- "LD"; .case <- "HD" }

.sheets <- grep("^Results_[0-9]+_Biomarkers$", readxl::excel_sheets(.res_file), value = TRUE)
.ks     <- sort(as.integer(sub("^Results_([0-9]+)_Biomarkers$", "\\1", .sheets)))
if (!length(.ks)) stop("No Results_<k>_Biomarkers sheets in the workbook.", call. = FALSE)

# best panel per size: the top row, sorted again by balanced accuracy then AUC
FIG_BEST <- lapply(.ks, function(k) {
  d <- as.data.frame(readxl::read_excel(.res_file, sheet = sprintf("Results_%d_Biomarkers", k)))
  d <- d[order(-d$Balanced_Accuracy, -d$AUC, na.last = TRUE), , drop = FALSE]
  list(k = k, markers = as.character(unlist(d[1, paste0("Biomarker", seq_len(k))])), row = d[1, , drop = FALSE])
})
names(FIG_BEST) <- as.character(.ks)
.ba  <- vapply(FIG_BEST, function(b) b$row$Balanced_Accuracy, numeric(1))
.auc <- vapply(FIG_BEST, function(b) b$row$AUC, numeric(1))
BEST_OVERALL_K <- .ks[order(-.ba, -.auc)[1]]


# ==========================================================================
# 4) DATA AND FOLDS (as in LogReg_Fast)
# ==========================================================================
# Rows with a missing confounder are left out before the folds are built,
# condL is a factor with the reference first, and the folds are read off a
# throwaway caret::train() after set.seed(seed): caret draws a private seed
# before building its resamples, so createMultiFolds() would not match.

.conf <- .wb_conf  # model covariates come from the selected workbook, not Script 1
.stopm <- character(0)
if (!identical(.pget("Script 1 registry key"), registry_key))
  .stopm <- c(.stopm, sprintf("registry key in this session (%s) differs from the workbook's (%s).",
                              registry_key, .pget("Script 1 registry key")))

FIG_DAT <- as.data.frame(data.dfPH)
# Match the common cohort used for ALL four models, even when BMI is not a predictor.
.required <- c("Age_DOS", "Sex", "BMI")
if (length(setdiff(.required, names(FIG_DAT)))) stop("Missing common-cohort columns: Age_DOS, Sex, BMI.")
.sx <- trimws(as.character(FIG_DAT$Sex)); .sx[is.na(FIG_DAT$Sex) | .sx == ""] <- NA_character_
FIG_DAT$Sex <- factor(.sx)
.ok <- complete.cases(FIG_DAT[, .required, drop = FALSE]) & !is.na(FIG_DAT$condL)
.source_rows <- which(.ok)
FIG_DAT <- droplevels(FIG_DAT[.ok, , drop = FALSE])
if (nlevels(FIG_DAT$Sex) != 2L) stop("Expected two observed Sex levels.")
# Use the exact saved folds when available; otherwise reproduce caret's seed procedure.
.fold_file <- sprintf("DMFT_Four_Models_Resampling_%s.rds", registry_key)
.saved_folds <- if (file.exists(.fold_file)) readRDS(.fold_file) else NULL
if (!is.null(.saved_folds) && (!identical(.saved_folds$source_rows, .source_rows) ||
    !identical(.saved_folds$sex_levels, levels(FIG_DAT$Sex)) ||
    !identical(.saved_folds$registry_key, registry_key) ||
    !identical(as.character(.saved_folds$outcome), as.character(FIG_DAT$condL))))
  stop("Saved resampling cohort/coding does not match this preprocessing session.")
FIG_DAT$condL <- factor(as.character(FIG_DAT$condL), levels = c(.ref, .case))
levels(FIG_DAT$condL) <- make.names(levels(FIG_DAT$condL))
.lev <- levels(FIG_DAT$condL)

if (!identical(as.integer(.pget("Total n")), nrow(FIG_DAT)))
  .stopm <- c(.stopm, sprintf("%d rows in data.dfPH after dropping missing confounders; the workbook has %s.",
                              nrow(FIG_DAT), .pget("Total n")))
.all_markers <- unique(unlist(lapply(FIG_BEST, `[[`, "markers")))
if (length(setdiff(.all_markers, colnames(FIG_DAT))))
  .stopm <- c(.stopm, paste("Marker(s) not in data.dfPH:", paste(setdiff(.all_markers, colnames(FIG_DAT)), collapse = ", ")))
if (length(.stopm))
  stop("\n\n===== SESSION DOES NOT MATCH THE WORKBOOK =====\n  ", paste(.stopm, collapse = "\n  "),
       "\n  Re-run Script 1 with the settings of the workbook's run.\n", call. = FALSE)

.custom_summary <- function(data, lev = NULL, model = NULL) {
  ts  <- twoClassSummary(data, lev, model)
  ds  <- defaultSummary(data, lev, model)
  out <- c(ts[["ROC"]], ts[["Spec"]], ts[["Sens"]], ds[["Accuracy"]], ds[["Kappa"]],
           mean(c(ts[["Spec"]], ts[["Sens"]])))
  names(out) <- c("ROC", "Sens", "Spec", "Accuracy", "Kappa", "BalancedAccuracy")
  out
}
.control <- trainControl(method = "repeatedcv", number = FIG_CV_FOLDS, repeats = FIG_CV_REPEATS,
                         summaryFunction = .custom_summary, classProbs = TRUE,
                         allowParallel = FALSE, savePredictions = "final")

.fold_index <- function(y) {
  d <- data.frame(condL = y, z = sin(seq_along(y)))
  set.seed(FIG_CV_SEED)
  fit <- suppressWarnings(caret::train(condL ~ z, data = d, method = "glm", family = "binomial",
                                       metric = "BalancedAccuracy", trControl = .control))
  fit$control$index
}

# Youden cut as in LogReg_Fast: highest Sensitivity + Specificity; among tied
# cut-offs (within half a step), the one closest to 0.5
.youden_best <- function(roc_obj) {
  m <- coords(roc_obj, "all",
              ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv"),
              transpose = FALSE)
  J    <- m$sensitivity + m$specificity
  tol  <- 0.5 / (length(roc_obj$cases) * length(roc_obj$controls))
  tied <- which(J >= max(J, na.rm = TRUE) - tol)
  m[tied[order(abs(m$threshold[tied] - 0.5), -m$threshold[tied])[1]], , drop = FALSE]
}

# Refit one panel through caret on the fixed folds. Results are cached in the
# session (dmft_fig_cache), so re-running the script to restyle a figure is quick.
.cache <- list() # fresh cache prevents reuse across datasets or sex coding changes
.refit <- function(markers) {
  preds <- c(markers, .conf)
  key <- paste(registry_key, nrow(FIG_DAT), FIG_CV_FOLDS, FIG_CV_REPEATS, FIG_CV_SEED, paste(preds, collapse = "+"),
               format(sum(vapply(FIG_DAT[, preds, drop = FALSE], function(z) sum(as.numeric(z)), numeric(1))),
                      digits = 15))
  if (!is.null(.cache[[key]])) return(.cache[[key]])
  if (is.null(.cache[[".folds"]]) || !identical(.cache[[".folds_key"]], paste(registry_key, nrow(FIG_DAT), FIG_CV_FOLDS, FIG_CV_REPEATS, FIG_CV_SEED))) {
    .cache[[".folds"]]     <<- if (!is.null(.saved_folds)) .saved_folds$index else .fold_index(FIG_DAT$condL)
    .cache[[".folds_key"]] <<- paste(registry_key, nrow(FIG_DAT), FIG_CV_FOLDS, FIG_CV_REPEATS, FIG_CV_SEED)
  }
  ctrl <- .control
  ctrl$index <- .cache[[".folds"]]
  sub <- FIG_DAT[, c(preds, "condL"), drop = FALSE]
  set.seed(FIG_CV_SEED)
  fit <- suppressWarnings(caret::train(as.formula(paste("condL ~", paste(preds, collapse = " + "))),
                                       data = sub, method = "glm", family = "binomial",
                                       metric = "BalancedAccuracy", trControl = ctrl))
  pred    <- fit$pred
  roc_obj <- pROC::roc(pred$obs, pred[[.lev[2]]], levels = .lev, direction = "<", quiet = TRUE)
  yb      <- .youden_best(roc_obj)
  out <- list(markers = markers, predictors = preds,
              pred = data.frame(rowIndex = pred$rowIndex, Resample = pred$Resample,
                                obs = pred$obs, p = pred[[.lev[2]]]),
              roc = roc_obj, youden = yb, auc = as.numeric(roc_obj$auc),
              ba = (yb$sensitivity + yb$specificity) / 2, final = fit$finalModel)
  .cache[[key]] <<- out
  out
}


# ==========================================================================
# 5) LABELS
# ==========================================================================
# Marker names become display labels, e.g.
#   Actinomyces_sp._HMT_169            -> *Actinomyces* sp. HMT 169
#   Moraxella_Moraxella_unclassified   -> unclassified *Moraxella*
#   Lachnospiraceae_G7_bacterium_HMT_163 -> Lachnospiraceae [G-7] bacterium HMT 163
# *...* marks italics (genus and species); family and higher ranks stay upright.

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

# words with their italic flag and whether a space precedes them
.md_words <- function(s) {
  chunks <- strsplit(s, "*", fixed = TRUE)[[1]]
  words <- list(); it <- FALSE; sp_next <- FALSE
  for (ch in chunks) {
    if (nzchar(ch)) {
      parts <- strsplit(ch, " ", fixed = TRUE)[[1]]
      first <- TRUE
      for (w in parts) {
        if (!nzchar(w)) { sp_next <- TRUE; next }
        sp <- if (first) sp_next else TRUE
        if (!length(words)) sp <- FALSE
        words[[length(words) + 1]] <- list(w = w, it = it, sp = sp)
        sp_next <- FALSE; first <- FALSE
      }
      if (endsWith(ch, " ")) sp_next <- TRUE
    }
    it <- !it
  }
  words
}

# plotmath: Greek letters as symbols, italic runs in italic()
.greek_expr <- function(txt) {
  g <- c("\u03b1" = "alpha", "\u03b2" = "beta", "\u03b3" = "gamma", "\u03ba" = "kappa")
  chars <- strsplit(txt, "")[[1]]
  if (!any(chars %in% names(g))) return(txt)
  pieces <- list(); buf <- ""
  for (ch in chars) {
    if (ch %in% names(g)) {
      if (nzchar(buf)) pieces[[length(pieces) + 1]] <- buf
      pieces[[length(pieces) + 1]] <- as.name(g[[ch]]); buf <- ""
    } else buf <- paste0(buf, ch)
  }
  if (nzchar(buf)) pieces[[length(pieces) + 1]] <- buf
  Reduce(function(a, b) call("*", a, b), pieces)
}
.words_expr <- function(words) {
  segs <- list()
  for (wd in words) {
    txt <- paste0(if (wd$sp) " " else "", wd$w)
    n <- length(segs)
    if (n && identical(segs[[n]]$it, wd$it)) segs[[n]]$txt <- paste0(segs[[n]]$txt, txt)
    else segs[[n + 1]] <- list(it = wd$it, txt = txt)
  }
  parts <- lapply(segs, function(sg) { e <- .greek_expr(sg$txt); if (sg$it) call("italic", e) else e })
  Reduce(function(a, b) call("*", a, b), parts)
}
.words_text <- function(words) paste0(vapply(words, function(wd) paste0(if (wd$sp) " " else "", wd$w), ""), collapse = "")

# label as an expression on one line, or on two lines split where the longer line is shortest
label_expr <- function(s, two_lines = FALSE) {
  words <- .md_words(s)
  if (!two_lines || length(words) < 2) return(as.expression(.words_expr(words)))
  cand <- which(vapply(words, `[[`, logical(1), "sp"))
  cand <- cand[cand > 1] - 1L
  if (!length(cand)) return(as.expression(.words_expr(words)))
  cost <- vapply(cand, function(j) max(nchar(.words_text(words[1:j])),
                                       nchar(trimws(.words_text(words[(j + 1):length(words)])))), numeric(1))
  j  <- cand[which.min(cost)]
  w2 <- words[(j + 1):length(words)]
  w2[[1]]$sp <- FALSE
  as.expression(call("atop", .words_expr(words[1:j]), .words_expr(w2)))
}

.conf_phrase <- function(conf) {
  if (!length(conf)) return("unadjusted")
  lab <- c(Age_DOS = "age", Teeth = "number of teeth", BMI = "BMI", Sex = "sex",
           Cigsmoke = "smoking", Diabetes = "diabetes")
  l <- ifelse(conf %in% names(lab), lab[conf], conf)
  paste("adjusted for", if (length(l) == 1) l else paste(paste(l[-length(l)], collapse = ", "), "and", l[length(l)]))
}
MODEL_TITLE <- paste("Logistic regression,", .conf_phrase(.conf))
.n_markers  <- function(k) sprintf("%d biomarker%s", k, if (k == 1) "" else "s")


# ==========================================================================
# 6) OUTPUT HELPERS
# ==========================================================================

FILE_TAG <- paste(.pget("Model configuration"), registry_key, sep = "_")
.written <- character(0)
save_figure <- function(stem, width, height, draw, pointsize = 12) {
  for (fmt in FIG_FORMATS) {
    f <- file.path(OUT_DIR, paste0(stem, ".", fmt))
    switch(fmt,
           png  = png(f, width = width, height = height, units = "in", res = PNG_DPI, pointsize = pointsize),
           pdf  = pdf(f, width = width, height = height, pointsize = pointsize, useDingbats = FALSE),
           tiff = tiff(f, width = width, height = height, units = "in", res = PNG_DPI,
                       pointsize = pointsize, compression = "lzw"),
           stop("FIG_FORMATS: use png, pdf or tiff."))
    tryCatch(draw(), finally = dev.off())
    .written <<- c(.written, f)
  }
}
write_table <- function(df, stem) {
  f <- file.path(OUT_DIR, paste0(stem, ".csv"))
  utils::write.csv(df, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  .written <<- c(.written, f)
}


# ==========================================================================
# 7) REFITS AND THE CHECK AGAINST THE WORKBOOK
# ==========================================================================

cat(sprintf("\n===== DMFT figures: %s =====\n", basename(.res_file)))
cat(sprintf("  n = %d (%s) | CV: %d-fold x %d repeats, seed %d | covariates: %s\n",
            nrow(FIG_DAT), paste(names(table(FIG_DAT$condL)), table(FIG_DAT$condL), collapse = ", "),
            FIG_CV_FOLDS, FIG_CV_REPEATS, FIG_CV_SEED, if (length(.conf)) paste(.conf, collapse = ", ") else "none"))
cat(sprintf("  Best overall panel: %s (balanced accuracy %.4f, AUC %.4f)\n", .n_markers(BEST_OVERALL_K),
            FIG_BEST[[as.character(BEST_OVERALL_K)]]$row$Balanced_Accuracy, FIG_BEST[[as.character(BEST_OVERALL_K)]]$row$AUC))

FIG_REFITS <- list(); FIG_REFIT_COV <- NULL
.ks_refit <- c(if (MAKE_ROC) .ks, BEST_OVERALL_K)
.ks_refit <- sort(unique(.ks_refit))
if (length(.ks_refit)) {
  cat(sprintf("  Refitting %d panel(s) on the workbook's folds%s...\n", length(.ks_refit),
              if (MAKE_ROC && ROC_COVARIATES_ONLY && length(.conf)) " and the covariates-only model" else ""))
  for (k in .ks_refit) {
    t0 <- Sys.time()
    FIG_REFITS[[as.character(k)]] <- .refit(FIG_BEST[[as.character(k)]]$markers)
    cat(sprintf("    %-12s AUC %.4f  balanced accuracy %.4f  (%.1fs)\n", .n_markers(k),
                FIG_REFITS[[as.character(k)]]$auc, FIG_REFITS[[as.character(k)]]$ba,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  if (MAKE_ROC && ROC_COVARIATES_ONLY && length(.conf)) {
    FIG_REFIT_COV <- .refit(character(0))
    cat(sprintf("    %-12s AUC %.4f  balanced accuracy %.4f\n", "covariates", FIG_REFIT_COV$auc, FIG_REFIT_COV$ba))
  }
  assign("dmft_fig_cache", .cache, envir = globalenv())

  FIG_CHECK <- do.call(rbind, lapply(names(FIG_REFITS), function(kk) {
    w <- FIG_BEST[[kk]]$row; r <- FIG_REFITS[[kk]]
    data.frame(Model = .n_markers(as.integer(kk)), Panel = paste(label_plain(vapply(r$markers, marker_label, "")), collapse = "; "),
               AUC_workbook = w$AUC, AUC_refit = r$auc,
               Balanced_Accuracy_workbook = w$Balanced_Accuracy, Balanced_Accuracy_refit = r$ba,
               Threshold_workbook = w$Threshold, Threshold_refit = r$youden$threshold,
               Sensitivity_refit = r$youden$sensitivity, Specificity_refit = r$youden$specificity,
               Accuracy_refit = r$youden$accuracy, PPV_refit = r$youden$precision, NPV_refit = r$youden$npv,
               stringsAsFactors = FALSE)
  }))
  .maxdiff <- max(abs(c(FIG_CHECK$AUC_workbook - FIG_CHECK$AUC_refit,
                        FIG_CHECK$Balanced_Accuracy_workbook - FIG_CHECK$Balanced_Accuracy_refit,
                        FIG_CHECK$Threshold_workbook - FIG_CHECK$Threshold_refit)), na.rm = TRUE)
  if (.maxdiff > 1e-6) {
    cat(sprintf("  !! Refits DIFFER from the workbook (largest difference %.2g); see the warning below.\n", .maxdiff))
    stop(sprintf(paste0("The refits differ from the workbook (largest difference %.2g). The data, folds or ",
                           "settings in this session are not those of the workbook's run; the figures show ",
                           "different refits. Export stopped to avoid using an unmatched model."), .maxdiff), call. = FALSE)
  } else {
    cat(sprintf("  Refits match the workbook (largest difference in AUC, balanced accuracy, threshold: %.1e)\n",
                .maxdiff))
  }
  if (!is.null(FIG_REFIT_COV)) {
    FIG_CHECK <- rbind(FIG_CHECK, data.frame(
      Model = "Covariates only", Panel = paste(vapply(.conf, marker_label, ""), collapse = "; "),
      AUC_workbook = NA, AUC_refit = FIG_REFIT_COV$auc, Balanced_Accuracy_workbook = NA,
      Balanced_Accuracy_refit = FIG_REFIT_COV$ba, Threshold_workbook = NA,
      Threshold_refit = FIG_REFIT_COV$youden$threshold, Sensitivity_refit = FIG_REFIT_COV$youden$sensitivity,
      Specificity_refit = FIG_REFIT_COV$youden$specificity, Accuracy_refit = FIG_REFIT_COV$youden$accuracy,
      PPV_refit = FIG_REFIT_COV$youden$precision, NPV_refit = FIG_REFIT_COV$youden$npv, stringsAsFactors = FALSE))
  }
  write_table(FIG_CHECK, paste0("Refit_check_", FILE_TAG))
}


# ==========================================================================
# 8) ROC CURVES OF THE BEST PANEL AT EACH SIZE
# ==========================================================================

draw_roc <- function() {
  op <- par(mar = c(4.6, 4.6, 3.2, 1.6), pty = "s")
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
       main = MODEL_TITLE, cex.main = 0.98)
  abline(0, 1, col = "grey65")

  if (!is.null(FIG_REFIT_COV))
    lines(1 - FIG_REFIT_COV$roc$specificities, FIG_REFIT_COV$roc$sensitivities, col = "grey45", lwd = 2, lty = 2)
  for (kk in names(FIG_REFITS)) {
    k <- as.integer(kk); r <- FIG_REFITS[[kk]]
    lines(1 - r$roc$specificities, r$roc$sensitivities, col = SIZE_COLS[k], lwd = 2)
  }
  if (ROC_MARK_YOUDEN) {
    for (kk in names(FIG_REFITS)) {
      k <- as.integer(kk); yb <- FIG_REFITS[[kk]]$youden
      points(1 - yb$specificity, yb$sensitivity, pch = SIZE_PCH[k], col = "white", cex = 1.9, lwd = 3)
      points(1 - yb$specificity, yb$sensitivity, pch = SIZE_PCH[k], col = SIZE_COLS[k], cex = 1.35, lwd = 1.5)
    }
  }
  ks  <- as.integer(names(FIG_REFITS))
  leg <- sprintf("%s (AUC = %.3f)", vapply(ks, .n_markers, ""), vapply(FIG_REFITS, `[[`, numeric(1), "auc"))
  col <- SIZE_COLS[ks]; lty <- rep(1, length(ks)); pch <- if (ROC_MARK_YOUDEN) SIZE_PCH[ks] else rep(NA, length(ks))
  if (!is.null(FIG_REFIT_COV)) {
    leg <- c(leg, sprintf("Covariates only (AUC = %.3f)", FIG_REFIT_COV$auc))
    col <- c(col, "grey45"); lty <- c(lty, 2); pch <- c(pch, NA)
  }
  legend("bottomright", legend = leg, col = col, lty = lty, lwd = 2, pch = pch, pt.cex = 1.2,
         bty = "n", inset = c(0.01, 0.06), cex = 0.9, seg.len = 2.6)
  if (ROC_MARK_YOUDEN)
    mtext("Symbols: Youden cut of each panel", side = 1, line = -1.3, adj = 0.98, cex = 0.7, col = "grey35")
}

if (MAKE_ROC) {
  save_figure(paste0("Figure_ROC_", FILE_TAG), 6.5, 6.5, draw_roc)
  cat("  ROC curves drawn.\n")
}


# ==========================================================================
# 9) RADAR PLOT OF THE BEST PANEL AT EACH SIZE
# ==========================================================================

.RADAR_METRICS <- c(Balanced_Accuracy = "Balanced accuracy", Sensitivity = "Sensitivity",
                    Specificity = "Specificity", Precision = "PPV", NPV = "NPV", AUC = "AUC")

RADAR_VALUES <- do.call(rbind, lapply(FIG_BEST, function(b)
  data.frame(Model = .n_markers(b$k), Panel = paste(label_plain(vapply(b$markers, marker_label, "")), collapse = "; "),
             b$row[, names(.RADAR_METRICS)], check.names = FALSE, stringsAsFactors = FALSE)))
rownames(RADAR_VALUES) <- NULL

# Area of each polygon as a fraction of the polygon with every metric at 1.
# With equal angles between the axes this is sum(r[i] * r[i + 1]) / m around
# the polygon, r being the radii: the metrics themselves ("full", centre = 0)
# or the radii as drawn ("plotted", centre = RADAR_MIN). The area depends on
# the order of the axes.
.radar_floor <- function(vals) if (is.null(RADAR_MIN)) max(0, floor(min(vals, na.rm = TRUE) * 10) / 10) else RADAR_MIN
.radar_area  <- function(r) sum(r * c(r[-1], r[1])) / length(r)
.rv      <- as.matrix(RADAR_VALUES[, names(.RADAR_METRICS)])
.lo_area <- switch(RADAR_AREA, full = 0, plotted = .radar_floor(.rv),
                   stop("RADAR_AREA must be \"full\" or \"plotted\"."))
RADAR_VALUES$Radar_area <- apply(.rv, 1, function(v) .radar_area(pmin(1, pmax(0, (v - .lo_area) / (1 - .lo_area)))))

# Base R (no fmsb: loading fmsb replaces pROC's plot.roc and print.roc in the session)
draw_radar <- function() {
  vals  <- as.matrix(RADAR_VALUES[, names(.RADAR_METRICS)])
  lo    <- .radar_floor(vals)
  ticks <- seq(lo, 1, length.out = max(1L, round((1 - lo) / 0.1)) + 1)
  m     <- ncol(vals)
  ang   <- pi / 2 - 2 * pi * (seq_len(m) - 1) / m          # first axis at 12 o'clock, clockwise
  rad   <- function(v) pmin(1, pmax(0, (v - lo) / (1 - lo)))
  ks    <- vapply(FIG_BEST, `[[`, numeric(1), "k")

  op <- par(mar = c(6, 1, 3.8, 1), xpd = NA)
  on.exit(par(op))
  plot.new()
  plot.window(xlim = c(-1.3, 1.3), ylim = c(-1.12, 1.12), asp = 1)
  for (t in ticks) polygon(rad(t) * cos(ang), rad(t) * sin(ang), border = "grey82", lwd = 0.8)
  segments(0, 0, cos(ang), sin(ang), col = "grey82", lwd = 0.8)
  text(0.025, rad(ticks), sprintf("%.2f", ticks), adj = c(0, -0.25), cex = 0.7, col = "grey35")

  ca <- cos(ang); sa <- sin(ang)
  for (j in seq_len(m))
    text(1.1 * ca[j], 1.1 * sa[j], .RADAR_METRICS[[j]], cex = 0.95,
         adj = c(if (ca[j] > 0.2) 0 else if (ca[j] < -0.2) 1 else 0.5,
                 if (sa[j] > 0.2) 0 else if (sa[j] < -0.2) 1 else 0.5))

  for (i in seq_along(ks)) {
    r <- rad(vals[i, ])
    polygon(r * ca, r * sa, border = SIZE_COLS[ks[i]], lwd = 2, col = NA)
    points(r * ca, r * sa, pch = SIZE_PCH[ks[i]], col = SIZE_COLS[ks[i]], cex = 1.1, lwd = 1.5)
  }
  title(main = "Selected panel at each size", cex.main = 1.1, line = 2.4)
  mtext(sprintf("%s; radial axis from %.2f (centre) to 1.00", sub("^Logistic regression, ", "", MODEL_TITLE), lo),
        side = 3, line = 1.1, cex = 0.8)
  legend(x = 0, y = -1.25, xjust = 0.5, yjust = 1,
         legend = sprintf("%s (area %.3f)", vapply(ks, .n_markers, ""), RADAR_VALUES$Radar_area),
         col = SIZE_COLS[ks], lwd = 2, pch = SIZE_PCH[ks], ncol = 2, bty = "n", cex = 0.9, seg.len = 2.4,
         x.intersp = 0.8, text.width = 0.95)
  mtext(if (RADAR_AREA == "full")
          "Area: polygon area with each metric on a 0-1 scale, as a fraction of the area with every metric at 1"
        else sprintf("Area: polygon as drawn (centre %.2f), as a fraction of the area with every metric at 1", lo),
        side = 1, line = 4.5, cex = 0.72, col = "grey30")
}

if (MAKE_RADAR) {
  save_figure(paste0("Figure_Radar_", FILE_TAG), 7, 7.4, draw_radar)
  write_table(RADAR_VALUES, paste0("Figure_Radar_values_", FILE_TAG))
  cat("  Radar plot drawn. Area per panel size:",
      paste(sprintf("%d: %.3f", vapply(FIG_BEST, `[[`, numeric(1), "k"), RADAR_VALUES$Radar_area), collapse = " | "), "\n")
}


# ==========================================================================
# 10) HEATMAP OF THE BEST OVERALL PANEL
# ==========================================================================

prepare_heatmap <- function() {
  r       <- FIG_REFITS[[as.character(BEST_OVERALL_K)]]
  markers <- r$markers
  covs    <- if (HEATMAP_COVARIATES) .conf else character(0)

  p_mean <- tapply(r$pred$p, r$pred$rowIndex, mean)
  p_all  <- rep(NA_real_, nrow(FIG_DAT)); p_all[as.integer(names(p_mean))] <- p_mean
  grp    <- FIG_DAT$condL
  is_case <- grp == .lev[2]

  cut <- if (is.numeric(HEATMAP_CUTOFF)) HEATMAP_CUTOFF else switch(HEATMAP_CUTOFF,
    pooled      = r$youden$threshold,
    participant = .youden_best(pROC::roc(grp, p_all, levels = .lev, direction = "<", quiet = TRUE))$threshold,
    stop("HEATMAP_CUTOFF must be \"pooled\", \"participant\" or a number."))
  pred_case <- p_all >= cut
  mis <- pred_case != is_case

  dmft <- if (exists("dmft_cols") && all(dmft_cols %in% colnames(FIG_DAT)))
    rowSums(sapply(FIG_DAT[, dmft_cols, drop = FALSE], function(z) as.numeric(as.character(z)))) else rep(NA_real_, nrow(FIG_DAT))

  key2 <- if (HEATMAP_ORDER == "dmft" && !all(is.na(dmft))) dmft else rep(0, nrow(FIG_DAT))
  ord  <- order(is_case, key2, p_all)

  cols  <- c(markers, covs)
  scale01 <- function(x) { x <- as.numeric(x); rg <- range(x, na.rm = TRUE)
                           if (!is.finite(diff(rg)) || diff(rg) == 0) rep(0.5, length(x)) else (x - rg[1]) / diff(rg) }
  S <- sapply(cols, function(v) scale01(FIG_DAT[[v]]))[ord, , drop = FALSE]

  n_ref <- sum(!is_case); n_case <- sum(is_case)
  sens  <- mean(pred_case[is_case]); spec <- mean(!pred_case[!is_case])
  list(r = r, markers = markers, covs = covs, cols = cols, S = S, ord = ord,
       dmft = dmft[ord], p = p_all[ord], mis = mis[ord], is_case = is_case[ord], cut = cut,
       n_ref = n_ref, n_case = n_case, sens = sens, spec = spec,
       tp = sum(pred_case[is_case]), tn = sum(!pred_case[!is_case]),
       id = if ("StudyID" %in% colnames(FIG_DAT)) as.character(FIG_DAT$StudyID)[ord] else as.character(ord))
}

draw_heatmap <- function(hm, W, Ht) {
  pal <- switch(HEATMAP_PALETTE,
                blue_orange = grDevices::colorRampPalette(c("#0072B2", "#E69F00"))(101),
                cm    = grDevices::cm.colors(101),
                blues = grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(101),
                stop("Unknown HEATMAP_PALETTE."))
  op <- par(mar = c(0, 0, 0, 0), xpd = NA)
  on.exit(par(op))
  plot.new()
  plot.window(xlim = c(0, W), ylim = c(0, Ht), xaxs = "i", yaxs = "i")

  n  <- nrow(hm$S); rh <- HEATMAP_ROW_IN; cw <- 0.56; gap <- 0.16
  x_dmft <- 0.42; x_p <- 0.86; x0 <- 1.12; x_lab0 <- 0.14
  nm <- length(hm$markers); nc <- length(hm$covs)
  x_left <- c(x0 + (seq_len(nm) - 1) * cw, if (nc) x0 + nm * cw + gap + (seq_len(nc) - 1) * cw)
  x_end  <- max(x_left) + cw
  y_top  <- Ht - 1.72
  row_t  <- y_top - (seq_len(n) - 1) * rh
  row_b  <- row_t - rh
  cex_row <- min(0.8, rh * 72 * 0.9 / par("ps"))

  # title
  text(x_lab0, Ht - 0.24, sprintf("%s  |  %s", MODEL_TITLE, .n_markers(BEST_OVERALL_K)),
       adj = c(0, 0.5), font = 2, cex = 1.05)
  text(x_lab0, Ht - 0.48, sprintf("Mean-OOF classification: sensitivity %.3f (%d/%d), specificity %.3f (%d/%d), cutoff %.3f",
                                  hm$sens, hm$tp, hm$n_case, hm$spec, hm$tn, hm$n_ref, hm$cut),
       adj = c(0, 0.5), cex = 0.9)

  # cells
  cidx <- round(hm$S * 100) + 1
  colm <- matrix(pal[cidx], nrow = n)
  colm[is.na(hm$S)] <- "grey85"
  xl <- rep(x_left, each = n); yb <- rep(row_b, times = length(x_left))
  rect(xl, yb, xl + cw, yb + rh, col = as.vector(colm), border = "white", lwd = 0.5)

  # Misclassified participants: outline every cell without changing its fill.
  if (any(hm$mis)) for (xx in x_left)
    rect(xx, row_b[hm$mis], xx + cw, row_t[hm$mis], col = NA, border = "black", lwd = 1.2)

  # row labels
  lab_col <- "black"; lab_font <- ifelse(hm$mis, 2, 1)
  ymid <- (row_t + row_b) / 2
  if (!all(is.na(hm$dmft)))
    text(x_dmft, ymid, ifelse(is.na(hm$dmft), "", format(hm$dmft, trim = TRUE)), cex = cex_row, col = lab_col, font = lab_font)
  text(x_p, ymid, sprintf("%.2f", hm$p), cex = cex_row, col = lab_col, font = lab_font)
  if (!all(is.na(hm$dmft))) text(x_dmft, y_top + 0.06, "DMFT", adj = c(0.5, 0), cex = 0.85)
  text(x_p, y_top + 0.06, sprintf("P(%s)", .lev[2]), adj = c(0.5, 0), cex = 0.85)

  # column labels, rotated
  for (j in seq_along(hm$cols)) {
    lab <- if (hm$cols[j] == "Sex") paste0("Sex: ", paste(levels(FIG_DAT$Sex), collapse = " / ")) else
      label_expr(marker_label(hm$cols[j]), two_lines = j <= nm)
    text(x_left[j] + 0.16, y_top + 0.06, lab, srt = 45, adj = c(0, 0), cex = 0.85)
  }

  # group line and brackets
  n_ref <- sum(!hm$is_case)
  y_sep <- row_b[n_ref]
  segments(x_lab0, y_sep, x_end, y_sep, lwd = 1.3)
  x_br <- x_end + 0.12
  segments(x_br, row_t[1] - 0.01, x_br, y_sep + 0.03, lwd = 1)
  segments(x_br, y_sep - 0.03, x_br, row_b[n] + 0.01, lwd = 1)
  text(x_br + 0.13, (row_t[1] + y_sep) / 2, sprintf("Low DMFT (%s), n = %d", .lev[1], n_ref), srt = 90, cex = 0.85)
  text(x_br + 0.13, (y_sep + row_b[n]) / 2, sprintf("High DMFT (%s), n = %d", .lev[2], n - n_ref), srt = 90, cex = 0.85)

  # legends
  y1 <- row_b[n] - 0.34; bw <- 1.5
  xs <- seq(x0, x0 + bw, length.out = 102)
  rect(xs[-102], y1 - 0.065, xs[-1], y1 + 0.065, col = pal, border = NA)
  rect(x0, y1 - 0.065, x0 + bw, y1 + 0.065, border = "grey40", lwd = 0.6)
  text(x0, y1 - 0.11, "0", adj = c(0.5, 1), cex = 0.75)
  text(x0 + bw, y1 - 0.11, "1", adj = c(0.5, 1), cex = 0.75)
  text(x0 + bw + 0.15, y1, "Scaled 0 (lowest) to 1 (highest) within each column", adj = c(0, 0.5), cex = 0.85)
  y2 <- y1 - 0.42
  rect(x0, y2 - 0.065, x0 + 0.3, y2 + 0.065, col = NA, border = "black", lwd = 1.2)
  text(x0 + 0.45, y2, "Misclassified at the cutoff: false positive above the line, false negative below",
       adj = c(0, 0.5), cex = 0.85)
}

if (MAKE_HEATMAP) {
  FIG_HM <- prepare_heatmap()
  .W  <- 1.12 + 0.56 * length(FIG_HM$markers) + (if (length(FIG_HM$covs)) 0.16 + 0.56 * length(FIG_HM$covs) else 0) + 0.62
  .Ht <- 1.72 + nrow(FIG_HM$S) * HEATMAP_ROW_IN + 0.95
  save_figure(paste0("Figure_Heatmap_", FILE_TAG), max(.W, 7.2), .Ht, function() draw_heatmap(FIG_HM, max(.W, 7.2), .Ht),
              pointsize = 10)

  write_table(data.frame(StudyID = FIG_HM$id, Group = ifelse(FIG_HM$is_case, .lev[2], .lev[1]), DMFT = FIG_HM$dmft,
                         P_HD_mean = round(FIG_HM$p, 4), Predicted = ifelse(FIG_HM$p >= FIG_HM$cut, .lev[2], .lev[1]),
                         Misclassified = FIG_HM$mis, check.names = FALSE),
              paste0("Figure_Heatmap_participants_", FILE_TAG))

  cat(sprintf("  Heatmap drawn: %s | cutoff %.4f | sensitivity %.3f (%d of %d) | specificity %.3f (%d of %d)\n",
              .n_markers(BEST_OVERALL_K), FIG_HM$cut, FIG_HM$sens, FIG_HM$tp, FIG_HM$n_case, FIG_HM$spec, FIG_HM$tn, FIG_HM$n_ref))
  cat(sprintf("  Misclassified participants: %d false positive(s), %d false negative(s)\n",
              sum(FIG_HM$mis & !FIG_HM$is_case), sum(FIG_HM$mis & FIG_HM$is_case)))
}



# Export a complete selected model independently of which figures are requested.
.fm <- FIG_REFITS[[as.character(BEST_OVERALL_K)]]$final
.markers <- FIG_BEST[[as.character(BEST_OVERALL_K)]]$markers
.cf <- summary(.fm)$coefficients
.st <- if (exists("stats_PH_vif")) stats_PH_vif else if (exists("stats_PH")) stats_PH else NULL
if (is.null(.st)) stop("LD reference means and SDs are required for the calculator export.")
.mu <- .st$Mean[match(.markers, .st$Column)]
.sd <- .st$SD[match(.markers, .st$Column)]
if (any(!is.finite(.mu)) || any(!is.finite(.sd)) || any(.sd <= 0))
  stop("Missing or invalid LD reference statistics for selected biomarkers.")
.beta <- coef(.fm)
if (any(!is.finite(.beta))) stop("Non-finite model coefficients; no calculator export.")
.raw <- .beta
.raw[.markers] <- .beta[.markers] / .sd
.raw["(Intercept)"] <- .beta["(Intercept)"] - sum(.beta[.markers] * .mu / .sd)
.input <- setNames(rep("value in original units", length(.raw)), names(.raw))
.input["(Intercept)"] <- ""
.input[.markers] <- "log2(raw value + 5e-07)"
.refs <- character(0)
for (v in .conf) if (is.factor(FIG_DAT[[v]])) {
  lev <- levels(FIG_DAT[[v]])
  for (ll in lev[-1]) .input[paste0(v, ll)] <- sprintf("indicator: %s = %s", v, ll)
  .refs[paste0("(Reference level: ", v, ")")] <- lev[1]
}
# Verify the raw-input equation against the fitted model for every participant.
.xx <- model.matrix(.fm)
for (j in seq_along(.markers)) .xx[, .markers[j]] <- .xx[, .markers[j]] * .sd[j] + .mu[j]
if (max(abs(drop(.xx %*% .raw) - predict(.fm, type = "link"))) > 1e-8)
  stop("Raw-input equation verification failed.")
FIG_COEFS <- data.frame(Term = names(.raw), Label = names(.raw),
  Coefficient_raw_input = unname(.raw), Input = unname(.input[names(.raw)]), stringsAsFactors = FALSE)
.meta <- data.frame(Term = c("(Threshold)", "(Delta)", "(Case group)", "(Reference group)", names(.refs)),
  Label = c("Classification threshold", "log2 offset", "high DMFT (DMFT > 13)", "low DMFT (DMFT <= 13)", unname(.refs)),
  Coefficient_raw_input = c(FIG_BEST[[as.character(BEST_OVERALL_K)]]$row$Threshold, 5e-7, NA, NA, rep(NA, length(.refs))),
  Input = "", stringsAsFactors = FALSE)
FIG_COEFS <- rbind(FIG_COEFS, .meta)
write_table(FIG_COEFS, paste0("Model_coefficients_best_panel_", FILE_TAG))
write_table(FIG_COEFS, "dmft_selected_model")
for (kk in names(FIG_REFITS)) write_table(FIG_REFITS[[kk]]$pred, paste0("OOF_predictions_k", kk, "_", FILE_TAG))
message("Load dmft_selected_model.csv in the Shiny app, or place it next to app.R and restart.")

if (MAKE_PERFORMANCE) {
  .perf <- do.call(rbind, lapply(FIG_BEST, function(x) data.frame(
    k = x$k, Balanced_Accuracy = x$row$Balanced_Accuracy, AUC = x$row$AUC,
    Lower = x$row$Balanced_Accuracy_CI_Lower, Upper = x$row$Balanced_Accuracy_CI_Upper)))
  write_table(.perf, paste0("Best_panel_performance_", FILE_TAG))
  save_figure(paste0("Figure_Performance_", FILE_TAG), 8, 4.2, function() {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    plot(.perf$k, 100 * .perf$Balanced_Accuracy, type = "b", pch = 16, col = "#0072B2",
         ylim = range(100 * c(.perf$Lower, .perf$Upper)), xlab = "Number of biomarkers",
         ylab = "Balanced accuracy (%)", main = "Classification performance")
    arrows(.perf$k, 100 * .perf$Lower, .perf$k, 100 * .perf$Upper,
           angle = 90, code = 3, length = 0.04, col = "#0072B2")
    plot(.perf$k, .perf$AUC, type = "b", pch = 16, col = "#E69F00",
         xlab = "Number of biomarkers", ylab = "AUC", main = "Discrimination")
  })
  message("Performance bars are 2.5th–97.5th repeat percentiles, not independent-validation confidence intervals.")
}

# ==========================================================================
# 11) SUMMARY
# ==========================================================================
cat("\n===== Files written =====\n")
for (.f in .written) cat("  ", .f, "\n")
cat("\nAll done.\n")
