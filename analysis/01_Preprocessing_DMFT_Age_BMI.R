# ============================================================
# SCRIPT 1 — DATA PREPROCESSING & FEATURE ENGINEERING
# ============================================================
# Builds the analysis data for one of three projects (PROJECT, section 1):
#
#   caries      no caries vs high caries experience (Cavitated_Caries), or low
#               vs high caries burden (Clinical_Caries), among subjects with
#               healthy periodontium or gingivitis
#   DMFT        low vs high DMFT (decayed + missing + filled teeth), split at
#               DMFT_CUTOFF
#   validation  the panels of Miller, Yan et al. 2025 (J Clin Periodontol,
#               doi:10.1111/jcpe.70049) tested on this dataset, in groups of
#               periodontitis cases against the same healthy controls
#               (PanelValidation), or the paper's own comparison (PaperDesign)
#   perio       P vs H, G vs H and G vs P, with periodontal status from the
#               original diagnosis, PISA or AdjustPRA
#
# What each project fixes (everything else is a setting in section 1):
#
#                candidate markers   TIMP-1 ratios  bacterial ratios  VIF     S100A8
#   caries       selected or all     no             always            toggle  kept
#   DMFT         selected or all     no             toggle            toggle  kept
#   validation   the paper's 31 *    yes            no                no      excluded *
#   perio        selected or all     no             toggle            toggle  kept
#   (* EXCLUDE_S100A8 = FALSE keeps S100A8: the paper's 33)
#
# condL is a factor with the reference group as the first level and the case
# as the second, e.g. H / G for G vs H and LC / HC for HC vs LC. Script 2 takes
# the second level as the case.
#
# Steps: read and merge -> taxon names -> TIMP-1 ratios (validation) ->
# KNN imputation -> candidate markers -> log2 -> outcome labels -> z-scores
# against the reference group -> bacterial ratio features -> VIF ->
# biomarker_registry.
#
# Inputs : three Excel files in DATA_DIR; script1_panels.R (marker lists)
# Outputs: data.dfPH (perio also data.dfGH, data.dfGP), log2 and z-scored;
#          data.df7_imputed (raw scale, imputed); validation_sets;
#          stats_*, biomarkers_*_vif; an entry in biomarker_registry.
#          Then run Script 2 (combinatorial logistic regression).
# ============================================================

library(readxl)
library(dplyr)
library(VIM)


# ============================================================
# 1) SETTINGS
# ============================================================
PROJECT <- "DMFT"              # "caries" | "DMFT" | "validation" | "perio"

# --- Confounders (all projects) ---------------------------------------------
# Any columns of the data, e.g. Age_DOS, Teeth, BMI, Sex, Age_DOS_cat / BMI_cat
# (need FACTORIZE_AGE_BMI = TRUE), Cigsmoke, Diabetes, T1Diabetes, T2Diabetes,
# DMFT. NULL = none. They keep their original units (not log2 or z-scored).
# Age and BMI are the study adjustment variables.
DMFT_MODEL_SETS <- list(
  Biomarkers_Age_BMI = list(biomarkers = TRUE, covariates = c("Age_DOS", "BMI"))
)
CONFOUNDERS <- c("Age_DOS", "BMI")

FACTORIZE_AGE_BMI <- TRUE            # adds Age_DOS_cat and BMI_cat
AGE_CUTOFF        <- 50
BMI_CUTOFF        <- 30

# --- caries -----------------------------------------------------------------
# Cavitated_Caries: CF (no caries) if <= CARIES_NOCARIES_MAX, CP (high caries
#   experience) if >= CARIES_HIGH_MIN; values in between are dropped.
# Clinical_Caries: LC if <= CARIES_BINARY_CUTOFF, HC if above (cutoff from
#   Youden's J against Cavitated_Caries).
# Bacterial ratio features are always on.
CARIES_SOURCE_COL         <- "Cavitated_Caries"     # or "Clinical_Caries"
CARIES_NOCARIES_MAX       <- 0
CARIES_HIGH_MIN           <- 2
CARIES_BINARY_CUTOFF      <- 5
CARIES_PERIODONTAL_FILTER <- c("Healthy periodontium", "Gingivitis")   # original diagnosis; NULL = everyone
CARIES_MARKERS            <- "selected"             # "selected" (script1_panels.R) | "all"
CARIES_VIF                <- TRUE

# --- DMFT -------------------------------------------------------------------
# DMFT = the sum of DMFT_COMPONENTS (NA if any of them is NA). LD (low DMFT) if
# DMFT <= DMFT_CUTOFF, HD (high DMFT) if above. DMFT_CUTOFF = NULL splits at the
# median of the analysed subjects (printed when the script runs).
DMFT_COMPONENTS         <- c("Missing_Teeth", "Restored_Teeth", "Teeth_wCarries")
DMFT_CUTOFF             <- NULL
DMFT_PERIODONTAL_FILTER <- NULL          # original diagnosis kept, e.g. c("Healthy periodontium", "Gingivitis"); NULL = everyone
DMFT_MARKERS            <- "all"         # "selected" (biomarkers_DMFT in script1_panels.R) | "all"
DMFT_RATIO_FEATURES     <- FALSE
DMFT_VIF                <- FALSE

# --- validation -------------------------------------------------------------
# Always: the paper's markers and TIMP-1 ratios, no bacterial ratios, no VIF.
# PanelValidation: groups of periodontitis cases, each against the same controls:
#   1) all periodontitis   2) periodontitis with DM_COL = Yes
#   3) periodontitis with type 1 or type 2 diabetes (VALIDATION_ANY_DIABETES_GROUP;
#      not added when DM_COL = "Diabetes")
#   4) periodontitis without diabetes
# PaperDesign: periodontitis with DM_COL = Yes vs healthy periodontium without
#   diabetes, as in the paper.
# "Without diabetes" always means no diabetes of either type.
VALIDATION_DESIGN             <- "PanelValidation"      # or "PaperDesign"
VALIDATION_CONTROLS           <- "Healthy_NoDiabetes"   # or "Healthy_All" (PanelValidation)
VALIDATION_ANY_DIABETES_GROUP <- TRUE
DM_COL                        <- "T2Diabetes"           # "T2Diabetes" | "T1Diabetes" | "Diabetes" (either type)
DM_PERIO_DEFINITION           <- "original"             # periodontal status from "original" | "PISA" | "AdjustPRA"
EXCLUDE_S100A8                <- TRUE                   # drops S100A8 and S100A8/TIMP-1 (33 -> 31 markers)

# --- perio ------------------------------------------------------------------
PERIO_DEFINITION     <- "original"   # "original" | "PISA" | "AdjustPRA"
PERIO_MARKERS        <- "selected"   # "selected" (script1_panels.R) | "all"
PERIO_RATIO_FEATURES <- TRUE
PERIO_VIF            <- FALSE

# --- Shared -----------------------------------------------------------------
TOP_N_RATIOS  <- 10                  # bacterial ratio features kept (largest case / reference fold change)
VIF_THRESHOLD <- 10
# PISA / AdjustPRA: Healthy below *_HEALTHY_MAX, Gingivitis up to *_MOD_MAX,
# Periodontitis above (perio, and validation when DM_PERIO_DEFINITION uses them)
PISA_HEALTHY_MAX      <- 4.765
PISA_MOD_MAX          <- 16.235
ADJUSTPRA_HEALTHY_MAX <- 5.50
ADJUSTPRA_MOD_MAX     <- 14.69

# --- Files ------------------------------------------------------------------
# The first of these folders that exists on this computer is used.
DATA_DIR_CHOICES <- c("C:/Users/qya231/OneDrive - University of Kentucky",
                      "C:/Users/qiyan/OneDrive - University of Kentucky")
# script1_panels.R is read from this script's folder (when run with Source);
# otherwise from the working directory.
SCRIPT_DIR  <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
PANELS_FILE <- file.path(SCRIPT_DIR, "script1_panels.R")


# ============================================================
# 2) DERIVED SETTINGS (set from PROJECT; do not edit)
# ============================================================
check_choice <- function(x, choices, name) {
  if (length(x) != 1 || is.na(x) || !x %in% choices)
    stop(sprintf("%s must be one of: %s", name, paste(sprintf("'%s'", choices), collapse = ", ")),
         call. = FALSE)
}
check_choice(PROJECT, c("caries", "DMFT", "validation", "perio"), "PROJECT")
if (PROJECT == "caries") {
  check_choice(CARIES_SOURCE_COL, c("Cavitated_Caries", "Clinical_Caries"), "CARIES_SOURCE_COL")
  check_choice(CARIES_MARKERS, c("selected", "all"), "CARIES_MARKERS")
}
if (PROJECT == "DMFT") {
  check_choice(DMFT_MARKERS, c("selected", "all"), "DMFT_MARKERS")
  if (!is.null(DMFT_CUTOFF) && (length(DMFT_CUTOFF) != 1 || !is.numeric(DMFT_CUTOFF)))
    stop("DMFT_CUTOFF must be one number, or NULL for the median.", call. = FALSE)
}
if (PROJECT == "validation") {
  check_choice(VALIDATION_DESIGN, c("PanelValidation", "PaperDesign"), "VALIDATION_DESIGN")
  check_choice(VALIDATION_CONTROLS, c("Healthy_NoDiabetes", "Healthy_All"), "VALIDATION_CONTROLS")
  check_choice(DM_COL, c("T2Diabetes", "T1Diabetes", "Diabetes"), "DM_COL")
  check_choice(DM_PERIO_DEFINITION, c("original", "PISA", "AdjustPRA"), "DM_PERIO_DEFINITION")
}
if (PROJECT == "perio") {
  check_choice(PERIO_DEFINITION, c("original", "PISA", "AdjustPRA"), "PERIO_DEFINITION")
  check_choice(PERIO_MARKERS, c("selected", "all"), "PERIO_MARKERS")
}

CONDL_SOURCE <- switch(PROJECT, caries = "Caries", DMFT = "DMFT", validation = VALIDATION_DESIGN,
                       perio = PERIO_DEFINITION)
USE_ALL_BIOMARKERS <- switch(PROJECT, caries = CARIES_MARKERS == "all", DMFT = DMFT_MARKERS == "all",
                             validation = FALSE, perio = PERIO_MARKERS == "all")
USE_PUBLICATION_BIOMARKERS <- PROJECT == "validation"
USE_TIMP_RATIOS            <- PROJECT == "validation"
USE_RATIO_FEATURES <- switch(PROJECT, caries = TRUE, DMFT = DMFT_RATIO_FEATURES, validation = FALSE,
                             perio = PERIO_RATIO_FEATURES)
USE_VIF            <- switch(PROJECT, caries = CARIES_VIF, DMFT = DMFT_VIF, validation = FALSE,
                             perio = PERIO_VIF)
EXCLUDE_S100A8     <- PROJECT == "validation" && EXCLUDE_S100A8
SINGLE_COMPARISON  <- PROJECT != "perio"        # every project but perio uses data.dfPH only
DM_COMPARISON      <- PROJECT == "validation"   # groups defined by diabetes status
TIMP_RATIO_SET     <- "publication"             # recorded in the registry
PUBLICATION_SET    <- "paper33"


# ============================================================
# 3) HELPER FUNCTIONS
# ============================================================

# "Bacteria;...;Genus;Species(100)" -> "Genus_Species"
shorten_bacteria_name_v2 <- function(name) {
  parts    <- strsplit(name, ";")[[1]]
  parts    <- parts[parts != ""]
  last_two <- tail(parts, 2)
  last_two <- gsub("\\(100\\)", "", last_two)
  paste(last_two, collapse = "_")
}

rename_bacteria_cols <- function(df, bact_cols, short_names) {
  current_names      <- colnames(df)
  idx                <- match(bact_cols, current_names)
  idx                <- idx[!is.na(idx)]
  current_names[idx] <- short_names[bact_cols %in% current_names]
  colnames(df)       <- current_names
  df
}

# Syntactic column names: brackets and hyphens dropped, "&" -> "and", spaces -> "_"
clean_name <- function(name) {
  name <- gsub("\\[", "", name);  name <- gsub("\\]", "", name)
  name <- gsub("\\(", "", name);  name <- gsub("\\)", "", name)
  name <- gsub("&", "and", name); name <- gsub(" ", "_", name)
  name <- gsub("-", "", name)
  name
}

# REDCap checkbox export ("Checked" / "Unchecked") -> factor No / Yes, reference
# "No". REDCap writes an unanswered checkbox as "Unchecked", so "No" means "not
# checked", not a confirmed no.
recode_checkbox <- function(x, var_name) {
  x_chr <- tolower(trimws(as.character(x)))
  x_chr[x_chr == ""] <- NA_character_
  yes <- c("checked",   "yes", "1")    # "yes" / "no" accepted, so re-running is harmless
  no  <- c("unchecked", "no",  "0")
  unexpected <- setdiff(unique(x_chr[!is.na(x_chr)]), c(yes, no))
  if (length(unexpected) > 0)
    warning(sprintf("%s: unexpected value(s) set to NA: %s",
                    var_name, paste(unexpected, collapse = ", ")))
  out <- rep(NA_character_, length(x_chr))
  out[x_chr %in% yes] <- "Yes"
  out[x_chr %in% no]  <- "No"
  factor(out, levels = c("No", "Yes"))
}

# PISA / AdjustPRA score -> periodontal class
recode_from_score <- function(x, low, high) {
  x   <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x <  low]             <- "Healthy periodontium"
  out[!is.na(x) & x >= low & x <= high] <- "Gingivitis"
  out[!is.na(x) & x >  high]            <- "Periodontitis"
  out
}

# z-scores against the reference group's mean and SD (centred only if SD is 0)
standardize_by_reference <- function(df, biomarkers, ref_group_label, group_col = "condL") {
  ref_data <- df[df[[group_col]] == ref_group_label, biomarkers, drop = FALSE]
  stats_df <- data.frame(
    Column = biomarkers,
    Mean   = sapply(ref_data, mean, na.rm = TRUE),
    SD     = sapply(ref_data, sd,   na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  df_std <- df
  for (col in biomarkers) {
    m <- stats_df$Mean[stats_df$Column == col]
    s <- stats_df$SD[stats_df$Column == col]
    df_std[[col]] <- if (!is.na(s) && s > 0) (df[[col]] - m) / s else df[[col]] - m
  }
  list(data = df_std, stats = stats_df)
}

# Bacterial ratio features: every bacterium with a higher raw mean in cases
# over every bacterium with a lower one; the top_N pairs by case / reference
# fold change of the mean ratio are added as log2(num / den) columns.
build_ratio_features <- function(data_df, raw_df, bacteria_bm,
                                 case_label, ref_label,
                                 top_N = 10, eps = 1e-9) {
  cat(sprintf("\n  Ratio features %s%s (case = %s, reference = %s)\n",
              case_label, ref_label, case_label, ref_label))

  raw_sub       <- raw_df[raw_df$StudyID %in% data_df$StudyID, c("StudyID", bacteria_bm), drop = FALSE]
  raw_sub$condL <- data_df$condL[match(raw_sub$StudyID, data_df$StudyID)]
  ref_rows      <- raw_sub[raw_sub$condL == ref_label,  , drop = FALSE]
  case_rows     <- raw_sub[raw_sub$condL == case_label, , drop = FALSE]

  mean_ref  <- colMeans(ref_rows[,  bacteria_bm, drop = FALSE], na.rm = TRUE)
  mean_case <- colMeans(case_rows[, bacteria_bm, drop = FALSE], na.rm = TRUE)
  direction_df <- data.frame(
    Biomarker     = bacteria_bm,
    Mean_ref      = mean_ref,
    Mean_case     = mean_case,
    Diff_case_ref = mean_case - mean_ref,
    Pct_Change    = (mean_case - mean_ref) / (mean_ref + eps) * 100,
    row.names     = NULL
  )
  numerators   <- direction_df[direction_df$Diff_case_ref > 0, ]
  denominators <- direction_df[direction_df$Diff_case_ref < 0, ]
  numerators   <- numerators[order(numerators$Diff_case_ref, decreasing = TRUE), ]
  denominators <- denominators[order(denominators$Diff_case_ref), ]

  if (nrow(numerators) == 0 || nrow(denominators) == 0) {
    cat("  No bacterium is higher (or none lower) in cases; no ratio features.\n")
    return(list(data_df = data_df, pair_results = NULL, ratio_cols = character(0)))
  }

  # all numerator x denominator pairs, numerators in the outer order
  pairs <- expand.grid(Denominator = denominators$Biomarker, Numerator = numerators$Biomarker,
                       stringsAsFactors = FALSE)
  mean_ratio <- function(rows) {
    mapply(function(num, den) mean((rows[[num]] + eps) / (rows[[den]] + eps), na.rm = TRUE),
           pairs$Numerator, pairs$Denominator, USE.NAMES = FALSE)
  }
  ratio_ref  <- mean_ratio(ref_rows)
  ratio_case <- mean_ratio(case_rows)
  pair_results <- data.frame(
    Numerator       = pairs$Numerator,
    Denominator     = pairs$Denominator,
    Mean_Ratio_ref  = ratio_ref,
    Mean_Ratio_case = ratio_case,
    Ratio_Fold      = ratio_case / (ratio_ref + eps),
    Log2_Fold       = log2((ratio_case + eps) / (ratio_ref + eps)),
    stringsAsFactors = FALSE
  )
  pair_results <- pair_results[order(pair_results$Ratio_Fold, decreasing = TRUE), ]
  rownames(pair_results) <- NULL

  top_pairs  <- head(pair_results, top_N)
  ratio_cols <- paste0("ratio__", top_pairs$Numerator, "__div__", top_pairs$Denominator)
  rows       <- match(data_df$StudyID, raw_df$StudyID)
  for (i in seq_along(ratio_cols)) {
    num_vals <- raw_df[[top_pairs$Numerator[i]]][rows]
    den_vals <- raw_df[[top_pairs$Denominator[i]]][rows]
    data_df[[ratio_cols[i]]] <- log2((num_vals + eps) / (den_vals + eps))
  }
  cat(sprintf("  Added %d log2 ratio features\n", length(ratio_cols)))
  list(data_df = data_df, pair_results = pair_results, ratio_cols = ratio_cols)
}

# Numeric, non-constant columns among bm
safe_cols <- function(df, bm) {
  X <- df[, intersect(bm, colnames(df)), drop = FALSE]
  X <- X[, sapply(X, is.numeric),                          drop = FALSE]
  X <- X[, sapply(X, function(z) sd(z, na.rm = TRUE) > 0), drop = FALSE]
  colnames(X)
}

# Backward VIF pruning: drop the column with the highest VIF until every VIF is
# below Threshold. VIF = 1 / (1 - R^2) from regressing each column on the
# others. VIFs are compared at 15 significant digits, so near-ties go to the
# earlier column.
VIF.fn <- function(in_frame, Threshold = 10, Trace = TRUE) {
  dat <- in_frame
  repeat {
    vars <- names(dat)
    vif  <- vapply(vars, function(v) {
      fit <- lm(as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = "+"))), data = dat)
      1 / (1 - summary(fit)$r.squared)
    }, numeric(1))
    vif   <- as.numeric(as.character(vif))
    worst <- which(vif == max(vif, na.rm = TRUE))[1]
    if (is.na(worst) || vif[worst] < Threshold) break
    if (Trace) cat(sprintf("    removed %s (VIF %.2f)\n", vars[worst], vif[worst]))
    dat <- dat[, -worst, drop = FALSE]
  }
  if (Trace) cat(sprintf("    %d kept, all VIF < %g\n", ncol(dat), Threshold))
  names(dat)
}

run_vif_prune <- function(df, bm, Threshold = 10, Trace = TRUE) {
  X <- df[, safe_cols(df, bm), drop = FALSE]
  if (ncol(X) < 2) return(colnames(X))
  VIF.fn(X, Threshold = Threshold, Trace = Trace)
}

# condL as a factor: reference level first, case second. Script 2's
# factor(condL) keeps this order; on plain text it would sort alphabetically.
as_outcome <- function(d, ref, case) {
  if (is.null(d)) return(NULL)
  d$condL <- factor(d$condL, levels = c(ref, case))
  if (anyNA(d$condL)) stop(sprintf("condL has values other than %s and %s.", ref, case))
  d
}

# Selected-marker list -> cleaned names present in the data (NULL stays NULL)
select_markers <- function(markers, label) {
  if (is.null(markers)) return(NULL)
  cleaned <- sapply(markers, clean_name, USE.NAMES = FALSE)
  missing <- setdiff(cleaned, colnames(data.df7_imputed))
  if (length(missing) > 0)
    cat(sprintf("    %s: not in the data, skipped: %s\n", label, paste(missing, collapse = ", ")))
  cleaned[cleaned %in% colnames(data.df7_imputed)]
}

# The paper's table labels -> column names; stops on an unknown label
translate_paper <- function(x, what) {
  unknown <- setdiff(x, names(paper_label_map))
  if (length(unknown) > 0)
    stop(what, ": label(s) not in paper_label_map: ", paste(unknown, collapse = ", "),
         " (script1_panels.R).")
  unname(paper_label_map[x])
}

# Why one of the paper's markers can or cannot be used
marker_status <- function(cur) {
  ifelse(cur %in% excluded_markers, "excluded (EXCLUDE_S100A8)",
         ifelse(!cur %in% publication_markers, "not in the paper's marker set",
                ifelse(cur %in% colnames(data.df7_imputed), "ok", "not in data")))
}

# On the log2 scale x / TIMP1 = x - TIMP1 (up to the pseudo-count), so a glm
# panel holding x, TIMP1 and x / TIMP1 together is not identifiable.
report_triples <- function(bm, label) {
  if (is.null(bm) || !"TIMP1" %in% bm) return(invisible(NULL))
  nums <- sub("_div_TIMP1$", "", intersect(timp_ratios_full, bm))
  trip <- nums[nums %in% bm]
  if (length(trip) > 0)
    cat(sprintf(paste0("  NOTE %s: TIMP1, %s and their TIMP-1 ratios are all candidates;\n",
                       "       a glm panel holding x, TIMP1 and x/TIMP1 together is not identifiable.\n"),
                label, paste(trip, collapse = ", ")))
}


# ============================================================
# 4) READ AND MERGE DATA
# ============================================================
DATA_DIR <- DATA_DIR_CHOICES[dir.exists(DATA_DIR_CHOICES)][1]
if (is.na(DATA_DIR)) stop("None of the folders in DATA_DIR_CHOICES exists on this computer.")
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/")
if (USE_ALL_BIOMARKERS && !USE_PUBLICATION_BIOMARKERS) {
  # The all-marker DMFT analysis needs no external selected-panel definitions.
  publication_markers_all <- character(0)
  prior_panels_reported <- data.frame()
} else {
  if (!file.exists(PANELS_FILE)) stop("Marker lists not found: ", PANELS_FILE, ". Set PANELS_FILE in section 1.")
  source(PANELS_FILE, local = TRUE)
}

PCR_DentalBehvaior_2025_12_17_1359 <- read_excel(
  file.path(DATA_DIR, "PCR_DentalBehvaior_2025-12-17_1359.xlsx"),
  sheet = "Trimmed to 118 Subjects", skip = 1)
MERGED_PCR_RESULTS_16JAN2026 <- read_excel(
  file.path(DATA_DIR, "MERGED_PCR_RESULTS_16JAN2026.xlsx"),
  sheet = "11Cytos & OTUs-Sorted2-166", skip = 1)
MedHx_2025_12_15_1210 <- read_excel(
  file.path(DATA_DIR, "MedHx_2025-12-15_1210.xlsx"),
  sheet = "Trimmed", skip = 1)
setwd(file.path(DATA_DIR, "Shared with Qi"))

data.df1_raw    <- MERGED_PCR_RESULTS_16JAN2026
cols_to_numeric <- c("Age_DOS", "BMI", "Teeth",
                     "IL6", "MMP9", "S100A8", "MMP8", "IL1B_IL1F2", "CCL3_MIP1a",
                     "BAFF", "Resistin", "IFNa", "TIMP1", "PGE2",
                     "Per_PlaqueIndex", "PISA", "AdjustPRA", "Clinical_Caries", "Cavitated_Caries",
                     grep("Bacteria", colnames(data.df1_raw), value = TRUE))
cols_to_numeric <- intersect(unique(cols_to_numeric), colnames(data.df1_raw))
data.df1_typed  <- data.df1_raw %>% mutate(across(all_of(cols_to_numeric), as.numeric))

data.df3_merged <- data.df1_typed %>%
  left_join(MedHx_2025_12_15_1210,              by = "StudyID") %>%
  left_join(PCR_DentalBehvaior_2025_12_17_1359, by = "StudyID")

# Analysed sample: subjects 47 and 68 removed, then the first 115 rows. Rows are
# kept by position, so this depends on the row order of the PCR sheet.
data.df4_filtered <- data.df3_merged[!data.df3_merged$Subject %in% c("47", "68"), ]
data.df4_filtered <- data.df4_filtered[1:115, ]
if (anyNA(data.df4_filtered$StudyID))
  warning("Some of the first 115 rows have no StudyID; check the row order of the PCR sheet.")


# ============================================================
# 5) TAXON NAMES, MARKER LISTS, CONFOUNDER NAMES
# ============================================================
# Taxon columns come from the PCR file only, so medical-history fields such as
# BacterialEndocarditis are never treated as bacteria.
bacteria_cols <- intersect(grep("^Bacteria", colnames(data.df1_raw), value = TRUE),
                           colnames(data.df4_filtered))
stray_cols <- setdiff(grep("^Bacteria", colnames(data.df4_filtered), value = TRUE), bacteria_cols)
if (length(stray_cols) > 0)
  cat("Not treated as bacteria:", paste(stray_cols, collapse = ", "), "\n")
cat("Bacteria columns:", length(bacteria_cols), "\n")

biomarkers_short_v2 <- sapply(bacteria_cols, shorten_bacteria_name_v2, USE.NAMES = FALSE)
data.df5_renamed    <- rename_bacteria_cols(data.df4_filtered, bacteria_cols, biomarkers_short_v2)

cytokines_full   <- c("IL6", "MMP9", "S100A8", "MMP8", "IL1B_IL1F2", "CCL3_MIP1a",
                      "BAFF", "Resistin", "IFNa", "TIMP1", "PGE2")
S100A8_MARKERS   <- c("S100A8", "S100A8_div_TIMP1")
excluded_markers <- if (EXCLUDE_S100A8) S100A8_MARKERS else character(0)
cytokines_full   <- setdiff(cytokines_full, excluded_markers)
if (length(excluded_markers) > 0)
  cat(sprintf("Excluded from every candidate list: %s\n", paste(excluded_markers, collapse = ", ")))
biomarkers_full  <- c(biomarkers_short_v2, cytokines_full)

colnames(data.df5_renamed) <- sapply(colnames(data.df5_renamed), clean_name, USE.NAMES = FALSE)
biomarkers_full <- sapply(biomarkers_full, clean_name, USE.NAMES = FALSE)
cytokines_full  <- sapply(cytokines_full,  clean_name, USE.NAMES = FALSE)

confounders <- if (length(CONFOUNDERS) == 0) character(0) else
  sapply(CONFOUNDERS, clean_name, USE.NAMES = FALSE)

biomarkers_full <- biomarkers_full[biomarkers_full %in% colnames(data.df5_renamed)]
cytokines_full  <- cytokines_full[cytokines_full   %in% colnames(data.df5_renamed)]

is_num <- vapply(biomarkers_full, function(v) is.numeric(data.df5_renamed[[v]]), logical(1))
if (!all(is_num))
  stop("Non-numeric biomarker column(s): ", paste(biomarkers_full[!is_num], collapse = ", "))

data.df6_cleaned <- data.df5_renamed


# ============================================================
# 6) THE PAPER'S MARKERS AND TIMP-1 RATIOS (validation)
# ============================================================
# The paper's 33 candidates (script1_panels.R), minus the excluded S100A8 markers
publication_markers <- publication_markers_all[!publication_markers_all %in% excluded_markers]

# Cytokine / TIMP-1 ratios, built from raw concentrations before imputation,
# log2 and z-scoring; from here on they are treated like measured markers. A
# ratio is NA when either value is NA or TIMP-1 is 0, and KNN fills the cell.
timp_ratio_numerators <- setdiff(c("BAFF", "IFNa", "PGE2", "IL1B_IL1F2", "IL6", "CCL3_MIP1a",
                                   "MMP8", "MMP9", "S100A8"), excluded_markers)
timp_ratios_full <- character(0)

if (USE_TIMP_RATIOS) {
  missing_cyt <- setdiff(c(timp_ratio_numerators, "TIMP1"), colnames(data.df6_cleaned))
  if (length(missing_cyt) > 0)
    stop("TIMP-1 ratios need column(s) not found: ", paste(missing_cyt, collapse = ", "))

  timp1_raw <- data.df6_cleaned$TIMP1
  n_zero    <- sum(timp1_raw == 0, na.rm = TRUE)
  if (n_zero > 0)
    warning(sprintf("TIMP1 is 0 for %d subject(s); their TIMP-1 ratios are set to NA and imputed.", n_zero))

  cat(sprintf("\n===== TIMP-1 ratios (raw scale): %d =====\n", length(timp_ratio_numerators)))
  cat(sprintf("  TIMP1: NA = %d | zero = %d\n", sum(is.na(timp1_raw)), n_zero))
  for (num in timp_ratio_numerators) {
    col_name <- paste0(num, "_div_TIMP1")
    r        <- data.df6_cleaned[[num]] / timp1_raw
    r[!is.finite(r)] <- NA_real_
    data.df6_cleaned[[col_name]] <- r
    timp_ratios_full <- c(timp_ratios_full, col_name)
    cat(sprintf("  %-22s NA (to be imputed) = %d\n", col_name, sum(is.na(r))))
  }
  biomarkers_full <- c(biomarkers_full, timp_ratios_full)
}


# ============================================================
# 7) SMOKING, DIABETES AND DMFT
# ============================================================
# Cigsmoke, T1Diabetes and T2Diabetes -> No / Yes factors; Diabetes = type 1
# OR type 2: Yes if either box is checked, No only when both are unchecked.
checkbox_vars <- c("Cigsmoke", "T1Diabetes", "T2Diabetes")
missing_cb    <- setdiff(checkbox_vars, colnames(data.df6_cleaned))
if (length(missing_cb) > 0)
  stop("Column(s) not found after the merge: ", paste(missing_cb, collapse = ", "))

for (v in checkbox_vars)
  data.df6_cleaned[[v]] <- recode_checkbox(data.df6_cleaned[[v]], v)

is_t1 <- data.df6_cleaned$T1Diabetes == "Yes"
is_t2 <- data.df6_cleaned$T2Diabetes == "Yes"
if ("Diabetes" %in% colnames(data.df6_cleaned) && !is.factor(data.df6_cleaned$Diabetes))
  warning("An existing 'Diabetes' column is overwritten by the type 1 / type 2 combination.")
data.df6_cleaned$Diabetes <- factor(ifelse(is_t1 | is_t2, "Yes", "No"), levels = c("No", "Yes"))

n_both <- sum(is_t1 & is_t2, na.rm = TRUE)
if (n_both > 0)
  warning(sprintf("%d subject(s) have BOTH type 1 and type 2 checked (counted once in Diabetes).", n_both))

cat("\n===== Smoking and diabetes (all subjects) =====\n")
for (v in c(checkbox_vars, "Diabetes")) {
  x <- data.df6_cleaned[[v]]
  cat(sprintf("  %-11s No = %3d | Yes = %3d | NA = %d\n", v,
              sum(x == "No", na.rm = TRUE), sum(x == "Yes", na.rm = TRUE), sum(is.na(x))))
}

# DMFT = sum of DMFT_COMPONENTS, NA when any of them is NA. Built for every
# project when the columns exist; the DMFT project needs them.
dmft_cols    <- sapply(DMFT_COMPONENTS, clean_name, USE.NAMES = FALSE)
dmft_missing <- setdiff(dmft_cols, colnames(data.df6_cleaned))
if (length(dmft_missing) == 0) {
  if ("DMFT" %in% colnames(data.df6_cleaned))
    warning("The data already has a DMFT column; it is replaced by the sum of DMFT_COMPONENTS.")
  dmft_parts <- do.call(cbind, lapply(dmft_cols, function(v) suppressWarnings(as.numeric(data.df6_cleaned[[v]]))))
  data.df6_cleaned$DMFT <- rowSums(dmft_parts)
  cat(sprintf("\n  DMFT = %s: NA = %d | %s\n", paste(dmft_cols, collapse = " + "),
              sum(is.na(data.df6_cleaned$DMFT)),
              paste(names(summary(data.df6_cleaned$DMFT))[1:6],
                    round(summary(data.df6_cleaned$DMFT)[1:6], 1), sep = " ", collapse = ", ")))
} else if (PROJECT == "DMFT") {
  similar <- grep("teeth|caries|carries|missing|restored|filled|decay", colnames(data.df6_cleaned),
                  ignore.case = TRUE, value = TRUE)
  stop("DMFT component(s) not found: ", paste(dmft_missing, collapse = ", "),
       ". Similar columns in the data: ", if (length(similar)) paste(similar, collapse = ", ") else "none",
       ". Set DMFT_COMPONENTS in section 1.")
} else {
  cat(sprintf("\n  DMFT not built (column(s) not found: %s)\n", paste(dmft_missing, collapse = ", ")))
}


# ============================================================
# 8) COLUMN ORDER
# ============================================================
# Demographics and clinical columns, cytokines, TIMP-1 ratios, bacteria (sorted), rest
key_cols <- c("StudyID", "Subject", "condL", "Age_DOS", "BMI", "Teeth", "Sex",
              "Cigsmoke", "Diabetes", "T1Diabetes", "T2Diabetes",
              "Per_4mm", "MeanPD", "Per_BOP", "MeanCAL", "Per_AlvBoneLoss",
              "AdjustPRA", "PRA", "Per_PlaqueIndex", "PISA", "Clinical_Caries", "Cavitated_Caries", "DMFT")
key_cols <- intersect(key_cols, colnames(data.df6_cleaned))

cytokines_present   <- intersect(cytokines_full,      colnames(data.df6_cleaned))
timp_ratios_present <- intersect(timp_ratios_full,    colnames(data.df6_cleaned))
bacteria_present    <- intersect(biomarkers_short_v2, colnames(data.df6_cleaned))
bacteria_sorted     <- sort(bacteria_present)

all_special    <- c(key_cols, cytokines_present, timp_ratios_present, bacteria_sorted)
remaining_cols <- setdiff(colnames(data.df6_cleaned), all_special)
data.df6_cleaned <- data.df6_cleaned[, c(key_cols, cytokines_present, timp_ratios_present,
                                         bacteria_sorted, remaining_cols)]


# ============================================================
# 9) KNN IMPUTATION (k = 5, raw scale)
# ============================================================
data.df7_imputed <- data.df6_cleaned
if (USE_PUBLICATION_BIOMARKERS) {
  # The paper's markers are imputed from each other only; the rest separately,
  # so data.df7_imputed is complete for other scripts.
  knn_pub   <- intersect(unname(publication_markers), biomarkers_full)
  knn_other <- setdiff(biomarkers_full, knn_pub)
  cat(sprintf("\n===== KNN imputation: the paper's %d markers among themselves, %d others separately =====\n",
              length(knn_pub), length(knn_other)))
  data.df7_imputed[, knn_pub] <-
    kNN(data.df6_cleaned[, knn_pub, drop = FALSE], k = 5, imp_var = FALSE)
  if (length(knn_other) > 0)
    data.df7_imputed[, knn_other] <-
      kNN(data.df6_cleaned[, knn_other, drop = FALSE], k = 5, imp_var = FALSE)
} else {
  data.df7_imputed[, biomarkers_full] <-
    kNN(data.df6_cleaned[, biomarkers_full, drop = FALSE], k = 5, imp_var = FALSE)
}


# ============================================================
# 10) AGE / BMI CATEGORIES AND CONFOUNDERS
# ============================================================
if (FACTORIZE_AGE_BMI) {
  if ("Age_DOS" %in% colnames(data.df7_imputed)) {
    age_lo <- paste0("<",  AGE_CUTOFF)
    age_hi <- paste0(">=", AGE_CUTOFF)
    data.df7_imputed$Age_DOS_cat <- factor(
      ifelse(data.df7_imputed$Age_DOS >= AGE_CUTOFF, age_hi, age_lo), levels = c(age_lo, age_hi))
  }
  if ("BMI" %in% colnames(data.df7_imputed)) {
    bmi_lo <- paste0("<",  BMI_CUTOFF)
    bmi_hi <- paste0(">=", BMI_CUTOFF)
    data.df7_imputed$BMI_cat <- factor(
      ifelse(data.df7_imputed$BMI >= BMI_CUTOFF, bmi_hi, bmi_lo), levels = c(bmi_lo, bmi_hi))
  }
}

# Preserve source sex codes/labels, but fit categorical effects rather than a
# numerical slope through the codes. Missing/blank values stay missing.
if (PROJECT == "DMFT" && "Sex" %in% names(data.df7_imputed)) {
  sx <- trimws(as.character(data.df7_imputed$Sex))
  sx[is.na(data.df7_imputed$Sex) | sx == ""] <- NA_character_
  data.df7_imputed$Sex <- factor(sx)
  cat("\nSex coding retained from source (verify these labels):\n")
  print(table(data.df7_imputed$Sex, useNA = "ifany"))
}
missing_conf <- setdiff(confounders, colnames(data.df7_imputed))
if (length(missing_conf) > 0)
  stop("Confounder(s) not found in the data: ", paste(missing_conf, collapse = ", "),
       if (!FACTORIZE_AGE_BMI && any(missing_conf %in% c("Age_DOS_cat", "BMI_cat")))
         " (Age_DOS_cat / BMI_cat need FACTORIZE_AGE_BMI = TRUE)" else "",
       ". Check CONFOUNDERS in section 1.")

# Variables that define the groups cannot also be adjusted for: diabetes in
# validation, DMFT and its components in the DMFT project
outcome_vars <- if (DM_COMPARISON) c("Diabetes", "T1Diabetes", "T2Diabetes") else
  if (PROJECT == "DMFT") c("DMFT", dmft_cols) else character(0)
dropped_conf <- intersect(confounders, outcome_vars)
if (length(dropped_conf) > 0) {
  warning("Removed from confounders (they define the groups): ", paste(dropped_conf, collapse = ", "))
  confounders <- setdiff(confounders, dropped_conf)
}
if (PROJECT == "DMFT" && "Teeth" %in% confounders)
  cat("  NOTE: Teeth (teeth present) mirrors Missing_Teeth, part of DMFT; adjusting for it",
      "removes most of the M in DMFT.\n")

cat("\n===== Confounders =====\n")
if (length(confounders) == 0) {
  cat("  (none)\n")
} else {
  for (cf in confounders) {
    cat(sprintf("  %s (%s)\n", cf, class(data.df7_imputed[[cf]])[1]))
    if (is.factor(data.df7_imputed[[cf]]))
      print(table(data.df7_imputed[[cf]], useNA = "ifany"))
    if (any(is.na(data.df7_imputed[[cf]])))
      cat(sprintf("    !! %d NA(s): those rows will be dropped by glm / rpart\n",
                  sum(is.na(data.df7_imputed[[cf]]))))
  }
}


# ============================================================
# 11) THE PAPER'S PANELS (validation)
# ============================================================
# publication_biomarkers: the paper's markers usable in this dataset.
# prior_panels: the entries of Tables 2-4 whose markers are all usable, as
# column names; prior_panel_check says why any other entry is not.
if (PROJECT == "validation") {
  pub_status <- marker_status(unname(publication_markers))
  publication_biomarkers <- unname(publication_markers)[pub_status == "ok"]

  cat(sprintf("\n===== The paper's markers: %d of %d usable%s =====\n",
              length(publication_biomarkers), length(publication_markers),
              if (EXCLUDE_S100A8) " (33 minus S100A8 and S100A8/TIMP-1)" else ""))
  for (s in setdiff(unique(pub_status), "ok"))
    cat(sprintf("  %s: %s\n", s,
                paste(sprintf("%s (%s)", names(publication_markers)[pub_status == s],
                              publication_markers[pub_status == s]), collapse = ", ")))
  if (length(publication_biomarkers) == 0) stop("None of the paper's markers is usable.")
  if (length(publication_biomarkers) < length(publication_markers))
    warning(sprintf("Only %d of %d of the paper's markers are usable; see the list above.",
                    length(publication_biomarkers), length(publication_markers)))

  prior_panels      <- list()
  prior_panel_check <- data.frame(Panel = character(0), k = integer(0), Paper = character(0),
                                  Current = character(0), Status = character(0),
                                  stringsAsFactors = FALSE)
  for (pn in names(prior_panels_paper)) {
    lab <- prior_panels_paper[[pn]]
    cur <- translate_paper(lab, paste0("prior_panels_paper$", pn))
    st  <- marker_status(cur)
    bad <- st != "ok"
    prior_panel_check <- rbind(prior_panel_check, data.frame(
      Panel   = pn, k = length(lab),
      Paper   = paste(lab, collapse = " + "),
      Current = paste(cur, collapse = " + "),
      Status  = if (any(bad)) paste(sprintf("%s: %s", lab[bad], st[bad]), collapse = "; ") else "ready",
      stringsAsFactors = FALSE))
    if (!any(bad)) prior_panels[[pn]] <- cur
  }
  cat(sprintf("  Tables 2-4: %d of %d entries ready to score\n", length(prior_panels), length(prior_panels_paper)))
  for (i in which(prior_panel_check$Status != "ready"))
    cat(sprintf("    NOT READY %-8s %s -> %s\n", prior_panel_check$Panel[i],
                prior_panel_check$Paper[i], prior_panel_check$Status[i]))
  cat(sprintf("  Table 4 CIs of accuracy: %d printed, %d usable (the 6-marker CIs look shifted one row)\n",
              sum(!is.na(prior_panels_reported$Accuracy_CI_Lower)),
              sum(prior_panels_reported$Accuracy_CI_usable)))
} else {
  publication_biomarkers <- character(0)
  prior_panels           <- list()
  prior_panel_check      <- data.frame(Panel = character(0), k = integer(0), Paper = character(0),
                                       Current = character(0), Status = character(0),
                                       stringsAsFactors = FALSE)
}


# ============================================================
# 12) CANDIDATE MARKERS
# ============================================================
cat(sprintf("\n===== Candidate markers (%s, %s) =====\n", PROJECT, CONDL_SOURCE))

if (USE_PUBLICATION_BIOMARKERS) {
  biomarkers_for_PH <- publication_biomarkers
  biomarkers_for_GH <- NULL
  biomarkers_for_GP <- NULL
  cat(sprintf("  The paper's markers: %d\n", length(biomarkers_for_PH)))
  report_triples(biomarkers_for_PH, "PH")
} else if (USE_ALL_BIOMARKERS) {
  biomarkers_for_PH <- biomarkers_full
  biomarkers_for_GH <- if (SINGLE_COMPARISON) NULL else biomarkers_full
  biomarkers_for_GP <- if (SINGLE_COMPARISON) NULL else biomarkers_full
  cat(sprintf("  All biomarkers: %d\n", length(biomarkers_full)))
} else {
  lists <- switch(PROJECT,
                  caries = list(PH = selected_markers$Caries[[CARIES_SOURCE_COL]]),
                  DMFT   = selected_markers$DMFT,
                  perio  = selected_markers[[PERIO_DEFINITION]])
  if (is.null(lists$PH))
    stop("No selected markers for ", PROJECT, " yet: set ", PROJECT, "_MARKERS <- \"all\", or paste a ",
         "screened list into script1_panels.R (biomarkers_DMFT).")
  biomarkers_for_PH <- select_markers(lists$PH, "PH")
  biomarkers_for_GH <- select_markers(lists$GH, "GH")
  biomarkers_for_GP <- select_markers(lists$GP, "GP")
  cat(sprintf("  Selected markers (script1_panels.R): PH %d%s\n", length(biomarkers_for_PH),
              if (SINGLE_COMPARISON) "" else
                sprintf(" | GH %d | GP %d", length(biomarkers_for_GH), length(biomarkers_for_GP))))
}

overlap <- intersect(confounders, c(biomarkers_for_PH, biomarkers_for_GH, biomarkers_for_GP))
if (length(overlap) > 0)
  warning("Confounder(s) also among the candidate markers: ", paste(overlap, collapse = ", "),
          ". They will enter the model formula twice.")

# Markers that are log2-transformed (all comparisons together)
biomarkers_all_panels <- unique(c(biomarkers_for_PH, biomarkers_for_GH, biomarkers_for_GP))
biomarkers_all_panels <- biomarkers_all_panels[biomarkers_all_panels %in% colnames(data.df7_imputed)]

# The paper's panel markers that are not candidates are still log2-transformed
# and z-scored, so every entry of Tables 2-4 can be scored on data.dfPH.
prior_extra_PH <- if (PROJECT == "validation") setdiff(unique(unlist(prior_panels)), biomarkers_for_PH) else character(0)
if (length(prior_extra_PH) > 0) {
  cat(sprintf("  The paper's panel markers transformed but not candidates: %s\n",
              paste(prior_extra_PH, collapse = ", ")))
  biomarkers_all_panels <- unique(c(biomarkers_all_panels, prior_extra_PH))
}

# Cytokines and TIMP-1 ratios; kept out of the bacterial ratio features
nonbacterial_full <- c(cytokines_full, timp_ratios_full)
cytokines_PH <- intersect(nonbacterial_full, biomarkers_for_PH)
cytokines_GH <- if (!is.null(biomarkers_for_GH)) intersect(nonbacterial_full, biomarkers_for_GH) else character(0)
cytokines_GP <- if (!is.null(biomarkers_for_GP)) intersect(nonbacterial_full, biomarkers_for_GP) else character(0)


# ============================================================
# 13) LOG2 TRANSFORM (pseudo-count 5e-7)
# ============================================================
data.df8_log2 <- data.df7_imputed
data.df8_log2[, biomarkers_all_panels] <-
  data.df8_log2[, biomarkers_all_panels, drop = FALSE] %>%
  mutate(across(everything(), ~ . + (1e-06 / 2)))
data.df8_log2[, biomarkers_all_panels] <-
  lapply(data.df8_log2[, biomarkers_all_panels, drop = FALSE], log2)


# ============================================================
# 14) OUTCOME LABELS (condL)
# ============================================================
# Subjects outside every group get NA and are dropped at the end of this section.
cat(sprintf("\n===== Outcome: %s =====\n", CONDL_SOURCE))
CARIES_N_GAP_DROPPED <- NULL
DMFT_CUTOFF_USED     <- NULL

if (PROJECT == "perio") {

  if (PERIO_DEFINITION == "original") {
    cat("  Original diagnosis (condL) used as is.\n")
  } else {
    if (!PERIO_DEFINITION %in% colnames(data.df8_log2))
      stop(sprintf("Column '%s' not found.", PERIO_DEFINITION))
    score_vec <- data.df8_log2[[PERIO_DEFINITION]]
    low  <- if (PERIO_DEFINITION == "AdjustPRA") ADJUSTPRA_HEALTHY_MAX else PISA_HEALTHY_MAX
    high <- if (PERIO_DEFINITION == "AdjustPRA") ADJUSTPRA_MOD_MAX     else PISA_MOD_MAX
    data.df8_log2$condL_original <- data.df8_log2$condL
    data.df8_log2$condL          <- recode_from_score(score_vec, low, high)
    cat(sprintf("  %s: < %g Healthy | %g-%g Gingivitis | > %g Periodontitis\n",
                PERIO_DEFINITION, low, low, high, high))
    cat("  Score summary:\n"); print(summary(as.numeric(score_vec)))
  }

} else if (PROJECT == "caries") {

  if (!CARIES_SOURCE_COL %in% colnames(data.df8_log2))
    stop(sprintf("Column '%s' not found.", CARIES_SOURCE_COL))
  score_vec <- suppressWarnings(as.numeric(data.df8_log2[[CARIES_SOURCE_COL]]))
  data.df8_log2$condL_original <- data.df8_log2$condL

  if (CARIES_SOURCE_COL == "Clinical_Caries") {
    data.df8_log2$condL <- ifelse(is.na(score_vec), NA_character_,
                                  ifelse(score_vec <= CARIES_BINARY_CUTOFF, "Low_Caries_Burden",
                                         "High_Caries_Burden"))
    cat(sprintf("  Clinical_Caries: <= %g Low_Caries_Burden | > %g High_Caries_Burden\n",
                CARIES_BINARY_CUTOFF, CARIES_BINARY_CUTOFF))
  } else {
    data.df8_log2$condL <- ifelse(is.na(score_vec), NA_character_,
                                  ifelse(score_vec <= CARIES_NOCARIES_MAX, "Caries_Free",
                                         ifelse(score_vec >= CARIES_HIGH_MIN, "Caries_Present",
                                                NA_character_)))
    CARIES_N_GAP_DROPPED <- sum(!is.na(score_vec) &
                                  score_vec > CARIES_NOCARIES_MAX & score_vec < CARIES_HIGH_MIN)
    cat(sprintf("  Cavitated_Caries: <= %g Caries_Free | >= %g Caries_Present | in between dropped: %d\n",
                CARIES_NOCARIES_MAX, CARIES_HIGH_MIN, CARIES_N_GAP_DROPPED))
  }
  cat("  Score summary:\n"); print(summary(score_vec))
  cat("\n  Raw value counts:\n"); print(table(score_vec, useNA = "ifany"))

} else if (PROJECT == "DMFT") {

  dmft <- data.df8_log2$DMFT
  data.df8_log2$condL_original <- data.df8_log2$condL
  # the subjects analysed: DMFT known and, with a filter, the kept diagnoses
  in_pool <- !is.na(dmft) & (is.null(DMFT_PERIODONTAL_FILTER) |
                               data.df8_log2$condL_original %in% DMFT_PERIODONTAL_FILTER)
  DMFT_CUTOFF_USED <- if (is.null(DMFT_CUTOFF)) median(dmft[in_pool]) else DMFT_CUTOFF
  data.df8_log2$condL <- ifelse(is.na(dmft), NA_character_,
                                ifelse(dmft <= DMFT_CUTOFF_USED, "Low_DMFT", "High_DMFT"))
  cat(sprintf("  DMFT: <= %g Low_DMFT | > %g High_DMFT (cut-off: %s)\n", DMFT_CUTOFF_USED, DMFT_CUTOFF_USED,
              if (is.null(DMFT_CUTOFF)) "median of the analysed subjects" else "DMFT_CUTOFF"))
  cat("  DMFT summary (analysed subjects):\n"); print(summary(dmft[in_pool]))
  cat("\n  DMFT value counts (analysed subjects):\n"); print(table(dmft[in_pool]))

} else {   # validation

  data.df8_log2$condL_original <- data.df8_log2$condL

  if (DM_PERIO_DEFINITION == "original") {
    perio_class <- as.character(data.df8_log2$condL)
  } else {
    if (!DM_PERIO_DEFINITION %in% colnames(data.df8_log2))
      stop(sprintf("Column '%s' not found.", DM_PERIO_DEFINITION))
    low  <- if (DM_PERIO_DEFINITION == "AdjustPRA") ADJUSTPRA_HEALTHY_MAX else PISA_HEALTHY_MAX
    high <- if (DM_PERIO_DEFINITION == "AdjustPRA") ADJUSTPRA_MOD_MAX     else PISA_MOD_MAX
    perio_class <- recode_from_score(data.df8_log2[[DM_PERIO_DEFINITION]], low, high)
  }
  data.df8_log2$condL_perio <- perio_class

  # "With diabetes" follows DM_COL; "without" means no diabetes of either type.
  # Subjects with only the other type are in neither.
  dm_with    <- data.df8_log2[[DM_COL]] %in% "Yes"
  dm_without <- data.df8_log2$Diabetes  %in% "No"
  dm_unknown <- is.na(data.df8_log2[[DM_COL]]) | is.na(data.df8_log2$Diabetes)
  dm_status  <- ifelse(dm_with, "With", ifelse(dm_without, "Without",
                                               ifelse(dm_unknown, NA_character_, "Other type only")))

  if (CONDL_SOURCE == "PanelValidation") {

    # Every periodontitis case plus the controls. Perio_DM_group:
    #   "With" (DM_COL = Yes) | "Without" (no diabetes) |
    #   "Other type" (only the other type) | "Unknown" (status missing)
    in_case <- perio_class %in% "Periodontitis"
    in_ctrl <- perio_class %in% "Healthy periodontium"
    is_ctrl <- if (VALIDATION_CONTROLS == "Healthy_All") in_ctrl else in_ctrl & dm_without

    data.df8_log2$condL <- ifelse(in_case, "Periodontitis",
                                  ifelse(is_ctrl, "Healthy periodontium", NA_character_))
    data.df8_log2$Perio_DM_group <- ifelse(
      !in_case, NA_character_,
      ifelse(dm_with, "With", ifelse(dm_without, "Without",
                                     ifelse(dm_unknown, "Unknown", "Other type"))))

    cat(sprintf("  Cases   : every periodontitis case (%s)\n", DM_PERIO_DEFINITION))
    cat(sprintf("  Controls: healthy periodontium (%s), %s\n", DM_PERIO_DEFINITION,
                if (VALIDATION_CONTROLS == "Healthy_All") "with or without diabetes"
                else "no diabetes of either type (as in the paper)"))
    cat(sprintf(paste0("  Periodontitis: %d (with %s: %d | without diabetes: %d | other type only: %d | ",
                       "unknown: %d) | controls: %d\n"),
                sum(in_case), DM_COL, sum(in_case & dm_with), sum(in_case & dm_without),
                sum(in_case & dm_status %in% "Other type only"), sum(in_case & dm_unknown), sum(is_ctrl)))
    if (VALIDATION_CONTROLS == "Healthy_NoDiabetes")
      cat(sprintf("  Healthy periodontium left out: with diabetes %d | diabetes status unknown %d\n",
                  sum(in_ctrl & !dm_without & !dm_unknown), sum(in_ctrl & dm_unknown)))
    t1 <- data.df8_log2$T1Diabetes %in% "Yes"
    t2 <- data.df8_log2$T2Diabetes %in% "Yes"
    cat(sprintf(paste0("  Periodontitis with diabetes of either type: %d (type 2 only %d | type 1 only %d | ",
                       "both boxes checked %d)\n"),
                sum(in_case & (t1 | t2)), sum(in_case & t2 & !t1), sum(in_case & t1 & !t2),
                sum(in_case & t1 & t2)))

  } else {   # PaperDesign

    in_case <- perio_class %in% "Periodontitis"
    in_ctrl <- perio_class %in% "Healthy periodontium"
    is_case <- in_case & dm_with
    is_ctrl <- in_ctrl & dm_without

    data.df8_log2$condL <- ifelse(is_case, "Perio_Diabetes",
                                  ifelse(is_ctrl, "Healthy_NoDiabetes", NA_character_))

    cat(sprintf("  Cases   : periodontitis (%s) with %s = Yes: %d\n", DM_PERIO_DEFINITION, DM_COL, sum(is_case)))
    cat(sprintf("  Controls: healthy periodontium (%s), no diabetes of either type: %d\n",
                DM_PERIO_DEFINITION, sum(is_ctrl)))
    cat(sprintf(paste0("  Left out: periodontitis without diabetes %d | periodontitis with the other type only %d | ",
                       "healthy with diabetes %d | gingivitis or unclassified %d | unknown diabetes status %d\n"),
                sum(in_case & dm_without), sum(in_case & dm_status %in% "Other type only"),
                sum(in_ctrl & !dm_without & !dm_unknown), sum(!in_case & !in_ctrl),
                sum((in_case | in_ctrl) & dm_unknown)))
    if (sum((in_case | in_ctrl) & dm_unknown) > 0)
      warning(sprintf("%d subject(s) left out because diabetes status is missing.",
                      sum((in_case | in_ctrl) & dm_unknown)))
  }

  cat("\n  Periodontal class x diabetes status (all subjects):\n")
  print(table(Perio = perio_class, Diabetes = dm_status, useNA = "ifany"))
  cat("\n  Diabetes type among periodontitis cases:\n")
  print(table(T1Diabetes = data.df8_log2$T1Diabetes[in_case],
              T2Diabetes = data.df8_log2$T2Diabetes[in_case], useNA = "ifany"))
}

n_before      <- nrow(data.df8_log2)
data.df8_log2 <- data.df8_log2[!is.na(data.df8_log2$condL), ]
n_dropped     <- n_before - nrow(data.df8_log2)
if (n_dropped > 0) cat(sprintf("  Dropped %d row(s) with NA condL.\n", n_dropped))

cat(sprintf("\n  condL (n = %d):\n", nrow(data.df8_log2)))
print(table(data.df8_log2$condL, useNA = "ifany"))
if (CONDL_SOURCE != "original" && "condL_original" %in% colnames(data.df8_log2)) {
  cat("\n  Original diagnosis vs new condL:\n")
  print(table(Original = data.df8_log2$condL_original, New = data.df8_log2$condL, useNA = "ifany"))
}


# ============================================================
# 15) COMPARISON DATASETS AND SHORT LABELS
# ============================================================
# Case / reference labels:
#   perio      data.dfPH P / H, data.dfGH G / H, data.dfGP P / G
#   caries     data.dfPH CP / CF (Cavitated) or HC / LC (Clinical)
#   DMFT       data.dfPH HD / LD
#   validation data.dfPH P / H
# data.dfGH and data.dfGP are NULL except for perio.
if (PROJECT == "perio") {

  data.dfGH <- data.df8_log2 %>% filter(condL %in% c("Gingivitis", "Healthy periodontium")) %>%
    mutate(condL = ifelse(condL == "Gingivitis", "G", "H"))
  data.dfPH <- data.df8_log2 %>% filter(condL %in% c("Periodontitis", "Healthy periodontium")) %>%
    mutate(condL = ifelse(condL == "Periodontitis", "P", "H"))
  data.dfGP <- data.df8_log2 %>% filter(condL %in% c("Gingivitis", "Periodontitis")) %>%
    mutate(condL = ifelse(condL == "Gingivitis", "G", "P"))

} else if (PROJECT == "caries") {

  # Restricted by the original periodontal diagnosis, so periodontitis-associated
  # bacteria do not drive the caries comparison
  data.dfPH_all_caries <- data.df8_log2[!is.na(data.df8_log2$condL), ]
  if (!is.null(CARIES_PERIODONTAL_FILTER)) {
    data.dfPH <- data.dfPH_all_caries %>% filter(condL_original %in% CARIES_PERIODONTAL_FILTER)
    cat(sprintf("\n  Original diagnosis kept: %s (removed %d)\n",
                paste(CARIES_PERIODONTAL_FILTER, collapse = ", "),
                nrow(data.dfPH_all_caries) - nrow(data.dfPH)))
    print(table(data.dfPH$condL_original, useNA = "ifany"))
  } else {
    data.dfPH <- data.dfPH_all_caries
    cat("\n  No periodontal restriction.\n")
  }
  data.dfPH <- if (CARIES_SOURCE_COL == "Clinical_Caries") {
    data.dfPH %>% mutate(condL = ifelse(condL == "High_Caries_Burden", "HC", "LC"))
  } else {
    data.dfPH %>% mutate(condL = ifelse(condL == "Caries_Present", "CP", "CF"))
  }
  data.dfGH <- NULL
  data.dfGP <- NULL

} else if (PROJECT == "DMFT") {

  data.dfPH <- data.df8_log2
  if (!is.null(DMFT_PERIODONTAL_FILTER)) {
    data.dfPH <- data.dfPH %>% filter(condL_original %in% DMFT_PERIODONTAL_FILTER)
    cat(sprintf("\n  Original diagnosis kept: %s (removed %d)\n",
                paste(DMFT_PERIODONTAL_FILTER, collapse = ", "), nrow(data.df8_log2) - nrow(data.dfPH)))
  }
  data.dfPH <- data.dfPH %>% mutate(condL = ifelse(condL == "High_DMFT", "HD", "LD"))
  data.dfGH <- NULL
  data.dfGP <- NULL
  n_hd <- sum(data.dfPH$condL == "HD")
  n_ld <- sum(data.dfPH$condL == "LD")
  cat(sprintf("\n  DMFT: HD = %d | LD = %d\n", n_hd, n_ld))
  if (min(n_hd, n_ld) == 0)
    stop("One DMFT group is empty (HD = ", n_hd, ", LD = ", n_ld, "). Change DMFT_CUTOFF.")
  if (min(n_hd, n_ld) < 20)
    warning(sprintf("Smaller DMFT group has only %d subjects; 3-6 marker logistic models are unstable at this size.",
                    min(n_hd, n_ld)))

} else if (CONDL_SOURCE == "PaperDesign") {

  data.dfPH <- data.df8_log2
  data.dfGH <- NULL
  data.dfGP <- NULL
  n_case <- sum(data.dfPH$condL == "Perio_Diabetes")
  n_ctrl <- sum(data.dfPH$condL == "Healthy_NoDiabetes")
  cat(sprintf("\n  PaperDesign: P = %d, H = %d (the paper: 28 vs 28)\n", n_case, n_ctrl))
  if (min(n_case, n_ctrl) == 0)
    stop("One PaperDesign group is empty (P = ", n_case, ", H = ", n_ctrl,
         "). Check DM_COL and DM_PERIO_DEFINITION.")
  if (min(n_case, n_ctrl) < 20)
    warning(sprintf(paste0("Smaller PaperDesign group has only %d subjects. 5-fold CV puts ~%d of them ",
                           "in each test fold, and 3-6 marker logistic models are unstable at this size."),
                    min(n_case, n_ctrl), floor(min(n_case, n_ctrl) / 5)))
  data.dfPH <- data.dfPH %>% mutate(condL = ifelse(condL == "Perio_Diabetes", "P", "H"))

} else {   # PanelValidation: the groups are cut from data.dfPH in section 19

  data.dfPH <- data.df8_log2
  data.dfGH <- NULL
  data.dfGP <- NULL
  if (sum(data.dfPH$condL == "Periodontitis") == 0 || sum(data.dfPH$condL == "Healthy periodontium") == 0)
    stop("PanelValidation: no periodontitis cases or no controls. Check DM_PERIO_DEFINITION, ",
         "VALIDATION_CONTROLS and the diabetes columns.")
  data.dfPH <- data.dfPH %>% mutate(condL = ifelse(condL == "Periodontitis", "P", "H"))
}

PH_CASE_LABEL <- switch(PROJECT, caries = if (CARIES_SOURCE_COL == "Clinical_Caries") "HC" else "CP",
                        DMFT = "HD", "P")
PH_REF_LABEL  <- switch(PROJECT, caries = if (CARIES_SOURCE_COL == "Clinical_Caries") "LC" else "CF",
                        DMFT = "LD", "H")

# condL -> factor, reference level first (see as_outcome)
if (PROJECT == "perio") {
  data.dfPH <- as_outcome(data.dfPH, "H", "P")
  data.dfGH <- as_outcome(data.dfGH, "H", "G")
  data.dfGP <- as_outcome(data.dfGP, "G", "P")
} else {
  data.dfPH <- as_outcome(data.dfPH, PH_REF_LABEL, PH_CASE_LABEL)
}


# ============================================================
# 16) Z-SCORES (against the reference group)
# ============================================================
# Biomarkers only; confounders keep their units, so their coefficients are not
# on the biomarkers' scale.
if (SINGLE_COMPARISON) {
  result_PH <- standardize_by_reference(data.dfPH, c(biomarkers_for_PH, prior_extra_PH),
                                        ref_group_label = PH_REF_LABEL)
  data.dfPH <- result_PH$data
  stats_PH  <- result_PH$stats
  stats_GH  <- NULL
  stats_GP  <- NULL
} else {
  result_GH <- standardize_by_reference(data.dfGH, biomarkers_for_GH, ref_group_label = "H")
  data.dfGH <- result_GH$data
  stats_GH  <- result_GH$stats

  result_PH <- standardize_by_reference(data.dfPH, c(biomarkers_for_PH, prior_extra_PH),
                                        ref_group_label = "H")
  data.dfPH <- result_PH$data
  stats_PH  <- result_PH$stats

  result_GP <- standardize_by_reference(data.dfGP, biomarkers_for_GP, ref_group_label = "G")
  data.dfGP <- result_GP$data
  stats_GP  <- result_GP$stats
}
cat(sprintf("\n  z-scored against %s\n", if (SINGLE_COMPARISON) PH_REF_LABEL else "H (PH, GH) and G (GP)"))


# ============================================================
# 17) BACTERIAL RATIO FEATURES
# ============================================================
# Built from raw (imputed) values of the candidate bacteria only.
bacteria_biomarkers_PH <- setdiff(biomarkers_for_PH, cytokines_PH)
bacteria_biomarkers_GH <- if (!is.null(biomarkers_for_GH)) setdiff(biomarkers_for_GH, cytokines_GH) else character(0)
bacteria_biomarkers_GP <- if (!is.null(biomarkers_for_GP)) setdiff(biomarkers_for_GP, cytokines_GP) else character(0)

ratio_cols_PH <- character(0); pair_results_PH <- NULL
ratio_cols_GH <- character(0); pair_results_GH <- NULL
ratio_cols_GP <- character(0); pair_results_GP <- NULL

if (USE_RATIO_FEATURES) {
  cat(sprintf("\n===== Bacterial ratio features (top %d) =====\n", TOP_N_RATIOS))

  result_ratios_PH <- build_ratio_features(data.dfPH, data.df7_imputed, bacteria_biomarkers_PH,
                                           PH_CASE_LABEL, PH_REF_LABEL, TOP_N_RATIOS)
  data.dfPH       <- result_ratios_PH$data_df
  ratio_cols_PH   <- result_ratios_PH$ratio_cols
  pair_results_PH <- result_ratios_PH$pair_results

  if (!SINGLE_COMPARISON) {
    result_ratios_GH <- build_ratio_features(data.dfGH, data.df7_imputed, bacteria_biomarkers_GH,
                                             "G", "H", TOP_N_RATIOS)
    data.dfGH       <- result_ratios_GH$data_df
    ratio_cols_GH   <- result_ratios_GH$ratio_cols
    pair_results_GH <- result_ratios_GH$pair_results

    result_ratios_GP <- build_ratio_features(data.dfGP, data.df7_imputed, bacteria_biomarkers_GP,
                                             "P", "G", TOP_N_RATIOS)
    data.dfGP       <- result_ratios_GP$data_df
    ratio_cols_GP   <- result_ratios_GP$ratio_cols
    pair_results_GP <- result_ratios_GP$pair_results
  }
}


# ============================================================
# 18) VIF PRUNING
# ============================================================
# Biomarkers and ratio features only; confounder-biomarker collinearity is not screened.
biomarkers_PH_full <- c(biomarkers_for_PH, ratio_cols_PH)
biomarkers_GH_full <- if (SINGLE_COMPARISON) character(0) else c(biomarkers_for_GH, ratio_cols_GH)
biomarkers_GP_full <- if (SINGLE_COMPARISON) character(0) else c(biomarkers_for_GP, ratio_cols_GP)

if (USE_VIF) {
  cat(sprintf("\n===== VIF pruning (threshold %g) =====\n", VIF_THRESHOLD))
  cat("  PH\n")
  biomarkers_PH_vif <- run_vif_prune(data.dfPH, biomarkers_PH_full, Threshold = VIF_THRESHOLD)
  if (!SINGLE_COMPARISON) {
    cat("  GH\n")
    biomarkers_GH_vif <- run_vif_prune(data.dfGH, biomarkers_GH_full, Threshold = VIF_THRESHOLD)
    cat("  GP\n")
    biomarkers_GP_vif <- run_vif_prune(data.dfGP, biomarkers_GP_full, Threshold = VIF_THRESHOLD)
  } else {
    biomarkers_GH_vif <- character(0)
    biomarkers_GP_vif <- character(0)
  }
} else {
  biomarkers_PH_vif <- safe_cols(data.dfPH, biomarkers_PH_full)
  biomarkers_GH_vif <- if (!SINGLE_COMPARISON) safe_cols(data.dfGH, biomarkers_GH_full) else character(0)
  biomarkers_GP_vif <- if (!SINGLE_COMPARISON) safe_cols(data.dfGP, biomarkers_GP_full) else character(0)
}

stats_PH_vif <- stats_PH[stats_PH$Column %in% biomarkers_PH_vif, , drop = FALSE]
stats_GH_vif <- if (!is.null(stats_GH)) stats_GH[stats_GH$Column %in% biomarkers_GH_vif, , drop = FALSE] else NULL
stats_GP_vif <- if (!is.null(stats_GP)) stats_GP[stats_GP$Column %in% biomarkers_GP_vif, , drop = FALSE] else NULL


# ============================================================
# 19) VALIDATION GROUPS (PanelValidation)
# ============================================================
# Cut from data.dfPH, which is already z-scored to the controls. All groups
# share those controls, so one standardisation applies to all. condL is
# P (case) / H (control) in each.
validation_sets         <- NULL
validation_group_counts <- NULL
if (CONDL_SOURCE == "PanelValidation") {
  is_H <- data.dfPH$condL == "H"
  is_P <- data.dfPH$condL == "P"
  dm_type_label <- c(T2Diabetes = "type 2 diabetes", T1Diabetes = "type 1 diabetes",
                     Diabetes   = "diabetes (type 1 or 2)")[[DM_COL]]
  add_any_dm <- VALIDATION_ANY_DIABETES_GROUP && DM_COL != "Diabetes"
  validation_group_labels <- c(
    All_periodontitis               = "All periodontitis",
    Periodontitis_with_diabetes     = paste("Periodontitis with", dm_type_label),
    Periodontitis_with_any_diabetes = "Periodontitis with diabetes (type 1 or 2)",
    Periodontitis_without_diabetes  = "Periodontitis without diabetes"
  )
  validation_sets <- list(
    All_periodontitis               = data.dfPH[is_H | is_P, ],
    Periodontitis_with_diabetes     = data.dfPH[is_H | (is_P & data.dfPH$Perio_DM_group %in% "With"), ],
    Periodontitis_with_any_diabetes = if (add_any_dm) data.dfPH[is_H | (is_P & data.dfPH$Diabetes %in% "Yes"), ],
    Periodontitis_without_diabetes  = data.dfPH[is_H | (is_P & data.dfPH$Perio_DM_group %in% "Without"), ]
  )
  validation_sets <- validation_sets[!vapply(validation_sets, is.null, logical(1))]
  validation_group_counts <- data.frame(
    Group = names(validation_sets),
    Label = unname(validation_group_labels[names(validation_sets)]),
    n_P   = vapply(validation_sets, function(d) sum(d$condL == "P"), integer(1)),
    n_H   = vapply(validation_sets, function(d) sum(d$condL == "H"), integer(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  validation_group_counts$Total <- validation_group_counts$n_P + validation_group_counts$n_H

  cat(sprintf("\n===== Validation groups (controls: %s) =====\n", VALIDATION_CONTROLS))
  for (i in seq_len(nrow(validation_group_counts)))
    cat(sprintf("  %-42s P = %3d | H = %3d | total = %3d\n", validation_group_counts$Label[i],
                validation_group_counts$n_P[i], validation_group_counts$n_H[i], validation_group_counts$Total[i]))
  n_oth <- sum(is_P & data.dfPH$Perio_DM_group %in% "Other type")
  n_unk <- sum(is_P & data.dfPH$Perio_DM_group %in% "Unknown")
  if (n_oth > 0)
    cat(sprintf("  %d case(s) with only the other diabetes type are in 'All periodontitis'%s.\n",
                n_oth, if (add_any_dm) " and the type 1 or 2 group" else " only"))
  if (n_unk > 0)
    cat(sprintf("  %d case(s) with unknown diabetes status are in 'All periodontitis' only.\n", n_unk))
  if (add_any_dm && n_oth == 0)
    cat("  No case has only type 1 diabetes, so the type 1 or 2 group equals the type 2 group.\n")
  for (i in which(validation_group_counts$n_P < 20))
    warning(sprintf(paste0("%s has only %d cases: 5-fold CV puts ~%d in each test fold, and 3-6 marker ",
                           "logistic models are unstable at this size."),
                    validation_group_counts$Label[i], validation_group_counts$n_P[i],
                    floor(validation_group_counts$n_P[i] / 5)))
  if (any(validation_group_counts$n_P == 0))
    stop("A validation group has no cases: ",
         paste(validation_group_counts$Label[validation_group_counts$n_P == 0], collapse = ", "))
}


# ============================================================
# 20) SUMMARY
# ============================================================
cat("\n===== Summary =====\n")
cat(sprintf("  PROJECT = %s | condL from %s | subjects after labelling: %d\n",
            PROJECT, CONDL_SOURCE, nrow(data.df8_log2)))
if (PROJECT == "perio") {
  cat(sprintf("  PH: %d | GH: %d | GP: %d\n", nrow(data.dfPH), nrow(data.dfGH), nrow(data.dfGP)))
} else {
  cat(sprintf("  data.dfPH: %s = %d (case) | %s = %d (reference)\n",
              PH_CASE_LABEL, sum(data.dfPH$condL == PH_CASE_LABEL),
              PH_REF_LABEL,  sum(data.dfPH$condL == PH_REF_LABEL)))
}
if (PROJECT == "DMFT")
  cat(sprintf("  DMFT cut-off: %g (%s)%s\n", DMFT_CUTOFF_USED,
              if (is.null(DMFT_CUTOFF)) "median" else "DMFT_CUTOFF",
              if (is.null(DMFT_PERIODONTAL_FILTER)) "" else
                paste(" | original diagnosis kept:", paste(DMFT_PERIODONTAL_FILTER, collapse = ", "))))
if (PROJECT == "validation")
  cat(sprintf("  DM_COL = %s | DM_PERIO_DEFINITION = %s%s\n", DM_COL, DM_PERIO_DEFINITION,
              if (CONDL_SOURCE == "PanelValidation") paste(" | controls:", VALIDATION_CONTROLS) else ""))
cat(sprintf("  Confounders: %s\n", if (length(confounders) == 0) "none" else paste(confounders, collapse = ", ")))
cat(sprintf("  Markers: %s | TIMP-1 ratios: %d | bacterial ratios: %s | VIF: %s | S100A8: %s\n",
            if (USE_PUBLICATION_BIOMARKERS) "the paper's" else if (USE_ALL_BIOMARKERS) "all" else "selected",
            length(timp_ratios_full), if (USE_RATIO_FEATURES) "on" else "off",
            if (USE_VIF) sprintf("on (%g)", VIF_THRESHOLD) else "off",
            if (EXCLUDE_S100A8) "excluded" else "kept"))
cat(sprintf("  Final features: PH %d%s\n", length(biomarkers_PH_vif),
            if (SINGLE_COMPARISON) "" else
              sprintf(" | GH %d | GP %d", length(biomarkers_GH_vif), length(biomarkers_GP_vif))))


# ============================================================
# 21) BIOMARKER REGISTRY
# ============================================================
# One entry per setting combination; the key names every setting that changes
# the data, so different runs do not overwrite each other. "_panel" in a key
# means selected markers.
if (!exists("biomarker_registry")) biomarker_registry <- list()

caries_periodontal_key <- if (CONDL_SOURCE == "Caries" && !is.null(CARIES_PERIODONTAL_FILTER)) {
  paste0("_periodontal_", paste(gsub("[^A-Za-z0-9]+", "", CARIES_PERIODONTAL_FILTER), collapse = "Plus"))
} else ""

caries_cutoff_key <- if (CONDL_SOURCE == "Caries") {
  if (CARIES_SOURCE_COL == "Clinical_Caries")
    paste0("_", CARIES_SOURCE_COL, "_cutoff", CARIES_BINARY_CUTOFF)
  else
    paste0("_", CARIES_SOURCE_COL, "_no", CARIES_NOCARIES_MAX, "_high", CARIES_HIGH_MIN)
} else ""

dmft_key <- if (CONDL_SOURCE == "DMFT") {
  paste0("_cutoff", DMFT_CUTOFF_USED, if (is.null(DMFT_CUTOFF)) "median" else "",
         if (!is.null(DMFT_PERIODONTAL_FILTER))
           paste0("_periodontal_", paste(gsub("[^A-Za-z0-9]+", "", DMFT_PERIODONTAL_FILTER), collapse = "Plus"))
         else "")
} else ""

confounder_key <- if (length(confounders) == 0) {
  "_noConfounders"
} else {
  paste0("_conf_", paste(gsub("[^A-Za-z0-9]+", "", confounders), collapse = "_"))
}

dm_key <- if (CONDL_SOURCE == "PaperDesign") {
  paste0("_", DM_PERIO_DEFINITION, "_", DM_COL)
} else if (CONDL_SOURCE == "PanelValidation") {
  paste0("_", DM_PERIO_DEFINITION, "_", DM_COL, "_", VALIDATION_CONTROLS)
} else ""

registry_key <- paste0(
  CONDL_SOURCE,
  caries_cutoff_key,
  caries_periodontal_key,
  dmft_key,
  dm_key,
  confounder_key,
  if (USE_ALL_BIOMARKERS) "_allBiomarkers" else if (USE_PUBLICATION_BIOMARKERS)
    paste0("_publication", length(publication_markers)) else "_panel",
  if (EXCLUDE_S100A8) "_noS100A8" else "",
  if (!USE_TIMP_RATIOS) "" else paste0("_timpRatios", length(timp_ratios_full)),
  if (USE_RATIO_FEATURES) "_ratios" else "_noRatios",
  if (USE_VIF) "_vif" else "_noVif"
)

registry_key <- paste0(registry_key, "_exhaustive12_beam3")

biomarker_registry[[registry_key]] <- list(
  project              = PROJECT,
  condl_source         = CONDL_SOURCE,
  caries_source_col    = if (CONDL_SOURCE == "Caries") CARIES_SOURCE_COL    else NULL,
  caries_binary_cutoff = if (CONDL_SOURCE == "Caries") CARIES_BINARY_CUTOFF else NULL,
  caries_nocaries_max  = if (CONDL_SOURCE == "Caries") CARIES_NOCARIES_MAX  else NULL,
  caries_high_min      = if (CONDL_SOURCE == "Caries") CARIES_HIGH_MIN      else NULL,
  caries_n_gap_dropped = CARIES_N_GAP_DROPPED,
  caries_periodontal_filter = if (CONDL_SOURCE == "Caries") CARIES_PERIODONTAL_FILTER else NULL,
  caries_short_labels  = if (CONDL_SOURCE == "Caries") {
    if (CARIES_SOURCE_COL == "Clinical_Caries") "LC / HC" else "CF / CP"
  } else NULL,
  caries_cutoff_method = if (CONDL_SOURCE == "Caries") {
    if (CARIES_SOURCE_COL == "Clinical_Caries") "Youden_J_vs_Caviated_gold_standard"
    else "no_caries_vs_high_caries_experience"
  } else NULL,
  dmft_components      = if (CONDL_SOURCE == "DMFT") dmft_cols else NULL,
  dmft_cutoff          = DMFT_CUTOFF_USED,
  dmft_cutoff_rule     = if (CONDL_SOURCE == "DMFT") { if (is.null(DMFT_CUTOFF)) "median" else "fixed" } else NULL,
  dmft_periodontal_filter = if (CONDL_SOURCE == "DMFT") DMFT_PERIODONTAL_FILTER else NULL,
  confounders          = confounders,
  dm_perio_definition  = if (DM_COMPARISON) DM_PERIO_DEFINITION else NULL,
  dm_col               = if (DM_COMPARISON) DM_COL else NULL,
  publication_set      = PUBLICATION_SET,
  exclude_s100a8       = EXCLUDE_S100A8,
  validation_controls  = if (CONDL_SOURCE == "PanelValidation") VALIDATION_CONTROLS else NULL,
  validation_group_counts = validation_group_counts,
  excluded_markers     = excluded_markers,
  ph_labels            = c(case = PH_CASE_LABEL, reference = PH_REF_LABEL),
  use_publication_biomarkers = USE_PUBLICATION_BIOMARKERS,
  publication_biomarkers     = publication_biomarkers,
  prior_panels               = prior_panels,
  prior_panel_check          = prior_panel_check,
  prior_panels_reported      = prior_panels_reported,
  use_timp_ratios      = USE_TIMP_RATIOS,
  timp_ratio_set       = TIMP_RATIO_SET,
  timp_ratios          = timp_ratios_full,
  use_all_biomarkers   = USE_ALL_BIOMARKERS,
  use_ratio_features   = USE_RATIO_FEATURES,
  use_vif              = USE_VIF,
  vif_threshold        = VIF_THRESHOLD,
  registry_key         = registry_key,
  timestamp            = Sys.time(),
  PH = list(absolute   = biomarkers_PH_vif, relative = NULL,
            stats      = stats_PH_vif,      ratio_cols = ratio_cols_PH),
  GH = list(absolute   = biomarkers_GH_vif, relative = NULL,
            stats      = stats_GH_vif,      ratio_cols = ratio_cols_GH),
  GP = list(absolute   = biomarkers_GP_vif, relative = NULL,
            stats      = stats_GP_vif,      ratio_cols = ratio_cols_GP)
)

cat(sprintf("\n===== biomarker_registry: '%s' =====\n", registry_key))
cat(sprintf("  All keys: %s\n", paste(names(biomarker_registry), collapse = ", ")))
cat("\nScript 1 complete.\n")

# Configuration for the age/BMI regression script;
# regression uses participants complete for Age_DOS and BMI; Sex does not determine eligibility.
if (PROJECT == "DMFT") {
  biomarker_registry[[registry_key]]$model_sets <- DMFT_MODEL_SETS
  cat("\nAge/BMI preprocessing ready. Run 02_LogReg_DMFT_Age_BMI.R next.\n")
}
