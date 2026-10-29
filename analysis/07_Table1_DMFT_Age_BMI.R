# ============================================
# TABLE 1: Demographic and Clinical Characteristics -- DMFT project
# Low DMFT (LD) vs High DMFT (HD), Script 1 run with PROJECT = "DMFT"
# Continuous: mean (SD) and median (Q1, Q3); Wilcoxon (exact = FALSE) or
#             Welch t-test, chosen per variable (CONT_TEST below)
# Categorical: n (%); chi-square or Fisher's exact, chosen from expected counts
#
# Built from the PH / GH / GP Table 1 code (v2): same test helpers, smoker and
# diabetes rows, a footnote naming the test used for every row, fold
# differences, and .tex / .docx export. DMFT changes are marked "### DMFT".
#   - one comparison: data.dfPH, where Script 1 stores LD / HD
#   - rows follow the manuscript's Table 1: age, sex, race/ethnicity, BMI,
#     BMI category, number of teeth, periodontal diagnosis, (smoker,
#     diabetes), DMFT and its components
#   - DMFT is computed from Script 1's dmft_cols and checked against the LD / HD
#     split; it defines the groups, so it is shown but not tested
#   - console block with the numbers the manuscript asks for (N, cutoff,
#     DMFT distribution, overlap with cavitated lesions)
#
# Run after Script 1, in the same session:
#   source("script1_preprocessing.R")   # PROJECT = "DMFT"
#   source("07_Table1_DMFT_Age_BMI.R")
# ============================================

### DMFT: Hmisc is no longer needed (labels go through tbl_summary(label = ))
.need_pkgs <- c("dplyr", "gtsummary", "gt", "flextable", "officer")
.miss_pkgs <- .need_pkgs[!vapply(.need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkgs) > 0)
  stop("Missing package(s): ", paste(.miss_pkgs, collapse = ", "),
       "\n  install.packages(c(", paste0('"', .miss_pkgs, '"', collapse = ", "), "))", call. = FALSE)
suppressMessages({
  library(dplyr); library(gtsummary); library(gt); library(flextable); library(officer)
})


# ============================================
# SETTINGS
# ============================================
### DMFT: continuous test per variable.
###   "auto"   = Welch t-test when Shapiro-Wilk does not reject normality in
###              either group (p >= NORMALITY_ALPHA), otherwise Wilcoxon
###              rank-sum (exact = FALSE). The caries Table 1 footnote (BMI by
###              Welch t-test; age and teeth by Wilcoxon) follows this pattern.
###   "wilcox" = Wilcoxon for every continuous row, as in the PH / GH / GP tables.
CONT_TEST       <- "auto"
NORMALITY_ALPHA <- 0.05
# Force a test for named variables, e.g. c(BMI = "t.test", Teeth = "wilcox.test")
CONT_TEST_FORCE <- c()
# DMFT components are counts with many zeros: Wilcoxon for them whatever CONT_TEST says
COMPONENTS_WILCOXON <- TRUE

# One line per statistic for every continuous row
CONT_STATS <- c("{mean} ({sd})", "{median} ({p25}, {p75})")

INCLUDE_SMOKING_DIABETES <- FALSE    # cigarette smoker and diabetes rows (Script 1 section 7.2)
INCLUDE_PERIO_MEASURES   <- FALSE   # BOP, PD >= 4 mm, mean PD, CAL, PISA, PRA, plaque, bone loss
RACE_UNKNOWN_AS_MISSING  <- TRUE    # "Not reported or unknown" counted in the Missing row
BMI_CAT_CUTOFF <- if (exists("BMI_CUTOFF")) BMI_CUTOFF else 30

### DMFT: data codes -> table labels (codes as in the PH / GH / GP Table 1 code).
### CHECK against the codebook: the caries manuscript, from the same study,
### reports the opposite sex split (72.2% male) and Hispanic participants, which
### this coding does not give for the DMFT cohort.
SEX_CODES  <- c("1" = "Female", "2" = "Male")
RACE_CODES <- c("1" = "White",
                "2" = "Black or African American",
                "3" = "Hispanic or Latino",
                "4" = "Asian",
                "5" = "Not reported or unknown",
                "6" = "American Indian",
                "7" = "Native Hawaiian or Pacific Islander")

# DMFT component columns: Script 1's dmft_cols. Set by hand only if your
# Script 1 does not create dmft_cols, e.g. c("Carious", "Missing", "Restored").
DMFT_COMPONENT_COLS <- if (exists("dmft_cols")) dmft_cols else NULL


# ============================================
# Prerequisites from Script 1
# ============================================
.msgs <- character(0)
if (!exists("data.dfPH")) {
  .msgs <- "data.dfPH not found: run Script 1 with PROJECT = \"DMFT\" in this session first."
} else {
  if (exists("PROJECT") && !identical(PROJECT, "DMFT"))
    .msgs <- c(.msgs, sprintf("PROJECT = \"%s\"; this table is for PROJECT = \"DMFT\".", PROJECT))
  .grp <- sort(unique(as.character(data.dfPH$condL)))
  if (!identical(.grp, c("HD", "LD")))
    .msgs <- c(.msgs, sprintf("condL in data.dfPH should be LD and HD; found: %s.",
                              paste(.grp, collapse = ", ")))
  if (is.null(DMFT_COMPONENT_COLS)) {
    .msgs <- c(.msgs, "dmft_cols not found: set DMFT_COMPONENT_COLS to the DMFT component columns.")
  } else if (length(setdiff(DMFT_COMPONENT_COLS, colnames(data.dfPH))) > 0) {
    .msgs <- c(.msgs, paste("DMFT component column(s) not in data.dfPH:",
                            paste(setdiff(DMFT_COMPONENT_COLS, colnames(data.dfPH)), collapse = ", ")))
  }
  .need_cols <- setdiff(c("Age_DOS", "Sex", "Race", "BMI", "Teeth"), colnames(data.dfPH))
  if (length(.need_cols) > 0)
    .msgs <- c(.msgs, paste("Column(s) not in data.dfPH:", paste(.need_cols, collapse = ", ")))
}
if (length(.msgs) > 0)
  stop("\n\n===== PREREQUISITES NOT MET =====\n  ", paste(.msgs, collapse = "\n  "), "\n", call. = FALSE)


# ============================================
# Helper: select chi-square vs Fisher's exact (unchanged from v2)
# ============================================
determine_cat_test <- function(data, var, group) {
  cont_table      <- table(data[[var]], data[[group]])
  expected_counts <- suppressWarnings(chisq.test(cont_table))$expected
  if (sum(expected_counts < 5) > (0.2 * length(expected_counts)) ||
      any(expected_counts < 1)) {
    return("fisher.test")
  } else {
    return("chisq.test")
  }
}

# ============================================
### DMFT: helper: select Welch t-test vs Wilcoxon (see CONT_TEST)
# ============================================
determine_cont_test <- function(data, var, group, always_wilcox = character(0)) {
  if (var %in% names(CONT_TEST_FORCE)) return(unname(CONT_TEST_FORCE[[var]]))
  if (identical(CONT_TEST, "wilcox") || var %in% always_wilcox) return("wilcox.test")
  normal <- vapply(split(data[[var]], data[[group]]), function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 3) return(FALSE)
    p <- tryCatch(shapiro.test(x)$p.value, error = function(e) 0)
    p >= NORMALITY_ALPHA
  }, logical(1))
  if (all(normal)) "t.test" else "wilcox.test"
}

# ============================================
# Helper: fold difference (group1 median / group2 median) (unchanged from v2)
# ============================================
calculate_fold_diff <- function(data, variable, group_var, group1, group2) {
  m1 <- median(data[[variable]][data[[group_var]] == group1], na.rm = TRUE)
  m2 <- median(data[[variable]][data[[group_var]] == group2], na.rm = TRUE)
  if (!is.na(m1) && !is.na(m2) && m2 > 0) sprintf("%.2f", m1 / m2) else NA_character_
}

# Test names for the footnote. gtsummary's "chisq.test" keeps R's Yates
# continuity correction, which applies to 2 x 2 tables.
test_label <- function(test, var, data) {
  switch(test,
         t.test      = "Welch two-sample t-test",
         wilcox.test = "Wilcoxon rank-sum test",
         fisher.test = "Fisher's exact test",
         chisq.test  = if (length(unique(na.omit(as.character(data[[var]])))) == 2)
                         "Pearson chi-square test with Yates' continuity correction"
                       else "Pearson chi-square test",
         test)
}

# modify_footnote() became modify_footnote_header() in gtsummary 2.1
set_pvalue_footnote <- function(tbl, text) {
  if ("modify_footnote_header" %in% getNamespaceExports("gtsummary")) {
    gtsummary::modify_footnote_header(tbl, footnote = text, columns = "p.value", replace = TRUE)
  } else {
    gtsummary::modify_footnote(tbl, p.value = text)
  }
}


# ============================================
# Core: build Table 1 (LD vs HD)
# ============================================
build_table1_dmft <- function(dataset) {

  cat("\n", strrep("=", 50), "\n")
  cat("  Building Table 1: DMFT (Low DMFT vs High DMFT)\n")
  cat(strrep("=", 50), "\n")

  df <- as.data.frame(dataset)
  if (!all(c("Age_DOS", "BMI", "condL") %in% names(df))) stop("Run age/BMI preprocessing first.")
  keep <- complete.cases(df[c("Age_DOS", "BMI", "condL")])
  message("Table 1 age/BMI cohort: ", sum(keep), "; excluded ", sum(!keep))
  df <- droplevels(df[keep, , drop = FALSE])
  lev_LD <- "Low DMFT (LD)"
  lev_HD <- "High DMFT (HD)"

  # ------ recode condL ------
  ### DMFT: LD first, as in the manuscript (Overall | LD | HD | p-value)
  df$condL <- factor(as.character(df$condL), levels = c("LD", "HD"), labels = c(lev_LD, lev_HD))

  # ------ DMFT and its components ------
  ### DMFT: computed from Script 1's component columns
  comp_cols <- DMFT_COMPONENT_COLS
  comp_num  <- lapply(df[comp_cols], function(z) suppressWarnings(as.numeric(as.character(z))))
  bad_txt   <- mapply(function(raw, num) sum(!is.na(raw) & is.na(num)), df[comp_cols], comp_num)
  if (any(bad_txt > 0))
    warning("Non-numeric DMFT component values set to NA: ",
            paste(names(bad_txt)[bad_txt > 0], bad_txt[bad_txt > 0], sep = " = ", collapse = ", "))
  df$DMFT_total <- rowSums(as.data.frame(comp_num))          # NA if any component is NA

  if ("DMFT" %in% colnames(df)) {                            # a DMFT column stored by Script 1
    n_diff <- sum(abs(as.numeric(df$DMFT) - df$DMFT_total) > 1e-8, na.rm = TRUE)
    if (n_diff > 0)
      warning(sprintf("DMFT column differs from the sum of %s for %d subject(s); the sum is used.",
                      paste(comp_cols, collapse = " + "), n_diff))
  }

  cutoff <- if (exists("DMFT_CUTOFF_USED")) DMFT_CUTOFF_USED else NA_real_
  cut_rule <- if (!exists("DMFT_CUTOFF")) "rule unknown" else
    if (is.null(DMFT_CUTOFF)) "median of the analysed subjects" else "fixed"
  if (!is.na(cutoff)) {
    expect_HD <- df$DMFT_total > cutoff
    n_mis <- sum(expect_HD != (df$condL == lev_HD), na.rm = TRUE)
    if (n_mis > 0)
      warning(sprintf(paste0("%d subject(s) are in a different LD / HD group than DMFT <= %g implies. ",
                             "Check that DMFT_COMPONENT_COLS are the columns Script 1 used."),
                      n_mis, cutoff))
  }

  # component labels: D / M / F when each column's role is clear from its name
  role_of <- function(nm) {
    n   <- tolower(nm)
    tok <- strsplit(n, "[^a-z]+")[[1]]                  # "D_teeth" -> "d", "teeth"
    if (grepl("decay|cari|carr", n) || any(tok %in% c("d", "dt"))) return("D")   # "carr": Teeth_wCarries
    if (grepl("miss|extract", n)  || any(tok %in% c("m", "mt"))) return("M")
    if (grepl("fill|restor", n)   || any(tok %in% c("f", "ft"))) return("F")
    NA_character_
  }
  roles     <- vapply(comp_cols, role_of, character(1))
  use_roles <- !anyNA(roles) && !anyDuplicated(roles)
  role_lab  <- c(D = "Decayed teeth (D)", M = "Missing teeth (M)", F = "Filled teeth (F)")
  ord       <- if (use_roles) order(match(roles, c("D", "M", "F"))) else seq_along(comp_cols)
  comp_vars <- character(0)
  comp_labs <- list()
  for (i in ord) {
    v <- if (use_roles) paste0("DMFT_", roles[[i]]) else paste0("DMFT_part", i)
    df[[v]] <- comp_num[[i]]
    comp_vars <- c(comp_vars, v)
    comp_labs[[v]] <- if (use_roles) role_lab[[roles[[i]]]] else comp_cols[i]
  }
  cat(sprintf("\nDMFT = %s%s\n", paste(comp_cols, collapse = " + "),
              if (use_roles) "" else "  (component roles not recognised; column names used as labels)"))

  # ------ demographic factors ------
  ### DMFT: codes from SEX_CODES / RACE_CODES (settings); text labels pass through
  sx <- as.character(df$Sex)
  sx <- ifelse(sx %in% names(SEX_CODES), unname(SEX_CODES[sx]), sx)
  df$Sex <- factor(sx, levels = c("Male", "Female"))     ### DMFT: Male first, as in the manuscript

  rc <- as.character(df$Race)
  rc <- ifelse(rc %in% names(RACE_CODES), unname(RACE_CODES[rc]), rc)
  race_alias <- c("White/Caucasian" = "White",                      # v2 text labels
                  "Not Reported/Unknown" = "Not reported or unknown",
                  "Native Hawaiian/Pacific Islander" = "Native Hawaiian or Pacific Islander")
  rc <- ifelse(rc %in% names(race_alias), unname(race_alias[rc]), rc)
  if (RACE_UNKNOWN_AS_MISSING) rc[rc %in% RACE_CODES[["5"]]] <- NA
  df$Race <- droplevels(factor(rc, levels = unname(RACE_CODES)))

  for (v in c("Sex", "Race")) {                          # values that matched no label
    lost <- unique(as.character(dataset[[v]])[!is.na(dataset[[v]]) & is.na(df[[v]])])
    lost <- setdiff(lost, if (v == "Race" && RACE_UNKNOWN_AS_MISSING) c("5", RACE_CODES[["5"]]))
    if (length(lost) > 0)
      warning(sprintf("%s value(s) not in %s_CODES, shown as Missing: %s", v, toupper(v),
                      paste(lost, collapse = ", ")))
  }

  ### DMFT: BMI category, as in the caries Table 1
  df$BMI <- as.numeric(df$BMI)
  bmi_levels <- c(paste0("<", BMI_CAT_CUTOFF), paste0("\u2265", BMI_CAT_CUTOFF))
  df$BMI_cat <- factor(ifelse(is.na(df$BMI), NA,
                              ifelse(df$BMI >= BMI_CAT_CUTOFF, bmi_levels[2], bmi_levels[1])),
                       levels = bmi_levels)
  df$Age_DOS <- as.numeric(df$Age_DOS)
  df$Teeth   <- as.numeric(df$Teeth)

  ### DMFT: periodontal diagnosis (Script 1 keeps the original classification
  ### in condL_original when condL is recoded)
  perio_col <- intersect(c("condL_original", "condL_perio"), colnames(df))[1]
  has_perio <- !is.na(perio_col)
  if (has_perio) {
    pd <- as.character(df[[perio_col]])
    pd_levels <- c("Healthy periodontium", "Gingivitis", "Periodontitis")
    df$Perio_dx <- if (any(pd %in% pd_levels)) droplevels(factor(pd, levels = pd_levels)) else factor(pd)
  } else {
    cat("\nNo condL_original column: the periodontal diagnosis row is left out.\n")
  }

  # ------ smoking / diabetes (v2) ------
  yn_vars <- c("Cigsmoke", "Diabetes")
  if (INCLUDE_SMOKING_DIABETES) {
    if (!all(yn_vars %in% colnames(df)) ||
        !all(vapply(df[yn_vars], function(z) identical(levels(z), c("No", "Yes")), logical(1))))
      stop("Cigsmoke / Diabetes must be No/Yes factors. Re-run Script 1 (section 7.2) first, ",
           "or set INCLUDE_SMOKING_DIABETES <- FALSE.")
  }

  # ------ variables, in manuscript order ------
  perio_meas <- c(Per_BOP = "Bleeding on Probing (%)",
                  Per_4mm = "Sites with PD \u{2265}4mm (%)",
                  MeanPD = "Mean Pocket Depth (mm)",
                  MeanCAL = "Mean Clinical Attachment Loss (mm)",
                  PISA = "PISA",
                  AdjustPRA = "Adjusted PRA",
                  Per_PlaqueIndex = "Plaque Index (%)",
                  Per_AlvBoneLoss = "Alveolar Bone Loss (%)")
  perio_meas <- if (INCLUDE_PERIO_MEASURES) perio_meas[intersect(names(perio_meas), colnames(df))] else c()
  for (v in names(perio_meas)) df[[v]] <- as.numeric(df[[v]])

  lab <- c(list(Age_DOS  = "Age (years)",
                Sex      = "Sex",
                Race     = "Race/ethnicity",
                BMI      = "Body mass index (kg/m\u00b2)",
                BMI_cat  = "BMI category",
                Teeth    = "Number of teeth"),
           if (has_perio) list(Perio_dx = "Periodontal diagnosis"),
           if (INCLUDE_SMOKING_DIABETES) list(Cigsmoke = "Cigarette smoker",
                                              Diabetes = "Diabetes (type 1 or 2)"),
           as.list(perio_meas),
           list(DMFT_total = "DMFT"),
           comp_labs)
  table_vars <- names(lab)
  cont_vars  <- intersect(c("Age_DOS", "BMI", "Teeth", names(perio_meas), "DMFT_total", comp_vars),
                          table_vars)
  dich_vars  <- intersect(yn_vars, table_vars)
  cat_vars   <- setdiff(table_vars, c(cont_vars, dich_vars))

  # ------ checks printed to the console ------
  cat("\nRace distribution:\n")
  print(table(df$Race, df$condL, useNA = "ifany"))
  if (has_perio) {
    cat("\nPeriodontal diagnosis (", perio_col, "):\n", sep = "")
    print(table(df$Perio_dx, df$condL, useNA = "ifany"))
  }
  if (INCLUDE_SMOKING_DIABETES) {
    yn_check <- intersect(c("Cigsmoke", "Diabetes", "T1Diabetes", "T2Diabetes"), colnames(df))
    cat("\nSmoking / diabetes, n 'Yes' by group:\n")
    print(sapply(split(df[yn_check], df$condL),
                 function(d) vapply(d, function(z) sum(z == "Yes", na.rm = TRUE), numeric(1))))
  }
  n_na <- colSums(is.na(df[table_vars]))
  if (any(n_na > 0))
    cat("\nMissing values (shown as 'Missing' rows; excluded from % and tests):",
        paste(names(n_na)[n_na > 0], n_na[n_na > 0], sep = " = ", collapse = ", "), "\n")
  if (exists("confounders") && length(confounders) > 0) {
    conf_in <- intersect(confounders, colnames(df))
    n_cc <- sum(!complete.cases(df[conf_in]))
    if (n_cc > 0)
      cat(sprintf(paste0("Note: %d subject(s) have a missing confounder (%s); adjusted panel searches ",
                         "leave them out, but they are in this table.\n"), n_cc, paste(conf_in, collapse = ", ")))
  }

  # ------ select tests ------
  ### DMFT: DMFT defines the groups and is not tested; neither is a row with
  ### fewer than two observed categories
  n_obs_levels <- function(v) length(unique(na.omit(as.character(df[[v]]))))
  tested_cont <- setdiff(cont_vars, "DMFT_total")
  tested_cat  <- Filter(function(v) n_obs_levels(v) >= 2, c(cat_vars, dich_vars))
  untested    <- setdiff(table_vars, c(tested_cont, tested_cat))

  cont_tests <- vapply(tested_cont, function(v)
    determine_cont_test(df, v, "condL", if (COMPONENTS_WILCOXON) comp_vars else character(0)),
    character(1))
  cat_tests  <- vapply(tested_cat,  function(v) determine_cat_test(df, v, "condL"),  character(1))
  all_tests  <- c(cont_tests, cat_tests)

  cat("\nTests selected:\n")
  for (v in names(all_tests))
    cat(sprintf("  %-28s %s\n", lab[[v]], test_label(all_tests[[v]], v, df)))
  for (v in untested)
    cat(sprintf("  %-28s not tested\n", lab[[v]]))

  wilcox_vars <- names(all_tests)[all_tests == "wilcox.test"]
  t_vars      <- names(all_tests)[all_tests == "t.test"]
  test_args   <- c(setNames(rep(list(list(exact = FALSE)), length(wilcox_vars)), wilcox_vars),
                   setNames(rep(list(list(var.equal = FALSE)), length(t_vars)), t_vars))

  # footnote: each test with the rows it was used for
  used_labels <- vapply(names(all_tests), function(v) test_label(all_tests[[v]], v, df), character(1))
  by_test     <- split(names(used_labels), factor(used_labels, levels = unique(used_labels)))
  footnote_txt <- paste0(
    "Statistical tests: ",
    paste(sprintf("%s: %s", names(by_test),
                  vapply(by_test, function(vs) paste(unlist(lab[vs]), collapse = ", "), character(1))),
          collapse = "; "),
    ". DMFT defines the groups and is not tested.",
    if (length(setdiff(untested, "DMFT_total")) > 0)
      paste0(" Not tested (one category observed): ",
             paste(unlist(lab[setdiff(untested, "DMFT_total")]), collapse = ", "), "."))

  # ------ build gtsummary table ------
  tbl <- df %>%
    select(condL, all_of(table_vars)) %>%
    tbl_summary(
      by           = condL,
      label        = lab,
      type         = c(setNames(rep(list("continuous2"), length(cont_vars)), cont_vars),
                       setNames(rep(list("dichotomous"), length(dich_vars)), dich_vars),
                       setNames(rep(list("categorical"), length(cat_vars)),  cat_vars)),
      value        = setNames(rep(list("Yes"), length(dich_vars)), dich_vars),
      statistic    = list(all_continuous2() ~ CONT_STATS,
                          all_categorical() ~ "{n} ({p}%)"),
      digits       = list(all_continuous2() ~ 2,
                          all_categorical() ~ c(0, 1)),
      missing      = "ifany",
      missing_text = "Missing"
    ) %>%
    add_p(
      test       = as.list(all_tests),
      test.args  = test_args,
      include    = -any_of(untested),
      pvalue_fun = function(x) style_pvalue(x, digits = 3)
    ) %>%
    add_overall(col_label = "**Overall**  \n(N = {style_number(N)})") %>%
    modify_header(label = "**Characteristic**",
                  all_stat_cols(stat_0 = FALSE) ~ "**{level}**  \n(N = {style_number(n)})",
                  p.value = "**p-value**") %>%
    modify_spanning_header(all_stat_cols(stat_0 = FALSE) ~ "**DMFT group**") %>%
    bold_labels() %>%
    modify_table_styling(columns = "p.value",
                         rows = variable == "DMFT_total" & row_type == "label",
                         missing_symbol = "\u2014") %>%
    set_pvalue_footnote(footnote_txt)

  # ------ fold differences (HD / LD medians) ------
  fold_df <- data.frame(
    Variable = cont_vars,
    Label    = unlist(lab[cont_vars]),
    FoldDiff = vapply(cont_vars, function(v)
      calculate_fold_diff(df, v, "condL", lev_HD, lev_LD), character(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  colnames(fold_df)[3] <- "Fold_HD_over_LD"
  cat("\nFold differences (HD median / LD median; NA when the LD median is 0):\n")
  print(fold_df, row.names = FALSE)

  # ------ numbers for the manuscript ------
  ### DMFT: values for the highlighted placeholders in Methods 2.1 and Results 3.1
  n_grp <- table(df$condL)
  d     <- df$DMFT_total
  q     <- quantile(d, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  filt  <- if (exists("DMFT_PERIODONTAL_FILTER") && !is.null(DMFT_PERIODONTAL_FILTER))
             paste(DMFT_PERIODONTAL_FILTER, collapse = ", ") else "none"
  cat("\n===== Numbers for the manuscript =====\n")
  cat(sprintf("  N = %d | LD = %d | HD = %d\n", nrow(df), n_grp[[lev_LD]], n_grp[[lev_HD]]))
  if (is.na(cutoff)) {
    cat("  Cutoff: DMFT_CUTOFF_USED not found in this session\n")
  } else {
    cat(sprintf("  Cutoff: DMFT <= %s -> LD, > %s -> HD (%s)\n", format(cutoff), format(cutoff), cut_rule))
  }
  cat(sprintf("  Periodontal filter: %s\n", filt))
  cat(sprintf("  DMFT: median %.1f (Q1 %.1f, Q3 %.1f); mean %.2f (SD %.2f); range %g-%g; missing %d\n",
              q[2], q[1], q[3], mean(d, na.rm = TRUE), sd(d, na.rm = TRUE),
              min(d, na.rm = TRUE), max(d, na.rm = TRUE), sum(is.na(d))))
  if (!is.na(cutoff))
    cat(sprintf("  Subjects with DMFT equal to the cutoff (all in LD): %d\n", sum(d == cutoff, na.rm = TRUE)))
  for (v in comp_vars)
    cat(sprintf("  %-18s HD subjects with 0: %d of %d\n", comp_labs[[v]],
                sum(df[[v]][df$condL == lev_HD] == 0, na.rm = TRUE), n_grp[[lev_HD]]))
  if ("Cavitated_Caries" %in% colnames(df)) {
    cav <- cut(suppressWarnings(as.numeric(df$Cavitated_Caries)), breaks = c(-Inf, 0, 2, Inf),
               labels = c("0", "1-2", ">=3"))
    cat("\n  Cavitated lesions by DMFT group (>=3 = the caries manuscript's HC definition):\n")
    print(table(DMFT_group = df$condL, Cavitated_lesions = cav, useNA = "ifany"))
  }

  # ------ export ------
  filt_tag  <- if (filt == "none") "allPerio" else
                 paste(gsub("[^A-Za-z0-9]+", "", DMFT_PERIODONTAL_FILTER), collapse = "Plus")
  stem      <- sprintf("table1_DMFT_HD_vs_LD_cut%s_%s",
                       if (is.na(cutoff)) "Unknown" else format(cutoff), filt_tag)
  tex_file  <- paste0(stem, ".tex")
  docx_file <- paste0(stem, ".docx")

  # gt drops the header line breaks from the LaTeX file and says so; that note is muffled
  withCallingHandlers(
    gt::gtsave(as_gt(tbl), filename = tex_file, path = getwd()),
    warning = function(w) if (grepl("HTML tags found", conditionMessage(w))) invokeRestart("muffleWarning"))
  ### DMFT: 10 pt on a landscape Letter page; at 10 pt the table is ~8 in wide,
  ### more than a portrait page holds
  ft <- as_flex_table(tbl) %>%
    flextable::fontsize(size = 10, part = "all") %>%
    flextable::padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
    flextable::autofit()
  save_as_docx(ft, path = docx_file,
               pr_section = officer::prop_section(
                 page_size = officer::page_size(width = 8.5, height = 11, orient = "landscape")))

  cat("\nExported:\n  -", tex_file, "\n  -", docx_file, "\n")

  invisible(list(table = tbl, fold_diff = fold_df, data = df))
}


# ============================================
# Run
# ============================================
results_DMFT <- build_table1_dmft(data.dfPH)

cat("\n", strrep("=", 50), "\n  PREVIEW: DMFT (LD vs HD)\n", strrep("=", 50), "\n")
print(results_DMFT$table)

cat("\nAll done. Two files written to:", getwd(), "\n")
