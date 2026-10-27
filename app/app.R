# ==========================================================================
# DMFT SALIVARY BIOMARKER CALCULATOR - WEB EDITION (patient prediction)
# ==========================================================================
# Estimates an individual's probability of high DMFT (DMFT > 13) from six
# absolute raw salivary biomarker concentrations, age, number of teeth and
# BMI, with the published six-marker logistic regression model or with a
# model saved from the desktop app ("Download model"). The page shows the
# equation with its coefficients, the calculation for the values entered and
# a speedometer gauge of the estimate.
#
# This is the patient-prediction module of the desktop app (app.R) without
# the data-analysis module, so that it runs entirely in a web browser through
# Shinylive (webR). It needs only shiny and base R. The model, helper and
# display functions are copied unchanged from the desktop app.
#
# Run locally:  shiny::runApp("app")
# Build site:   source("build.R")      (runs shinylive::export("app", "docs"))
# ==========================================================================

suppressPackageStartupMessages(library(shiny))

MINUS <- "−"
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


# ==========================================================================
# 1) THE PUBLISHED MODEL
# ==========================================================================
# Raw-input form: logit p(HD) = intercept + sum_j coef_j * log2(x_j + delta)
#                               + coef_age * Age + coef_teeth * Teeth + coef_bmi * BMI
# All coefficients are the Coefficient_raw_input column of
# dmft_model_coefficients.csv (the six-marker model refitted to all 111
# participants). The log2 transform is applied to the raw value; the LD-group
# centring and scaling used in fitting are absorbed into the intercept and
# slopes.

PUBLISHED_MODEL <- list(
  source      = "published",
  name        = "Published DMFT model",
  description = "Six salivary biomarkers with age, number of teeth and BMI (logistic regression).",
  case_label  = "high DMFT (DMFT > 13)",
  ref_label   = "low DMFT (DMFT ≤ 13)",
  case_short  = "HD",
  ref_short   = "LD",
  transform   = "log2",
  delta       = 5e-7,
  intercept   = -6.64840612175518,
  markers = data.frame(
    key   = c("Neisseriaceae_unclassified_Neisseriaceae_unclassified", "Actinomyces_sp._HMT_169",
              "Kingella_oralis", "Actinomyces_gerencseriae", "Moraxella_Moraxella_unclassified",
              "Bacteroides_heparinolyticus"),
    abbr  = c("NEI", "A169", "KO", "AG", "MOR", "BH"),
    label = c("unclassified Neisseriaceae", "*Actinomyces* sp. HMT 169", "*Kingella oralis*",
              "*Actinomyces gerencseriae*", "unclassified *Moraxella*", "*Bacteroides heparinolyticus*"),
    coef  = c(-0.048272380750755, 0.0843432793202229, 0.0609430607758506, 0.0964328703096935,
              0.106107184433687, -0.803513146783222),
    stringsAsFactors = FALSE),
  covariates = data.frame(
    var   = c("Age_DOS", "Teeth", "BMI"),
    term  = c("Age_DOS", "Teeth", "BMI"),
    abbr  = c("Age", "Teeth", "BMI"),
    label = c("Age", "Number of teeth present", "Body mass index"),
    unit  = c("years", "", "kg/m²"),
    type  = "numeric",
    level = NA_character_,
    coef  = c(0.115676884794072, -0.4641425581805, 0.0371782506733584),
    stringsAsFactors = FALSE),
  cov_levels  = list(),
  threshold   = 0.4878687616,
  cohort      = paste("Fitted in 111 adults (56 with DMFT ≤ 13 and 55 with DMFT > 13) recruited at the",
                      "University of Kentucky College of Dentistry between November 2020 and December 2025."),
  performance = c("Balanced accuracy" = 0.8993, "Sensitivity" = 0.9402, "Specificity" = 0.8584, "AUC" = 0.8953),
  performance_note = "Five-fold cross-validation repeated 100 times in the development cohort; not externally validated.",
  caution     = paste("In the development data, elevated *Bacteroides heparinolyticus* values occurred in only five",
                      "participants, all with low DMFT, so a detectable value lowers the estimate sharply.")
)

model_ready <- function(m) {
  !is.null(m) && is.finite(m$intercept) && all(is.finite(m$markers$coef)) &&
    (nrow(m$covariates) == 0 || all(is.finite(m$covariates$coef)))
}

# Illustrative values for the "Fill example values" button (not a real participant):
# NEI, A169, KO, AG, MOR, BH, then age, number of teeth and BMI
EXAMPLE <- list(markers = c(4.6, 18.2, 13.5, 1.4, 0, 0), covariates = c(40, 28, 27))


# ==========================================================================
# 2) HELPERS, MODEL FILES AND PREDICTION (as in the desktop app)
# ==========================================================================

# Number for display: 4 significant digits, Unicode minus, x 10^n for tiny or huge values
fmt_num <- function(x, sig = 4, html = TRUE, keep0 = FALSE) {
  if (length(x) != 1 || is.na(x)) return("NA")
  if (!is.finite(x)) return(if (x > 0) "∞" else paste0(MINUS, "∞"))
  if (x == 0) return("0")
  ax <- abs(x)
  s <- if (ax >= 1e-3 && ax < 1e7) {   # keep0: keep trailing zeros (13.00 rather than 13)
    formatC(signif(x, sig), format = "fg", digits = sig, big.mark = ",", flag = if (keep0) "#" else "")
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
# *text* -> <i>text</i>, everything else escaped
md_html <- function(x) gsub("\\*([^*]+)\\*", "<i>\\1</i>", htmltools::htmlEscape(x))
# Readable label from a column name
pretty_name <- function(x) gsub("_", " ", x, fixed = TRUE)
# "a", "a and b", "a, b and c"
and_list <- function(x) if (length(x) <= 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])

# Model loaded from a coefficient table saved by the desktop app ("Download model")
model_from_table <- function(tab, name) {
  need <- c("Term", "Coefficient_raw_input", "Input")
  if (!all(need %in% names(tab))) stop("A model file needs the columns Term, Coefficient_raw_input and Input.")
  tab$Term <- trimws(as.character(tab$Term)); tab$Input <- trimws(as.character(tab$Input %||% ""))
  val <- function(t) { v <- suppressWarnings(as.numeric(tab$Coefficient_raw_input[tab$Term == t])); if (length(v)) v[1] else NA_real_ }
  lab <- function(t, d) { v <- tab$Label[tab$Term == t]; if (length(v) && !is.na(v[1]) && nzchar(v[1])) v[1] else d }
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
       covariates = covs, cov_levels = cov_levels, threshold = if (is.finite(thr)) thr else 0.5,
       cohort = "Loaded from a saved model file.", performance = NULL,
       performance_note = if (is.finite(thr)) "" else "No threshold in the file; 0.5 is used.", caution = NULL)
}

# Linear predictor and probability, with each term of the sum
predict_model <- function(m, marker_vals, cov_vals) {
  tx <- if (m$transform == "log2") log2(marker_vals + m$delta) else marker_vals
  terms <- data.frame(Term = m$markers$label, Abbr = m$markers$abbr, Input = marker_vals, Transformed = tx,
                      Coefficient = m$markers$coef, Contribution = m$markers$coef * tx, kind = "marker",
                      stringsAsFactors = FALSE)
  if (nrow(m$covariates)) {
    cv <- m$covariates
    xin <- vapply(seq_len(nrow(cv)), function(i) {
      v <- cov_vals[[cv$var[i]]]
      if (cv$type[i] == "factor") as.numeric(identical(as.character(v), cv$level[i])) else as.numeric(v)
    }, numeric(1))
    terms <- rbind(terms, data.frame(Term = cv$label, Abbr = cv$abbr, Input = xin, Transformed = xin,
                                     Coefficient = cv$coef, Contribution = cv$coef * xin, kind = cv$type,
                                     stringsAsFactors = FALSE))
  }
  eta <- m$intercept + sum(terms$Contribution)
  list(terms = terms, eta = eta, p = 1 / (1 + exp(-eta)))
}


# ==========================================================================
# 3) DISPLAYS: EQUATION AND GAUGE (as in the desktop app)
# ==========================================================================

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
        parts <- c(parts, sprintf("<b>%s</b> is 1 if %s is %s and 0 otherwise", esc(cv$abbr[i]), esc(pretty_name(cv$var[i])),
                                  esc(cv$level[i])))
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


# ==========================================================================
# 4) USER INTERFACE
# ==========================================================================

CSS <- "
:root { --accent:#0b7a75; --accent-dark:#075e5a; --ink:#0b0b0b; --ink2:#52514e; --card:#f6f7f8; --line:#dde1e5; }
body { color: var(--ink); background:#ffffff;
       font-family: system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }
.topbar { display:flex; align-items:baseline; flex-wrap:wrap; gap:2px 14px; background:#ffffff;
          border-bottom:1px solid var(--line); margin:0 -15px 18px; padding:14px 15px 12px; }
.topbar .brand { color: var(--accent-dark); font-weight:600; font-size:20px; }
.topbar .brand-sub { color: var(--ink2); font-size:14px; }
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
.calc-wrap { overflow-x:auto; }
.calc-table { width:100%; min-width:560px; border-collapse:collapse; background:#fff; font-size:14px; }
.calc-table th, .calc-table td { border-bottom:1px solid var(--line); padding:6px 8px; text-align:right; }
.calc-table th:first-child, .calc-table td:first-child { text-align:left; }
.calc-table tr.total td { font-weight:600; border-top:2px solid #b5b3ac; }
.notice { background:#fdecec; border:1px solid #f3b4b4; border-radius:8px; padding:10px 12px; margin-bottom:12px; }
.card .shiny-input-container:not(.shiny-input-container-inline) { width:100% !important; max-width:100%; }
.footer-note { color: var(--ink2); font-size:12px; margin:8px 0; }
.footer-note:first-child { margin-top:24px; }
"

ui <- fluidPage(
  title = "Salivary biomarker DMFT calculator",
  tags$head(tags$style(HTML(CSS))),
  div(class = "topbar",
      span(class = "brand", "Salivary biomarker DMFT calculator"),
      span(class = "brand-sub", "Patient prediction")),
  fluidRow(column(12,
    p(class = "lead-muted", "Estimate an individual's probability of high DMFT (DMFT > 13) from six absolute raw",
      "salivary biomarker concentrations, age, number of teeth and body mass index."),
    radioButtons("model_choice", NULL, inline = TRUE,
                 choices = c("Published DMFT model" = "published", "Load a saved model" = "loaded")),
    conditionalPanel("input.model_choice == 'loaded'",
      fileInput("model_file", "Saved model file (.csv saved with \"Download model\" in the desktop app)", accept = ".csv")),
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
          actionButton("calc", "Calculate probability", class = "btn-accent"),
          actionButton("clear", "Clear", class = "btn-soft"),
          uiOutput("example_btn", inline = TRUE)))),
    column(7,
      div(class = "card",
        h4(textOutput("gauge_title", inline = TRUE)),
        plotOutput("gauge", height = "250px"),
        uiOutput("hero_ui")),
      div(class = "card", h4("What the probability represents"), uiOutput("represents_ui")))
  ),
  fluidRow(column(12, div(class = "card",
    h4("How the probability is calculated"),
    uiOutput("equation_ui"),
    uiOutput("calc_ui")))),
  fluidRow(column(12,
    p(class = "footer-note",
      "Research use only. Estimates come from a cross-sectional development cohort and have not been externally validated;",
      "they are not a diagnosis and do not replace a dental examination."),
    p(class = "footer-note",
      "This calculator runs entirely in your web browser: the values you enter are not sent to or stored on any server.")))
)


# ==========================================================================
# 5) SERVER
# ==========================================================================

server <- function(input, output, session) {

  rv <- reactiveValues(published = PUBLISHED_MODEL, loaded = NULL, last_pred = NULL)

  observeEvent(input$model_file, {
    tab <- tryCatch(utils::read.csv(input$model_file$datapath, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) e)
    m <- if (inherits(tab, "error")) tab else tryCatch(model_from_table(tab, input$model_file$name), error = function(e) e)
    if (inherits(m, "error")) { showNotification(conditionMessage(m), type = "error", duration = 8); return() }
    rv$loaded <- m
  })
  current_model <- reactive(switch(input$model_choice %||% "published", published = rv$published, loaded = rv$loaded))
  observeEvent(current_model(), rv$last_pred <- NULL, ignoreNULL = FALSE)

  output$model_header <- renderUI({
    m <- current_model()
    if (is.null(m)) return(div(class = "notice", "Choose a saved model file."))
    p(class = "lead-muted", HTML(paste0("<b>", htmltools::htmlEscape(m$name), "</b>. ", md_html(m$description %||% ""))))
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
          selectInput(paste0("cv_", i), pretty_name(v), choices = lv)
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

  # Example values (published model only): illustrative inputs, not a real participant
  output$example_btn <- renderUI({
    if (identical(input$model_choice %||% "published", "published"))
      actionButton("example", "Fill example values", class = "btn-soft",
                   title = "Illustrative values, not from a real participant")
  })
  observeEvent(input$example, {
    m <- current_model(); req(m, identical(m$source, "published"))
    for (j in seq_len(nrow(m$markers))) updateNumericInput(session, paste0("mk_", j), value = EXAMPLE$markers[j])
    for (i in seq_len(nrow(m$covariates))) updateNumericInput(session, paste0("cv_", i), value = EXAMPLE$covariates[i])
    rv$last_pred <- NULL
    showNotification("Illustrative values filled in (not a real participant). Click \"Calculate probability\".",
                     type = "message", duration = 6)
  })

  observeEvent(input$clear, {
    m <- current_model(); req(m)
    for (j in seq_len(nrow(m$markers))) updateNumericInput(session, paste0("mk_", j), value = NA)
    vars <- unique(m$covariates$var)
    for (i in seq_along(vars)) if (m$covariates$type[m$covariates$var == vars[i]][1] == "numeric") updateNumericInput(session, paste0("cv_", i), value = NA)
    rv$last_pred <- NULL
  })

  observeEvent(input$calc, {
    m <- current_model(); req(m)
    if (!model_ready(m)) { showNotification("This model's coefficients are incomplete, so no probability can be calculated.", type = "error"); return() }
    x <- vapply(seq_len(nrow(m$markers)), function(j) as.numeric(input[[paste0("mk_", j)]] %||% NA), numeric(1))
    if (anyNA(x)) { showNotification("Enter a value for every biomarker (0 if not detected).", type = "error"); return() }
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
    if (identical(m$source, "published")) {
      j <- match("Bacteroides_heparinolyticus", m$markers$key)
      if (!is.na(j) && x[j] > 0)
        warn <- c(warn, paste("*Bacteroides heparinolyticus* was detected. In the development data, elevated values occurred",
                              "only in low-DMFT participants, so this value pulls the estimate down sharply; interpret with caution."))
    }
    out <- predict_model(m, x, cv); out$warn <- warn
    rv$last_pred <- out
  })

  output$gauge_title <- renderText({
    m <- current_model(); if (is.null(m)) "Estimated probability" else paste("Estimated probability of", m$case_label)
  })
  output$gauge <- renderPlot({
    m <- current_model(); pr <- rv$last_pred
    # narrow screens (phones): smaller labels so that they fit beside the dial
    w <- session$clientData$output_gauge_width %||% 600
    op <- graphics::par(cex = min(1, max(0.6, w / 520))); on.exit(graphics::par(op))
    draw_gauge(if (is.null(pr)) NA else pr$p, if (is.null(m)) NA else m$threshold,
               lo_lab = "Lower", hi_lab = "Higher")
  }, res = 110)
  output$hero_ui <- renderUI({
    m <- current_model(); pr <- rv$last_pred
    if (is.null(m)) return(NULL)
    if (is.null(pr)) return(p(class = "hero-sub", "Enter the values and click Calculate probability."))
    above <- pr$p >= m$threshold
    tagList(
      div(class = "hero", fmt_prob(pr$p)),
      div(class = "hero-sub", sprintf("estimated probability of %s versus %s", m$case_label, m$ref_label)),
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
      if (!is.null(perf)) p(sprintf("Cross-validated performance at the %s threshold: %s.", fmt_pct(m$threshold), perf),
                            span(class = "help-note", m$performance_note)),
      if (!is.null(m$caution)) p(HTML(md_html(m$caution))),
      p(class = "help-note", "Research estimate only: calibration and performance in new participants have not been established,",
        "and the probability is not a diagnosis."))
  })
  output$equation_ui <- renderUI({
    m <- current_model(); req(m)
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
        sprintf("log<sub>2</sub>(%s + &delta;) = %s", fmt_num(t$Input[i], 7), fmt_num(t$Transformed[i])) else fmt_num(t$Input[i], 7)
      nm <- if (grepl(t$Abbr[i], t$Term[i], fixed = TRUE)) md_html(t$Term[i])
            else sprintf("%s (%s)", md_html(t$Term[i]), htmltools::htmlEscape(t$Abbr[i]))
      tags$tr(tags$td(HTML(nm)), tags$td(HTML(tx)),
              tags$td(HTML(fmt_num(t$Coefficient[i]))), tags$td(HTML(fmt_num(t$Contribution[i], keep0 = TRUE))))
    })
    tagList(
      h4("Calculation for the values entered", style = "margin-top:18px;"),
      div(class = "calc-wrap", tags$table(class = "calc-table",
        tags$tr(tags$th("Term"), tags$th("Value used"), tags$th("Coefficient"), tags$th("Coefficient × value")),
        tags$tr(tags$td("Intercept"), tags$td(""), tags$td(""), tags$td(HTML(fmt_num(m$intercept)))),
        rows,
        tags$tr(class = "total", tags$td(HTML("Linear predictor &eta; (sum)")), tags$td(""), tags$td(""), tags$td(HTML(fmt_num(pr$eta, keep0 = TRUE)))))),
      div(class = "equation", style = "margin-top:12px;",
          HTML(sprintf("<i>p</i> = 1 / [1 + exp(%s%s)] = %s, i.e. <b>%s</b>", MINUS,
                       if (pr$eta < 0) sprintf("(%s)", fmt_num(pr$eta)) else fmt_num(pr$eta),
                       if (pr$p > 0.999) sprintf("1 %s %s", MINUS, fmt_num(1 - pr$p, 3)) else fmt_num(pr$p, 4),
                       fmt_prob(pr$p)))))
  })
}

if (!isTRUE(getOption("dmft.app.test"))) shinyApp(ui, server)
