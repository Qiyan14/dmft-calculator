# Base-R checks of the calculator; synthetic values below are fixtures, never study results.
# Rscript tests/test_selected_model.R [path/to/app.R] [path/to/04b_Export_DMFT_Equation_Single_Biomarker.V9.R]
args <- commandArgs(trailingOnly=TRUE)
app_path <- if(length(args)) args[1] else 'app/app.R'
eq_path <- if(length(args)>1) args[2] else 'analysis/04b_Export_DMFT_Equation_Single_Biomarker.V9.R'
load_functions <- function(path) {
 e <- new.env(parent=globalenv())
 for(ex in parse(path)) if(is.call(ex) && identical(ex[[1]],as.name('<-'))) {
  nm <- as.character(ex[[2]])
  if ((is.call(ex[[3]]) && identical(ex[[3]][[1]],as.name('function'))) ||
      nm %in% c('PUBLISHED_MODEL','EXAMPLE','MINUS','HARD_EXHAUSTIVE_MAX','APP_VERSION')) eval(ex,e)
 }
 e
}
e <- load_functions(app_path)
stopifnot(is.character(e$APP_VERSION), nzchar(e$APP_VERSION))

# ---- bundled study model: one biomarker (Bacteroides heparinolyticus) + age -------------
m0 <- e$PUBLISHED_MODEL
stopifnot(e$model_ready(m0), nrow(m0$markers) == 1L, identical(m0$markers$key, 'Bacteroides_heparinolyticus'),
          identical(m0$covariates$var, 'Age_DOS'), all(m0$covariates$type == 'numeric'),
          m0$transform == 'log2', m0$delta == 5e-7, abs(m0$threshold - 0.3834346393) < 1e-12)
# Published equation, written out independently of the app (from 04b's dmft_age_model.csv)
published_eta <- function(x, age) {
 L <- log2(x + 5e-7)
 -21.4861170069685 - 0.842187455616679*L[1] + 0.10224870924323*age
}
# The biomarker slope equals log(OR per doubling) in the results workbook (Results_1_Biomarkers, top row)
or_doubling <- 0.4307672079
stopifnot(max(abs(m0$markers$coef - log(or_doubling))) < 1e-9)
set.seed(20261030); worst <- 0
for (i in 1:10000) {
 x <- ifelse(runif(1) < 0.3, 0, exp(runif(1, log(1e-6), log(1e3))))   # about 30% zeros (not detected)
 age <- runif(1, 18, 100)
 pr <- e$predict_model(m0, x, list(Age_DOS = age))
 worst <- max(worst, abs(pr$eta - published_eta(x, age)))
 stopifnot(abs(pr$p - plogis(published_eta(x, age))) < 1e-12)
}
stopifnot(worst < 1e-13)
cat(sprintf('PASS: bundled model matches the published equation for 10,000 random inputs (largest difference %.1e).\n', worst))
# Saved and reloaded model gives identical predictions
tab0 <- e$model_to_table(m0)
m1 <- e$model_from_table(tab0, 'selected.csv'); m2 <- e$selected_study_model(tab0)
x <- 3.5; cv <- list(Age_DOS = 52)
stopifnot(identical(e$predict_model(m1, x, cv)$p, e$predict_model(m0, x, cv)$p),
          identical(e$predict_model(m2, x, cv)$p, e$predict_model(m0, x, cv)$p))
# Desktop loader: the matching export is accepted, anything else is refused and the bundle kept
dd <- tempfile(); dir.create(dd)
write.csv(tab0, file.path(dd, 'dmft_age_model.csv'), row.names = FALSE, na = '')
r <- e$load_published_coefficients(m0, dd)
stopifnot(grepl('matches', r$msg), abs(e$predict_model(r$model, x, cv)$p - e$predict_model(m0, x, cv)$p) < 1e-12,
          identical(r$model$performance, m0$performance))   # CSV keeps 15 significant digits
# Example values shown by "Fill example values" and in the manuscript (not a participant)
ex <- e$predict_model(m0, e$EXAMPLE$markers, list(Age_DOS = e$EXAMPLE$covariates[1]))
stopifnot(length(e$EXAMPLE$covariates) == 1L, abs(ex$p - plogis(published_eta(e$EXAMPLE$markers, 40))) < 1e-12,
          round(100 * ex$p, 1) == 55.8,
          ex$p >= m0$threshold)
cat(sprintf('PASS: save/reload, desktop export check and example (%.1f%%, above the %.1f%% threshold).\n',
            100 * ex$p, 100 * m0$threshold))
# The age/BMI comparison chart was removed: BMI is continuous, not a category
stopifnot(!any(c('age_bmi_curve_data', 'draw_age_bmi_curves', 'age_bmi_variables') %in% ls(e)))

# ---- panel search: all pairs, then beam ------------------------------------------------
job <- list(p=472L,s=list(search_mode='beam',exhaustive_limit=10),parents=matrix(c(1L,2L),2))
pairs <- e$level_candidates(job,2L)
stopifnot(pairs$mode=='exhaustive',ncol(pairs$cand)==111156L,identical(pairs$cand[,ncol(pairs$cand)],c(471L,472L)))
triples <- e$level_candidates(job,3L)
stopifnot(triples$mode=='beam',ncol(triples$cand)==470L,all(triples$cand[1,]==1),all(triples$cand[2,]==2))
cat('PASS: all 111,156 pairs evaluated; beam extensions start at size 3.\n')

# ---- export script 04b -> app import (synthetic fixture) -------------------------------
# 04b reads a results workbook of script 02 and the session that script 01 leaves (data.dfPH,
# stats_PH, registry_key); here both are synthetic, and the workbook is written with writexl.
if (file.exists(eq_path) && requireNamespace('readxl', quietly=TRUE) && requireNamespace('writexl', quietly=TRUE)) {
eq_abs <- normalizePath(eq_path)
set.seed(643); n <- 120L
outcome <- factor(rep(c('LD','HD'),each=n/2),levels=c('LD','HD'))
raw <- exp(rnorm(n, ifelse(outcome=='HD', 1, 0), 2)); raw[sample(n, 15)] <- 0   # some not detected
age <- round(runif(n, 20, 70) + 8*(outcome=='HD'))
L <- log2(raw+5e-7); mu <- mean(L[outcome=='LD']); sdv <- sd(L[outcome=='LD'])
key <- 'test_DMFT_cutoff13.9_conf_AgeDOS_allBiomarkers_noRatios_noVif_exhaustive12_beam3'
ev <- new.env(parent=globalenv())
ev$data.dfPH <- data.frame(condL=outcome, Marker1=(L-mu)/sdv, Age_DOS=age)
ev$stats_PH <- data.frame(Column='Marker1', Mean=mu, SD=sdv)
ev$registry_key <- key
f <- glm(condL ~ Marker1 + Age_DOS, data=ev$data.dfPH, family=binomial())
pars <- c('Script 1 registry key'=key, 'Model configuration'='Biomarkers_Age', 'Confounders'='Age_DOS',
          'Reference group'='LD', 'Case group (positive)'='HD',
          'DMFT cut-off'='<= 13.9 -> LD (low DMFT) | > 13.9 -> HD (high DMFT) (fixed)',
          'Cross-validation'="5-fold CV x 5 repeats, seed 777 (caret's own folds)",
          'Total n'=n, 'n (LD)'=n/2, 'n (HD)'=n/2,
          'Cohort policy'='Participants complete for Age_DOS (every covariate of every set)')
res1 <- data.frame(Biomarker1='Marker1', Balanced_Accuracy=0.8, AUC=0.85, Threshold=0.47, Sensitivity=0.8,
                   Specificity=0.8, OR_perSD_Biomarker1=exp(coef(f)[['Marker1']]),
                   OR_perDoubling_Biomarker1=exp(coef(f)[['Marker1']]/sdv))
wd <- tempfile(); dir.create(wd); old_wd <- setwd(wd)
writexl::write_xlsx(list(Analysis_Parameters=data.frame(Parameter=names(pars), Value=unname(pars)),
                         Results_1_Biomarkers=res1), sprintf('Results_LogReg_HD_vs_LD_Biomarkers_Age_%s.xlsx', key))
invisible(capture.output(sys.source(eq_abs, envir=ev)))
setwd(old_wd)
out <- file.path(wd, 'DMFT_single_biomarker_equation_Biomarkers_Age')
stopifnot(ev$MAX_LINK_DIFF < 1e-8, all(ev$CHECKS$Relative_difference < 1e-9))
tab <- read.csv(file.path(out,'dmft_age_model.csv'), stringsAsFactors=FALSE)
m <- e$selected_study_model(tab)
for(i in 1:n) stopifnot(abs(e$predict_model(m, raw[i], list(Age_DOS=age[i]))$p - fitted(f)[[i]]) < 1e-10)
# A different model placed next to the desktop app is refused; the bundled model stays in use
r <- e$load_published_coefficients(m0, out)
stopifnot(grepl('was not used', r$msg), identical(r$model$markers$key, m0$markers$key), e$model_ready(r$model))
m2<-e$selected_study_model(e$model_to_table(m));stopifnot(e$model_ready(m2))
bad<-tab;bad$Label[bad$Term=='(Model configuration)']<-'Biomarkers_Age_Sex'
stopifnot(inherits(try(e$selected_study_model(bad),silent=TRUE),'try-error'))
bad<-tab[tab$Term!='(Threshold)',];stopifnot(inherits(try(e$selected_study_model(bad),silent=TRUE),'try-error'))
cat('PASS: 04b export, raw-input equivalence for 120 synthetic participants, import and stale-model rejection.\n')
} else cat('SKIPPED: 04b export check (needs', eq_path, 'and the readxl and writexl packages).\n')
cat('All checks passed.\n')
