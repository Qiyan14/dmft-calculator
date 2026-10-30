# Base-R checks of the calculator; synthetic values below are fixtures, never study results.
# Rscript tests/test_selected_model.R [path/to/app.R] [path/to/04_Export_DMFT_Equation_Age_BMI.R]
args <- commandArgs(trailingOnly=TRUE)
app_path <- if(length(args)) args[1] else 'app/app.R'
eq_path <- if(length(args)>1) args[2] else 'analysis/04_Export_DMFT_Equation_Age_BMI.R'
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

# ---- bundled selected model: six biomarkers + age + BMI -------------------------------
m0 <- e$PUBLISHED_MODEL
stopifnot(e$model_ready(m0), nrow(m0$markers) == 6L, identical(m0$covariates$var, c('Age_DOS','BMI')),
          all(m0$covariates$type == 'numeric'), m0$transform == 'log2', m0$delta == 5e-7,
          abs(m0$threshold - 0.4303368621) < 1e-12)
# Published equation, written out independently of the app (from 04's DMFT_equation.txt)
published_eta <- function(x, age, bmi) {
 L <- log2(x + 5e-7)
 -25.593207966142202 + 0.054716075620806985*L[1] + 0.93738465006011396*L[2] - 0.89967475391325979*L[3] -
   0.91717962460446523*L[4] + 0.86395256497572881*L[5] - 0.93066377384919741*L[6] +
   0.12236268144283166*age + 0.045571695632458567*bmi
}
# The biomarker slopes equal log(OR per doubling) in the results workbook (Best_By_Size, 6 biomarkers)
or_doubling <- c(1.0562406796, 2.5532949187, 0.4067019164, 0.3996446009, 2.3725197238, 0.3942919028)
stopifnot(max(abs(m0$markers$coef - log(or_doubling))) < 1e-9)
set.seed(20261030); worst <- 0
for (i in 1:10000) {
 x <- ifelse(runif(6) < 0.3, 0, exp(runif(6, log(1e-6), log(1e3))))   # about 30% zeros (not detected)
 age <- runif(1, 18, 100); bmi <- runif(1, 12, 70)
 pr <- e$predict_model(m0, x, list(Age_DOS = age, BMI = bmi))
 worst <- max(worst, abs(pr$eta - published_eta(x, age, bmi)))
 stopifnot(abs(pr$p - plogis(published_eta(x, age, bmi))) < 1e-12)
}
stopifnot(worst < 1e-13)
cat(sprintf('PASS: bundled model matches the published equation for 10,000 random inputs (largest difference %.1e).\n', worst))
# Saved and reloaded model gives identical predictions
tab0 <- e$model_to_table(m0)
m1 <- e$model_from_table(tab0, 'selected.csv'); m2 <- e$selected_age_bmi_model(tab0)
x <- c(0.2, 0, 3.5, 0, 1e-4, 12); cv <- list(Age_DOS = 52, BMI = 31.4)
stopifnot(identical(e$predict_model(m1, x, cv)$p, e$predict_model(m0, x, cv)$p),
          identical(e$predict_model(m2, x, cv)$p, e$predict_model(m0, x, cv)$p))
# Desktop loader: the matching export is accepted, anything else is refused and the bundle kept
dd <- tempfile(); dir.create(dd)
write.csv(tab0, file.path(dd, 'dmft_age_bmi_model.csv'), row.names = FALSE, na = '')
r <- e$load_published_coefficients(m0, dd)
stopifnot(grepl('matches', r$msg), abs(e$predict_model(r$model, x, cv)$p - e$predict_model(m0, x, cv)$p) < 1e-12,
          identical(r$model$performance, m0$performance))   # CSV keeps 15 significant digits
# Example values shown by "Fill example values" and in the manuscript (not a participant)
ex <- e$predict_model(m0, e$EXAMPLE$markers, list(Age_DOS = e$EXAMPLE$covariates[1], BMI = e$EXAMPLE$covariates[2]))
stopifnot(abs(ex$p - plogis(published_eta(e$EXAMPLE$markers, 40, 25))) < 1e-12, round(100 * ex$p, 1) == 73.8,
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

# ---- export script -> app import (synthetic fixture) -----------------------------------
if (file.exists(eq_path)) {
q <- load_functions(eq_path)
key <- 'test_DMFT_Age_BMI_exhaustive12_beam3'
set.seed(643); n <- 120L; markers <- paste0('Marker',1:6)
raw <- matrix(exp(rnorm(n*6)),n); colnames(raw)<-markers
outcome <- factor(rep(c('LD','HD'),each=n/2),levels=c('LD','HD'))
x <- log2(raw+5e-7);mu <- colMeans(x[outcome=='LD',]);sd <- apply(x[outcome=='LD',],2,sd)
z <- sweep(sweep(x,2,mu),2,sd,'/')
d <- data.frame(z,Age_DOS=runif(n,20,70),BMI=runif(n,18,40),condL=outcome)
ref <- data.frame(Column=markers,Mean=mu,SD=sd)
# Deliberately make a THREE-marker panel best: the export script hardcodes no panel.
tables <- setNames(lapply(1:6,function(k) {
 row <- as.data.frame(setNames(as.list(markers[1:k]),paste0('Biomarker',1:k)))
 row$Balanced_Accuracy <- if(k==3) .9 else .7
 row$AUC <- .85;row$Threshold <- .47;row
}),paste0('Results_',1:6,'_Biomarkers'))
f <- glm(condL ~ Marker1+Marker2+Marker3+Age_DOS+BMI,data=d,family=binomial())
for(j in 1:3) tables$Results_3_Biomarkers[[paste0('OR_perDoubling_Biomarker',j)]] <- exp(coef(f)[markers[j]]/sd[j])
pars <- c('Script 1 registry key'=key,'Reference group'='LD','Case group (positive)'='HD',
          'Model configuration'='Biomarkers_Age_BMI','Confounders'='Age_DOS, BMI',
          '  k = 2 search'='exhaustive','n (LD)'='60','n (HD)'='60','Total n'='120')
fakefile <- tempfile(fileext='.xlsx');file.create(fakefile)
wb <- list(path=fakefile,parameters=data.frame(Parameter=names(pars),Value=unname(pars)),tables=tables)
q$read_dmft_workbook <- function(path) wb
out <- tempfile();dir.create(out)
d <- rbind(d,d[1,]);d$BMI[n+1]<-NA
b <- q$export_dmft_equation(d,ref,key,output_dir=out,results_file=fakefile)
stopifnot(length(b$source_rows)==n,b$selected_result$k==3L,b$max_link_difference<1e-8)
tab <- read.csv(file.path(out,'dmft_age_bmi_model.csv'),stringsAsFactors=FALSE)
m <- e$selected_age_bmi_model(tab)
for(i in 1:n)stopifnot(abs(e$predict_model(m,raw[i,1:3],list(Age_DOS=d$Age_DOS[i],BMI=d$BMI[i]))$p-predict(f,type='response')[i])<1e-10)
# A different model placed next to the desktop app is refused; the bundled model stays in use
r <- e$load_published_coefficients(m0,out)
stopifnot(grepl('was not used', r$msg), identical(r$model$markers$key, m0$markers$key), e$model_ready(r$model))
m2<-e$selected_age_bmi_model(e$model_to_table(m));stopifnot(e$model_ready(m2))
bad<-tab;bad$Label[bad$Term=='(Model configuration)']<-'Biomarkers_Age_Sex'
stopifnot(inherits(try(e$selected_age_bmi_model(bad),silent=TRUE),'try-error'))
bad<-tab[tab$Term!='(Threshold)',];stopifnot(inherits(try(e$selected_age_bmi_model(bad),silent=TRUE),'try-error'))
wb$parameters$Value[trimws(wb$parameters$Parameter)=='k = 2 search']<-'beam'
stopifnot(inherits(try(q$select_dmft_panel(list(wb),key),silent=TRUE),'try-error'))
cat('PASS: export script, raw-input equivalence for 120 synthetic participants, import and stale-model rejection.\n')
}
cat('All checks passed.\n')
