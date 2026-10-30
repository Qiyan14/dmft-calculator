# Base-R tests; synthetic values below are fixtures, never a clinical model.
# Rscript test_age_bmi.R [path/to/app.R] [path/to/04_Export_DMFT_Equation_Age_BMI.R]
args <- commandArgs(trailingOnly=TRUE)
app_path <- if(length(args)) args[1] else '06_Shiny_DMFT_Age_BMI.R'
eq_path <- if(length(args)>1) args[2] else '04_Export_DMFT_Equation_Age_BMI.R'
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
stopifnot(!e$model_ready(e$PUBLISHED_MODEL))
job <- list(p=472L,s=list(search_mode='beam',exhaustive_limit=10),parents=matrix(c(1L,2L),2))
pairs <- e$level_candidates(job,2L)
stopifnot(pairs$mode=='exhaustive',ncol(pairs$cand)==111156L,identical(pairs$cand[,ncol(pairs$cand)],c(471L,472L)))
triples <- e$level_candidates(job,3L)
stopifnot(triples$mode=='beam',ncol(triples$cand)==470L,all(triples$cand[1,]==1),all(triples$cand[2,]==2))
cat('PASS: all 111,156 pairs evaluated; beam extensions start at size 3.\n')
q <- load_functions(eq_path)
key <- 'test_DMFT_Age_BMI_exhaustive12_beam3'
set.seed(643); n <- 120L; markers <- paste0('Marker',1:6)
raw <- matrix(exp(rnorm(n*6)),n); colnames(raw)<-markers
outcome <- factor(rep(c('LD','HD'),each=n/2),levels=c('LD','HD'))
x <- log2(raw+5e-7);mu <- colMeans(x[outcome=='LD',]);sd <- apply(x[outcome=='LD',],2,sd)
z <- sweep(sweep(x,2,mu),2,sd,'/')
d <- data.frame(z,Age_DOS=runif(n,20,70),BMI=runif(n,18,40),condL=outcome)
ref <- data.frame(Column=markers,Mean=mu,SD=sd)
# Deliberately make a THREE-marker panel best: no hardcoded six-marker model.
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
# No Sex column is supplied or required.
b <- q$export_dmft_equation(d,ref,key,output_dir=out,results_file=fakefile)
stopifnot(length(b$source_rows)==n,b$selected_result$k==3L,b$max_link_difference<1e-8)
tab <- read.csv(file.path(out,'dmft_age_bmi_model.csv'),stringsAsFactors=FALSE)
m <- e$selected_age_bmi_model(tab)
for(i in 1:n)stopifnot(abs(e$predict_model(m,raw[i,1:3],list(Age_DOS=d$Age_DOS[i],BMI=d$BMI[i]))$p-predict(f,type='response')[i])<1e-10)
stopifnot(e$model_ready(e$load_published_coefficients(e$PUBLISHED_MODEL,out)$model))
m2<-e$selected_age_bmi_model(e$model_to_table(m));stopifnot(e$model_ready(m2))
cat('PASS: best panel from workbook, age/BMI cohort, raw-input refit equivalence, export/import and 120 predictions.\n')
cv <- list(Age_DOS=40,BMI=25);v<-raw[1,1:3]
curves<-e$age_bmi_curve_data(m,v,cv,18:100)
stopifnot(nrow(curves)==166,all(unique(curves$bmi)==c(25,30)))
stopifnot(abs(curves$risk[curves$age==40 & curves$bmi==25]-e$predict_model(m,v,cv)$p)<1e-12)
stopifnot(max(abs(qlogis(curves$risk[curves$bmi==30])-qlogis(curves$risk[curves$bmi==25])-5*coef(f)['BMI']))<1e-9)
bad<-tab;bad$Label[bad$Term=='(Model configuration)']<-'Biomarkers_Age_Sex'
stopifnot(inherits(try(e$selected_age_bmi_model(bad),silent=TRUE),'try-error'))
bad<-tab[tab$Term!='(Threshold)',];stopifnot(inherits(try(e$selected_age_bmi_model(bad),silent=TRUE),'try-error'))
wb$parameters$Value[trimws(wb$parameters$Parameter)=='k = 2 search']<-'beam'
stopifnot(inherits(try(q$select_dmft_panel(list(wb),key),silent=TRUE),'try-error'))
cat('PASS: age/BMI curves, stale model/search rejection and missing-threshold rejection.\nAll checks passed.\n')
plot_file<-tempfile(fileext='.pdf');pdf(plot_file)
e$draw_age_bmi_curves(curves,40,25,m$threshold);dev.off()
stopifnot(file.info(plot_file)$size>1000);unlink(plot_file)
cat('PASS: comparison plot renders.\n')
