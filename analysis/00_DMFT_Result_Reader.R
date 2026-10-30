# Read either a detailed search workbook or its Best_By_Size/Best_Overall summary.
# Summary files contain no coefficients, CV settings or participant predictions.
dmft_summary_tables <- function(by_size, overall, path) {
  need <- c('Model','N','LD','HD','N_Biomarkers','Covariates','Panel','Balanced_Accuracy','AUC','Threshold')
  if (!all(need %in% names(by_size)) || !all(need %in% names(overall))) stop('Incomplete DMFT summary sheets.')
  if(nrow(overall)!=1L || nrow(by_size)!=6L || !setequal(by_size$N_Biomarkers,1:6)) stop('Expected one best panel at each size 1–6 and one overall winner.')
  for(d in list(by_size,overall)) {
    if(any(d$Model!='Biomarkers_Age_BMI') || any(vapply(strsplit(d$Covariates,',',fixed=TRUE),function(x)!identical(trimws(x),c('Age_DOS','BMI')),logical(1)))) stop('Expected age/BMI adjustment in every summary row.')
    if(anyNA(d[need]) || any(!is.finite(as.matrix(d[c('N','LD','HD','N_Biomarkers','Balanced_Accuracy','AUC','Threshold')])))) stop('Missing or nonfinite summary values.')
    if(any(d$N != d$LD+d$HD) || any(d$LD<1 | d$HD<1) || any(d$Balanced_Accuracy<0 | d$Balanced_Accuracy>1 | d$AUC<0 | d$AUC>1 | d$Threshold<=0 | d$Threshold>=1)) stop('Invalid summary counts or metrics.')
  }
  if(length(unique(by_size$N))!=1 || length(unique(by_size$LD))!=1 || length(unique(by_size$HD))!=1)stop('Panels use different cohorts.')
  win <- by_size[order(-by_size$Balanced_Accuracy,-by_size$AUC,by_size$N_Biomarkers)[1],]
  if(!identical(as.character(win$Panel),as.character(overall$Panel)) || win$N_Biomarkers!=overall$N_Biomarkers ||
     any(abs(as.numeric(win[1,c('N','LD','HD','Balanced_Accuracy','AUC','Threshold')])-as.numeric(overall[1,c('N','LD','HD','Balanced_Accuracy','AUC','Threshold')]))>1e-9))stop('Best_Overall disagrees with Best_By_Size.')
  key <- sub('\\.xlsx$','',sub('^DMFT_Age_BMI_Summary_','',basename(path)))
  if(!startsWith(basename(path),'DMFT_Age_BMI_Summary_') || !grepl('_conf_AgeDOS_BMI_.*_exhaustive12_beam3$',key))stop('Keep the original summary filename: it identifies the preprocessing registry.')
  tabs <- setNames(lapply(seq_len(nrow(by_size)),function(i){
    r<-by_size[i,,drop=FALSE];k<-as.integer(r$N_Biomarkers)
    markers<-trimws(strsplit(r$Panel,';',fixed=TRUE)[[1]])
    if(length(markers)!=k || any(!nzchar(markers)) || anyDuplicated(markers))stop('Invalid panel in summary row ',i+1)
    for(j in seq_len(k))r[[paste0('Biomarker',j)]]<-markers[j]
    r$.excel_row<-i+1L;r
  }),paste0('Results_',by_size$N_Biomarkers,'_Biomarkers'))
  pars<-c('Script 1 registry key'=key,'Model configuration'='Biomarkers_Age_BMI','Confounders'='Age_DOS, BMI',
          'Reference group'='LD','Case group (positive)'='HD','k = 2 search'='exhaustive',
          'Total n'=as.character(win$N),'n (LD)'=as.character(win$LD),'n (HD)'=as.character(win$HD))
  list(path=path,summary=TRUE,parameters=data.frame(Parameter=names(pars),Value=unname(pars)),
       tables=tabs,best_by_size=by_size,best_overall=overall,
       registry_key=key,source_sheets=setNames(rep('Best_By_Size',length(tabs)),names(tabs)))
}
dmft_read_results <- function(path) {
  if(!requireNamespace('readxl',quietly=TRUE))stop("Install readxl: install.packages('readxl')")
  sh<-readxl::excel_sheets(path)
  if(all(c('Best_By_Size','Best_Overall') %in% sh))return(dmft_summary_tables(
    as.data.frame(readxl::read_excel(path,sheet='Best_By_Size')),
    as.data.frame(readxl::read_excel(path,sheet='Best_Overall')),path))
  if(!'Analysis_Parameters' %in% sh)stop('Neither a DMFT summary nor a detailed search workbook: ',path)
  sheets<-grep('^Results_[1-6]_Biomarkers$',sh,value=TRUE)
  list(path=path,summary=FALSE,parameters=as.data.frame(readxl::read_excel(path,sheet='Analysis_Parameters',col_types='text')),
       tables=setNames(lapply(sheets,function(s)as.data.frame(readxl::read_excel(path,sheet=s))),sheets))
}
dmft_saved_cv <- function(saved,n) {
  ix<-saved$index;nm<-names(ix)
  if(!length(ix)||length(nm)!=length(ix)||any(!grepl('^Fold[0-9]+\\.Rep[0-9]+$',nm))||anyDuplicated(nm))stop('Invalid saved resampling names.')
  reps<-sub('^.*\\.','',nm);folds<-table(reps)
  if(length(unique(as.integer(folds)))!=1L)stop('Unequal numbers of folds per repeat.')
  for(v in ix)if(!is.numeric(v)||anyNA(v)||any(v!=as.integer(v))||anyDuplicated(v)||any(v<1|v>n)||!length(v)||length(v)>=n)stop('Invalid training indices.')
  for(r in unique(reps)) {
    tests<-unlist(lapply(ix[reps==r],function(v)setdiff(seq_len(n),v)),use.names=FALSE)
    if(!identical(sort(as.integer(tests)),seq_len(n)))stop('Each participant must be held out once per repeat.')
  }
  list(folds=unname(as.integer(folds[1])),repeats=length(folds))
}
