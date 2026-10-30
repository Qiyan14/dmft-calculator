# Held-out permutation importance for the selected age/BMI model.
# Run after 03_Figures_DMFT_Age_BMI.R (which runs this script by default).
# Within each saved held-out fold, shuffle one predictor at a time, without
# refitting or reselecting the threshold/panel. Pool predictions across repeats.
# Primary measure: baseline BA minus shuffled BA, in percentage points.
# Whiskers: 2.5th–97.5th percentiles across shuffles, NOT confidence intervals.
# Negative importance is retained. Correlated predictors can share/substitute
# information; importance is model reliance, not causality or a p-value.
VI_PERMUTATIONS <- 20L
VI_SEED <- 20261030L

vi_auc <- function(p,y) {
 n1<-sum(y==1);n0<-sum(y==0)
 (sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)
}
vi_score <- function(p,y,threshold) c(BA=(mean(p[y==1]>=threshold)+mean(p[y==0]<threshold))/2,AUC=vi_auc(p,y))
dmft_permutation_importance <- function(data,markers,covariates,index,threshold,
                                       B=VI_PERMUTATIONS,seed=VI_SEED) {
 vars<-c(markers,covariates)
 if(B<2L || B!=as.integer(B))stop('Use at least two permutations.')
 if(!all(vapply(data[vars],is.numeric,logical(1))) || any(!is.finite(as.matrix(data[vars]))))stop('Importance requires finite numeric inputs.')
 if(!identical(levels(data$condL),c('LD','HD')))stop('Expected LD then HD outcome levels.')
 dmft_saved_cv(list(index=index),nrow(data))
 y<-as.integer(data$condL=='HD')
 total<-sum(vapply(index,function(tr)nrow(data)-length(tr),integer(1)))
 baseline<-numeric(total);obs<-integer(total)
 shuffled<-setNames(lapply(vars,function(v)matrix(NA_real_,total,B)),vars)
 formula<-reformulate(vars,response='condL');at<-0L;set.seed(seed)
 for(fold in names(index)) {
  train<-index[[fold]];test<-setdiff(seq_len(nrow(data)),train)
  fit<-suppressWarnings(glm(formula,data=data[train,,drop=FALSE],family=binomial(),na.action=na.fail))
  beta<-coef(fit)
  if(!isTRUE(fit$converged)||any(!is.finite(beta)))stop('Nonconverged or unidentified fold fit: ',fold)
  X<-model.matrix(delete.response(terms(fit)),data[test,,drop=FALSE]);eta<-drop(X%*%beta)
  ix<-at+seq_along(test);at<-at+length(test)
  baseline[ix]<-plogis(eta);obs[ix]<-y[test]
  for(v in vars)for(b in seq_len(B)) {
   perm<-sample.int(length(test))
   shuffled[[v]][ix,b]<-plogis(eta+beta[v]*(X[perm,v]-X[,v]))
  }
 }
 base<-vi_score(baseline,obs,threshold)
 draws<-do.call(rbind,lapply(vars,function(v){
   met<-vapply(seq_len(B),function(b)vi_score(shuffled[[v]][,b],obs,threshold),numeric(2))
   data.frame(Variable=v,Permutation=seq_len(B),BA_drop_pp=100*(base['BA']-met['BA',]),AUC_drop=base['AUC']-met['AUC',],row.names=NULL)
 }))
 tab<-do.call(rbind,lapply(vars,function(v){
  d<-draws[draws$Variable==v,];q<-quantile(d$BA_drop_pp,c(.025,.975),names=FALSE)
  data.frame(Variable=v,Type=if(v%in%markers)'Biomarker' else 'Adjustment variable',
             Mean_BA_drop_pp=mean(d$BA_drop_pp),Shuffle_P2.5=q[1],Shuffle_P97.5=q[2],
             Mean_AUC_drop=mean(d$AUC_drop),Baseline_BA=unname(base['BA']),Baseline_AUC=unname(base['AUC']),
             Threshold=threshold,Permutations=B)
 }))
 tab<-tab[order(-tab$Mean_BA_drop_pp),];rownames(tab)<-NULL
 list(summary=tab,draws=draws,baseline=base)
}
draw_dmft_importance <- function(tab,labels=tab$Variable) {
 ord<-order(tab$Mean_BA_drop_pp);d<-tab[ord,];labs<-labels[ord];yy<-seq_len(nrow(d))
 labs<-vapply(labs,function(s)paste(strwrap(s,width=39),collapse='\n'),character(1))
 old<-par(mar=c(5,18,3,1));on.exit(par(old))
 xr<-range(c(0,d$Mean_BA_drop_pp,d$Shuffle_P2.5,d$Shuffle_P97.5));pad<-max(diff(xr)*.08,.1)
 plot(NA,xlim=xr+c(-pad,pad),ylim=c(.5,nrow(d)+.5),xlab='Decrease in balanced accuracy (percentage points)',
      ylab='',yaxt='n',bty='n',main='Variable importance of the selected DMFT panel')
 abline(v=0,col='grey65',lty=2)
 cols<-ifelse(d$Type=='Biomarker','#0072B2','#D55E00')
 segments(d$Shuffle_P2.5,yy,d$Shuffle_P97.5,yy,col=cols,lwd=2)
 points(d$Mean_BA_drop_pp,yy,pch=19,col=cols,cex=1.2)
 axis(2,at=yy,labels=labs,las=1,tick=FALSE,cex.axis=.8)
 legend('bottomright',legend=c('Biomarker','Age / BMI'),col=c('#0072B2','#D55E00'),pch=19,bty='n',cex=.85)
 mtext('Whiskers: shuffle percentiles; not confidence intervals',side=1,line=3.7,cex=.8)
}
if(!isTRUE(getOption('dmft.importance.skip_run'))) {
 if(!exists('FIG_VERIFIED') || !exists('FIG_REFITS') || !exists('FIG_BEST'))stop('Run 03_Figures_DMFT_Age_BMI.R first.')
 if(!identical(FIG_VERIFIED$registry_key,registry_key) || !identical(FIG_VERIFIED$data,FIG_DAT))stop('The figure data have changed since validation.')
 best<-FIG_BEST[[as.character(BEST_OVERALL_K)]]
 cat(sprintf('\nPermutation importance: %d variables, %d shuffles, %d saved folds\n',length(best$markers)+length(.conf),VI_PERMUTATIONS,length(FIG_VERIFIED$index)))
 VI_RESULT<-dmft_permutation_importance(FIG_DAT,best$markers,.conf,FIG_VERIFIED$index,best$row$Threshold)
 if(max(abs(VI_RESULT$baseline-c(BA=best$row$Balanced_Accuracy,AUC=best$row$AUC)))>1e-6)
   stop('Importance baseline does not reproduce the selected workbook metrics; no importance output written.')
 VI_RESULT$summary$Label<-label_plain(vapply(VI_RESULT$summary$Variable,marker_label,''))
 write_table(VI_RESULT$summary,paste0('Variable_importance_',FILE_TAG))
 write_table(VI_RESULT$draws,paste0('Variable_importance_shuffles_',FILE_TAG))
 save_figure(paste0('Figure_Variable_Importance_',FILE_TAG),10,7,
             function()draw_dmft_importance(VI_RESULT$summary,VI_RESULT$summary$Label))
 writeLines(c('Method: held-out permutation importance; pooled repeated-CV balanced accuracy at the saved threshold.',
  'One variable was shuffled within each held-out fold. The fitted panel and threshold stayed fixed.',
  'Points show mean decreases in balanced accuracy; whiskers show shuffle percentiles, not confidence intervals.',
  'Negative values are retained. Correlation can redistribute importance. This analysis is conditional on the selected development model.'),
  file.path(OUT_DIR,paste0('Variable_importance_caption_',FILE_TAG,'.txt')))
}
