## We want to look at Gini coefficient for each feature in the RF models instead
## Load each model, get difference in Gini between models
## Kiley Graim
# March 2017

iter1 <- 'trainedViews_PD.0325901/trainedViews_Iteration1/'
views <- unlist(read.table('rf_names.txt', sep='\t', header=F, stringsAsFactors=F))
iter2 <- 'trainedViews_PD.0325901/trainedViews_Iteration3/'

view.diffs <- matrix(NA, nrow=length(views), ncol=5)
colnames(view.diffs) <- c('mean', 'max', 'min', 'sum', 'nfeatures')
rownames(view.diffs) <- views
view.diffs <- as.data.frame(view.diffs)

for(view in views) {

  ## Load that view for each iteration
  dat.i1 <- readRDS( paste0(iter1,view) )
  view1 <- dat.i1$model$importance[,'MeanDecreaseGini',drop=F]

  ## Load that view for each iteration
  dat.i3 <- readRDS( paste0(iter2,view) )
  view2 <- dat.i3$model$importance[,'MeanDecreaseGini',drop=F]

  ## Calculate differences in feature weights
  view.diff <- view1[,1] - view2[,1]
  view.diff <- abs(view.diff)
  view.diffs[view,'mean'] <- mean(view.diff)
  view.diffs[view,'max'] <- max(view.diff)
  view.diffs[view,'min'] <- min(view.diff)
  view.diffs[view,'sum'] <- sum(view.diff)
  view.diffs[view,'nfeatures'] <- length(view.diff)

  ## Store the differences, for if we want to plot them later
  write.table(sort(signif(view.diff,digits=5),decreasing=T), file=paste(unlist(strsplit(view, '[.]'))[1],'differences','txt',sep='.'), quote=F, row.names=T, col.names=F)
  
  #write.table(dat$model$importance[,'MeanDecreaseGini',drop=F], file=paste0(view,'iteration3','Gini.txt'), sep='\t', col.names=T, row.names=T, quote=F)

}

view.diffs <- view.diffs[with(view.diffs, order(-mean,-sum,-max)),]
write.table(view.diffs, file='view_summary_Gini_differences.tab', sep='\t', col.names=T, row.names=T, quote=F)



dat <- readRDS('trainedViews_PD.0325901/trainedViews_Iteration1/rf_DTPathwayOncogenic.RDS')
dat$model$importance[,'MeanDecreaseGini',drop=F]



