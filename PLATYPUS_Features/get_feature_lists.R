## Kiley Graim
## Feb 2017
## Grab lists of features from 2 PLATYPUS iterations 


iter1 <- 'trainedViews_PD.0325901/trainedViews_Iteration1/'
iter2 <- 'trainedViews_PD.0325901/trainedViews_Iteration3/'

views <- unlist(read.table('view_names.txt', sep='\t', header=F, stringsAsFactors=F))

view.diffs <- matrix(NA, nrow=length(views), ncol=5)
colnames(view.diffs) <- c('mean', 'max', 'min', 'sum', 'nfeatures')
rownames(view.diffs) <- views
view.diffs <- as.data.frame(view.diffs)

## For each view
for(view in views) {

  ## Load that view for each iteration
  temp <- readRDS( paste0(iter1,view) )
  view1 <- as(temp, 'matrix')

  temp <- readRDS( paste0(iter2,view) )
  view2 <- as(temp, 'matrix')

  ## Calculate difference in weight and store
  view.diff <- view1[,1] - view2[,1]
  view.diff <- abs(view.diff)
  view.diffs[view,'mean'] <- mean(view.diff)
  view.diffs[view,'max'] <- max(view.diff)
  view.diffs[view,'min'] <- min(view.diff)
  view.diffs[view,'sum'] <- sum(view.diff)
  view.diffs[view,'nfeatures'] <- length(view.diff)

  ## Store the differences, for if we want to plot them later
  write.table(sort(signif(view.diff,digits=5),decreasing=T), file=paste(unlist(strsplit(view, '[.]'))[1],'differences','txt',sep='.'), quote=F, row.names=T, col.names=F)
}

view.diffs <- view.diffs[with(view.diffs, order(-mean,-sum,-max)),]
write.table(view.diffs, file='view_summary_differences.tab', sep='\t', col.names=T, row.names=T, quote=F)
