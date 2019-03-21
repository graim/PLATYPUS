##################################################################################################
### Given expanded output from CV tests, generate plots showing AUC change & labels learned
######################################################################################
##  Arguments: 
##    fn.labs             = filename of the outcome labels
##    folder              = output directory name 
##    plot.single.cvfolds = whether or not to plot lines for each cv fold, in addition to the genearl plots
################################################################

plot.cv <- function(fn.labs, folder, plot.single.cvfolds=FALSE) {

  ## Check arguments exist
  if(!file.exists(fn.labs)) { stop(paste("ERROR: Labels file does not exist:",fn.labs)) }
  if(!file.exists(folder)) { stop(paste("ERROR: Folder does not exist:",folder)) }

  ## Don't plot each cv fold
  #plot.single.cvfolds <- FALSE
  
  ## Set up some information & fn.cv.resultss
  fn.cv.results <- file.path(folder,'perf_platypus_expanded.tab') 

  ## Gracefully quit if the requisite output file is missing    
  if(!file.exists(fn.cv.results)){ stop('ERROR: missing PLATYPUS CV performance file'); flush.console() }
  accuracy.platypus.iterations <- read.table(fn.cv.results, sep='\t',header=T)
  
  ## Load the labels, remove unlabeled samples
  labs <- read.table(fn.labs, sep='\t',header=T, row.names=1, check.names=F, stringsAsFactors = FALSE)
  labs <- labs[which(!(is.na(labs[,1]))),,drop=F]
  compound <- colnames(labs)[1] # Grab the name of the prediction task from the column name
  no.views <- length(grep('balanced.accuracy.view.', colnames(accuracy.platypus.iterations))) # Have this value for all views, so pull num views from that

  ## Set output plot name
  pdf.name <- paste('cv_performance',compound, no.views, 'views.pdf', sep='_')
  pdf.name <- file.path(folder,pdf.name)
        
  cv.folds <- max(accuracy.platypus.iterations$cv.fold)
        
  ## calculate no. correct predictions from accuracy and coverage for all.agree and majority.agree
  no.correctly.predicted.balanced.all <- accuracy.platypus.iterations$coverage.all.agree * accuracy.platypus.iterations$balanced.accuracy.all.agree
  accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.all)
        
  no.correctly.predicted.balanced.majority <- accuracy.platypus.iterations$coverage.majority.agree * accuracy.platypus.iterations$balanced.accuracy.majority.agree
  accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.majority)
        
  ## extend accuracy data table to maximal number of iterations
  no.iterations <- max(accuracy.platypus.iterations$iteration)
      
  data.extended <- c()  
  for(cv in 1:cv.folds){
    iterations.in.fold <- max(accuracy.platypus.iterations[accuracy.platypus.iterations$cv.fold==cv,"iteration"])
        
    accuracy.platypus.iterations.cvfold <- accuracy.platypus.iterations[which(accuracy.platypus.iterations$cv.fold == cv),]
    # fill in the empty iterations with the result from the last iteration
    if(iterations.in.fold < no.iterations){
      for(i in (iterations.in.fold + 1):no.iterations){
        accuracy.platypus.iterations.cvfold <- rbind(accuracy.platypus.iterations.cvfold,accuracy.platypus.iterations.cvfold[dim(accuracy.platypus.iterations.cvfold)[[1]],])
      } # end for no.iterations
    } # end if
    data.extended <- rbind(data.extended,accuracy.platypus.iterations.cvfold)
  } # end for cv.folds
      
      
  # accuracy of MVL
  platypus.acc.avg.all <- c()
  platypus.acc.avg.majority <- c()
  platypus.acc.sd.all <- c()
  platypus.acc.sd.majority <- c()
  for(i in 1:no.iterations){
    acc.vec.all <- c()
    acc.vec.majority <- c()
    for(cv in 1:cv.folds){
      acc.vec.all <- c(acc.vec.all,data.extended[which(data.extended$cv.fold==cv),"balanced.accuracy.all.agree"][i])
      acc.vec.majority <- c(acc.vec.majority,data.extended[which(data.extended$cv.fold==cv),"balanced.accuracy.majority.agree"][i])
    } # end for cv.folds
    platypus.acc.avg.all <- c(platypus.acc.avg.all, mean(acc.vec.all))
    platypus.acc.sd.all <- c(platypus.acc.sd.all, sd(acc.vec.all))
    platypus.acc.avg.majority <- c(platypus.acc.avg.majority, mean(acc.vec.majority))
    platypus.acc.sd.majority <- c(platypus.acc.sd.majority, sd(acc.vec.majority))
  } # end for no.iterations
      
  ## accuracy of views
  views.acc.avg <- c()
  views.acc.sd <- c()
  for(view.i in 1:no.views){
  #for(view.i in 1:length(grep("balanced.accuracy.view.", colnames(data.extended)))){ #TODO
    view.acc.avg <- c()
    view.acc.sd <- c()
    for(i in 1:no.iterations){
      acc.vec <- c()
      for(cv in 1:cv.folds){
        acc.vec <- c(acc.vec,data.extended[which(data.extended$cv.fold==cv),paste0("balanced.accuracy.view.",view.i)][i]) 
      } # end for cv.folds
      view.acc.avg <- c(view.acc.avg, mean(acc.vec))
      view.acc.sd <- c(view.acc.sd, sd(acc.vec))
    } # end for no.iterations
    views.acc.avg <- cbind(views.acc.avg,view.acc.avg)
    views.acc.sd <- cbind(views.acc.sd,view.acc.sd)
  } # end for view.
      
  ## prediction coverage
  platypus.cov.avg.all <- c()
  platypus.cov.avg.majority <- c()
  platypus.cov.sd.all <- c()
  platypus.cov.sd.majority <- c()
  for(i in 1:no.iterations){
    cov.vec.all <- c()
    cov.vec.majority <- c()
    for(cv in 1:cv.folds){
      cov.vec.all <- c(cov.vec.all,data.extended[which(data.extended$cv.fold==cv),"coverage.all.agree"][i])
      cov.vec.majority <- c(cov.vec.majority,data.extended[which(data.extended$cv.fold==cv),"coverage.majority.agree"][i])
    } # end for cv.folds
    platypus.cov.avg.all <- c(platypus.cov.avg.all, mean(cov.vec.all))
    platypus.cov.sd.all <- c(platypus.cov.sd.all, sd(cov.vec.all))
    platypus.cov.avg.majority <- c(platypus.cov.avg.majority, mean(cov.vec.majority))
    platypus.cov.sd.majority <- c(platypus.cov.sd.majority, sd(cov.vec.majority))
  } # end for no.iterations
      
  ## inclusion of unlabelled data in training
  platypus.inclusion.avg <- c()
  platypus.inclusion.sd <- c()
  for(i in 1:no.iterations){
    inclusion.vec <- c()
    for(cv in 1:cv.folds){
      inclusion.vec <- c(inclusion.vec,data.extended[which(data.extended$cv.fold==cv),"no.ids.labelled"][i])
    } # end for cv.folds
    platypus.inclusion.avg <- c(platypus.inclusion.avg, mean(inclusion.vec))
    platypus.inclusion.sd <- c(platypus.inclusion.sd, sd(inclusion.vec))
  } # end for no.iterations
      
      
  # thresholds of MVL labelling
  ll.upper.avg <- c()
  ll.upper.sd <- c()
  for(i in 1:no.iterations){
    upper.vec <- c()
    for(cv in 1:cv.folds){
      upper.vec <- c(upper.vec,data.extended[which(data.extended$cv.fold==cv),"weighting.threshold.upper"][i])
    } # end for cv.folds
    ll.upper.avg <- c(ll.upper.avg, mean(upper.vec))
    ll.upper.sd <- c(ll.upper.sd, sd(upper.vec))
  } # end for no.iterations
      
  ll.lower.avg <- c()
  ll.lower.sd <- c()
  for(i in 1:no.iterations){
    lower.vec <- c()
    for(cv in 1:cv.folds){
      lower.vec <- c(lower.vec,data.extended[which(data.extended$cv.fold==cv),"weighting.threshold.lower"][i])
      #lower.vec <- c(lower.vec,data.extended[which(data.extended[,"cv.fold"] == cv),"weighting.threshold.lower"][i])
    } # end for cv.folds
    ll.lower.avg <- c(ll.lower.avg, mean(lower.vec))
    ll.lower.sd <- c(ll.lower.sd, sd(lower.vec))
  } # end for no.iterations
      
  #############
  # PLOTTING
  #############
  
  ## Create the pdf where all plots are saved
  pdf(file=pdf.name,width=8,height=7)
      
  green.transparent <- rgb(col2rgb("green")[1,1],col2rgb("green")[2,1],col2rgb("green")[3,1],100,maxColorValue=255)
  par(mar=c(4, 4, 0, 0) + 0.1,par(xaxs='i',yaxs='i')) 
  layout(as.matrix(rbind(1,2,3)), widths=c(1,1,1), heights=c(1,1,1.25))
  par(mar=c(1, 4, 0.6, 0.1) + 0.1)

  ## prediction accuracy plot
  plot(1:no.iterations,platypus.acc.avg.all,ylim=c(0.5,1),xlim=c(1,no.iterations),type='l',lwd=2,col="red3",xlab="",ylab="CV Accuracy",xaxt='n')
  grid()

  if(plot.single.cvfolds){
    for(cv in 1:cv.folds){
      lines(1:no.iterations,data.extended[which(data.extended$cv.fold==cv),"balanced.accuracy.all.agree"],type='l',col="tomato")
    } # end for cv.folds
  } # end if
      
  upper <- platypus.acc.avg.majority + platypus.acc.sd.majority
  lower <- platypus.acc.avg.majority - platypus.acc.sd.majority
  lightblue.transparent <- rgb(col2rgb("lightblue")[1,1],col2rgb("lightblue")[2,1],col2rgb("lightblue")[3,1],100,maxColorValue=255)
      
  polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)), col = lightblue.transparent, border = NA)
      
  if(plot.single.cvfolds){
    for(cv in 1:cv.folds){
      lines(1:no.iterations,data.extended[which(data.extended$cv.fold==cv),"balanced.accuracy.majority.agree"],type='l',col="lightblue")
    } # end for
  } # end if
     
  #red polygon again
  upper <- platypus.acc.avg.all + platypus.acc.sd.all
  lower <- platypus.acc.avg.all - platypus.acc.sd.all
  tomato.transparent <- rgb(col2rgb("tomato")[1,1],col2rgb("tomato")[2,1],col2rgb("tomato")[3,1],100,maxColorValue=255)
  polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)), col = tomato.transparent, border = NA)
      
  #single views
  for(view.i in 1:dim(views.acc.avg)[[2]]){
    lines(1:no.iterations,views.acc.avg[,view.i],pch=20,col=green.transparent,lwd=2)
  } # end for

  # lines again
  lines(1:no.iterations,platypus.acc.avg.majority,col="blue",lwd=2)
  lines(1:no.iterations, platypus.acc.avg.all, col = 'red3',lwd=2)
      
      
  legend("bottomright",legend=c("all agree","majority agrees","single views"), lty=c(1,1,1),lwd=c(2,2,2), col=c("red3","blue","green"), bty='n')
      
  ## inclusion of unlabelled data plot
  par(mar=c(1, 4, 0.6, 0.1) + 0.1)
  plot(1:no.iterations,platypus.inclusion.avg,ylim=c(0,max(platypus.inclusion.avg)),type='l',lwd=2,col="black",xlab="Iterations",ylab="# Labelled Samples",xaxt='n')
  grid()
      
  upper <- platypus.inclusion.avg + platypus.inclusion.sd
  lower <- platypus.inclusion.avg - platypus.inclusion.sd
      
  polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)), col = "grey", border = NA)
  lines(1:no.iterations,platypus.inclusion.avg,col="black",lwd=2)
      
  if(plot.single.cvfolds){
    for(cv in 1:cv.folds){
      lines(1:no.iterations,data.extended[which(data.extended$cv.fold==cv),"no.ids.labelled"],type='l',col="grey")
      lines(1:no.iterations,platypus.inclusion.avg,ylim=c(0,max(platypus.inclusion.avg)),type='l',lwd=2,col="black") # Re-plot black line so it's on top
    } # end for
  } # end if
      
      
  ### Threshold plot
  par(mar=c(4, 4, 0.6, 0.1) + 0.1)
  magenta.transparent <- rgb(col2rgb("magenta")[1,1],col2rgb("magenta")[2,1],col2rgb("magenta")[3,1],100,maxColorValue=255)
      
  plot(1:no.iterations,ll.upper.avg,ylim=c(0,max(ll.upper.avg+ll.upper.sd)),type='l',lwd=2,col="magenta3",xlab="Iteration",ylab="Voting Sum")
  grid()
      
  upper <- ll.upper.avg + ll.upper.sd
  lower <- ll.upper.avg - ll.upper.sd
      
  polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)), col = magenta.transparent, border = NA)
  lines(1:no.iterations,ll.upper.avg,col="red3",lwd=2)
      
  if(plot.single.cvfolds){
    for(cv in 1:cv.folds){
      lines(1:no.iterations,data.extended[which(data.extended$cv.fold==cv),"weighting.threshold.upper"],type='l',col="tomato")
    } # end for
    lines(1:no.iterations, ll.upper.avg, col = 'red3',lwd=2)
  } # end if
      
  #lower
  upper <- ll.lower.avg + ll.lower.sd
  lower <- ll.lower.avg - ll.lower.sd
  lightblue.transparent <- rgb(col2rgb("lightblue")[1,1],col2rgb("lightblue")[2,1],col2rgb("lightblue")[3,1],100,maxColorValue=255)
  polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)), col = green.transparent, border = NA)
      
  if(plot.single.cvfolds){
    for(cv in 1:cv.folds){
      lines(1:no.iterations,data.extended[which(data.extended$cv.fold==cv),"weighting.threshold.lower"],type='l',col="lightblue")
    } # end for
  } # end if
  lines(1:no.iterations,ll.lower.avg,type='l',lwd=2,col="green3")
  legend("left",legend=c("min. voting for favored label","max. voting for contrary label"), lty=c(1,1),lwd=c(2,2), col=c("magenta3","green3"), bty='n')

  dev.off()

  print(paste('Finished! Success. Plot saved as:', pdf.name));flush.console()

} # end plot.cv function
