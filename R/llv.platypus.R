#' Label learning validation
#'
#' Similar to cross-fold validation, label learning validation for platypus is used to help identify the number of iterations to run when training a platypus model, so that label learning is most effective.
#' @param view.list List of view objects 
#' @param fn.labs File containing outcome labels
#' @param llv.folds number of folds for label learning validation (similar to cross validation folds), default=10
#' @param n.iters Maximal number of iterations for each platypus run, default=100
#' @param majority.threshold.percent Percent agreement required to learn a sample's class label, default=100
#' @param nfolds Number of cross-validation folds
#' @param expanded.output Expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
#' @param updating Updating the accuracies of the single views in each iteration, default=FALSE
#' @param ignore.label Label class to ignore, if any. Defaults to 'intermediate'
#' @param parallel Whether or not to run in parallel mode TODO remove? numcores enables anyway
#' @param classcol.labs Column containing the task label data. Default 1.
#' @param output.folder Name of the folder where output is stored.
#' @return A list containing fold.accuracy, labelling.matrix,labelling.matrices.views
#' @keywords platypus
#' @export 
llv.platypus <- function(view.list,fn.labs,llv.folds=10,n.iters=100,majority.threshold.percent=100,nfolds=10,expanded.output=TRUE,updating=FALSE,ignore.label='intermediate',parallel=FALSE,classcol.labs=1,output.folder=NA) {
  
  ## Set debug flag on/off for testing - currently we don't use this
  #flag.debug <- TRUE
  flag.debug <- FALSE
  if(flag.debug) { print('Debug is on') }

  print(paste('Ignoring labels for samples labelled:',ignore.label))

  ## Create output directory if it doesn't already exist TODO untested
  if(!dir.exists(output.folder)) { dir.create(output.folder) }
 
  ## Load the label data
  labs <- load.label.data(fn.labs,classcol.labs)
  
  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)
  
  ## TODO: should the real unlabelled data be excluded from the label learning validation? If so, it has to be deleted from the inputfiles and new matrix-files and parameter-files have to be saved 
  ###(maybe work around by marking the real unlabelled data as 'testing' like in the cv, so that they are not used by the platypus-function)
  
  ## Sort the labels in k subsets
  # Randomly create fold number for each label (balance folds to equal size)    TODO: shall we balance folds for equal class instances (based on the overall occurrence of each class in the data set)?
  # TODO: give the option to give user-defined folds as input
  # TODO: check to make sure that each sample has data in at least 1 view before assigning a fold
  fold.vec <- c(rep(seq(llv.folds), each=floor(dim(labs)[[1]]/llv.folds)),sample(seq(llv.folds), dim(labs)[[1]]-llv.folds*(floor(dim(labs)[[1]]/llv.folds)), replace=FALSE))
  fold.vec <- sample(fold.vec)
  labs$fold <- fold.vec
    
  
  ###################
  # function to do one llv-run for one fold
  ###################
  do.one.llvfold <- function(k) {
    
    print( paste('Processing fold:',k) )
    
    ## Define reduced label file
    # The cell lines with a label to be ignored (usually the intermediate label), have to be kept in order to be excluded from the feature matrices of the views
    labs.reduced <- labs[which(labs$fold != k | labs[,classcol.labs] == ignore.label),,drop=F]
    fn.labels.reduced <- tempfile()
    utils::write.table(labs.reduced, file= fn.labels.reduced, sep="\t",row.names=T, col.names=T, quote=FALSE)
    
    # get platypus result list
    platypus.result <- platypus(view.list=view.list, fn.labs=fn.labels.reduced, i=n.iters, m=majority.threshold.percent,e=expanded.output,u=updating, b=ignore.label,nfolds=nfolds)
    #platypus.result <- platypus(view.list=view.list, fn.labs=fn.labels.reduced, i=n.iters, m=majority.threshold.percent,expanded.output=expanded.output,updating=updating, ignore.label=ignore.label)
    
    ## Calculate final ll-performance
    # get rid of unused iterations
    labelling.matrix <- platypus.result$labelling.matrix[,!apply(is.na(platypus.result$labelling.matrix),2,all),drop=F]
    
    perf.llvfold <- get.labelling.performance(labelling.matrix[,dim(labelling.matrix)[[2]],drop=F],labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),classcol.labs,drop=F],unique.labels)
    perf.llvfold <- c(k,perf.llvfold[1,-1])
    names(perf.llvfold)[1] <- "llv.fold"
    
    ## If expanded output flag is on - calculate performance in each iteration
    if(expanded.output){
      perf.iterations <- get.labelling.performance(labelling.matrix,labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),classcol.labs,drop=F],unique.labels)
      llv.fold <- rep(k,dim(perf.iterations)[[1]])
      
      iteration.information <- platypus.result$iteration.information
#      if(weighting){
        weighting.threshold.upper <- c()
        weighting.threshold.lower <- c()
        weighting.threshold <- c()
        no.ids.labelled <- c()
        no.ids.unlabelled <- c()
        no.new.labelled <- c()
        no.ids.left.unlabelled <- c()

        for(i in seq(length(iteration.information))){
          weighting.threshold.upper <- c(weighting.threshold.upper,iteration.information[[i]]$weighting.threshold.upper)
          weighting.threshold.lower <- c(weighting.threshold.lower,iteration.information[[i]]$weighting.threshold.lower)
          weighting.threshold <- c(weighting.threshold,iteration.information[[i]]$weighting.threshold)
          no.ids.labelled <- c(no.ids.labelled,iteration.information[[i]]$no.ids.labelled)
          no.ids.unlabelled <- c(no.ids.unlabelled,iteration.information[[i]]$no.ids.unlabelled)
          no.new.labelled <- c(no.new.labelled,iteration.information[[i]]$no.new.labelled)
          no.ids.left.unlabelled <- c(no.ids.left.unlabelled,iteration.information[[i]]$no.ids.left.unlabelled)
        }
        
        perf.iterations <- cbind(llv.fold,weighting.threshold.upper,weighting.threshold.lower,weighting.threshold,no.ids.labelled,no.ids.unlabelled,no.new.labelled,no.ids.left.unlabelled,perf.iterations)

#      } else {
#        iterations.woChange <- c()
#        majority <- c()
#        majority.missingData <- c()
#        majority.threshold <- c()
#        no.ids.labelled <- c()
#        no.ids.unlabelled <- c()
#        no.new.labelled <- c()
#        no.ids.left.unlabelled <- c()
#
#        for(i in 1:length(iteration.information)){
#          iterations.woChange <- c(iterations.woChange,iteration.information[[i]]$iterations.woChange)
#          majority <- c(majority,iteration.information[[i]]$majority)
#          majority.missingData <- c(majority.missingData,iteration.information[[i]]$majority.missingData)
#          majority.threshold <- c(majority.threshold,iteration.information[[i]]$majority.threshold)
#          no.ids.labelled <- c(no.ids.labelled,iteration.information[[i]]$no.ids.labelled)
#          no.ids.unlabelled <- c(no.ids.unlabelled,iteration.information[[i]]$no.ids.unlabelled)
#          no.new.labelled <- c(no.new.labelled,iteration.information[[i]]$no.new.labelled)
#          no.ids.left.unlabelled <- c(no.ids.left.unlabelled,iteration.information[[i]]$no.ids.left.unlabelled)
#        }
#        
#        perf.iterations <- cbind(llv.fold,iterations.woChange,majority,majority.missingData,majority.threshold,no.ids.labelled,no.ids.unlabelled,no.new.labelled,no.ids.left.unlabelled,perf.iterations)
#
#      } # end if(weighting)
      
    } # end if(expanded.output)

    return.list <- list(perf.llvfold=perf.llvfold)
    if(expanded.output){
      return.list <- list(perf.llvfold=perf.llvfold,perf.iterations=perf.iterations,labelling.matrix=platypus.result$labelling.matrix,labelling.matrices.views=platypus.result$labelling.matrices.views)
    }
    return(return.list) 
  } # end function do.one.llvfold
 
  if(parallel) { llv.result.list <- parallel::mclapply(seq(llv.folds), do.one.llvfold) }
          else { llv.result.list <- lapply(seq(llv.folds), do.one.llvfold) }
  
  ## Get output together
  # Collect accuracy over llv-iterations
  accuracy.llvfolds <- c()
  for(k in seq(llv.folds)){
    accuracy.llvfolds <- rbind(accuracy.llvfolds,llv.result.list[[k]]$perf.llvfold)
  }

  if(expanded.output){
    ## Collect accuracy over platypus-iterations and labelling information
    accuracy.platypus.iterations <- c()
    for(k in seq(llv.folds)){
      accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations,llv.result.list[[k]]$perf.iterations)
    }
    
    if(!is.na(output.folder)){
      utils::write.table(accuracy.llvfolds, file= paste0(output.folder,"/perf_llv.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
      utils::write.table(accuracy.platypus.iterations, file= paste0(output.folder,"/perf_llv_expanded.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
    }
    
    labelling.matrix.llvlist <- list()
    labelling.matrices.views.llvlist <- list()
    for(k in seq(llv.folds)){
      labelling.matrix.llvlist[[k]] <- llv.result.list[[k]]$labelling.matrix
      labelling.matrices.views.llvlist[[k]] <- llv.result.list[[k]]$labelling.matrices.views
    }
    
    save(labelling.matrix.llvlist,file =paste0(output.folder,"/labelling.matrix.llvlist.Rdata") )
    save(labelling.matrices.views.llvlist,file =paste0(output.folder,"/labelling.matrices.views.llvlist.Rdata") )
  }

  # TODO: this should return an LLV object
  return( list(fold.accuracy=accuracy.llvfolds, labelling.matrix=labelling.matrix.llvlist,labelling.matrices.views=labelling.matrices.views.llvlist) )
}
