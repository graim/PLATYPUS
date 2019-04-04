#' pseudocode:
#'     for (v in views)
#'       load v
#'     load labels
#'
#'     divide labelled data in folds
#'     for(f in folds):
#'       unlabelled.data <- fold f
#'       delete f from labels
#'       while( not converged on labels OR no changes anymore )
#'         for( v in views )
#'           train v on all sets but f
#'           try to predict fold f
#'         add newly found labels to training data
#'         calculate performance for added labels
#'       calculate overall performance for a fold
#'     average performance over folds


#' Label learning validation
#'
#' Similar to cross-fold validation, label learning validation for platypus is used to help identify the number of iterations to run when training a platypus model, so that label learning is most effective.
#' @param view.list List of view objects 
#' @param fn.labs File containing outcome labels
#' @param k number of folds for label learning validation (similar to cross validation folds), eg. 5, default=10
#' @param i Maximal number of iterations for each platypus run, default=100
#' @param m Percent agreement required to learn a sample's class label, default=100
#' @param u Updating the accuracies of the single views in each iteration, default=FALSE
#' @param e Expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
#' @param b Label class to ignore, if any. Defaults to 'intermediate'
#' @param o Name of the folder where output is stored.
#' @param p The number of cores to use. Enables parallelization.
#' @return A list containing fold.accuracy, labelling.matrix,labelling.matrices.views
#' @keywords platypus
#' @export 
#' @examples
#' TODO show how to generate config.files and fn.labs
#' llv.platypus(config.files,fn.labs)
#' llv.platypus(view.list=view.list,fn.labs=fn.labs,no.iterations=5,majority.threshold.percent=75,output.folder='platypus_output')
llv.platypus <- function(view.list,fn.labs,llv.folds=10,no.iterations=100,majority.threshold.percent=100,expanded.output=TRUE,updating=FALSE,ignore.label='intermediate',parallel=FALSE,num.cores=25,classcol.labs=1,output.folder=NA) {
  
  ## Set debug flag on/off for testing - currently we don't use this
  #flag.debug <- TRUE
  flag.debug <- FALSE
  if(flag.debug) { print('Debug is on');flush.console() }

  ## TODO: Don't load libraries this way :)
  require(foreach)
  require(methods)
  if(parallel) {
    require(doParallel)
    cl <- makeCluster(num.cores,outfile="")
    registerDoParallel(cl, cores = num.cores)
  }

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
    write.table(labs.reduced, file= fn.labels.reduced, sep="\t",row.names=T, col.names=T, quote=FALSE)
    
    # get platypus result list
    platypus.result <- platypus(view.list=view.list, fn.labs=fn.labels.reduced, i=no.iterations, m=majority.threshold.percent,expanded.output=expanded.output,updating=updating, ignore.label=ignore.label)
    
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
  
  if (parallel) {
    llv.result.list <- foreach(k=seq(llv.folds), .export=c("platypus","drop.features","ElasticNet","RandomForest",
        "load.parameterfile","load.data","load.data.ElasticNet","load.data.RandomForest",
        "load.label.data","get.unique.labels","view.train","view.train.ElasticNet","view.train.RandomForest","view.predict","view.predict.ElasticNet",
        "view.predict.RandomForest","platypus.predict","update.accuracies","update.accuracy","update.accuracy.ElasticNet","update.accuracy.RandomForest",
        "calculate.accuracy","get.majority.weighting","get.new.labels.majorityWeighted","normalize.accuracies",
        "normalize.accuracy.log","calculate.performance","calculate.performance.view","get.labelling.performance"),
        .verbose=TRUE, .packages=c("glmnet","randomForest")) %dopar%
      do.one.llvfold(k = k)
  } else {
    llv.result.list <- lapply(seq(llv.folds), do.one.llvfold)
    #llv.result.list <- foreach(k=seq(llv.folds)) %do% do.one.llvfold(k = k)  # Old but saving it for now
  }  
  
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
      write.table(accuracy.llvfolds, file= paste0(output.folder,"/perf_llv.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
      write.table(accuracy.platypus.iterations, file= paste0(output.folder,"/perf_llv_expanded.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
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
