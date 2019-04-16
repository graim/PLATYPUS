# Pseudocode:
#     get platypus parameters and view information
#     load labelled data
#     divide labelled data in k subsets
#     for fold.i in 1:k
#       hold out subset fold.i of labelled data
#       train platypus.R on all-but-fold.i labelled (and unlabelled) data with the given views   # unlabelled data is optional, but making use of unlabelled data is one of the biggest advantages of platypus
#       predict fold.i subset on trained platypus views



#' k-fold cross validation for platypus
#'
#' @param view.list List of view objects 
#' @param fn.labs File containing outcome labels
#' @param classcol.labs Which column from the labels file to use for learning
#' @param cv.folds number of folds for label learning validation (similar to cross validation folds), eg. 5, default=10
#' @param n.iters Maximal number of iterations for each platypus run, default=100
#' @param majority.threshold.percent Percent agreement required to learn a sample's class label, default=100
#' @param expanded.output Expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
#' @param updating Updating the accuracies of the single views in each iteration, default=FALSE
#' @param ignore.label Label class to ignore, if any. Defaults to 'intermediate'
#' @param parallel Whether or not to run in parallel mode.
#' @param num.cores The number of cores to use. Enables parallelization.
#' @param output.folder Name of the folder where output is stored.
#' @return A list containing fold.accuracy, labelling.matrix,labelling.matrices.views
#' @keywords platypus
#' @export
cv.platypus <- function(view.list,fn.labs,classcol.labs=1,cv.folds=10,n.iters=100,majority.threshold.percent=100,expanded.output=FALSE,updating=FALSE,ignore.label='intermediate',parallel=FALSE,num.cores=25,output.folder=NA) {
 
  ## Set debug flag on/off for testing
  flag.debug <- TRUE
  #flag.debug <- FALSE
  if(flag.debug) { print('Debug is on') }


  ## Load libraries, install if not already installed
  ## TODO: Move this into the package installation, then just load libraries normally
 require(foreach)
  require(methods)

#  if(!require(foreach)) {
#    install.packages('foreach')
#   library(foreach) 
#  }
#  requireNamespace(methods) # TODO untested # TODO untested
  
  # set parallel background if parallel flag is set
#  if(parallel){
    ## TODO: Move this into the package installation, then just load libraries normally
 #   if(!require(doParallel)) {
 #     install.packages('doParallel')
  #    library(doParallel)
  #  }
#    cl <- makeCluster(num.cores,outfile="")
#    registerDoParallel(cl, cores = num.cores)
#  }

  ## Create output directory if it doesn't already exist
  if(!dir.exists(output.folder)) { dir.create(output.folder) }  

  ## Load the label data
  if(flag.debug) { print('Loading labels') }
  labs <- load.label.data(fn.labs,classcol.labs)
  
  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)

  ## Sort the labels in k subsets
  # TODO: check to make sure that each sample has data in at least 1 view before assigning a fold
  if(is.numeric(cv.folds)){
    print(paste('Number of cv folds provided:',cv.folds))
    fold.vec <- c(rep(seq(cv.folds), each=floor(dim(labs)[[1]]/cv.folds)),sample(seq(cv.folds), dim(labs)[[1]]-cv.folds*(floor(dim(labs)[[1]]/cv.folds)), replace=FALSE))
    fold.vec <- sample(fold.vec)
    labs$fold <- fold.vec
  } else if(is.vector(cv.folds) & length(cv.folds)==nrow(labs)) { 
    print('Using provided cv fold assignments')
    labs$fold <- cv.folds
  } else {
    print('ERROR: cv.folds must be a vector or integer')
    quit(save='no',status=1)
  }
  
  ###################
  # function to do a cv run for one fold
  ###################
  do.one.cvfold <- function(k) {
    
    print( paste('Processing fold:',k) ) 

    ## Define reduced label file
    # mark labels of the hold out data set with 'testing', so that they can also be excluded from the unlabelled data in the platypus.R 
    labs.reduced <- labs
    labs.reduced[which(labs$fold == k),classcol.labs] <- 'testing'
    fn.labels.reduced <- tempfile()
    write.table(labs.reduced, file= fn.labels.reduced, sep="\t",row.names=T, col.names=T, quote=FALSE)
#    print(apply(labs.reduced,2,table)) 
#    print(dim(labs))
#    print(dim(labs.reduced))
    if(flag.debug) { print(table(labs.reduced)) }
 
    # get platypus result list
    platypus.result <- platypus(view.list=view.list, fn.labs=fn.labels.reduced, i=n.iters, m=majority.threshold.percent,expanded.output=expanded.output,updating=updating,ignore.label=ignore.label)
 
    # TODO: THIS NEXT CODE BLOCK IS THE BROKEN BIT FOR SSC EXAMPLE
    # Error in na.fail.default(list(labels = c(`11031.p1` = 1L, `11033.p1` = 2L,  : 
    #   missing values in object
    ## Predict hold-out data subset with platypus result
    test.ids <- rownames(labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),,drop=F])
    predictions <- platypus.predict(platypus.result$final.views, platypus.result$weighting.threshold, test.ids,unique.labels,labs[,classcol.labs,drop=FALSE])

    if(flag.debug) { print('Held out data subset predictions made successfully') }
    if(flag.debug) { print(head(predictions)) }
    
    
    ## Calculate final platypus performance for each cv-fold
    perf.all.agree <- calculate.performance(predictions[,c("final","category.all")],labs[,classcol.labs,drop=F],unique.labels)   
    if(flag.debug) { print('Calculated all.agree') }
    perf.majority.agree <- calculate.performance(predictions[,c("final","category.majority")],labs[,classcol.labs,drop=F],unique.labels)
    if(flag.debug) { print('Calculated majority.agree') }
    
    ## Calculate final performance values for each view
    perf.views <- list()
    accuracy.views <- c()
    balanced.accuracy.views <- c()
    for(view.i in seq(length(platypus.result$final.views))){
      perf.view <- calculate.performance.view(predictions[,view.i,drop=F],labs[,classcol.labs,drop=F],unique.labels)
      perf.views[[view.i]] <- perf.view
      accuracy.views <- c(accuracy.views,perf.view$accuracy)
      balanced.accuracy.views <- c(balanced.accuracy.views,perf.view$balanced.accuracy)
    }
    
    accuracy.cvfolds <- c(k,perf.all.agree$accuracy,perf.all.agree$balanced.accuracy,perf.all.agree$coverage
                          ,perf.majority.agree$accuracy,perf.majority.agree$balanced.accuracy,perf.majority.agree$coverage
                          ,accuracy.views, balanced.accuracy.views)
    performance.cvfolds <- list(perf.all.agree=perf.all.agree,perf.majority.agree=perf.majority.agree,perf.views=perf.views)
    
    ## Calculate performance in each iteration, if expanded output flag is on
    if(expanded.output){
      iteration.information <- platypus.result$iteration.information
      performance.iterations <- list()
      accuracy.platypus.iterations <- c()
      for(i in seq(length(iteration.information))){
        ## Predict hold-out data subset with platypus view from iteration
        view.list <- iteration.information[[i]]$view.list
        test.ids <- rownames(labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),,drop=F])
        predictions <- platypus.predict(view.list, platypus.result$weighting.threshold, test.ids, unique.labels,labs[,classcol.labs,drop=FALSE]) 
        
        ## Calculate performance of the iteration platypus and views
        if(flag.debug) { print('Caluculating iteration performances') }
        perf.all.agree <- calculate.performance(predictions[,c("final","category.all")],labs[,classcol.labs,drop=F],unique.labels)   
        perf.majority.agree <- calculate.performance(predictions[,c("final","category.majority")],labs[,classcol.labs,drop=F],unique.labels)
        
        perf.views <- list()
        accuracy.views <- c()
        balanced.accuracy.views <- c()
        weighting.views <- c()
        weighting.norm.views <- c()
        for(view.i in seq(length(view.list))){
          perf.view <- calculate.performance.view(predictions[,view.i,drop=F],labs[,classcol.labs,drop=F],unique.labels)
          perf.views[[view.i]] <- perf.view
          accuracy.views <- c(accuracy.views,perf.view$accuracy)
          balanced.accuracy.views <- c(balanced.accuracy.views,perf.view$balanced.accuracy)
          
          weighting.views <- c(weighting.views,view.list[[view.i]]$acc)
          weighting.norm.views <- c(weighting.norm.views,view.list[[view.i]]$acc.norm)
        }
        
        accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations
                                         ,c(k,i
                                            ,iteration.information[[i]]$weighting.threshold.upper,iteration.information[[i]]$weighting.threshold.lower,iteration.information[[i]]$weighting.threshold
                                            ,weighting.views,weighting.norm.views
                                            ,iteration.information[[i]]$no.ids.labelled,iteration.information[[i]]$no.ids.unlabelled
                                            ,iteration.information[[i]]$no.new.labelled, iteration.information[[i]]$no.ids.left.unlabelled
                                            ,perf.all.agree$accuracy,perf.all.agree$balanced.accuracy,perf.all.agree$coverage
                                            ,perf.majority.agree$accuracy,perf.majority.agree$balanced.accuracy,perf.majority.agree$coverage
                                            ,accuracy.views,balanced.accuracy.views
                                         )
        )
        
        performance.iterations[[i]] <- list(perf.all.agree=perf.all.agree,perf.majority.agree=perf.majority.agree,perf.views=perf.views) 
      }
    } # end if expanded.output
    return.list <- list(accuracy.cvfolds=accuracy.cvfolds,performance.cvfolds=performance.cvfolds)
    if(expanded.output){
      return.list <- list(accuracy.cvfolds=accuracy.cvfolds,performance.cvfolds=performance.cvfolds,accuracy.platypus.iterations=accuracy.platypus.iterations,performance.iterations=performance.iterations
                          ,labelling.matrix=platypus.result$labelling.matrix,labelling.matrices.views=platypus.result$labelling.matrices.views)
    } 
    if(flag.debug) { print('leaving do.one.cvfold') }
    return(return.list)
  } # end do.one.cvfold fxn

  if(parallel) { cv.result.list <- parallel::mclapply(seq(cv.folds), do.one.cvfold) }
          else { cv.result.list <- lapply(seq(cv.folds), do.one.cvfold) }

  if(flag.debug) { print('out of cv folds loop') }

  ## Collect accuracy over cv-iterations
  accuracy.cvfolds <- c()
  performance.cvfolds <- list()
  for(k in seq(cv.folds)){
    accuracy.cvfolds <- rbind(accuracy.cvfolds,cv.result.list[[k]]$accuracy.cvfolds)
    performance.cvfolds[[k]] <- cv.result.list[[k]]$performance.cvfolds
  }
  
  colnames(accuracy.cvfolds) <- c("cv.fold","accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree","accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                  ,paste0("accuracy.view.",seq(length(view.list))),paste0("balanced.accuracy.view.",seq(length(view.list))))

  if(flag.debug) { print('collected cv iteration performances') }

  if(!is.na(output.folder)) {
    write.table(accuracy.cvfolds,file=file.path(output.folder,"perf_platypus.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
    save(performance.cvfolds,file=file.path(output.folder,"performance.cvfolds.Rdata") )
  }


  if(flag.debug) { print('Stored results to file') }

  if(expanded.output){
    if(flag.debug) { print('Entering expanded.output loop') }
    ## Collect accuracy over platypus-iterations
    accuracy.platypus.iterations <- c()
    performance.iterations <- list()
    for(k in seq(cv.folds)){
      accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations,cv.result.list[[k]]$accuracy.platypus.iterations)
      performance.iterations[[k]] <- cv.result.list[[k]]$performance.iterations
    }
    
      colnames(accuracy.platypus.iterations) <- c("cv.fold","iteration","weighting.threshold.upper","weighting.threshold.lower","weighting.threshold"
                                             ,paste0("weighting.view.",seq(length(view.list))),paste0("weighting.norm.view.",seq(length(view.list)))
                                             ,"no.ids.labelled","no.ids.unlabelled"
                                             ,"no.new.labelled", "no.ids.left.unlabelled"
                                             ,"accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree"
                                             ,"accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                             ,paste0("accuracy.view.",seq(length(view.list))),paste0("balanced.accuracy.view.",seq(length(view.list))))

    if(flag.debug) { print('Renamed the results') }
    if(!is.na(output.folder)) { 
      write.table(accuracy.platypus.iterations, file= file.path(output.folder,"perf_platypus_expanded.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
      save(performance.iterations,file =file.path(output.folder,"performance.iterations.Rdata") )
    }
    if(flag.debug) { print('Renamed the results') }
    
    labelling.matrix.cvlist <- list()
    labelling.matrices.views.cvlist <- list()
    for(k in seq(cv.folds)){
      labelling.matrix.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrix
      labelling.matrices.views.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrices.views
    }
    if(flag.debug) { print('Stored the labelling matrix') }
   
    if(!is.na(output.folder)) { 
      save(labelling.matrix.cvlist,file =file.path(output.folder,"labelling.matrix.cvlist.Rdata") )
      save(labelling.matrices.views.cvlist,file =file.path(output.folder,"labelling.matrices.views.cvlist.Rdata") )
    }
    if(flag.debug) { print('Wrote to files the labelling matrix') }
    
  } # end if(expanded.output)

  return( list(fold.accuracy=accuracy.cvfolds, labelling.matrix=labelling.matrix.cvlist, labelling.matrices.views=labelling.matrices.views.cvlist) )

} # end cv.platypus

