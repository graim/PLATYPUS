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
#' @param fn.views List of view files
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
#' cv.platypus(fn.views=config.files,fn.labs=fn.labs,no.iterations=5,majority.threshold.percent=75,output.folder='platypus_output')
#' cv.platypus(config.files,fn.labs)
cv.platypus <- function(fn.views,fn.labs,classcol.labs=1,cv.folds=10,no.iterations=100,majority.threshold.percent=100,expanded.output=FALSE,updating=FALSE,ignore.label='intermediate',parallel=FALSE,num.cores=25,output.folder=NA) {
 
  ## Set debug flag on/off for testing
  flag.debug <- TRUE
  #flag.debug <- FALSE
  if(flag.debug) { print('Debug is on');flush.console() }


  ## Load libraries, install if not already installed
  ## TODO: Move this into the package installation, then just load libraries normally
  if(!require(foreach)) {
    install.packages('foreach')
   library(foreach) 
  }
  if(!require(methods)) { # TODO: where is this used????
    install.packages('methods')
    library(methods)
  }
  
  # set parallel background if parallel flag is set
  if(parallel){
    ## TODO: Move this into the package installation, then just load libraries normally
    if(!require(doParallel)) {
      install.packages('doParallel')
      library(doParallel)
    }
    cl <- makeCluster(num.cores,outfile="")
    registerDoParallel(cl, cores = num.cores)
  }

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
    print(paste('Number of cv folds provided:',cv.folds));flush.console()
    fold.vec <- c(rep(seq(cv.folds), each=floor(dim(labs)[[1]]/cv.folds)),sample(seq(cv.folds), dim(labs)[[1]]-cv.folds*(floor(dim(labs)[[1]]/cv.folds)), replace=FALSE))
    fold.vec <- sample(fold.vec)
    labs$fold <- fold.vec
  } else if(is.vector(cv.folds) & length(cv.folds)==nrow(labs)) { 
    print('Using provided cv fold assignments');flush.console()
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
#    print(apply(labs.reduced,2,table));flush.console() 
#    print(dim(labs));flush.console()
#    print(dim(labs.reduced));flush.console()
    if(flag.debug) { print(table(labs.reduced));flush.console() }
 
    # get platypus result list
    platypus.result <- platypus(fn.views=fn.views, fn.labs=fn.labels.reduced, i=no.iterations, m=majority.threshold.percent,expanded.output=expanded.output,updating=updating,ignore.label=ignore.label)
 
    # TODO: THIS NEXT CODE BLOCK IS THE BROKEN BIT FOR SSC EXAMPLE
    # Error in na.fail.default(list(labels = c(`11031.p1` = 1L, `11033.p1` = 2L,  : 
    #   missing values in object
    ## Predict hold-out data subset with platypus result
    test.ids <- rownames(labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),,drop=F])
    predictions <- platypus.predict(platypus.result$final.views, platypus.result$weighting.threshold, test.ids,unique.labels,labs[,classcol.labs,drop=FALSE])

    if(flag.debug) { print('Held out data subset predictions made successfully');flush.console() }
    if(flag.debug) { print(head(predictions));flush.console() }
    
    
    ## Calculate final platypus performance for each cv-fold
    perf.all.agree <- calculate.performance(predictions[,c("final","category.all")],labs[,classcol.labs,drop=F],unique.labels)   
    if(flag.debug) { print('Calculated all.agree');flush.console() }
    perf.majority.agree <- calculate.performance(predictions[,c("final","category.majority")],labs[,classcol.labs,drop=F],unique.labels)
    if(flag.debug) { print('Calculated majority.agree');flush.console() }
    
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
        if(flag.debug) { print('Caluculating iteration performances');flush.console() }
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
    if(flag.debug) { print('leaving do.one.cvfold');flush.console() }
    return(return.list)
  } # end do.one.cvfold fxn


  
  if (parallel) {
    print("working parallel")
    # TODO: I've been removing functions but haven't been very careful about making sure they're removed from these lists
    cv.result.list <- foreach(k=seq(cv.folds), .export=c("platypus", "drop.features" ,"ElasticNet" ,"RandomForest"
      ,"load.parameterfile" ,"load.data" ,"load.data.ElasticNet" ,"load.data.RandomForest"
      ,"load.label.data" ,"get.unique.labels" ,"view.train" ,"view.train.ElasticNet" ,"view.train.RandomForest" ,"view.predict" ,"view.predict.ElasticNet"
      ,"view.predict.RandomForest" ,"platypus.predict" ,"update.accuracies" ,"update.accuracy" ,"update.accuracy.ElasticNet" ,"update.accuracy.RandomForest"
      ,"calculate.accuracy" ,"get.majority.weighting" ,"get.new.labels.majorityWeighted" ,"normalize.accuracies"
      ,"normalize.accuracy.log" ,"calculate.performance" ,"calculate.performance.view" ,"get.labelling.performance")
      , .verbose=TRUE
      , .packages=c("glmnet","randomForest")) %dopar% do.one.cvfold(k = k)
  } else {
    print("working non-parallel")
    cv.result.list <- lapply(seq(cv.folds), do.one.cvfold)
    #cv.result.list <- foreach(k=seq(cv.folds)) %do% do.one.cvfold(k = k) # Old but saving it for now
  }

  if(flag.debug) { print('out of cv folds loop');flush.console() }

  ## Collect accuracy over cv-iterations
  accuracy.cvfolds <- c()
  performance.cvfolds <- list()
  for(k in seq(cv.folds)){
    accuracy.cvfolds <- rbind(accuracy.cvfolds,cv.result.list[[k]]$accuracy.cvfolds)
    performance.cvfolds[[k]] <- cv.result.list[[k]]$performance.cvfolds
  }
  
  colnames(accuracy.cvfolds) <- c("cv.fold","accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree","accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                  ,paste0("accuracy.view.",seq(length(fn.views))),paste0("balanced.accuracy.view.",seq(length(fn.views))))

  if(flag.debug) { print('collected cv iteration performances');flush.console() }

  if(!is.na(output.folder)) {
    write.table(accuracy.cvfolds,file=file.path(output.folder,"perf_platypus.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
    save(performance.cvfolds,file=file.path(output.folder,"performance.cvfolds.Rdata") )
  }


  if(flag.debug) { print('Stored results to file');flush.console() }

  if(expanded.output){
    if(flag.debug) { print('Entering expanded.output loop');flush.console() }
    ## Collect accuracy over platypus-iterations
    accuracy.platypus.iterations <- c()
    performance.iterations <- list()
    for(k in seq(cv.folds)){
      accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations,cv.result.list[[k]]$accuracy.platypus.iterations)
      performance.iterations[[k]] <- cv.result.list[[k]]$performance.iterations
    }
    
      colnames(accuracy.platypus.iterations) <- c("cv.fold","iteration","weighting.threshold.upper","weighting.threshold.lower","weighting.threshold"
                                             ,paste0("weighting.view.",seq(length(fn.views))),paste0("weighting.norm.view.",seq(length(fn.views)))
                                             ,"no.ids.labelled","no.ids.unlabelled"
                                             ,"no.new.labelled", "no.ids.left.unlabelled"
                                             ,"accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree"
                                             ,"accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                             ,paste0("accuracy.view.",seq(length(fn.views))),paste0("balanced.accuracy.view.",seq(length(fn.views))))

    if(flag.debug) { print('Renamed the results');flush.console() }
    if(!is.na(output.folder)) { 
      write.table(accuracy.platypus.iterations, file= file.path(output.folder,"perf_platypus_expanded.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
      save(performance.iterations,file =file.path(output.folder,"performance.iterations.Rdata") )
    }
    if(flag.debug) { print('Renamed the results');flush.console() }
    
    labelling.matrix.cvlist <- list()
    labelling.matrices.views.cvlist <- list()
    for(k in seq(cv.folds)){
      labelling.matrix.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrix
      labelling.matrices.views.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrices.views
    }
    if(flag.debug) { print('Stored the labelling matrix');flush.console() }
   
    if(!is.na(output.folder)) { 
      save(labelling.matrix.cvlist,file =file.path(output.folder,"labelling.matrix.cvlist.Rdata") )
      save(labelling.matrices.views.cvlist,file =file.path(output.folder,"labelling.matrices.views.cvlist.Rdata") )
    }
    if(flag.debug) { print('Wrote to files the labelling matrix');flush.console() }
    
  } # end if(expanded.output)

  return( list(fold.accuracy=accuracy.cvfolds, labelling.matrix=labelling.matrix.cvlist, labelling.matrices.views=labelling.matrices.views.cvlist) )

} # end cv.platypus

