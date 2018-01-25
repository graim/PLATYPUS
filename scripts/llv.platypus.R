## Functions for the MVL framework - label learning validation
##
## Created: July 2015, Verena Friedl
##    Updated Sept 2015, Verena Friedl
##    Updated: April 2016, Kiley Graim
##    Updated: Feb 2017, Kiley Graim
## Last Updated: Jan 2018, Kiley Graim

# pseudo-code:
#
#     for (v in views)
#       load v
#     load labels
#
#     divide labelled data in folds
#     for(f in folds):
#       unlabelled.data <- fold f
#       delete f from labels
#       while( not converged on labels OR no changes anymore )
#         for( v in views )
#           train v on all sets but f
#           try to predict fold f
#         add newly found labels to training data
#         calculate performance for added labels
#       calculate overall performance for a fold
#     average performance over folds
#       



################################################################################
###  Main Function  #########################################################
################################################################################

# call llv.platypus.R
# first pass the filepath for the labs file and the column name or number of the class (if not given, first column is default)
# pass the filepath of a parameter-file for each view
# [OPTIONS] <filename_labs> [-c <class_col>] for each view(<filename_viewfile>) 
# OPTIONS:
# -k <number of folds for label learning validation (similar to cross validation folds), eg. 5, default=10>
# -i <maximal number of iterations for each platypus run, eg. 100, default=100>
# -m <majority threshold in percent, eg. 75, default=100>
# -w flag for weighting the preditions by accuracy, default=FALSE
# -u flag for updating the accuracies of the single views in each iteration, default=FALSE
# -e flag for expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
# -b <class_name> flag for excluding cell lines that fall into class 'class_name' for the binary drug response definition, default='intermediate'
# -o <output folder>: folder to save output to, default=~/
# -p <num of cores>: give the number of cores to use and turn on parallelization




###################
# main function
###################
llv.platypus <- function(fn.views,fn.labs,llv.folds=10,no.iterations=100,majority.threshold.percent=100,expanded.output=TRUE,weighting=FALSE,updating=FALSE,ignore.label='intermediate',parallel=FALSE,num.cores=25,classcol.labs=1,output.folder=NA) {
  
#  source(paste0(Sys.getenv("HOME"),'/MVL/scripts/platypus.R'))
#  source(paste0(Sys.getenv("HOME"),'/MVL/scripts/platypus.basicFunctions.R'))
  require(foreach)
  require(methods)
  
  if(parallel) {
    require(doParallel)
    cl <- makeCluster(num.cores,outfile="")
    registerDoParallel(cl, cores = num.cores)
  }
 
  ## Load the label data
  labs <- load.label.data(fn.labs,classcol.labs)
  
  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)
  
  ## TODO: should the real unlabelled data be excluded from the label learning validation? If so, it has to be deleted from the inputfiles and new matrix-files and parameter-files have to be saved 
  ###(maybe work around by marking the real unlabelled data as 'testing' like in the cv, so that they are not used by the platypus-function)
  
  ## Sort the labels in k subsets
  # Randomly create fold number for each label (balance folds to equal size)    TODO: shall we balance folds for equal class instances (based on the overall occurrence of each class in the data set)?
  # TODO: give the option to give user-defined folds as input
  fold.vec <- c(rep(1:llv.folds, each=floor(dim(labs)[[1]]/llv.folds)),sample(1:llv.folds, dim(labs)[[1]]-llv.folds*(floor(dim(labs)[[1]]/llv.folds)), replace=FALSE))
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
    platypus.result <- platypus(fn.views=fn.views, fn.labs=fn.labels.reduced, i=no.iterations, m=majority.threshold.percent,expanded.output=expanded.output,updating=updating)
    
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
      if(weighting){
        weighting.threshold.upper <- c()
        weighting.threshold.lower <- c()
        weighting.threshold <- c()
        no.ids.labelled <- c()
        no.ids.unlabelled <- c()
        no.new.labelled <- c()
        no.ids.left.unlabelled <- c()

        for(i in 1:length(iteration.information)){
          weighting.threshold.upper <- c(weighting.threshold.upper,iteration.information[[i]]$weighting.threshold.upper)
          weighting.threshold.lower <- c(weighting.threshold.lower,iteration.information[[i]]$weighting.threshold.lower)
          weighting.threshold <- c(weighting.threshold,iteration.information[[i]]$weighting.threshold)
          no.ids.labelled <- c(no.ids.labelled,iteration.information[[i]]$no.ids.labelled)
          no.ids.unlabelled <- c(no.ids.unlabelled,iteration.information[[i]]$no.ids.unlabelled)
          no.new.labelled <- c(no.new.labelled,iteration.information[[i]]$no.new.labelled)
          no.ids.left.unlabelled <- c(no.ids.left.unlabelled,iteration.information[[i]]$no.ids.left.unlabelled)
        }
        
        perf.iterations <- cbind(llv.fold,weighting.threshold.upper,weighting.threshold.lower,weighting.threshold,no.ids.labelled,no.ids.unlabelled,no.new.labelled,no.ids.left.unlabelled,perf.iterations)

      } else {
        iterations.woChange <- c()
        majority <- c()
        majority.missingData <- c()
        majority.threshold <- c()
        no.ids.labelled <- c()
        no.ids.unlabelled <- c()
        no.new.labelled <- c()
        no.ids.left.unlabelled <- c()

        for(i in 1:length(iteration.information)){
          iterations.woChange <- c(iterations.woChange,iteration.information[[i]]$iterations.woChange)
          majority <- c(majority,iteration.information[[i]]$majority)
          majority.missingData <- c(majority.missingData,iteration.information[[i]]$majority.missingData)
          majority.threshold <- c(majority.threshold,iteration.information[[i]]$majority.threshold)
          no.ids.labelled <- c(no.ids.labelled,iteration.information[[i]]$no.ids.labelled)
          no.ids.unlabelled <- c(no.ids.unlabelled,iteration.information[[i]]$no.ids.unlabelled)
          no.new.labelled <- c(no.new.labelled,iteration.information[[i]]$no.new.labelled)
          no.ids.left.unlabelled <- c(no.ids.left.unlabelled,iteration.information[[i]]$no.ids.left.unlabelled)
        }
        
        perf.iterations <- cbind(llv.fold,iterations.woChange,majority,majority.missingData,majority.threshold,no.ids.labelled,no.ids.unlabelled,no.new.labelled,no.ids.left.unlabelled,perf.iterations)

      } # end if(weighting)
      
    } # end if(expanded.output)

    return.list <- list(perf.llvfold=perf.llvfold)
    if(expanded.output){
      return.list <- list(perf.llvfold=perf.llvfold,perf.iterations=perf.iterations,labelling.matrix=platypus.result$labelling.matrix,labelling.matrices.views=platypus.result$labelling.matrices.views)
    }
    return(return.list) 
  } # end function do.one.llvfold
  
  if (parallel) {
    llv.result.list <- foreach(k=1:llv.folds, .export=c("platypus","drop.features","ElasticNet","RandomForest","setAlpha","setMeasure","setMtry",
        "setNtree","setDrop","setDropTo","setAcc","setAccNorm","load.parameterfile","load.data","load.data.ElasticNet","load.data.RandomForest",
        "load.label.data","get.unique.labels","addX","view.train","view.train.ElasticNet","view.train.RandomForest","view.predict","view.predict.ElasticNet",
        "view.predict.RandomForest","platypus.predict","update.accuracies","update.accuracy","update.accuracy.ElasticNet","update.accuracy.RandomForest",
        "calculate.accuracy","get.majority.counting","get.majority.weighting","get.new.labels.majorityCount","get.new.labels.majorityWeighted","normalize.accuracies",
        "normalize.accuracy.linear","normalize.accuracy.log","calculate.performance","calculate.performance.view","get.labelling.performance"),
        .verbose=TRUE, .packages=c("glmnet","randomForest")) %dopar%
      do.one.llvfold(k = k)
  } else {
    llv.result.list <- foreach(k=1:llv.folds) %do%
      do.one.llvfold(k = k)
  }  
  
  ## Get output together
  # Collect accuracy over llv-iterations
  accuracy.llvfolds <- c()
  for(k in 1:llv.folds){
    accuracy.llvfolds <- rbind(accuracy.llvfolds,llv.result.list[[k]]$perf.llvfold)
  }

  if(expanded.output){
    ## Collect accuracy over platypus-iterations and labelling information
    accuracy.platypus.iterations <- c()
    for(k in 1:llv.folds){
      accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations,llv.result.list[[k]]$perf.iterations)
    }
    
    if(!is.na(output.folder)){
      write.table(accuracy.llvfolds, file= paste0(output.folder,"/perf_llv.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
      write.table(accuracy.platypus.iterations, file= paste0(output.folder,"/perf_llv_expanded.tab"), sep="\t",row.names=F, col.names=T, quote=FALSE)
    }
    
    labelling.matrix.llvlist <- list()
    labelling.matrices.views.llvlist <- list()
    for(k in 1:llv.folds){
      labelling.matrix.llvlist[[k]] <- llv.result.list[[k]]$labelling.matrix
      labelling.matrices.views.llvlist[[k]] <- llv.result.list[[k]]$labelling.matrices.views
    }
    
    save(labelling.matrix.llvlist,file =paste0(output.folder,"labelling.matrix.llvlist.Rdata") )
    save(labelling.matrices.views.llvlist,file =paste0(output.folder,"labelling.matrices.views.llvlist.Rdata") )
  }

  # TODO: this should return an LLV object
  return( list(fold.accuracy=accuracy.llvfolds, labelling.matrix=labelling.matrix.llvlist,labelling.matrices.views=labelling.matrices.views.llvlist) )
}
