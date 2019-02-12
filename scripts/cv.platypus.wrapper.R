#!/usr/bin/env Rscript
## k-fold cross validation of platypus
##
## Created: July 2015, Verena Friedl
##         Updated: September 2015, Verena Friedl
##    Last Updated: Jan 2018, Kiley Graim


# pseudo-code:
#
#     get platypus parameters and view information
#     load labelled data    
#     divide labelled data in k subsets
#     for fold.i in 1:k
#       hold out subset fold.i of labelled data
#       train platypus.R on all-but-fold.i labelled (and unlabelled) data with the given views   # unlabelled data is optional, but making use of unlabelled data is one of the biggest advantages of platypus
#       predict fold.i subset on trained platypus views
#
#



################################################################################
###  Main Function  #########################################################
################################################################################


# call cv.platypus.R
# first pass the filepath for the labs file and the column name or number of the class (if not given, first column is default)
# pass the filepath of a parameter-file for each view
# [OPTIONS] <filename_labs> [-c <class_col>] for each view(<filename_viewfile>) 
# OPTIONS:
# -k <number of folds for cross validation, eg. 5, default=10>
# -i <number of iterations, eg. 100, default=100>
# -m <majority threshold in percent, eg. 75, default=100>
# -w flag for weighting the preditions by accuracy, default=FALSE
# -u flag for updating the accuracies of the single views in each iteration, default=FALSE
# -e flag for expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
# -b <class_name> flag for excluding cell lines that fall into class 'class_name' for the binary drug response definition, default='intermediate'
# -o <output folder>: folder to save output to, default=~/
# -p <num of cores>: give the number of cores to use and turn on parallelization

## Example: cv.platypus.R -k 5 -i 10 -m 100 -w -o platypus_output -u -c  TODO



cv.platypus <- function(argv) {
 
#/Users/kgraim/Documents/MVL/PLATYPUS/scripts 
# TODO: Don't hard code this
  #source(paste0(Sys.getenv("HOME"),'/Documents/MVL/scripts/platypus.R'))
  #source(paste0(Sys.getenv("HOME"),'/Documents/MVL/scripts/platypus.basicFunctions.R'))
  source('platypus.R')
  source('platypus.basicFunctions.R')
  require(foreach)
  require(methods)
  
  
  ## Read in the command line arguments
  # set number of folds for cross validation
  cv.folds <- 10
  if("-k" %in% argv){
    cv.folds <- as.numeric(as.character(argv[which(argv == "-k")+1]))
    argv <- argv[-c(which(argv == "-k"),which(argv == "-k")+1)]
  }
  
  # set maximal number of iterations for each platypus run, if given
  no.iterations <- NA
  if("-i" %in% argv){
    no.iterations <- as.numeric(as.character(argv[which(argv == "-i")+1]))
    argv <- argv[-c(which(argv == "-i"),which(argv == "-i")+1)]
  }
  
  # set majority threshold, if given
  majority.threshold.percent <- NA
  if("-m" %in% argv){
    majority.threshold.percent <- as.numeric(as.character(argv[which(argv == "-m")+1]))
    argv <- argv[-c(which(argv == "-m"),which(argv == "-m")+1)]
  }
  
  # set expanded output flag, if given
  expanded.output <- FALSE
  if("-e" %in% argv){
    expanded.output <- TRUE
    argv <- argv[-c(which(argv == "-e"))]
  }
  
  # set weighting flag, if given
  weighting = FALSE
  if("-w" %in% argv){
    weighting <- TRUE
    argv <- argv[-c(which(argv == "-w"))]
  }
  
  # set updating flag, if given
  updating = FALSE
  if("-u" %in% argv){
    updating <- TRUE
    argv <- argv[-c(which(argv == "-u"))]
  }
  
  # set class to ignore from training and predicting
  ignore.label <- "intermediate"
  if("-b" %in% argv){
    ignore.label <- as.character(argv[which(argv == "-b")+1])
    argv <- argv[-c(which(argv == "-b"),which(argv == "-b")+1)]
  }
  
  # set parallel background if parallel flag is set
  parallel <- FALSE
  num.cores <- 25
  if("-p" %in% argv){
    parallel <- TRUE
    num.cores <- as.numeric(as.character(argv[which(argv == "-p")+1])) 
    
    require(doParallel)
    
    cl <- makeCluster(num.cores,outfile="")
    registerDoParallel(cl, cores = num.cores)

    argv <- argv[-c(which(argv == "-p"),which(argv == "-p")+1)]
  }
  
  # set output folder, if given
  output.folder <- paste0(Sys.getenv("HOME"),'/')
  if("-o" %in% argv){
    output.folder <- as.character(argv[which(argv == "-o")+1])
    argv <- argv[-c(which(argv == "-o"),which(argv == "-o")+1)]
  }
  
  # retrieve filenames
  fn.labs <- argv[1]
  classcol.labs <- 1
  if ("-c" %in% argv) {
    col <- argv[which(argv == "-c")+1]
    if (is.na(as.numeric(as.character(col)))){
      classcol.labs <- col
    } else{
      classcol.labs <- as.numeric(as.character(col))
    }
    argv <- argv[-c(which(argv == "-c"),which(argv == "-c")+1)]
  }
  
  fn.views <- argv[2:length(argv)]

  ## Load the label data
  labs <- load.label.data(fn.labs,classcol.labs)
  
  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)
  
  ## Sort the labels in k subsets
  # Randomly create fold number for each label (balance folds to equal size)    TODO: shall we balance folds for equal class instances (based on the overall occurrence of each class in the data set)?
  # TODO: give the option to give user-defined folds as input
  fold.vec <- c(rep(1:cv.folds, each=floor(dim(labs)[[1]]/cv.folds)),sample(1:cv.folds, dim(labs)[[1]]-cv.folds*(floor(dim(labs)[[1]]/cv.folds)), replace=FALSE))
  fold.vec <- sample(fold.vec)
  labs$fold <- fold.vec
  
  ###################
  # function to do a cv run for one fold
  ###################
  do.one.cvfold <- function(k) {
    
    print(k)
    
    ## Define reduced label file
    # mark labels of the hold out data set with 'testing', so that they can also be excluded from the unlabelled data in the platypus.R
    labs.reduced <- labs
    labs.reduced[which(labs$fold == k),classcol.labs] <- 'testing'
    fn.labels.reduced <- tempfile()
    write.table(labs.reduced
                , file= fn.labels.reduced
                , sep="\t",row.names=T, col.names=T, quote=FALSE)
    
    ## Call platypus 
    # build call vector
    argv.for.platypus <- c()
    
    if(!(is.na(no.iterations))){
      argv.for.platypus <-c(argv.for.platypus,"-i",no.iterations)
    }
    
    if(!(is.na(majority.threshold.percent))){
      argv.for.platypus <-c(argv.for.platypus,"-m",majority.threshold.percent)
    }
    
    if(expanded.output){
      argv.for.platypus <-c(argv.for.platypus,"-e")
    }
    
    if(weighting){
      argv.for.platypus <-c(argv.for.platypus,"-w")
    }
    
    if(updating){
      argv.for.platypus <-c(argv.for.platypus,"-u")
    }
    
    argv.for.platypus <- c(argv.for.platypus,"-b",ignore.label)
    
    argv.for.platypus <- c(argv.for.platypus,fn.labels.reduced,"-c",classcol.labs)
    
    argv.for.platypus <- c(argv.for.platypus,fn.views)
    
    # get platypus result list
    platypus.result <- platypus(argv.for.platypus)
    
    ## Predict hold-out data subset with platypus result
    test.ids <- rownames(labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),,drop=F])
    if(weighting){
      predictions <- platypus.predict(platypus.result$final.views, platypus.result$weighting.threshold, test.ids,weighting,unique.labels)
    } else {
      predictions <- platypus.predict(platypus.result$final.views, platypus.result$majority.threshold, test.ids,weighting,unique.labels)
    }
    
    
    ## Calculate final platypus performance for each cv-fold
    perf.all.agree <- calculate.performance(predictions[,c("final","category.all")],labs[,classcol.labs,drop=F],unique.labels)   
    perf.majority.agree <- calculate.performance(predictions[,c("final","category.majority")],labs[,classcol.labs,drop=F],unique.labels)
    
    ## Calculate final performance values for each view
    perf.views <- list()
    accuracy.views <- c()
    balanced.accuracy.views <- c()
    for(view.i in 1:length(platypus.result$final.views)){
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
      for(i in 1:length(iteration.information)){
        ## Predict hold-out data subset with platypus view from iteration
        view.list <- iteration.information[[i]]$view.list
        test.ids <- rownames(labs[which(labs$fold == k & labs[,classcol.labs] != ignore.label),,drop=F])
        
        if(weighting){
          predictions <- platypus.predict(view.list, platypus.result$weighting.threshold, test.ids, weighting,unique.labels)
        } else {
          predictions <- platypus.predict(view.list, platypus.result$majority.threshold, test.ids, weighting,unique.labels)
        }
        
        
        ## Calculate performance of the iteration platypus and views
        perf.all.agree <- calculate.performance(predictions[,c("final","category.all")],labs[,classcol.labs,drop=F],unique.labels)   
        perf.majority.agree <- calculate.performance(predictions[,c("final","category.majority")],labs[,classcol.labs,drop=F],unique.labels)
        
        perf.views <- list()
        accuracy.views <- c()
        balanced.accuracy.views <- c()
        weighting.views <- c()
        weighting.norm.views <- c()
        for(view.i in 1:length(view.list)){
          perf.view <- calculate.performance.view(predictions[,view.i,drop=F],labs[,classcol.labs,drop=F],unique.labels)
          perf.views[[view.i]] <- perf.view
          accuracy.views <- c(accuracy.views,perf.view$accuracy)
          balanced.accuracy.views <- c(balanced.accuracy.views,perf.view$balanced.accuracy)
          
          weighting.views <- c(weighting.views,view.list[[view.i]]$acc)
          weighting.norm.views <- c(weighting.norm.views,view.list[[view.i]]$acc.norm)
        }
        
        if(weighting){
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
        } else {
          accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations
                                           ,c(k,i,iteration.information[[i]]$iterations.woChange
                                              ,iteration.information[[i]]$majority,iteration.information[[i]]$majority.missingData,iteration.information[[i]]$majority.threshold
                                              ,iteration.information[[i]]$no.ids.labelled,iteration.information[[i]]$no.ids.unlabelled
                                              ,iteration.information[[i]]$no.new.labelled, iteration.information[[i]]$no.ids.left.unlabelled
                                              ,perf.all.agree$accuracy,perf.all.agree$balanced.accuracy,perf.all.agree$coverage
                                              ,perf.majority.agree$accuracy,perf.majority.agree$balanced.accuracy,perf.majority.agree$coverage
                                              ,accuracy.views,balanced.accuracy.views
                                           )
          )
        }
        
        
        performance.iterations[[i]] <- list(perf.all.agree=perf.all.agree,perf.majority.agree=perf.majority.agree,perf.views=perf.views) 
      }
    }
    return.list <- list(accuracy.cvfolds=accuracy.cvfolds,performance.cvfolds=performance.cvfolds)
    if(expanded.output){
      return.list <- list(accuracy.cvfolds=accuracy.cvfolds,performance.cvfolds=performance.cvfolds,accuracy.platypus.iterations=accuracy.platypus.iterations,performance.iterations=performance.iterations
                          ,labelling.matrix=platypus.result$labelling.matrix,labelling.matrices.views=platypus.result$labelling.matrices.views)
    }
    return(return.list)
  }


  
  if (parallel) {
    print("working parallel")
    cv.result.list <- foreach(k=1:cv.folds
                              , .export=c("platypus"
                                          ,"drop.features"
                                          ,"ElasticNet"
                                          ,"RandomForest"
                                          ,"setAlpha"
                                          ,"setMeasure"
                                          ,"setMtry"
                                          ,"setNtree"
                                          ,"setDrop"
                                          ,"setDropTo"
                                          ,"setAcc"
                                          ,"setAccNorm"
                                          ,"load.parameterfile"
                                          ,"load.data"
                                          ,"load.data.ElasticNet"
                                          ,"load.data.RandomForest"
                                          ,"load.label.data"
                                          ,"get.unique.labels"
                                          ,"addX"
                                          ,"view.train"
                                          ,"view.train.ElasticNet"
                                          ,"view.train.RandomForest"
                                          ,"view.predict"
                                          ,"view.predict.ElasticNet"
                                          ,"view.predict.RandomForest"
                                          ,"platypus.predict"
                                          ,"update.accuracies"
                                          ,"update.accuracy"
                                          ,"update.accuracy.ElasticNet"
                                          ,"update.accuracy.RandomForest"
                                          ,"calculate.accuracy"
                                          ,"get.majority.counting"
                                          ,"get.majority.weighting"
                                          ,"get.new.labels.majorityCount"
                                          ,"get.new.labels.majorityWeighted"
                                          ,"normalize.accuracies"
                                          ,"normalize.accuracy.linear"
                                          ,"normalize.accuracy.log"
                                          ,"calculate.performance"
                                          ,"calculate.performance.view"
                                          ,"get.labelling.performance"
                                          )
                              , .verbose=TRUE
                              , .packages=c("glmnet","randomForest")) %dopar%
      do.one.cvfold(k = k)
  } else {
    print("working non-parallel")
    cv.result.list <- foreach(k=1:cv.folds) %do%
      do.one.cvfold(k = k)
  }

  ## Collect accuracy over cv-iterations
  accuracy.cvfolds <- c()
  performance.cvfolds <- list()
  for(k in 1:cv.folds){
    accuracy.cvfolds <- rbind(accuracy.cvfolds,cv.result.list[[k]]$accuracy.cvfolds)
    performance.cvfolds[[k]] <- cv.result.list[[k]]$performance.cvfolds
  }
  
  colnames(accuracy.cvfolds) <- c("cv.fold","accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree","accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                  ,paste0("accuracy.view.",1:length(fn.views)),paste0("balanced.accuracy.view.",1:length(fn.views)))

  # TODO: better define output and output file/mode
  write.table(accuracy.cvfolds
             , file= paste0(output.folder,"/perf_platypus.tab")
             , sep="\t",row.names=F, col.names=T, quote=FALSE)

  save(performance.cvfolds,file =paste0(output.folder,"performance.cvfolds.Rdata") )

  if(expanded.output){
    ## Collect accuracy over platypus-iterations
    accuracy.platypus.iterations <- c()
    performance.iterations <- list()
    for(k in 1:cv.folds){
      accuracy.platypus.iterations <- rbind(accuracy.platypus.iterations,cv.result.list[[k]]$accuracy.platypus.iterations)
      performance.iterations[[k]] <- cv.result.list[[k]]$performance.iterations
    }
    
    if(weighting){
      colnames(accuracy.platypus.iterations) <- c("cv.fold","iteration","weighting.threshold.upper","weighting.threshold.lower","weighting.threshold"
                                             ,paste0("weighting.view.",1:length(fn.views)),paste0("weighting.norm.view.",1:length(fn.views))
                                             ,"no.ids.labelled","no.ids.unlabelled"
                                             ,"no.new.labelled", "no.ids.left.unlabelled"
                                             ,"accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree"
                                             ,"accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                             ,paste0("accuracy.view.",1:length(fn.views)),paste0("balanced.accuracy.view.",1:length(fn.views)))
    } else{
      colnames(accuracy.platypus.iterations) <- c("cv.fold","iteration","iterations.woChange","majority","majority.missingData","majority.threshold"
                                             ,"no.ids.labelled","no.ids.unlabelled"
                                             ,"no.new.labelled", "no.ids.left.unlabelled"
                                             ,"accuracy.all.agree","balanced.accuracy.all.agree","coverage.all.agree"
                                             ,"accuracy.majority.agree","balanced.accuracy.majority.agree","coverage.majority.agree"
                                             ,paste0("accuracy.view.",1:length(fn.views)),paste0("balanced.accuracy.view.",1:length(fn.views)))
    }
    
    write.table(accuracy.platypus.iterations
                , file= paste0(output.folder,"/perf_platypus_expanded.tab")
                , sep="\t",row.names=F, col.names=T, quote=FALSE)
    
    save(performance.iterations,file =paste0(output.folder,"performance.iterations.Rdata") )
    
    
    labelling.matrix.cvlist <- list()
    labelling.matrices.views.cvlist <- list()
    for(k in 1:cv.folds){
      labelling.matrix.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrix
      labelling.matrices.views.cvlist[[k]] <- cv.result.list[[k]]$labelling.matrices.views
    }
    
    save(labelling.matrix.cvlist,file =paste0(output.folder,"labelling.matrix.cvlist.Rdata") )
    save(labelling.matrices.views.cvlist,file =paste0(output.folder,"labelling.matrices.views.cvlist.Rdata") )
    
  }

}

cv.platypus(commandArgs(TRUE))
