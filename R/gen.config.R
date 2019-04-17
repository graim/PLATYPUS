## Fxn to write an elastic net config file
write.config.en <- function(x, v, task, fn.config='config_TEST.txt') {
  write(paste('data.fn',v,sep='\t'), file=fn.config, append=FALSE)
  write(paste('type','en', sep='\t'), file=fn.config, append=TRUE)
  write(paste('alpha',x[['alpha']], sep='\t'), file=fn.config, append=TRUE)
  #write(paste('alpha',x$alpha, sep='\t'), file=fn.config, append=TRUE) # TODO 
  write(paste('measure','auc', sep='\t'), file=fn.config, append=TRUE)
  write(paste('drop','FALSE', sep='\t'), file=fn.config, append=TRUE) # TODO: this tells if data should be dropped to x features. 
  write(paste('acc',x$accuracy, sep='\t'), file=fn.config, append=TRUE)
  write(paste('taskname',task, sep='\t'), file=fn.config, append=TRUE) # TODO: added for multitask learning but not yet used. Need to add to svm/rf versions too
}

## TODO: We do not have SVM support in PLATYPUS yet!!! 
## Fxn to write a SVM config file
write.config.svm <- function(x, fn.config='config_TEST.txt') {
  write(paste('data.fn',v,sep='\t'), file=fn.config, append=FALSE)
  write(paste('type','svm', sep='\t'), file=fn.config, append=TRUE)
  write(paste('C',x$C, sep='\t'), file=fn.config, append=TRUE)
  write(paste('acc',x$Accuracy, sep='\t'), file=fn.config, append=TRUE)
  write(paste('accSD',x$SDAccuracy, sep='\t'), file=fn.config, append=TRUE)
}

## Fxn to write a random forest config file
## TODO: if no fn.config is provided, should it return a view object? Could then have platypus read in either config files or view objects - it's an easy check 
write.config.rf <- function(x,v,fn.config='config_TEST.txt') {
  write(paste('data.fn',v,sep='\t'), file=fn.config, append=FALSE)
  write(paste('type','rf', sep='\t'), file=fn.config, append=TRUE)
  write(paste('mtry',x$mtry, sep='\t'), file=fn.config, append=TRUE)
  write(paste('ntree',x$ntree, sep='\t'), file=fn.config, append=TRUE)
  write(paste('acc',x$accuracy, sep='\t'), file=fn.config, append=TRUE)
}

## TODO: Do we really need all these input options???

#' Generate configuration files for platypus
#'
#' @param view.names List of files containing view feature data
#' @param fn.tasks File containing all task labels, one column per task
#' @param config.loc Where the config files should be stored
#' @param model.type Type of classifier to use (select from en, rf, svm)
#' @param delim Delimiter for the task file
#' @param delim.v Delimiter for the view data file
#' @param n.iters Number of iterations to run
#' @param ignore.label Label to ignore in the task file (default 'intermediate')
#' @param nfolds Number of folds to use
#' @param mtry For random forest models, what mtry value(s) to use
#' @param ntree For random forest models, number of trees to build
#'
#' @return List of config filenames, for use in platypus
#'
#' @export
gen.config <- function(view.names, fn.tasks, config.loc='config', model.type=c('en','rf','svm'), delim=',', delim.v='\t', n.iters=10, ignore.label='intermediate', nfolds=10, mtry=NA, ntree=c(500,1000,1500,2000)) {

  ## For each task - load the task
  ##   For each view - load the view
  ##     Find optimal parameters for view/task pair
  ##     generate config file
  ##     add config filename to return list
  ## return list of config filenames

  ## Make sure model type is in our current list of options
  model.type=match.arg(model.type)  

  ## Set up options
  alpha.seq <- seq(0.1, 0.9, 0.1)

  ## Load file with tasks (one task per column, NA/blank values for missing labels)
  tasks <- read.table(fn.tasks,sep=delim, header=TRUE, row.names=1,check.names=FALSE,stringsAsFactors=FALSE)

  ## Store list of config filenames for returning
  fns.config <- list()

  ## Main loop
  for( v in view.names ) {
    print(paste('View',v))

    ## Load the view data
    X <- read.table(v,sep=delim.v, header=TRUE, row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
    X <- X[complete.cases(X),]

    for( task in colnames(tasks) ) {
      print(paste('Task',task))

      ## Use the current task labels - for multiview learning there's just the 1 task
      y <- tasks[,colnames(tasks)==task]
      names(y) <- rownames(tasks)

      ## Set up filename for this config
      ## TODO: this takes the view name from the filename. Might not want to do that
      fn.config <- switch(model.type,
        en = file.path(config.loc, paste0('config_en_',task,'_',unlist(strsplit(basename(v),'.', fixed=TRUE))[1],'.txt')),
        rf = file.path(config.loc, paste0('config_rf_',task,'_',unlist(strsplit(basename(v),'.', fixed=TRUE))[1],'.txt')),
        svm= file.path(config.loc, paste0('config_svm_',task,'_',unlist(strsplit(basename(v),'.', fixed=TRUE))[1],'.txt'))
      )
      print( paste('Generating config ',fn.config) )
      fns.config <- c(fns.config, fn.config)

      ## Parameter sweep based on task type
      if(model.type=='en') {
        res <- single.elasticNet.predictor( X, y, alpha = alpha.seq, iterations = n.iters, nfolds=nfolds)
        write.config.en(res, v, task, fn.config=fn.config)

      } else if(model.type=='rf') {
        if(is.na(mtry)) { mtry <- seq(ceiling(sqrt(ncol(X)))) } # If mtry not provided, use this default
        res <- single.randomForest.predictor(X, y, mtry=seq(ceiling(sqrt(ncol(X)))) )
        write.config.rf(res, v, fn.config=fn.config)

      } else if(model.type=='svm') {
        res <- single.svm.predictor( X, y )
        write.config.svm(res, fn.config=fn.config)
      }
      rm(res) 
      
    } # end for tasks
  } # end for views

  print('Finished, success!')
  return(unlist(fns.config)) # TODO: added unlist() because was returning a nested list. Haven't rerun to make sure is correct
}
