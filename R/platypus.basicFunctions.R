## Basic Functions for the MVL framework

## Subset to k most variable features
## TODO: Do we really need to provide this??
drop.features <- function(dat, k) {
  if(ncol(dat) > k ) {
    vars <- apply(dat, 2, stats::var)
    vars <- sort(vars, decreasing=TRUE)
    thresh <- vars[k]
    dat <- dat[,names(which(vars>=thresh))]
    write( paste('Shrinking to', k, 'most variable features...\n', ncol(dat)), stdout() )
    if(ncol(dat) > k) {
      write( paste('NOTE: Many variables have the same variance, so', ncol(dat), 'variables were kept\n'), stdout() )
    }
  } else {
    write(paste('WARNING: Data matrix already has', k, 'features\n'),stdout() )
  }
  return(dat)
}

## TODO Make a platypus object class, as well as cv.platypus and llv.platypus
# fn.labs, view.list, ignore.label='intermediate', i=100, m=100, u=FALSE, e=FALSE,updating=FALSE,expanded.output=FALSE
# Careful- we already define platypus() as a function for mvl
#Platypus <- function(labs, views, iters=NA) {
#  if (is.na(iters)) iters <- 10
#  if (!is.character(labs)) stop("labs must be character")
#  if (!is.list(views)) stop("views must be a list")
#  if (!all(is.character(sapply(views,class)))) stop('views must be a list of strings')
#  if(!is.numeric(iters)) stop ('iters must be numeric')
#  structure(list(views=views,labs=labs,iters=iters), class = "foo")
#}
#platypus('moo', list('oink','baa'))


## create a class for each view-type
ElasticNet <- function(param.file="",data.matrix=c(),data.fn="",alpha=0.9,measure="auc",drop=TRUE,drop.to=5000,model=c(),acc=0.5, acc.norm=0, family='binomial') {
  me <- list(
    param.file = param.file,
    data.matrix = data.matrix,
    data.fn = data.fn,
    alpha = alpha,
    measure = measure,
    family = family,
    drop = drop,
    drop.to = drop.to,
    model = model,
    acc = acc,
    acc.norm = acc.norm
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"ElasticNet")
  return(me)
}

RandomForest <- function(param.file="",data.matrix=c(),data.fn="",mtry="sqrt",ntree=500,drop=TRUE,drop.to=5000,model=c(),acc=0.5,acc.norm=0) {
  me <- list(
    param.file = param.file,
    data.matrix = data.matrix,
    data.fn = data.fn,   
    mtry = mtry,  
    ntree = ntree,
    drop = drop,
    drop.to = drop.to,
    model = model,
    acc = acc,
    acc.norm = acc.norm
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"RandomForest")
  return(me)
}

## TODO: I created this but didn't actually update it or change content from RandomForest object
SupportVectorMachine <- function(param.file="",data.matrix=c(),data.fn="",model=c(),acc=0.5,accSD=0) {
  me <- list(
    param.file = param.file,
    data.matrix = data.matrix,
    data.fn = data.fn,
#    type = type, 
#    C = C,
    model = model,
    acc = acc,
    accSD = accSD
  )

  ## Set the name for the class
  class(me) <- append(class(me),"SupportVectorMachine")
  return(me)
}

## Read in the parameter file for a view and return a view object
## TODO: Add in option for view name
## TODO: add in SVM models
#' Create a view from a configuration file
#'
#' @param filename File containing config information
#' @param delim Optional delimeter for the file containing view data (not the config file!)
#' @return View object
#' @keywords platypus
#' @export
load.parameterfile <- function(filename, delim='\t'){
 
  # Load parameters file 
  param.table <- utils::read.table(filename, sep='\t',header=FALSE, row.names=1)
  
  # Check the type and create a view object, set type-specific parameters
  type <- param.table["type",1]
  if(type == "en"){
    view.object <- ElasticNet(param.file=filename,data.fn=toString(param.table["data.fn",1]))
    if("alpha" %in% rownames(param.table)){ view.object$alpha <- as.numeric(as.character(param.table["alpha",1])) }
    if("measure" %in% rownames(param.table)){ view.object$measure <- as.character(param.table["measure",1]) }
  } else if(type == 'rf'){
    view.object <- RandomForest(param.file=filename,data.fn=toString(param.table["data.fn",1]))
    if("mtry" %in% rownames(param.table)){ view.object$mtry <- as.numeric(as.character(param.table["mtry",1])) }
    if("ntree" %in% rownames(param.table)){ view.object$ntree <- as.numeric(as.character(param.table["ntree",1])) }
  } else if(type == 'svm'){
    stop("Sorry, SVMs are not working yet!")
    ## TODO: implement svm types :)
  } else { 
    stop(paste("Could not find valid type specification in parameter file:",filename,"; Valid types are: rf, en, svm."))
  }

  # Parameters common to all view types
  if("drop" %in% rownames(param.table)){ view.object$drop <- as.logical(param.table["drop",1]) } # TODO: Remove this option
  if("drop.to" %in% rownames(param.table)){ view.object$drop.to <- as.numeric(as.character(param.table["drop.to",1])) } # TODO: Remove this option
  if("acc" %in% rownames(param.table)){ view.object$acc <- as.numeric(as.character(param.table["acc",1])) }
  if("data.fn" %in% rownames(param.table)){ view.object$data.matrix <- data.matrix( utils::read.table(view.object$data.fn, sep=delim, header=TRUE, row.names=1, check.names=FALSE) ) } # TODO: Do we really want to force everything to be numeric???

  return(view.object)
}


## load the feature matrix from file to a matrix object
## TODO: Do we really need these fxns??
load.data <- function(view.object){
  UseMethod("load.data",view.object)
}
load.data.ElasticNet <- function(view.object,delim='\t'){
  # read the matrix from the file
  mat <- data.matrix( utils::read.table(view.object$data.fn, sep=delim, header=TRUE, row.names=1, check.names=FALSE) )
  print(paste('data loaded for Elastic Net view',view.object$data.fn))
  
  # drop features
  if(view.object$drop){
    mat <- drop.features(mat, view.object$drop.to)
  }
  
  # set the matrix
  view.object$data.matrix <- mat
  
  return(view.object)
}
load.data.RandomForest <- function(view.object,delim='\t'){
  # read the matrix from the file
  mat <- data.matrix( utils::read.table(view.object$data.fn, sep=delim,header=TRUE, row.names=1, check.names=FALSE) )
  print(paste('data loaded for Random Forest view',view.object$data.fn))
  
  # drop features
  if(view.object$drop){
    mat <- drop.features(mat, view.object$drop.to)
  }
  
  # set the matrix
  view.object$data.matrix <- mat
  
  # set mtry parameter to sqrt(#features) if default
  if(view.object$mtry == "sqrt"){
    view.object$mtry <- floor(sqrt(ncol(mat)))
    #view.object <- setMtry(view.object,floor(sqrt(ncol(mat))))
  }
  return(view.object)
}
load.data.SupportVectorMachine <- function(view.object,delim='\t'){
  # read the matrix from the file
  mat <- data.matrix( utils::read.table(view.object$data.fn, sep=delim,header=TRUE, row.names=1, check.names=FALSE) )
  print(paste('data loaded for SVM view',view.object$data.fn))

  # drop features
  if(view.object$drop){
    mat <- drop.features(mat, view.object$drop.to)
  }

  # set the matrix
  view.object$data.matrix <- mat

  return(view.object)
}

## load the label file
## TODO: I updated this to return the full labels matrix, but each fxn will now need to iterate over each column
load.label.data.old <- function(fn.labs,delim='\t'){
  # read file
  labs <- utils::read.table(fn.labs, sep=delim,header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE) 
  # take out 'NA' values - only in cases where ALL labels are NA for a given sample
  labs <- labs[!apply(labs,1,function(x){any(is.na(x))}),,drop=FALSE]
  return(labs)
}
## load the label file
## TODO: Keeping for now so I can run old platypus - retire it once we've added mtl
load.label.data <- function(fn.labs,classcol.labs, delim='\t'){
  # read file
  labs <- utils::read.table(fn.labs, sep=delim,header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE) 
  # take out 'NA' values
  labs <- labs[which(!(is.na(labs[,classcol.labs]))),,drop=FALSE]
  return(labs)
}

## retrieve the two labels for training/predicting
get.unique.labels <- function(label.vec,ignore.label){

  #print(paste('ignoring',ignore.label))
  ## Get the the two labels
  unique.labels <- unique(label.vec)
  
  # exclude a label if given by ignore.label, e.g. 'intermediate'
  if(ignore.label %in% unique.labels){
    unique.labels <- unique.labels[-c(which(unique.labels == ignore.label))]
  }
  # exclude label 'testing' - this is introduced by cv.platypus() to mark the test set
  if('testing' %in% unique.labels){
    unique.labels <- unique.labels[-c(which(unique.labels == 'testing'))]
  }

  unique.labels <- as.vector(sort(unique.labels))
  #print(unique.labels)
  
  # MVL classification for 2 labels
  if(length(unique.labels) != 2){
    print(paste0("ATTENTION: ",length(unique.labels)," labels for a binary task. Results might not be correct!"))
  }
  return(unique.labels)
}

## Train one view, given a specific model type
view.train <- function( labels, view.object ) {
  UseMethod("view.train",view.object)
}
view.train.ElasticNet <- function(labels, view.object, nfolds=10 ){
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(rownames(labels),rownames(view.object$data.matrix))

#    library(glmnet)
  # Reduce number of folds if fewer than 10 samples per fold
  if(length(ids)/nfolds < 10) { nfolds <- floor(length(ids)/10); print(paste('Using',nfolds,'because not enough samples')) }
  view.object$model <- cv.glmnet( view.object$data.matrix[ids,]
                                  ,as.factor(labels[ids,1]), family=view.object$family
                                  ,type.measure=view.object$measure
                                  ,alpha=view.object$alpha 
                                  ,nfolds=nfolds
                                  ,keep=TRUE)
  return(view.object)
}
view.train.RandomForest <- function( labels, view.object  ){
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(rownames(labels),rownames(view.object$data.matrix))

  # Train model
  #  library(randomForest)
  view.object$model <- randomForest( view.object$data.matrix[ids,], 
                                     as.factor(labels[ids,1]), 
                                     mtry=view.object$mtry, 
                                     ntree=view.object$ntree )
  return(view.object)
}
view.train.SupportVectorMachine <- function( labels, view.object  ){
  ## TODO: UNTESTED
  # Intersect IDs in labels and in the feature data
  ids <- intersect(rownames(labels),rownames(view.object$data.matrix))

#  print( summary( labels[ids,1] ) ) # TODO

  # Train model
  #if(!require(e1071)) {
  #  install.packages('e1071')
    #library(e1071)
  #}
  view.object$model <- e1071::svm(labels[ids,1] ~ ., data=view.object$data.matrix[ids,], 
                                       kernel=view.object$kernel,  
                                       cost=view.object$cost, 
                                       gamma=view.object$gamma)

  return(view.object)
}



## Take a trained view and new labels, return predictions
view.predict <- function(ids.unlabelled, view.object) {
  UseMethod("view.predict",view.object)
}
view.predict.ElasticNet <- function(ids.unlabelled, view.object) {
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(ids.unlabelled,rownames(view.object$data.matrix))
  if(length(ids) > 0){
    return( glmnet::predict.cv.glmnet(view.object$model, view.object$data.matrix[ids,,drop=FALSE], type='class', s='lambda.min') )  # TODO: use lambdaMin or default s="lambda.1se" ?  # TODO: changed to predict.cv.glmnet but not tested, from predict()
  }
  
  return(as.character(c()))
}
view.predict.RandomForest <- function(ids.unlabelled, view.object) {
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(ids.unlabelled,rownames(view.object$data.matrix))
  if(length(ids) > 0){
    return( as.character( predict(view.object$model, view.object$data.matrix[ids,,drop=FALSE], type='response')  ) )
    #return( as.character( randomForest::predict.randomForest(view.object$model, view.object$data.matrix[ids,,drop=FALSE], type='response')  ) )
  }
  
  return(as.character(c()))
}
view.predict.SupportVectorMachine <- function(ids.unlabelled, view.object) {
  ## TODO: UNTESTED
  # Intersect IDs in labels and in the feature data
  ids <- intersect(ids.unlabelled,rownames(view.object$data.matrix))
  if(length(ids) > 0){
    return( as.character( e1071::predict.svm(view.object$model, view.object$data.matrix[ids,,drop=FALSE])  ) )
  }

  return(as.character(c()))
}


### PLATYPUS has 2 options for label learning- either using a stacked model or our homebrew optimization function (the ensemble approach)
## TODO: This function should also be something the user can use to make new predictions, given a set of views?
#' Given a trained platypus model, make predictions. Generally used internally by platypus
#'
#' @param view.list List of views to use for predictions
#' @param majority Confidence threshold for predictions
#' @param test.ids Samples to be predicted
#' @param unique.labels Class labels
#' @param labels Labels for known samples
#' @param join.fxn Function used to make group votes
#'
#' @return Matrix of predictions (One per view, overall, etc)
#'
#' @export
platypus.predict <- function(view.list, majority, test.ids,unique.labels,labels,join.fxn='ensemble'){
  if(join.fxn=='ensemble') { return(platypus.predict.ensemble(view.list, majority, test.ids,unique.labels)) }
  else if(join.fxn=='stacked') { return(platypus.predict.stacked(view.list, majority, test.ids,unique.labels,labels)) }
  else {} # Only can handle ensemble/stacked learning for now TODO add error message & gracefully quit

}

## 20180709 - Kiley updates to add in stacked learning
## Take a platypus result and predict new labels with it
platypus.predict.stacked <- function(view.list, majority, test.ids,unique.labels,labels){

  ## Debug flag can be manually activated, for testing purposes 
  #flag.debug <- TRUE
  flag.debug <- FALSE 
  if(flag.debug) { print('Debug is on') }

  ## Return values:
  #  Test.ids is the list of sample IDs we are testing.
  #  predictions is the data frame of samples (rows) by views label predictions (eg. sensitive/non-sensitive)
  #  final is a list of outcome label predictions - NOTE: some of these have NA values! Names are sample IDs
  #  category.all is a list of agreed/not.agreed labels. Names are sample IDs
  #  category.majority is also a list of agreed/not.agreed labels. Names are sample IDs
  
  if(flag.debug) { print('platypus.predict') }
  #  predictions is the data frame of samples (rows) by views label predictions (eg. sensitive/non-sensitive)
  ## Make per-view predictions
  predictions <- matrix(data=NA, nrow=length(test.ids), ncol=length(view.list),dimnames=list(test.ids, paste0("view.",seq(length(view.list))))) # TODO: Update 'view.' to be view names
  for(view.i in seq(length(view.list))){
    ids <- intersect(test.ids,rownames(view.list[[view.i]]$data.matrix))
    if(length(ids) > 0){
      predictions[ids,view.i] <- view.predict(ids,view.list[[view.i]])
    }
  } 
  predictions[is.na(predictions)] <- 'not.predicted' # TEMP removing NA values so model predictions work
  if(flag.debug) { print('Created predictions matrix') }

  # samples fall in 3 categories if view predictions are combined in an ensemble
  # agreed: there was enough agreement between the single views to make an ensemble prediction
  # not.agreed: not enough agreement between the single views to make an ensemble prediction
  # not.assessed: not enough data was available in the single views to ever reach an ensemble prediction under the given requirements

  # category.all is a list of agreed/not.agreed labels. Names are sample IDs
  category.all <- rep(NA,length(test.ids))#dim(predictions)[[1]])
  names(category.all) <- test.ids # rownames(predictions)
  category.all[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x))>=length(view.list))),,drop=FALSE])] <- "not.agreed"
  category.all[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x))<length(view.list))),,drop=FALSE])] <- "not.assessed"
  category.all[rownames(predictions[which(apply(predictions,1,function(x) table(x)[1]==length(view.list))),,drop=FALSE])] <- "agreed"
  category.majority <- category.all # TODO: TEMP. This used to be different and part of weighting, not sure going to use in stacked version 
  if(flag.debug) { print('Created category.all and category.majority') }

  ## Build the stacked model, get predictions of samples w/greater or equal to threshodl agreement levels
  #if(!require(randomForest)) {
  #  install.packages('randomForest')
    #library(randomForest)
  #}
  #require(randomForest) 
  #for.model.predictions <- cbind(labels=labels[ids,1], predictions[ids,]) # Combine labels and predictions from each view to train the stacked model
  for.model.predictions <- cbind(labels=labels[rownames(predictions),1], predictions) # Combine labels and predictions from each view to train the stacked model
  #print(paste('Class for.model.predictions',class(for.model.predictions)))
  if(flag.debug) { print(unique(for.model.predictions[,1])) } # Print the unique labels 
  preds.model       <- randomForest( labels ~ ., for.model.predictions, ntree=100,replace=TRUE)
  if(flag.debug) { print('Make stacked model predictions') }
  preds.stckd       <- predict(preds.model, predictions, predict.all=TRUE)
  #preds.stckd       <- randomForest::predict.randomForest(preds.model, predictions, predict.all=TRUE)
  if(flag.debug) { print('Make preds.stckd') }
  preds.stacked.all <- preds.stckd$individual
  if(flag.debug) { print('Make preds.stckd.indiv') }

  if(flag.debug) { print('Stacked model predictions made') }
  tbl.preds.stckd           <- matrix(NA, nrow=nrow(preds.stacked.all), ncol=2)
  rownames(tbl.preds.stckd) <- rownames(preds.stacked.all)
  colnames(tbl.preds.stckd) <- unique.labels
  tbl.preds.stckd[,1]       <-   apply(preds.stacked.all, 1, function(x) {sum(x==unique.labels[1])} ) 
  tbl.preds.stckd[,2]       <-   apply(preds.stacked.all, 1, function(x) {sum(x==unique.labels[2])} ) 
  tbl.preds.stckd           <- tbl.preds.stckd/ncol(preds.stacked.all)
 
  #category.majority <- predict(preds.model, predictions, predict.all=TRUE)$aggregate #preds.lbls #majority.res.list$category.majority #TODO: updated the typo, not tested yet 
  #category.majority <- preds.stckd$aggregate #preds.lbls #majority.res.list$category.majority #TODO: updated the typo, not tested yet 

  if(flag.debug) { print('Created tbl.preds.stckd') }
  #  final is a list of outcome label predictions - NOTE: some of these have NA values! Names are sample IDs
  final <- rep(NA, length(test.ids))
  names(final) <- test.ids
  final[rownames(preds.stckd)] <- preds.stckd$aggregate # tbl.preds.stckd[apply(tbl.preds.stckd,1,max)>=majority,]
  if(flag.debug) { print('created final') }

  # Reset the NA values from before
  predictions[predictions=='not.predicted'] <- NA # Put these back to NA values
  if(flag.debug) { 
    print(unique(predictions))
    print('Re-Updated predictions matrix')
  }

  ## Debug only
  if(flag.debug) { 
    print('Dimensions of test.ids, predictions, final, category.all, category.majority:')
    print(length(test.ids))
    print(dim(predictions))
    print(length(final)) # TODO: Final in this version is a list, in stacked version is matrix. update!
    print(length(category.all))
    print(length(category.majority))

    print('TEST.IDS:')
    print(utils::head(test.ids))
    print('PREDICTIONS:')
    print(utils::head(predictions))
    print('FINAL:')
    print(utils::head(final))
    print('CATEGORY.ALL:')
    print(utils::head(category.all))
    print('CATEGORY.MAJORITY:')
    print(utils::head(category.majority))
  }

  predictions <- cbind(predictions,final,category.all,category.majority)

  if(flag.debug) { print('leaving platypus.predict') }
  return(predictions)

}

## Old platypus predict function - uses ensemble
platypus.predict.ensemble <- function(view.list, majority, test.ids,unique.labels){
  flag.debug <- FALSE

  predictions <- matrix(data=NA, nrow=length(test.ids), ncol=length(view.list),dimnames=list(test.ids, paste0("view.",seq(length(view.list))))) # TODO: Update 'view.' to be view names
  for(view.i in seq(length(view.list))){
    ids <- intersect(test.ids,rownames(view.list[[view.i]]$data.matrix))
    if(length(ids) > 0){
      predictions[ids,view.i] <- view.predict(ids,view.list[[view.i]]) 
    }
  }

  # samples fall in 3 categories if view predictions are combined in an ensemble
  # agreed: there was enough agreement between the single views to make an ensemble prediction
  # not.agreed: not enough agreement between the single views to make an ensemble prediction
  # not.assessed: not enough data was available in the single views to ever reach an ensemble prediction under the given requirements

  category.all <- rep(NA,dim(predictions)[[1]])
  names(category.all) <- rownames(predictions)
  category.all[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x))>=length(view.list))),,drop=FALSE])] <- "not.agreed"
  category.all[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x))<length(view.list))),,drop=FALSE])] <- "not.assessed"
  category.all[rownames(predictions[which(apply(predictions,1,function(x) table(x)[1]==length(view.list))),,drop=FALSE])] <- "agreed"

  majority.res.list <- get.majority.weighting(view.list,majority,predictions,unique.labels)
  category.majority <- majority.res.list$category.majority
  final <- majority.res.list$final

  ## Debug only
  if(flag.debug) {
    print('Dimensions of test.ids, predictions, final, category.all, category.majority:')
    print(length(test.ids))
    print(dim(predictions))
    print(length(final)) # TODO: Final in this version is a list, in stacked version is matrix. update!
    print(length(category.all))
    print(length(category.majority))
  
    print('TEST.IDS:')
    print(utils::head(test.ids))
    print('PREDICTIONS:')
    print(utils::head(predictions))
    print('FINAL:')
    print(utils::head(final))
    print('CATEGORY.ALL:')
    print(utils::head(category.all))
    print('CATEGORY.MAJORITY:')
    print(utils::head(category.majority))
  }

  predictions <- cbind(predictions,final,category.all,category.majority)
  return(predictions)

}

## Adjusting view accuracy in each platypus iteration to ensure correct weighting
## Calculate the current accuracy of each view by looking at the originally known labels predicted during training of the view (by elastic net cv or random forest oob)
## The views are trained on the known and newly learned labels, to get the change in accuracy by adding new labels
update.accuracies <- function(view.list,known.labels){
  for(view.i in seq(length(view.list))){
    view <- view.list[[view.i]]
    view <- update.accuracy(view,known.labels)
    view.list[[view.i]] <- view
  print( paste('Updating view accuracies', view) )
  }
  
  return(view.list)
}
update.accuracy <- function(view, known.labels){
  UseMethod("update.accuracy",view)
}
## method predicts on samples also used in training to avoid manual cross-validation
## assessment of prediction accuracy of the view is therefore highly overestimated
## TODO: Make sure these update methods are correct
update.accuracy.ElasticNet <- function(view, known.labels){
  
  predictions <- view.predict(rownames(known.labels),view)
  label.tab <- merge(known.labels, as.data.frame(predictions), by.x="row.names", by.y="row.names", all=TRUE, sort=TRUE) 
  
  acc.list <- calculate.accuracy(label.tab,known.labels)
  view$acc <- acc.list$b.acc
  view$acc.norm <- normalize.accuracy.log(acc.list$b.acc)
  
  return(view)
}
## used oob predictions of random forest - correct estimation of model accuracy
update.accuracy.RandomForest <- function(view, known.labels){
  oob.predictions <- view$model$predicted
  label.tab <- merge(known.labels, as.data.frame(oob.predictions), by.x="row.names", by.y="row.names", all=TRUE, sort=TRUE) 
  
  acc.list <- calculate.accuracy(label.tab,known.labels)
  view$acc <- acc.list$b.acc
  view$acc.norm <- normalize.accuracy.log(acc.list$b.acc)
  
  return(view)
}
## TODO: update.accuracy.svm function needs creation

## TODO: Double check this is correct. Why is it so complicated???
calculate.accuracy <- function(label.tab,known.labels){
  correct <- c()
  total <- c()
  for(label in as.character(unique(known.labels[,1]))){
    c <- length(which(label.tab[,2] == label & label.tab[,2] == label.tab[,3]))
    correct <- c(correct,c)
    t <- length(which(label.tab[,2] == label))
    total <- c(total,t)
  }
  
  acc <- sum(correct) / sum(total)
  b.acc <- ( ( correct[1] / total[1] ) + (correct[2] / total[2] ) )/2
  return(list(b.acc=b.acc,acc=acc))
}

## get ensemble predictions for a majority setup - weighting votes by accuracy
# TODO: Update this to work with a stacked model version
get.majority.weighting <- function(view.list,majority,predictions,unique.labels){
  
  category.majority <- rep(NA,dim(predictions)[[1]])
  names(category.majority) <- rownames(predictions)
  final <- rep(NA,dim(predictions)[[1]])
  names(final) <- rownames(predictions)
  
  for(i in seq(length(category.majority))){
    p1 <- 0
    p2 <- 0
    for(view.i in seq(length(view.list))){
      if(!(is.na(predictions[i,view.i]))){
        if(predictions[i,view.i] == unique.labels[1]){
          p1 <- p1 + view.list[[view.i]]$acc.norm
        }
        if(predictions[i,view.i] == unique.labels[2]){
          p2 <- p2 + view.list[[view.i]]$acc.norm
        }
      }
    }
    if(p1 + p2 < majority){
      category.majority[i] <- "not.assessed"
    } else {
      if(p1 < majority & p2 < majority){
        category.majority[i] <- "not.agreed"
      }
    }
    if(p1 > p2 & p1 > majority){
      category.majority[i] <- "agreed"
      final[i] <- unique.labels[1]
    }
    if(p2 > p1 & p2 > majority){
      category.majority[i] <- "agreed"
      final[i] <- unique.labels[2]
    }
    
  }
  
  return(list(category.majority=category.majority,final=final))
}

# weigh the predictions by accuracy of the single views
# TODO: Update this to work with a stacked model version
get.new.labels.majorityWeighted <- function(predictions,view.list,unique.labels){

  # get normalized accuracies of the views
  acc.norm <- c()
  for(view.i in seq(length(view.list))){
    acc.norm <- c(acc.norm,view.list[[view.i]]$acc.norm)
  }
  
  # sum up the accuracies for both labels
  p1.vec <- apply(predictions,1,function(x) 
    sum(acc.norm[which(x == unique.labels[1])])
  )
  p2.vec <- apply(predictions,1,function(x) 
    sum(acc.norm[which(x == unique.labels[2])])
  )
  
  # retrieve the maximal reached sum of accuracies and get all the predictions reaching it
  max.pred.value <- max(c(p1.vec,p2.vec))
  max.preds <- cbind(predictions, p1.vec, p2.vec)
  
  max.preds <- max.preds[which(max.preds[,"p1.vec"] >= max.pred.value | max.preds[,"p2.vec"] >= max.pred.value),,drop=FALSE]
  
  # from those, retrieve the predictions with the minimal value for a contrary prediction
  min.pred.value <- min(c(max.preds[,"p1.vec"],max.preds[,"p2.vec"]))
  
  min.max.preds <- max.preds[which(max.preds[,"p1.vec"] <= min.pred.value | max.preds[,"p2.vec"] <= min.pred.value),,drop=FALSE]
  
  # return them as the new labels
  new.labelled <- min.max.preds[,1,drop=FALSE]
  new.labelled[which(min.max.preds[,"p1.vec"] > min.max.preds[,"p2.vec"]),1] <- unique.labels[1]
  new.labelled[which(min.max.preds[,"p2.vec"] > min.max.preds[,"p1.vec"]),1] <- unique.labels[2]
  

  return(list(new.labelled=new.labelled, weighting.threshold.upper=max.pred.value, weighting.threshold.lower=min.pred.value))
  
}

## normalize accuracies from [0.5,1] to [0,1] and log-transform them
## TODO: Do I need this? Really only call it 1 time
## TODO: I changed this, make sure it actually works
normalize.accuracies <- function(view.list){
  for(view.i in seq(length(view.list))){
    view.list[[view.i]]$acc.norm <- normalize.accuracy.log(view.list[[view.i]]$acc)
  }
  return(view.list)
}

# normalize to range [0,1] and log-transform
normalize.accuracy.log <- function(acc){
  if(is.na(acc)) {
    stop(paste0("Non-numeric accuracy value")) #TODO: TEMP
  }
  if(acc >= 0 & acc <= 1){
    if(acc <= 0.5) { 
      acc <- 0.51 # Don't want negative values
    }
    if(acc == 1){
      acc <- 0.99  ## avoids to get Inf-values
    }
    acc.norm <- (acc - 0.5) / (1-0.5)
    acc.norm <- -log(1-acc.norm)
  } else{
    stop(paste0("Accuracy value not between 0.5 and 1"))
  }
  return(acc.norm)
}

## calculate performance values for a prediction ( used in cv.platypus )
## TODO: update this function
calculate.performance <- function(predictions,labels,unique.labels){
  labels <- as.character(labels[rownames(predictions),])
  correct <- rep(0,dim(predictions)[[1]])
  names(correct) <- rownames(predictions)
  correct[which(labels == predictions[,1])] <- 1
  
  predictions <- cbind(predictions,labels, correct)

  predictions.assessed <- predictions[which(predictions[,2] != "not.assessed"),,drop=FALSE]
  predictions.agreed <- predictions[which(predictions[,2] == "agreed"),,drop=FALSE]
  
  acc <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1),,drop=FALSE])[[1]] / dim(predictions.agreed)[[1]]
  cov <- dim(predictions.agreed)[[1]] / dim(predictions.assessed)[[1]]
  
  # calculate balanced accuracy
  if(dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]] != 0){
    sens <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1 & predictions.agreed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]] / dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]]
  } else{
    sens <- 0
  }
  if(dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[2]),,drop=FALSE])[[1]] != 0){
    spec <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1 & predictions.agreed[,"labels",drop=FALSE] == unique.labels[2]),,drop=FALSE])[[1]] / dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[2]),,drop=FALSE])[[1]]
  } else{
    spec <- 0
  }
  
  bal.acc <- ( sens + spec ) / 2
  
  no.correct <- length(which(predictions.agreed[,"correct"] == 1))
  correct.classes <- c()
  for(class in unique.labels){
    correct.classes <- c(correct.classes,length(which(predictions.agreed[,"correct"] == 1 & predictions.agreed[,1] == class)))
  }
  names(correct.classes) <- unique.labels
  no.incorrect <- length(which(predictions.agreed[,"correct"] == 0))
  incorrect.classes <- c()
  for(class in unique.labels){
    incorrect.classes <- c(incorrect.classes,length(which(predictions.agreed[,"correct"] == 0 & predictions.agreed[,1] == class)))
  }
  names(incorrect.classes) <- unique.labels
  no.not.agreed <- dim(predictions.assessed)[[1]] - dim(predictions.agreed)[[1]]
  miss <- dim(predictions)[[1]] - dim(predictions.assessed)[[1]]
  
  return(list(accuracy=acc, coverage=cov, sensitivity=sens, specificity=spec, balanced.accuracy=bal.acc, correct=no.correct, correct.classes = correct.classes
              , incorrect=no.incorrect, incorrect.classes=incorrect.classes, not.agreed=no.not.agreed, missingData=miss))
}

## calculate performance values for each view
calculate.performance.view <- function(predictions,labels,unique.labels){
  labels <- as.character(labels[rownames(predictions),])
  correct <- rep(0,dim(predictions)[[1]])
  names(correct) <- rownames(predictions)
  correct[which(labels == predictions[,1])] <- 1
  
  predictions <- cbind(predictions,labels, correct)
  
  predictions.assessed <- predictions[which(!(is.na(predictions[,1]))),]
  
  acc <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1),,drop=FALSE])[[1]] / dim(predictions.assessed)[[1]]
  # calculate balanced accuracy
  if(dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]] != 0){
    sens <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1 & predictions.assessed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]] / dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[1]),,drop=FALSE])[[1]]
  } else{
    sens <- 0
  }
  if(dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[2]),,drop=FALSE])[[1]] != 0){
    spec <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1 & predictions.assessed[,"labels"] == unique.labels[2]),,drop=FALSE])[[1]] / dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[2]),,drop=FALSE])[[1]]
  } else{
    spec <- 0
  }
  
  bal.acc <- ( sens + spec ) / 2
  
  no.correct <- length(which(predictions.assessed[,"correct"] == 1))
  correct.classes <- c()
  for(class in unique.labels){
    correct.classes <- c(correct.classes,length(which(predictions.assessed[,"correct"] == 1 & predictions.assessed[,1] == class)))
  }
  names(correct.classes) <- unique.labels
  no.incorrect <- length(which(predictions.assessed[,"correct"] == 0))
  incorrect.classes <- c()
  for(class in unique.labels){
    incorrect.classes <- c(incorrect.classes,length(which(predictions.assessed[,"correct"] == 0 & predictions.assessed[,1] == class)))
  }
  names(incorrect.classes) <- unique.labels
  miss <- dim(predictions)[[1]] - dim(predictions.assessed)[[1]]
  
  return(list(accuracy=acc, sensitivity=sens, specificity=spec, balanced.accuracy=bal.acc, correct=no.correct, correct.classes = correct.classes
              , incorrect=no.incorrect, incorrect.classes=incorrect.classes, missingData=miss))
  
  
}


## calculate performance values for platypus label learning validation (llv.platypus)
get.labelling.performance <- function(labelling.matrix,labels,unique.labels){
  
  ids <- intersect(rownames(labels),rownames(labelling.matrix))
  labelling.matrix <- labelling.matrix[ids,,drop=FALSE]
  labels <- labels[ids,,drop=FALSE]
  
  result.table <- c()
  for(i in seq(dim(labelling.matrix)[[2]])){
    vec <- labelling.matrix[,i] == labels[,1]
    acc <- length(vec[!(is.na(vec)) & vec==TRUE])/length(vec[!(is.na(vec))])
    cov <- length(vec[!(is.na(vec))]) / length(vec)
    
    unique.labels <- unique(labels[,1])
    acc.labels <- c()
    cov.labels <- c()
    for(l in seq(length(unique.labels))){
      pred.pos <- rownames(labelling.matrix[which(labelling.matrix[,i] == unique.labels[l]),,drop=FALSE])
      
      is.pos <- rownames(labels[which(labels[,1] == unique.labels[l]),,drop=FALSE])
      
      pred.na <- rownames(labelling.matrix[is.na(labelling.matrix[,i]),,drop=FALSE])
      
      if(length(is.pos[!(is.pos %in% pred.na)]) != 0){
        acc.label <- length(intersect(pred.pos,is.pos))/ length(is.pos[!(is.pos %in% pred.na)])
      } else{
        acc.label <- 0
      }
      
      acc.labels <- c(acc.labels,acc.label )
      cov.labels <- c(cov.labels, 1-(length(intersect(pred.na, is.pos)) / length(is.pos) ))
    }
    
    bal.acc <- mean(acc.labels)
    
    result.table <- rbind(result.table, c(i,acc,bal.acc,cov,acc.labels,cov.labels))
    colnames(result.table) <- c("iteration","accuracy","balanced.accuracy","coverage",paste0(unique.labels,"_accuracy"),paste0(unique.labels,"_coverage"))
  }
  return(result.table)
}




