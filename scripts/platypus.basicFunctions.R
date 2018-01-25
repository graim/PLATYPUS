## Basic Functions for the MVL framework
##
## Created: Dec 2014, Kiley Graim
##    Updated: September 2015, Verena Friedl
##    Updated: Feb 2017, Kiley Graim
##    Last Updated: Jan 2018, Kiley Graim


## Subset to k most variable features
drop.features <- function(dat, k) {
  if(ncol(dat) > k ) {
    vars <- apply(dat, 2, var)
    vars <- sort(vars, decreasing=T)
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


## create a class for each view-type
ElasticNet <- function(param.file="",data.matrix=c(),data.fn="",alpha=0.9,measure="mse",drop=TRUE,drop.to=5000,model=c(),acc=0.5, acc.norm=0, family='binomial') {
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

## Set methods for the view objects
setAlpha <- function(view.object, newValue){
  view.object$alpha <- newValue
  return(view.object)
}

setMeasure <- function(view.object, newValue){
  view.object$measure <- newValue
  return(view.object)
}

setMtry <- function(view.object, newValue){
  view.object$mtry <- newValue
  return(view.object)
}

setNtree <- function(view.object, newValue){
  view.object$ntree <- newValue
  return(view.object)
}

setDrop <- function(view.object, newValue){
  view.object$drop <- newValue
  return(view.object)
}

setDropTo <- function(view.object, newValue){
  view.object$drop.to <- newValue
  return(view.object)
}

setAcc <- function(view.object, newValue){
  view.object$acc <- newValue
  return(view.object)
}

setAccNorm <- function(view.object, newValue){
  view.object$acc.norm <- newValue
  return(view.object)
}

## Read in the parameter file for a view and return a view object
load.parameterfile <- function(filename){
  
  param.table <- read.table(filename, sep='\t',header=F, row.names=1)
  #print( paste('Loaded parameter file', filename) )
  
  # check the type and create a view object
  type <- param.table["type",1]
  if(type == "en"){
    view.object <- ElasticNet(param.file=filename,data.fn=toString(param.table["data.fn",1]))
    if("alpha" %in% rownames(param.table)){
      view.object <- setAlpha(view.object,as.numeric(as.character(param.table["alpha",1])))
    }
    if("measure" %in% rownames(param.table)){
      view.object <- setMeasure(view.object,as.character(param.table["measure",1]))       
    }
    if("drop" %in% rownames(param.table)){
      view.object <- setDrop(view.object,as.logical(param.table["drop",1]))
    }
    if("drop.to" %in% rownames(param.table)){
      view.object <- setDropTo(view.object,as.numeric(as.character(param.table["drop.to",1])))
    }
    if("acc" %in% rownames(param.table)){
      view.object <- setAcc(view.object,as.numeric(as.character(param.table["acc",1])))
    }
    return(view.object)
  }
  if(type == 'rf'){
    view.object <- RandomForest(param.file=filename,data.fn=toString(param.table["data.fn",1]))
    if("mtry" %in% rownames(param.table)){
      view.object <- setMtry(view.object,as.numeric(as.character(param.table["mtry",1])))
    }
    if("ntree" %in% rownames(param.table)){
      view.object <- setNtree(view.object,as.numeric(as.character(param.table["ntree",1])))
    }
    if("drop" %in% rownames(param.table)){
      view.object <- setDrop(view.object,as.logical(param.table["drop",1]))
    }
    if("drop.to" %in% rownames(param.table)){
      view.object <- setDropTo(view.object,as.numeric(as.character(param.table["drop.to",1])))
    }
    if("acc" %in% rownames(param.table)){
      view.object <- setAcc(view.object,as.numeric(as.character(param.table["acc",1])))
    }
    
    return(view.object)
  }
  # if nothing was returned so far, send error message
  stop(paste0("Could not find type specified 'en' or 'rf' in parameter file ",filename))
}


## load the feature matrix from file to a matrix object
load.data <- function(view.object){
  UseMethod("load.data",view.object)
}
load.data.ElasticNet <- function(view.object){
  
  # read the matrix from the file
  mat <- data.matrix( read.table(view.object$data.fn, sep='\t',header=T, row.names=1, check.names=F) )
  print(paste('data loaded for Elastic Net view',view.object$data.fn))
  # TODO: It's not always tab-delimited!!!
  
  # drop features
  if(view.object$drop){
    mat <- drop.features(mat, view.object$drop.to)
  }
  
  # add an 'X' to all rownames starting with a number (To stay consistent, bc R adds 'X' to matrix columns starting with a number)
  #rownames(mat) <- addX(rownames(mat))
  
  # set the matrix
  view.object$data.matrix <- mat
  
  return(view.object)
}
load.data.RandomForest <- function(view.object){
  # read the matrix from the file
  mat <- data.matrix( read.table(view.object$data.fn, sep='\t',header=T, row.names=1, check.names=F) )
  print(paste('data loaded for Random Forest view',view.object$data.fn))
  ## TODO: It's not always tab-delimited!!!
  
  # drop features
  if(view.object$drop){
    mat <- drop.features(mat, view.object$drop.to)
  }
  
  # add an 'X' to all rownames starting with a number (To stay consistent, bc R adds 'X' to matrix columns starting with a number)
  rownames(mat) <- addX(rownames(mat))
  
  # set the matrix
  view.object$data.matrix <- mat
  
  # set mtry parameter to sqrt(#features) if default
  if(view.object$mtry == "sqrt"){
    view.object <- setMtry(view.object,floor(sqrt(ncol(mat))))
  }
  return(view.object)
}

## load the label file
load.label.data <- function(fn.labs,classcol.labs){
  # read file
  labs <- read.table(fn.labs, sep='\t',header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE)
  # take out 'NA' values
  labs <- labs[which(!(is.na(labs[,classcol.labs]))),,drop=F]
  # add X to rows starting with number
#  rownames(labs) <- addX(rownames(labs))
  
  return(labs)
}
## retrieve the two labels for training/predicting
get.unique.labels <- function(label.vec,ignore.label){
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
  
  # MVL classification for 2 labels
  if(length(unique.labels) != 2){
    print(paste0("ATTENTION: ",length(unique.labels)," labels for a binary task. Results might not be correct!"))
  }
  return(unique.labels)
}

## find strings starting with a number and add an 'X' in the front
# this is done to ensure consistence over all rownames and colnames for R matrices (R adds 'X' to colnames starting with a number)
addX <- function(vector){
  vector[grep("^[0-9].",vector)] <- paste0("X",vector[grep("^[0-9].",vector)])
  return(vector)
}


## Train one view, given a specific model type
view.train <- function( labels, view.object ) {
  UseMethod("view.train",view.object)
}
view.train.ElasticNet <- function(labels, view.object ){
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(rownames(labels),rownames(view.object$data.matrix))

#  print( summary( labels[ids,1] ) ) # TODO
  
  # Train model
  require(glmnet)
  #print(paste(levels(as.factor(labels[ids,1])),table(as.factor(labels[ids,1])))) # TODO: TEMP PRINT STATEMENT
  view.object$model <- cv.glmnet( view.object$data.matrix[ids,],
                                  as.factor(labels[ids,1]), family=view.object$family, 
                                  type.measure=view.object$measure, 
                                  alpha=view.object$alpha 
                                  , keep=T)  #defaults to nfold=10
  return(view.object)
}
view.train.RandomForest <- function( labels, view.object  ){
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(rownames(labels),rownames(view.object$data.matrix))

#  print( summary( labels[ids,1] ) ) # TODO

  # Train model
  require(randomForest)
  view.object$model <- randomForest( view.object$data.matrix[ids,], 
                                     as.factor(labels[ids,1]), 
                                     mtry=view.object$mtry, 
                                     ntree=view.object$ntree )
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
    return( predict(view.object$model, view.object$data.matrix[ids,,drop=F], type='class', s='lambda.min') )  # TODO: use lambdaMin or default s="lambda.1se" ? 
  }
  
  return(as.character(c()))
}
view.predict.RandomForest <- function(ids.unlabelled, view.object) {
  
  # Intersect IDs in labels and in the feature data
  ids <- intersect(ids.unlabelled,rownames(view.object$data.matrix))
  if(length(ids) > 0){
    return( as.character( predict(view.object$model, view.object$data.matrix[ids,,drop=F], type='response')  ) )
  }
  
  return(as.character(c()))
}


## Take a platypus result and predict new labels with it
platypus.predict <- function(view.list, majority, test.ids,weighting,unique.labels){

  predictions <- matrix(data=NA, nrow=length(test.ids), ncol=length(view.list),dimnames=list(test.ids, paste0("view.",1:length(view.list))))
  for(view.i in 1:length(view.list)){
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

  if(weighting){
    majority.res.list <- get.majority.weighting(view.list,majority,predictions,unique.labels)
    category.majority <- majority.res.list$category.majority
    final <- majority.res.list$final
  } else {
    majority.res.list <- get.majority.counting(majority,predictions)
    category.majority <- majority.res.list$categroy.majority
    final <- majority.res.list$final
  }

  predictions <- cbind(predictions,final,category.all,category.majority)
  return(predictions)
  
}

## Adjusting view accuracy in each platypus iteration to ensure correct weighting
## Calculate the current accuracy of each view by looking at the originally known labels predicted during training of the view (by elastic net cv or random forest oob)
## The views are trained on the known and newly learned labels, to get the change in accuracy by adding new labels
update.accuracies <- function(view.list,known.labels){
  for(view.i in 1:length(view.list)){
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
update.accuracy.ElasticNet <- function(view, known.labels){
  
  preditions <- view.predict(rownames(known.labels),view)
  label.tab <- merge(known.labels, as.data.frame(preditions), by.x="row.names", by.y="row.names", all=T, sort=T) 
  
  acc.list <- calculate.accuracy(label.tab,known.labels)
  view <- setAcc(view,acc.list$b.acc)
  view <- setAccNorm(view,normalize.accuracy.log(acc.list$b.acc))
  
  return(view)
}
## used oob predictions of random forest - correct estimation of model accuracy
update.accuracy.RandomForest <- function(view, known.labels){
  oob.predictions <- view$model$predicted
  label.tab <- merge(known.labels, as.data.frame(oob.predictions), by.x="row.names", by.y="row.names", all=T, sort=T) 
  
  acc.list <- calculate.accuracy(label.tab,known.labels)
  view <- setAcc(view,acc.list$b.acc)
  view <- setAccNorm(view,normalize.accuracy.log(acc.list$b.acc))
  
  return(view)
}
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

## get ensemble predictions for a majority setup - counting votes
get.majority.counting <- function(majority,predictions){
  category.majority <- rep(NA,dim(predictions)[[1]])
  names(category.majority) <- rownames(predictions)
  category.majority[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x)) >= majority)),,drop=FALSE])] <- "not.agreed"
  category.majority[rownames(predictions[which(apply(predictions,1,function(x) sum(table(x)) < majority)),,drop=FALSE])] <- "not.assessed"
  category.majority[rownames(predictions[which(apply(predictions,1,function(x) sort(table(x),decreasing=T)[1] >= majority)),,drop=FALSE])] <- "agreed"
  
  final <- rep(NA,dim(predictions)[[1]])
  names(final) <- rownames(predictions)
  final[names(category.majority[which(category.majority == "agreed")])] <- apply(predictions[names(category.majority[which(category.majority == "agreed")]),],1,function(x) names(sort(table(x),decreasing=T))[1])
  
  return(list(category.majority=category.majority,final=final))
}
## get ensemble predictions for a majority setup - weighting votes by accuracy
get.majority.weighting <- function(view.list,majority,predictions,unique.labels){
  
  category.majority <- rep(NA,dim(predictions)[[1]])
  names(category.majority) <- rownames(predictions)
  final <- rep(NA,dim(predictions)[[1]])
  names(final) <- rownames(predictions)
  
  for(i in 1:length(category.majority)){
    p1 <- 0
    p2 <- 0
    for(view.i in 1:length(view.list)){
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

## take single predictions for unlabeled samples and return all where view agreement meets requirements
# count number of views agreeing in a predictions, check if thresholds are met
get.new.labels.majorityCount <- function(predictions,majority,majority.missingData){
  
  new.labelled <- c()
  
  if(majority != majority.missingData){
    missingData.predictions <- predictions[which(apply(predictions,1,function(x) sum(table(x)) < length(view.list))),,drop=FALSE]
    new.labelled <- missingData.predictions[which(apply(missingData.predictions,1,function(x) length(names(table(x))[table(x) >= majority.missingData])>=1)),,drop=FALSE]
  }
  
  new.labelled <- rbind(new.labelled,predictions[which(apply(predictions,1,function(x) length(names(table(x))[table(x) >= majority])>=1)),,drop=FALSE])
  
  new.labelled[,1] <- apply(new.labelled,1,function(x) names(sort(table(x),decreasing=TRUE)[1]))
  new.labelled <- new.labelled[,1,drop=FALSE]
  new.labelled <- new.labelled[unique(rownames(new.labelled)),,drop=FALSE]
  
  return(new.labelled)
}

# weigh the predictions by accuracy of the single views
get.new.labels.majorityWeighted <- function(predictions,view.list,unique.labels){

  # get normalized accuracies of the views
  acc.norm <- c()
  for(view.i in 1:length(view.list)){
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
  
  max.preds <- max.preds[which(max.preds[,"p1.vec"] >= max.pred.value | max.preds[,"p2.vec"] >= max.pred.value),,drop=F]
  
  # from those, retrieve the predictions with the minimal value for a contrary prediction
  min.pred.value <- min(c(max.preds[,"p1.vec"],max.preds[,"p2.vec"]))
  
  min.max.preds <- max.preds[which(max.preds[,"p1.vec"] <= min.pred.value | max.preds[,"p2.vec"] <= min.pred.value),,drop=F]
  
  # return them as the new labels
  new.labelled <- min.max.preds[,1,drop=F]
  new.labelled[which(min.max.preds[,"p1.vec"] > min.max.preds[,"p2.vec"]),1] <- unique.labels[1]
  new.labelled[which(min.max.preds[,"p2.vec"] > min.max.preds[,"p1.vec"]),1] <- unique.labels[2]
  

  return(list(new.labelled=new.labelled, weighting.threshold.upper=max.pred.value, weighting.threshold.lower=min.pred.value))
  
}


## normalize accuracies from [0.5,1] to [0,1] and log-transform them
normalize.accuracies <- function(view.list){
  
  for(view.i in 1:length(view.list)){
    view <- setAccNorm(view.list[[view.i]],normalize.accuracy.log(view.list[[view.i]]$acc))
    view.list[[view.i]] <- view
  }
  
  return(view.list)
}

# normalize to range [0,1]
normalize.accuracy.linear <- function(acc){
  if(acc >= 0.5 & acc <= 1){
    acc.norm <- (acc - 0.5) / (1-0.5)
  } else{
    stop(paste0("Accuracy value not between 0.5 and 1"))
  }
  return(acc.norm)
}
# normalize to range [0,1] and log-transform
normalize.accuracy.log <- function(acc){
  if(acc >= 0.5 & acc <= 1){
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
calculate.performance <- function(predictions,labels,unique.labels){
  labels <- as.character(labels[rownames(predictions),])
  correct <- rep(0,dim(predictions)[[1]])
  names(correct) <- rownames(predictions)
  correct[which(labels == predictions[,1])] <- 1
  
  predictions <- cbind(predictions,labels, correct)

  predictions.assessed <- predictions[which(predictions[,2] != "not.assessed"),,drop=F]
  predictions.agreed <- predictions[which(predictions[,2] == "agreed"),,drop=F]
  
  acc <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1),,drop=F])[[1]] / dim(predictions.agreed)[[1]]
  cov <- dim(predictions.agreed)[[1]] / dim(predictions.assessed)[[1]]
  
  # calculate balanced accuracy
  if(dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[1]),,drop=F])[[1]] != 0){
    sens <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1 & predictions.agreed[,"labels"] == unique.labels[1]),,drop=F])[[1]] / dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[1]),,drop=F])[[1]]
  } else{
    sens <- 0
  }
  if(dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[2]),,drop=F])[[1]] != 0){
    spec <- dim(predictions.agreed[which(predictions.agreed[,"correct"] == 1 & predictions.agreed[,"labels",drop=F] == unique.labels[2]),,drop=F])[[1]] / dim(predictions.agreed[which(predictions.agreed[,"labels"] == unique.labels[2]),,drop=F])[[1]]
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
  
  acc <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1),,drop=F])[[1]] / dim(predictions.assessed)[[1]]
  # calculate balanced accuracy
  if(dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[1]),,drop=F])[[1]] != 0){
    sens <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1 & predictions.assessed[,"labels"] == unique.labels[1]),,drop=F])[[1]] / dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[1]),,drop=F])[[1]]
  } else{
    sens <- 0
  }
  if(dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[2]),,drop=F])[[1]] != 0){
    spec <- dim(predictions.assessed[which(predictions.assessed[,"correct"] == 1 & predictions.assessed[,"labels"] == unique.labels[2]),,drop=F])[[1]] / dim(predictions.assessed[which(predictions.assessed[,"labels"] == unique.labels[2]),,drop=F])[[1]]
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
  labelling.matrix <- labelling.matrix[ids,,drop=F]
  labels <- labels[ids,,drop=F]
  
  result.table <- c()
  for(i in 1:dim(labelling.matrix)[[2]]){
    vec <- labelling.matrix[,i] == labels[,1]
    acc <- length(vec[!(is.na(vec)) & vec==TRUE])/length(vec[!(is.na(vec))])
    cov <- length(vec[!(is.na(vec))]) / length(vec)
    
    unique.labels <- unique(labels[,1])
    acc.labels <- c()
    cov.labels <- c()
    for(l in 1:length(unique.labels)){
      pred.pos <- rownames(labelling.matrix[which(labelling.matrix[,i] == unique.labels[l]),,drop=F])
      
      is.pos <- rownames(labels[which(labels[,1] == unique.labels[l]),,drop=F])
      
      pred.na <- rownames(labelling.matrix[is.na(labelling.matrix[,i]),,drop=F])
      
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




