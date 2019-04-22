### single support vector machine predictor #######################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
#fuck it, using caret and a generic function for it. 
#' Parameter sweeps of svm, random forest, and elastic net models for config generation.
#'
#' @param X Feature matrix (sample per row)
#' @param y Named vector of labels
#'
#' @return list containing optimal model summary
#'
#' @examples
#' X <- matrix(rnorm(10000), nrow=100)
#' y <- c(rep('MOO',50),rep('OINK',50))
#' names(y) <- paste0('Sample',seq(nrow(X)))
#' rownames(X) <- paste0('Sample',seq(nrow(X)))
#' y <- factor(y)
#' res <- single.predictor(X,y,model='rf')
#' res <- single.predictor(X,y,model='en')
#' res <- single.predictor(X,y,model='svm')
#'
#' @export
single.predictor <- function(X,y,model,ignore.label='intermediate') {

  ## switch to caret syntax
  model <- switch(model,
    en = 'glmnet',
    rf = 'rf',
    svm = 'svmRadialCost'
  )

  ## Remove unlabeled samples and those with missing data
  ## Drop to set of samples in both labels and features data
  y   <- y[y!=ignore.label] # remove label user wants to ignore
  y   <- y[stats::complete.cases(y)]
  X   <- X[stats::complete.cases(X),,drop=F]
  X   <- X[intersect(rownames(X),names(y)),,drop=F]
  y   <- y[intersect(rownames(X),names(y))]

  ## Make sure we're using factor labels
  if(!is.factor(y)) {
    message('Converting labels to factor')
    y <- as.factor(y)
  }
 
  ## Parameter sweep and return 'best' model
  res.model <- caret::train(x = t(X), y = y, method=model, preProcess = 'center', metric = 'Accuracy', trControl = trainControl(method = 'repeatedcv', number = 10))
  return(res.model$results[which(res.model$results$Accuracy == max(res.model$results$Accuracy))[1],])

}
