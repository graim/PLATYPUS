### single support vector machine predictor #######################################################
# generalized function get best parameters for a given task/view pair
###################################################################################################
#fuck it, using caret and a generic function for it. 
#' Parameter sweeps of svm, random forest, and elastic net models for config generation.
#'
#' @param X Feature matrix (sample per row)
#' @param y Named vector of labels
#' @param model Type of single predictor model to build. Choose from en (elastic net), rf (random forest), or svm (support vector machine)
#' @param ignore.label Label class to ignore. Will be removed from labels set before training.
#'
#' @return list containing optimal model summary
#'
#' @import caret
#'
#' @examples
#' X <- matrix(rnorm(10000), nrow=100)
#' y <- c(rep('MOO',50),rep('OINK',50))
#' names(y) <- paste0('Sample',seq(nrow(X)))
#' rownames(X) <- paste0('Sample',seq(nrow(X)))
#' colnames(X) <- paste0('Feature',seq(ncol(X)))
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
  y   <- y[stats::complete.cases(y)]
  y   <- y[y!=ignore.label] # remove label user wants to ignore
  X   <- X[stats::complete.cases(X),,drop=F]
  X   <- X[intersect(rownames(X),names(y)),,drop=F]
  y   <- y[intersect(rownames(X),names(y))]

  y   <- factor(y) # Remove the excess factor level ignore.label

#  ## Make sure we're using factor labels
#  if(!is.factor(y)) {
#    message('Converting labels to factor')
#    y <- as.factor(y)
#  }
 
  ## Parameter sweep and return 'best' model
  res.model <- caret::train(x = X, y = y, method=model, preProcess = 'center', metric = 'Accuracy', trControl = caret::trainControl(method = 'repeatedcv', number = 10))
  res <- res.model$results[which(res.model$results$Accuracy == max(res.model$results$Accuracy))[1],]
  res$Model <- model # For gen.config use
  return( res )
#  return(res.model$results[which(res.model$results$Accuracy == max(res.model$results$Accuracy))[1],])

}
