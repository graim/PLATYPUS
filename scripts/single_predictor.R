### Functions used by gen_config

### single_elasticNet_predictor.R ################################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
## Inputs:
#   X = data frame/matrix of features where each row is 1 sample
#   y = vector of named labels
## Optional Inputs:
#   alpha = values of alpha to test
#   iterations = number of iterations per alpha to run the test
single.elasticNet.predictor <- function(X,y,alpha = seq(0,1,0.1),iterations = 10,nfolds=10, measure='auc', viewName='', ignore.label='intermediate') {

  ## Remove unlabeled samples and those with missing data
  ## Drop to set of samples in both labels and features data
  #print('Dropping to overlapping samples')
  y   <- y[y!=ignore.label] # remove label user wants to ignore
  y   <- y[complete.cases(y)]
  X   <- X[complete.cases(X),,drop=F]
  X   <- X[intersect(rownames(X),names(y)),,drop=F]
  y   <- y[intersect(rownames(X),names(y))]

  ## TODO: when do I use lambda???
  res.all <- data.frame(alpha=rep(alpha, each=iterations), 'num.samples'=nrow(X), error=NA, accuracy=NA, lambda=NA)

  ## Calculate results for each alpha, repeated iterations times - TODO this is really really slow :)
  for( i in 1:nrow(res.all)) {
    ## Shrink folds if not enough data
    if(nrow(X)/nfolds < 10) {
      #print(table(y))
      nfolds <- floor(nrow(X)/10)
      print(paste('Too few samples, using',nfolds,'instead.'))
    }
    res <- cv.glmnet(as.matrix(X),y,family="binomial",type.measure="auc", alpha=res.all[i,'alpha'], nfolds=nfolds) 
    res.all[i,'error'] <- res$cvsd[ res$glmnet.fit$lambda == res$lambda.min ]
    res.all[i,'accuracy'] <- res$cvm[res$glmnet.fit$lambda == res$lambda.min]
    res.all[i,'lambda'] <- res$lambda.min
  } # end for each alpha

  res <- aggregate(res.all, by=list(res.all$alpha), FUN=median)[,c('alpha','accuracy')] # TODO
  return(res[res$accuracy==max(res$accuracy),])
  
} # end single.elasticNet.predictor

#                 alpha     error  std.error    lambda  accuracy num.cell.lines

### single random forest predictor ################################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
## Inputs:
#   X = data frame/matrix of features where each row is 1 sample
#   y = vector of named labels
## Optional Inputs:
#   mtry - num predictors sampled for splitting at each node
#   ntree - number of trees grown
## Output: For best model list of,
#   accuracy 
#   mtry 
#   ntree
single.randomForest.predictor <- function(X, y, mtry=NA, ntree=c(500,1000,1500,2000)){
  # Intersect IDs in labels and in the feature data
  ids <- intersect(names(y),rownames(X))
  X   <- X[ids,]
  y   <- as.factor(y[ids])

  # set up mtry if not specified by function call
  if(is.na(mtry)) { mtry <- seq(max(floor(ncol(X)/3), 1)) } # TODO suppress warnings - don't care if mtry was specified as a list

  ## Grid sweep parameters, make lots of models
  # Results matrix
  params <- expand.grid(mtry, ntree)
  colnames(params) <- c('mtry','ntree')

  # Iterate over each parameter set
  # Train models
  randomForest.model <- apply(params, 1, function(x) { randomForest(X,y,mtry=x[1],ntry=x[2]) } )

  ## Calculate acc for each parameter set
  # acc = (TP + TN) / (P + N)
  params$acc <- sapply(randomForest.model, function(x) {signif(sum(diag(table(x$predicted,y)))/sum(table(x$predicted,y)), digits=3)})

  #print(params)

  # Select best model, return that result
  params.best <- randomForest.model[[max(params$acc, index=T)]]
  return( list( accuracy=params$acc[max(params$acc, index=T)], mtry=params.best$mtry, ntree=params.best$ntree ) )
}

### single support vector machine predictor #######################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
#TODO
