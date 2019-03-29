

## Set up test suite
library(devtools)
devtools::use_testthat()

## Source the functions
source('../R/single.predictor.R')

## Create dummy data
X <- matrix(rnorm(1000), nrow=100)
y <- c(rep('A',50), rep('B',50))
names(y) <- paste0('sample',1:100)
rownames(X) <- names(y)

## Test each function - TODO make this a full test suite and move to package test folder?
single.elasticNet.predictor(X,y, iterations=2, nfolds=3)
single.randomForest.predictor(X,y)
