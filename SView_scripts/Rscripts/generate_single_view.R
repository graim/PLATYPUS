#!/usr/bin/env Rscript

##### #!/inside/home/khoulaha/bin/Rscript

# load libraries
library(doParallel);
library(getopt);
library(caret);

# general parameters
date <- Sys.Date();
rscp.loc <- '/Users/kgraim/Documents/MVL/PLATYPUS/SView_scripts/Rscripts/' # TODO: this should be automagically set

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
  # required parameters
  'sample.feature.matrix'  ,  'f',   1, "character",
  'drug.response'          ,  'd',   1, "character",
  'model'                  ,  'm',   1, "character",
  # optional parameters
  'output.dir'             ,  'o',   2, "character",
  'mtry'                   ,  't',  2, "character",
  'ntree'                  ,  'n',   2, "character",
  'alpha'                  ,  'a',  2, "character",
  'cost'                   ,  'c',  2, "character",
  'iterations'             ,  'i',   2, "integer",
  'num.cores'              ,   'u',  2, "integer",
  'parallel'               ,  'p',  2, "logical",
  'best.parameters'        ,  'b',  2, "character",
  'rerun'                  ,  'r',  2, "integer",
  'importance'             ,  'e',  2, "logical",
  "flag"                   ,  'l',  2, "character"
  ), byrow = TRUE, ncol = 4);
opt = getopt(spec);

# set defaults
if (is.null(opt$mtry    )) { opt$mtry = "seq(
  round(sqrt(ncol(sample.feature.matrix))) - round(0.05*ncol(sample.feature.matrix)),
  round(sqrt(ncol(sample.feature.matrix))) + round(0.05*ncol(sample.feature.matrix)), 
  by = 1
  )"
  }
if (is.null(opt$ntree    )) {opt$ntree = "c(500,1000,1500,2000)"}
if (is.null(opt$alpha    )) {opt$alpha = "seq(0.2,0.8,0.2)"}
if (is.null(opt$cost    )) {opt$cost = "c(0.1,1,10,100)"}
if (is.null(opt$iterations  )) {opt$iterations = 100}
if (is.null(opt$num.cores  )) {opt$num.cores = 3}
if (is.null(opt$parallel  )) {opt$parallel = TRUE}
if (is.null(opt$rerun    )) {opt$rerun = 1}
if (is.null(opt$importance  )) {opt$importance = FALSE}

####################################################################################################
# function to calculate median value across all iterations 
calculate.median.value <- function(results, value = 'accuracy') {
  # unlist data frames
  listVec <- lapply(results, c, recursive=TRUE);
  # generate tmp matrix
  tmp.matrix <- do.call(cbind, listVec);

  # calculate median of value specified over all iterations for all drugs
  tmp <- apply(
    tmp.matrix[grep(value, rownames(tmp.matrix)),],
    1,
    median
    );
  return(tmp);
  }

# read in sample feature matrix
sample.feature.matrix <- read.delim(
  opt$sample.feature.matrix,
  sep = '\t',
  header = TRUE
  );

# read in drug response
drug.response <- read.delim(
  opt$drug.response,
  sep = '\t', 
  header = TRUE
  );

if (!is.null(opt$best.parameters)) {
  # read in best parameters file
  best.parameters <- read.delim(
    opt$best.parameters,
    sep = '\t',
    header = TRUE
    );
} else {
  # if no best parameters file specified, set as NULL
  best.parameters <- NULL;
}

if (opt$model == 'randomForest') {

  # evaluate parameter arguments
  mtry <- eval(parse(text = opt$mtry));
  ntree <- eval(parse(text = opt$ntree));

  # ensure that mtry is not below 2
  mtry <- mtry[mtry >= 2];

  # source code
  source(paste0(rscp.loc,"single_randomForest_predictor.R"));
  if (opt$rerun == 1) {
    # run random Forest with specified mtry and ntree
    sample.feature.matrix.results <- single.randomForest.predictor(
      sample.feature.matrix = sample.feature.matrix,
      drug.response = drug.response,
      mtry = mtry,
      ntree = ntree,
      best.parameters = best.parameters,
      output.value = 'all',
      parallel = opt$parallel,
      num.cores = opt$num.cores
      );
    } else {
      tmp.sample.feature.matrix.results <- list();
      # re-run random forest multiple times and take median accuracy value over all iterations
      for (i in 1:opt$rerun) {
        tmp.sample.feature.matrix.results[[i]] <- single.randomForest.predictor(
          sample.feature.matrix = sample.feature.matrix,
          drug.response = drug.response,
          mtry = mtry,
          ntree = ntree,
          best.parameters = best.parameters,
          output.value = 'best',
          parallel = opt$parallel,
          num.cores = opt$num.cores
          );
        }

      # calculate median oob error, sensitivity, specificity and accuracy 
      sample.feature.matrix.results <- data.frame(
        mtry = tmp.sample.feature.matrix.results[[1]][,'mtry'],
        ntree = tmp.sample.feature.matrix.results[[1]][,'ntree'],
        oob.err = calculate.median.value(tmp.sample.feature.matrix.results, value = 'oob.err'),
        sensitivity = calculate.median.value(tmp.sample.feature.matrix.results, value = 'sensitivity'),
        specificity = calculate.median.value(tmp.sample.feature.matrix.results, value = 'specificity'),
        accuracy = calculate.median.value(tmp.sample.feature.matrix.results, value = 'accuracy'),
        num.cell.lines = tmp.sample.feature.matrix.results[[1]][,'num.cell.lines']
        );
      rownames(sample.feature.matrix.results) <- rownames(tmp.sample.feature.matrix.results[[1]]);
    }

} else if (opt$model == 'elasticNet') {

  # evaluate parameter arguments
  alpha <- eval(parse(text = opt$alpha));

  # source code
  source(paste0(rscp.loc,"single_elasticNet_predictor.R"));
  if (opt$rerun == 1) {
    # run elastic with specified alpha
    sample.feature.matrix.results <- single.elasticNet.predictor(
      sample.feature.matrix = sample.feature.matrix,
      drug.response = drug.response,
      alpha = alpha,
      best.parameters = best.parameters,
      iterations = opt$iterations,
      output.value = 'all',
      parallel = opt$parallel,
      num.cores = opt$num.cores
      );
    } else {
      tmp.sample.feature.matrix.results <- list();
      # re-run elastic net multiple times and take median error, std.error and lambda
      for (i in 1:opt$rerun) {
        tmp.sample.feature.matrix.results[[i]] <- single.elasticNet.predictor(
          sample.feature.matrix = sample.feature.matrix,
          drug.response = drug.response,
          alpha = alpha,
          best.parameters = best.parameters,
          iterations = opt$iterations,
          output.value = 'best',
          parallel = opt$parallel,
          num.cores = opt$num.cores
          );
        }
      # calculate median error, std error and lambda over all re runs
      sample.feature.matrix.results <- data.frame(
        alpha = tmp.sample.feature.matrix.results[[1]][,'alpha'],
        error = calculate.median.value(tmp.sample.feature.matrix.results, value = '^error'),
        std.error = calculate.median.value(tmp.sample.feature.matrix.results, value = 'std.error'),
        lambda = calculate.median.value(tmp.sample.feature.matrix.results, value = 'lambda'),
        accuracy = calculate.median.value(tmp.sample.feature.matrix.results, value = 'accuracy'),
        num.cell.lines = tmp.sample.feature.matrix.results[[1]][,'num.cell.lines']
        );
      rownames(sample.feature.matrix.results) <- rownames(tmp.sample.feature.matrix.results[[1]]);
      }

} else if (opt$model == 'svm') {

  # evaluate parameter arguments 
  cost <- eval(parse(text = opt$cost));

  # source code
  source(paste0(rscp.loc,"single_svm_predictor.R"));
  if (opt$rerun == 1) {
    # run svm with specified cost
    sample.feature.matrix.results <- single.svm.predictor(
      sample.feature.matrix = sample.feature.matrix,
      drug.response = drug.response,
      cost = cost,
      best.parameters = best.parameters,
      output.value = 'all',
      parallel = opt$parallel,
      num.cores = opt$num.cores
      );
    } else {
      tmp.sample.feature.matrix.results <- list();
      # re-run elastic net multiple times and take median accuracy, kappa, sd accuracy and sd kappa
      for (i in 1:opt$rerun) {
        # run svm with specified cost
        tmp.sample.feature.matrix.results[[i]] <- single.svm.predictor(
          sample.feature.matrix = sample.feature.matrix,
          drug.response = drug.response,
          cost = cost,
          best.parameters = best.parameters,
          output.value = 'best',
          parallel = opt$parallel,
          num.cores = opt$num.cores
          );
        }
      # create output matrix with median accuracy, kappa, sd accuracy and sd kappa over all re runs
      sample.feature.matrix.results <- data.frame(
        cost = tmp.sample.feature.matrix.results[[1]][,'C'],
        Accuracy = calculate.median.value(tmp.sample.feature.matrix.results, value = '^Accuracy'),
        Kappa = calculate.median.value(tmp.sample.feature.matrix.results, value = '^Kappa'),
        AccuracySD = calculate.median.value(tmp.sample.feature.matrix.results, value = 'SDAccuracy'),
        KappaSD = calculate.median.value(tmp.sample.feature.matrix.results, value = 'SDKappa'),
        num.cell.lines = tmp.sample.feature.matrix.results[[1]][,'num.cell.lines']
        );
      rownames(sample.feature.matrix.results) <- rownames(tmp.sample.feature.matrix.results[[1]]);
    }
  
} else {
  stop("Please specify valid model. Options are 'randomForest', 'elasticNet' or 'svm'...")
}

if (opt$rerun == 1) {
  # write best parameters results to file
  write.table(
    x = sample.feature.matrix.results$best.parameters,
    file = paste0(
      opt$output.dir,
      date,
      "_SampleFeatureMatrix_",
      opt$model,
      "_BestParameters",
      opt$flag,
      ".txt"
      ),
    sep = '\t',
    quote = FALSE
    );
  } else {
    # write best parameters results to file
    write.table(
      x = sample.feature.matrix.results,
      file = paste0(
        opt$output.dir,
        date,
        "_SampleFeatureMatrix_",
        opt$model,
        "_BestParametersReRun",
        opt$flag,
        ".txt"
        ),
      sep = '\t',
      quote = FALSE
      );
  }

if (opt$importance) {
  # write coefficients or importance to file
  if (opt$model == 'elasticNet') {
    # coefficients specified as sparse matrix therefore write out as sparse matrix
    write.table(
      sample.feature.matrix.results$coefficients,
      file = paste0(
        opt$output.dir,
        date,
        "_Coefficients_",
        opt$model,
        opt$flag,
        ".txt"
        ),
      sep = '\t',
      quote = FALSE
      );
  } else if (opt$model == 'randomForest') {
    # write importance values to file
    write.table(
      x = sample.feature.matrix.results$importance,
      file = paste0(
        opt$output.dir,
        date,
        "_Importance_",
        opt$model,
        opt$flag,
        ".txt"
        ),
      sep = '\t',
      quote = FALSE
      );
    }
  }
  
