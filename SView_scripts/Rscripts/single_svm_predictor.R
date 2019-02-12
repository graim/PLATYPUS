### single_svm_predictor.R ########################################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
single.svm.predictor <- function(
  sample.feature.matrix,
  drug.response,
  cost = c(0.1,1,10,100),
  output.value = 'best',
  best.parameters = NULL,
  is.parallel = FALSE,
  num.cores = 20
  ) {

  # set parallel background if parallel is set to true
  if (is.parallel) {
    cl <- makeCluster(num.cores);
    registerDoParallel(cl, cores = num.cores);
  }

  # check that sample feature matrix is not mising any values
  if (any(is.na(sample.feature.matrix))) {
    stop("NAs are present in sample feature matrix. Please ensure all features have values ...");
  }

  # keep only cell lines that are present in both sample feature matrix and drug response
  sample.feature.matrix   <- sample.feature.matrix[rownames(sample.feature.matrix) %in% rownames(drug.response),];
  drug.response       <- drug.response[rownames(drug.response) %in% rownames(sample.feature.matrix),];

  # order both matrices to ensure rows are the same
  sample.feature.matrix   <- sample.feature.matrix[order(rownames(sample.feature.matrix)),];
  drug.response       <- drug.response[order(rownames(drug.response)),];

  # prepare training scheme
  control <- trainControl(method = 'repeatedcv', number = 10);

  # predict drug sensitivity in each cell line based on single group of features
  svm.results <- list();
  best.parameters.test <- list();
  for (compound in colnames(drug.response)) {
    # keep only cells that have sensitive or non-sensitive response
    compound.drug.response       <- drug.response[!is.na(drug.response[,compound]),];
    compound.sample.feature.matrix   <- sample.feature.matrix[!is.na(drug.response[,compound]),];

    # if best parameters specified using cost from file
    if (!is.null(best.parameters)) {
      cost <- best.parameters[compound,'C'];
    }

    # find number of overlapping cell lines
    num.cell.lines <- nrow(compound.drug.response);

    # factor drug response 
    compound.drug.response[,compound] <- factor(compound.drug.response[,compound]);

    # train svm model
    svm.model <- train(
      x = compound.sample.feature.matrix,
      y = compound.drug.response[,compound],
      method = 'svmRadialCost',
      preProcess = 'center',
      metric = 'Accuracy',
      tuneGrid = expand.grid(.C=cost),
      trControl = control
      );

    # add results to list
    # if more than once cost parameter gives same accuracy then use the lowest cost parameter (aka the first output)
    svm.results[['all.results']][[compound]] <- svm.model;
    best.parameters.test[[compound]] <- cbind(
      svm.model$results[which(svm.model$results$Accuracy == max(svm.model$results$Accuracy))[1],],
      num.cell.lines
      );
    colnames(best.parameters.test[[compound]]) <- c("C","Accuracy","Kappa","SDAccuracy","SDKappa","num.cell.lines");
    }

  # convert best parameters to data frame
  best.parameters <- do.call('rbind', best.parameters.test);

  if (output.value == 'best') {
      # if output value is set to best, return only the lowest oob error with its associated parameters
      return(best.parameters);

    } else if (output.value == 'all') {

      # add best summary 
      svm.results$best.parameters <- best.parameters;

      # return all results
      return(svm.results);
      }
    }

