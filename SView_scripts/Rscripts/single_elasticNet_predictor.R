### single_elasticNet_predictor.R ################################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
###################################################################################################
single.elasticNet.predictor <- function(
  sample.feature.matrix,
  drug.response,
  alpha = seq(0.2,0.8,0.2),
  iterations = 100,
  best.parameters = NULL,
  output.value = 'best',
  parallel = FALSE,
  num.cores = 20,
  coef.output.dir = NULL
  ) {

  # set parallel background if parallel is set to true
  if (parallel) {
    cl <- makeCluster(num.cores);
    registerDoParallel(cl, cores = num.cores);
  }

  # create function to run using foreach function 
  calculate.optimal.elasticNet.parameters <- function(alpha.test,iteration,compound,foldid) {

    # run elastic net cross validation
    elastic.net.results <- cv.glmnet(
      x = as.matrix(compound.sample.feature.matrix),
      y = compound.drug.response[,compound], 
      family = 'binomial', 
      alpha = alpha.test, 
      foldid = foldid,
      type.measure = "mse"
      );

    # extract coefficients
    extracted.coefficients <- coef(elastic.net.results);

    # combine results into list
    combined.elastic.net.results <- list()
    combined.elastic.net.results[['glmnet.fit']] <- elastic.net.results;
    combined.elastic.net.results[['coefficients']] <- extracted.coefficients;

    # return results
    return(combined.elastic.net.results);
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

  # predict drug sensitivity in each cell line based on single group of features
  elasticNet.results <- list();
  for (compound in colnames(drug.response)) {

    # keep only cells that have sensitive or non-sensitive response
    compound.drug.response       <- drug.response[!is.na(drug.response[,compound]),];
    compound.sample.feature.matrix   <- sample.feature.matrix[!is.na(drug.response[,compound]),];

    # ensure that cell lines in drug response match the order of those in sample feature matrix
    if (!all(rownames(compound.drug.response) == rownames(compound.sample.feature.matrix))) {
      stop("Cell lines do not match between drug response and sample feature matrix. Please fix ...");
    }

    # find number of cell lines overlapping
    num.overlapping.cell.lines <- nrow(compound.drug.response);

    # factor drug response 
    compound.drug.response[,compound] <- factor(compound.drug.response[,compound]);

    # if already determined best parameters can use the alpha instead of re-testing alphas
    if (is.null(best.parameters)) {
      # create vector of all parameters to test
      # in this case, combine each alpha with an iteration number 
      parameters.to.test <- expand.grid(alpha,1:iterations);
      } else {
        # create vector of best parameters to test
        parameters.to.test <- expand.grid(best.parameters[compound,'alpha'],1:iterations);
        # set alpha as alpha in best parameter file
        alpha <- best.parameters[compound,'alpha'];
      }

    # assign a fold to each cell line
    # find the size of each fold 
    number.in.fold   <- floor(nrow(compound.drug.response)/10);
    # set vector to store fold number associated with each cell line
    cell.line.fold   <- rep(NA, nrow(compound.drug.response));
    # create vector to keep track of cell lines not assigned a fold yet
    cell.lines.left <- 1:nrow(compound.drug.response); 
    for (i in 1:9) {
      # sample cell lines for each fold
      rows.to.choose <- sample(cell.lines.left, size = number.in.fold);
      # assigned sampled cell lines fold number 
      cell.line.fold[rows.to.choose] <- i;
      # remove cell lines assigned a fold from vector of remaining cell lines
      cell.lines.left <- cell.lines.left[!cell.lines.left %in% rows.to.choose];
      }
    # assign remaining cell lines fold 10 (may be less than the remaining folds)
    cell.line.fold[is.na(cell.line.fold)] <- 10;

    if (parallel) {
      # test all combinations of parameters running each elastic net in parallel if parallel argument set as true
      elasticNet.model <- foreach(alpha.test = parameters.to.test[,1], iteration = parameters.to.test[,2], .packages = 'glmnet') %dopar%
        calculate.optimal.elasticNet.parameters(alpha.test = alpha.test, iteration = iteration, compound = compound, foldid = cell.line.fold);
      } else {
        # if parallelize argument not set as true
        elasticNet.model <- foreach(alpha.test = parameters.to.test[,1], iteration = parameters.to.test[,2], .packages = 'glmnet') %do%
          calculate.optimal.elasticNet.parameters(alpha.test = alpha.test, iteration = iteration, compound = compound, foldid = cell.line.fold);
      }

    # name each iteration indicating the iteration number and alpha tested
    names(elasticNet.model) <- apply(parameters.to.test, 1, paste, collapse = '_');

    # create matrix of parameters
    elastic.net.parameters <- lapply(
      elasticNet.model,
      function(elastic.net.results) {
        data.frame(
           error = elastic.net.results$glmnet.fit$cvm[elastic.net.results$glmnet.fit$lambda == elastic.net.results$glmnet.fit$lambda.min],
           std.error = elastic.net.results$glmnet.fit$cvsd[elastic.net.results$glmnet.fit$lambda == elastic.net.results$glmnet.fit$lambda.min],
           lambda = elastic.net.results$glmnet.fit$lambda.min
           );
        }
      );
    all.iterations.data <- do.call('rbind', elastic.net.parameters);
    # add alpha tested to matrix of parameters
    all.iterations.data$alpha <- rep(alpha, each = iterations);

    # median error and lambda for each alpha across each iteration
    median.values.per.alpha <- by(
      all.iterations.data,
      all.iterations.data$alpha,
      function(x) {
        # create new data frame with median values
        tmp.df <- data.frame(
          alpha = x$alpha[1],
          error = median(x$error),
          std.error = median(x$std.error),
          lambda = median(x$lambda),
          accuracy = 1 - median(x$error),
          num.cell.lines = num.overlapping.cell.lines
          );
        return(tmp.df);
        }
      );

    # add parameter results to elasticNet.results list by compound
    elasticNet.results[[compound]][['all.iterations']]   <- elasticNet.model;
    elasticNet.results[[compound]][['median.values']]   <- do.call('rbind', median.values.per.alpha);
    }

    # find best parameters, in this case alpha, for each compound
    best.parameters <- lapply(
      elasticNet.results,
      function(compound.results) {
        compound.results$median.values[compound.results$median.values$error == min(compound.results$median.values$error),][1,];
        }
      );
    # convert list to data frame
    best.parameters <- do.call('rbind', best.parameters);

    if (output.value == 'best') {
      # if output value is set to best, return only the lowest oob error with its associated parameters
      return(best.parameters);

    } else if (output.value == 'all') {
      # if output value is set to all return both the summary of best parameters and the output of each test
      all.elasticNet.results <- list();

      # add best summary 
      all.elasticNet.results$best.parameters <- best.parameters;

      # all elastic net results
      all.elasticNet.results$all.results <- elasticNet.results;

      ## find median coefficient for each compound and output data frame of non-zero features
      # find coefficients of features
      feature.coefficients <- list();
      for (compound in rownames(best.parameters)) {
        # extract coefficients of each run with best alpha
        with.coefficients <- lapply(
          elasticNet.results[[compound]][['all.iterations']][grep(best.parameters[compound,'alpha'], names(elasticNet.results[[compound]][['all.iterations']]))],
          '[[',
          "coefficients"
          );
        # merge into one data frame
        with.coefficients <- as.matrix(do.call(cbind, with.coefficients));

        # find median coefficient 
        feature.coefficients[[compound]] <- data.frame(apply(with.coefficients, 1, median));
        }
      # merge into one data frame
      feature.coefficients <- do.call(cbind, feature.coefficients);
      colnames(feature.coefficients) <- rownames(best.parameters);
      # keep only features with non zero coefficients
      keep <- apply(feature.coefficients, 1, function(x) { all(x == 0)});
      feature.coefficients <- feature.coefficients[!keep,];

      # add coefficients to all elasticNet results
      all.elasticNet.results$coefficients <- feature.coefficients;

      return(all.elasticNet.results);
      }
    }
