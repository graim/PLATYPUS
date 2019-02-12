### single_randomForest_predictor.R ###############################################################
# generalized function to take in feature by sample expression matrix and predict drug sensitivity
# requires to R packages: "randomForest" and "foreach"
###################################################################################################
single.randomForest.predictor <- function(sample.feature.matrix, drug.response, mtry = seq( round(sqrt(ncol(sample.feature.matrix))) - round(0.05*ncol(sample.feature.matrix)), round(sqrt(ncol(sample.feature.matrix))) + round(0.05*ncol(sample.feature.matrix)), by = 1 ), ntree = c(500,1000,1500,2000), best.parameters = NULL, output.value = 'best', is.parallel = FALSE, num.cores = 20) {

  # set is.parallel background if is.parallel is set to true
  if (is.parallel) {
    cl <- makeCluster(num.cores);
    registerDoParallel(cl, cores = num.cores);
  }

  # create function to calculate scores
  calculate.scores <- function(predicted, response) {
    # create data frame of results
    tmp.evaluation <- matrix(nrow = 1, ncol = 7);
    colnames(tmp.evaluation) <- c(
      'true.positives',
      'false.positives',
      'true.negatives',
      'false.negatives',
      'sensitivity',
      'specificity',
      'balanced.accuracy'
      );

    # find tp, fp, tn and fn of training data only
    #print(table(predicted))
    labs.ids <- levels(predicted[,1])
    tmp.evaluation[1,'true.positives']   <- sum(predicted[,1] == labs.ids[[1]] & response[,1] == labs.ids[[1]]);
    tmp.evaluation[1,'false.positives']  <- sum(predicted[,1] == labs.ids[[1]] & response[,1] == labs.ids[[2]]);
    tmp.evaluation[1,'true.negatives']   <- sum(predicted[,1] == labs.ids[[2]] & response[,1] == labs.ids[[2]]);
    tmp.evaluation[1,'false.negatives']  <- sum(predicted[,1] == labs.ids[[2]] & response[,1] == labs.ids[[1]]);

    #tmp.evaluation[1,'true.positives']   <- sum(predicted[,1] == 'ASD' & response[,1] == 'ASD');
    #tmp.evaluation[1,'false.positives'] <- sum(predicted[,1] == 'ASD' & response[,1] == 'HEALTHY');
    #tmp.evaluation[1,'true.negatives']   <- sum(predicted[,1] == 'HEALTHY' & response[,1] == 'HEALTHY');
    #tmp.evaluation[1,'false.negatives'] <- sum(predicted[,1] == 'HEALTHY' & response[,1] == 'ASD');

    # find sensitivity, specificity and balanced accuracy of training data only 
    tmp.evaluation[1,'sensitivity']     <- tmp.evaluation[1,'true.positives']/(tmp.evaluation[1,'true.positives'] + tmp.evaluation[1,'false.negatives']);
    tmp.evaluation[1,'specificity']     <- tmp.evaluation[1,'true.negatives']/(tmp.evaluation[1,'true.negatives'] + tmp.evaluation[1,'false.positives']);
    tmp.evaluation[1,'balanced.accuracy']   <- (tmp.evaluation[1,'sensitivity'] + tmp.evaluation[1,'specificity'])/2;

    return(tmp.evaluation)
    } # end calculate.scores

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
  randomForest.results <- list();
  for (compound in colnames(drug.response)) {
    print(paste('Testing outcome label:',compound));flush.console()

    # keep only cells that have sensitive or HEALTHY response
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

    # for each fold run random forest 
    randomForest.model   <- list();
    parameter.evaluation   <- list();

    if (is.null(best.parameters)) {
      # create vector all combinations of mtry and ntree to test
      parameters.to.test   <- expand.grid(mtry, ntree);
    } else {
      # create vector of best parameters to test
      parameters.to.test <- best.parameters[compound,c('mtry', 'ntree')];
    }

    # create vector of names of parameters testing
    parameter.names <- apply(parameters.to.test,1,function(x) {paste("mtry",x[1],"ntree",x[2],sep = ".")});

    if (is.parallel) {
      # test all combinations of parameters running each randomForest in parallel if is.parallel argument set as true
      randomForest.model <- foreach(mtry.test = parameters.to.test[,1], ntree.test = parameters.to.test[,2], .packages = 'randomForest') %dopar%
        randomForest(
          x = compound.sample.feature.matrix,
          y = compound.drug.response[,compound],
          mtry = mtry.test,
          ntree = ntree.test,
          importance = TRUE
          );
        } else {
          # if is.parallelize argument not set as true 
          randomForest.model <- foreach(mtry.test = parameters.to.test[,1], ntree.test = parameters.to.test[,2], .packages = 'randomForest') %do%
            randomForest(
              x = compound.sample.feature.matrix,
              y = compound.drug.response[,compound],
              mtry = mtry.test,
              ntree = ntree.test,
              importance = TRUE
              );
        }

    # evaluate each combination of parameters by extracting the oob and calculating sensitive, specificity and balanced accuracy
    parameter.evaluation <- lapply(
      randomForest.model,
      function(parameter.model) {
        tmp.evaluation <- calculate.scores(
          predicted = data.frame(parameter.model$predicted),
          response = data.frame(compound.drug.response[,compound])
          );
        tmp.evaluation <- cbind(
          tmp.evaluation,
          mean(parameter.model$confusion[,3]),
          num.overlapping.cell.lines
          );
        colnames(tmp.evaluation) <- c(colnames(tmp.evaluation)[1:7],'oob.err','num.cell.lines');
        return(tmp.evaluation);
        }
      );

    # assign names of parameters tested to each list
    names(randomForest.model) <- parameter.names;
    names(parameter.evaluation) <- parameter.names;
      
    # store information into list
    randomForest.results[[compound]][['randomForest.model']] <- randomForest.model;
    randomForest.results[[compound]][['parameter.evaluation']] <- parameter.evaluation;
  }

  # for each compound find the best parameters
  best.parameters <- lapply(
    randomForest.results,
    function(compound.results) {
      # extract all oob errors
      oob.err <- unlist(
        lapply(
          compound.results$parameter.evaluation,
          function(x) {
            x[,'oob.err'];
            }
          )
        );
      # find min oob.err
      min.oob.err <- min(oob.err, na.rm=TRUE); # KILEY EDITED 20180131
      # find parameter name associated with min oob.err
      tmp.par.name <- substr(
        names(oob.err)[oob.err == min.oob.err],
        start = 1,
        stop = nchar(names(oob.err)[oob.err == min.oob.err])-8
        );
      # if multiple parameters give lowest error keep first one as this will likely require fewer trees
      if (length(tmp.par.name) > 1) {
        tmp.par.name <- tmp.par.name[1];
        }
      print("Finished another ..."); # TODO
      # find sensitivity, specificity and balanced accuracy of best parameters
      evaluation.best.parameters <- compound.results$parameter.evaluation[tmp.par.name];
      # split str ra
      split.tmp.par.name <- unlist(strsplit(tmp.par.name, split = "\\."));
      # set results matrix
      best.parameters.oob.err <- data.frame(
        mtry = split.tmp.par.name[2],
        ntree = split.tmp.par.name[4],
        oob.err = min.oob.err,
        sensitivity = evaluation.best.parameters[[1]][,'sensitivity'],
        specificity = evaluation.best.parameters[[1]][,'specificity'],
        accuracy = evaluation.best.parameters[[1]][,'balanced.accuracy'],
        num.cell.lines = compound.results$parameter.evaluation[[1]][1,'num.cell.lines']
        );
      return(best.parameters.oob.err);
      }
    );
  # convert list into data frame
  best.parameters <- do.call('rbind', best.parameters);

  if (output.value == 'best') {
    # if output value is set to best, return only the lowest oob error with its associated parameters
    return(best.parameters);

  } else if (output.value == 'all') {
    # if output value is set to all return both the summary of best parameters and the output of each test
    all.randomForest.results <- list();

    # add best summary 
    all.randomForest.results$best.parameters <- best.parameters;

    # all all randomForest reuslts
    all.randomForest.results$all.results <- randomForest.results;

    # extract importance gini score for each feature from run with best parameters
    tmp.importance <- list();
    for (compound in rownames(best.parameters)) {
      best.parameter.name <- paste("mtry", best.parameters[compound,'mtry'], 'ntree', best.parameters[compound,'ntree'], sep = '.');
      # extract importance measures for run with best parameters
      tmp.importance[[compound]] <- randomForest.results[[compound]][['randomForest.model']][[best.parameter.name]][['importance']][,'MeanDecreaseGini'];
    }
    # reformat in matrix for each drug
    importance.values <- do.call('cbind', tmp.importance);
    # add drug names as column names
    colnames(importance.values) <- rownames(best.parameters);

    # order by overall gene importance and add to all results list
    all.randomForest.results$importance <- importance.values[order(rowSums(importance.values), decreasing = TRUE),];

    return(all.randomForest.results);
    }
  }
