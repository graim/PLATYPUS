## generate_feature_data_summary.R
## Created: July 2015, Katie Houlahan
##      Updated: March 2016, Kiley Graim
## Last Updated: Jan 2018, Kiley Graim
## 
## summarizes gene data by feature indicated in inputted file 
## inputted file must be of form feature by gene with binary input indicating if gene is associated with feature
## can calculate median or variance as specified in value argument 

generate.feature.data.summary <- function(gene.sample.data, feature.gene.matrix,type = 'expression', value = 'median',num.gene.threshold = 5) {
  if(value=='kurtosis') { require('e1071') } # For the kurtosis option

  # keep only genes that are found in both data and feature/gene matrix 
  feature.gene.matrix <- feature.gene.matrix[,colnames(feature.gene.matrix) %in% colnames(gene.sample.data)];

  # remove any feature terms that aren't associated with any of the remaining genes
  all.zero <- apply( feature.gene.matrix, 1, function(x) { all(x == 0) } );
  feature.gene.matrix <- feature.gene.matrix[!all.zero,];

  # ensure the number of genes associated with each feature exceeds threshold 
  # remove any features without enough genes supporting it 
  feature.gene.matrix <- feature.gene.matrix[rowSums(feature.gene.matrix) >= num.gene.threshold,];

  # set output matrix
  sample.feature.data <- matrix(nrow = nrow(gene.sample.data), ncol = nrow(feature.gene.matrix));
  rownames(sample.feature.data) <- rownames(gene.sample.data);
  colnames(sample.feature.data) <- rownames(feature.gene.matrix);

  # populate matrix with median gene data corresponding to each feature 
  # if no genes corresponding to that feature are present in data file return NA
  for (feature in 1:nrow(feature.gene.matrix)) {

    # find genes corresponding to feature
    feature.genes <- colnames(feature.gene.matrix)[which(feature.gene.matrix[feature,] == 1)];

    # find data of genes associated with each feature for each cell line
    if (type %in% c('expression', 'CNV')) {
	# TODO: max() is an actual function- should we rename this? Also, variance->var mapping should be added for convenience
      if(value == 'max' ) {
        # find max value of median, variance and kurtosis
        sample.feature.data[,feature] <- apply(gene.sample.data,1,
            function(sample) {
              tmp.median     <- median(sample[colnames(gene.sample.data) %in% feature.genes]);
              tmp.var     <- var(sample[colnames(gene.sample.data) %in% feature.genes]);
              tmp.kurtosis  <- kurtosis(sample[colnames(gene.sample.data) %in% feature.genes]);
              max(tmp.median, tmp.var, tmp.kurtosis);
              }
            );
      } else {
        # Use whatever statistic the user prefers (mean, median, variance, etc)
        try(sample.feature.data[,feature] <- apply(gene.sample.data, 1, function(x) {do.call(value, list(unlist(x[colnames(gene.sample.data) %in% feature.genes])))} ) )
      }
    } else if (type == 'mutation') {
      # find ratio of mutations in pathway compared to mutations in cell line 
      sample.feature.data[,feature] <- apply(gene.sample.data,1,
          function(sample) {
            sample.genes <- length(feature.genes);
            feature.sample.mutations <- sum(sample[colnames(gene.sample.data) %in% feature.genes]);
            feature.sample.mutations/sample.genes;
            }
          );

    } else {
      stop("Please specify valid type. Options are 'expression', 'CNV' or 'mutation'...")
    }
  } # end for
  return(sample.feature.data);
} # end fxn
