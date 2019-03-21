## The standard platypus function
##
## Created: Dec 2014, Kiley Graim
##    Updated: September 2015, Verena Friedl
##    Updated: March 2016, Kiley Graim
##    Updated: June 2016, Kiley Graim
##    Updated: Feb 2017, Kiley Graim
##    Updated: Jan 2018, Kiley Graim
##    Last Updated: Aug 2018, Kiley Graim 
## Check github for dates of latest updates

# pseudo-code:
#     for (v in views)
#       load v
#     load labels
#
#     while( not converged on labels OR no changes anymore )
#       for( v in views )
#         train v on labelled data
#         predict unlabelled data
#       check for instances where predictions agree -> new labels
#       add newly found labels to labelled data
#  





################################################################################
###  Main Function  #########################################################
################################################################################

# call platypus.R
# first pass the filepath for the labs file and the column name or number of the class (if not given, first column is default)
# pass the filepath of a parameter-file for each view
# [OPTIONS] <filename_labs> [-c <class_col>] for each view(<filename_viewfile>) 
# OPTIONS:
# -i <maximal number of iterations for each platypus run, eg. 100, default=100>
# -m <majority threshold in percent, eg. 75, default=100>
# -u flag for updating the accuracies of the single views in each iteration, default=FALSE
# -e flag for expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
# -b <class_name> flag for excluding cell lines that fall into class 'class_name' for the binary drug response definition, default='intermediate'


# fn.labs: labels file
# fn.views: list of parameter files
platypus <- function(fn.labs, fn.views, ignore.label='intermediate', i=100, m=100, u=FALSE, e=FALSE,updating=FALSE,expanded.output=FALSE) {

  ## Debug flag can be manually activated, for testing purposes 
  #flag.debug <- TRUE
  flag.debug <- FALSE 
  if(flag.debug) { print('Debug is on');flush.console() }

  ## Set more readable names
  majority.threshold.percent<-m
  no.iterations<-i
  classcol.labs <- 1 # The column containing the output labels. Often is the first column, not including rownames
  
  ## Create view from each parameter file and store in a list of views 
  view.list <- lapply(fn.views, load.parameterfile )

  ## Load the data for labs and each view
  labs <- load.label.data(fn.labs,classcol.labs)
  if(flag.debug) {print(table(labs[,classcol.labs]));flush.console()}

  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)
  view.list <- lapply(view.list, load.data)
  if(flag.debug) { print(lapply(view.list, function(x){length(intersect(labs,rownames(x)))} ));flush.console()  } 

  ## Take all IDs for each data type
  all.ids <- unique( unlist(lapply(view.list, function(x) {rownames(x$data.matrix)} )) )
  labs <- labs[rownames(labs) %in% all.ids,,drop=FALSE]  # labs without data in any view can not be used for training and/or prediction

  ## Sort out labelled and unlabelled IDs
  # filter all ids marked with "ignore.label"
  # if this method is called by cv.platypus, the hold-out fold is marked with 'testing' - has to be excluded from ids.labelled and ids.unlabelled!
  ids.labelled <- rownames(labs[!(is.na(labs[,classcol.labs])) & labs[,classcol.labs] != 'testing' & labs[,classcol.labs] != ignore.label,,drop=FALSE])
  labels <- labs[ids.labelled,classcol.labs,drop=F]
  labels[,classcol.labs] <- factor(labels[,classcol.labs])
  known.labels <- labels[,classcol.labs,drop=F]
     # important to take all IDs, which are in the feature data, not just the ones which have 'NA' as class label in labs
  ids.unlabelled <- all.ids[!(all.ids %in% ids.labelled) 
                            & !(all.ids %in% rownames(labs[labs[,classcol.labs] == 'testing',,drop=F])) 
                            & !(all.ids %in% rownames(labs[labs[,classcol.labs] == ignore.label,,drop=F]))]
  
  # summary of newly found labels
  new.labels <- c()
  
  # number of views, which have to agree on a prediction, to take it into the training data
  view.list <- normalize.accuracies(view.list)

  sum.acc <- 0
  for(view.i in seq(length(view.list))){
    sum.acc = sum.acc + view.list[[view.i]]$acc.norm
  }
    
  weighting.threshold <- majority.threshold.percent*sum.acc/100
  #starting values for all views agree
  weighting.threshold.upper <- sum.acc
  weighting.threshold.lower <- 0
    
  # counter
  iterations.woChange <- 0
  
  # collect information about each iteration for expanded output
  if(expanded.output){
    collect.iteration.lists <- list()
    
    unlabelled.matrix <- matrix(data=NA,nrow=length(ids.unlabelled),ncol=no.iterations,dimnames=list(ids.unlabelled))
    unlabelled.matrices.list <- list()
    for (view.i in seq(length(view.list))) {
      unlabelled.matrices.list[[view.i]] <- matrix(data=NA,nrow=length(ids.unlabelled),ncol=no.iterations,dimnames=list(ids.unlabelled))
    }
  }
  
  ## Repeat until all labels are 'known' OR we hit the iterations limit OR no new labels are added anymore
  for (i in seq(no.iterations)) {
	print(paste('Iteration',i))
    
    ## Quit if there aren't any unlabelled IDs left
    if (length(ids.unlabelled) <= 0) { 
      print('No new labels to learn, stopping label learning.')
      break 
    }

    ## Collect starting information for the iteration
    if(expanded.output){
      iteration.list <- list(iteration=i,weighting.threshold=weighting.threshold
                             ,no.ids.labelled=dim(labels)[[1]],no.ids.unlabelled=length(ids.unlabelled))
    }
    
    ## Train each view
    view.list <- lapply(view.list, function(x) { view.train(labels,x) } )
    
    ## Test on the unknown labels
    predictions <- matrix(data=NA, nrow=length(ids.unlabelled), ncol=length(view.list),dimnames=list(ids.unlabelled, seq(length(view.list))))
    for(view.i in seq(length(view.list))){
      ids <- intersect(ids.unlabelled,rownames(view.list[[view.i]]$data.matrix))
      predictions[ids,view.i] <- view.predict(ids.unlabelled,view.list[[view.i]]) 
    }


    ## Add unknown labels to known, where view agreement meets requirements
    if(updating){
        ## update accuracies
        view.list <- update.accuracies(view.list, known.labels)
        ## update weighting.threshold
        sum.acc <- 0
        for(view.i in seq(length(view.list))){
          sum.acc = sum.acc + view.list[[view.i]]$acc.norm
        }
        weighting.threshold <- majority.threshold.percent*sum.acc/100
    }
      
    new.labelled.list <- get.new.labels.majorityWeighted(predictions,view.list,unique.labels)
    
    # Quit, if the maximal prediction value reached is lower than the given threshold
    if(new.labelled.list$weighting.threshold.upper < weighting.threshold){
      print('No longer learning new labels, stopping label learning.')
      break
    } else{
      new.labelled <- new.labelled.list$new.labelled
      weighting.threshold.upper <- new.labelled.list$weighting.threshold.upper
      weighting.threshold.lower <- new.labelled.list$weighting.threshold.lower
    }
    
    colnames(new.labelled) <- colnames(labels)
    
    ## Add new labels
    labels <- rbind(labels,new.labelled)
    new.labels <- rbind(new.labels, new.labelled)
    if(flag.debug) {print( paste('Number labeled samples', length(ids.labelled)) ) }

    ## TODO: If only 1 class in labels list, quit with a warning UNTESTED
    if(now(table(labels))!=2) { 
      print('No longer 2 classes, stopping iterations.')
      break
    }

    ## Remove newly found labels from unlabelled IDs
    ids.unlabelled <- ids.unlabelled[!(ids.unlabelled %in% rownames(labels))]
    print( paste('Remaining unlabeled IDs:', length(ids.unlabelled)) )


    ## Collect information about the iterations for expanded output
    if(expanded.output){
      iteration.list$weighting.threshold.upper <- weighting.threshold.upper
      iteration.list$weighting.threshold.lower <- weighting.threshold.lower
      iteration.list$no.new.labelled <- length(new.labelled)
      iteration.list$no.ids.left.unlabelled <- length(ids.unlabelled)
      iteration.list$view.list <- view.list
      collect.iteration.lists[[i]] <- iteration.list
      
      unlabelled.matrix[rownames(new.labels),i] <- new.labels
      for (view.i in seq(length(view.list))) {
        unlabelled.matrices.list[[view.i]][rownames(predictions),i] <- predictions[,view.i]
      }
    }
  } # end for no.iterations

  ## Collect information for the returned result
  final.result.list <- list(final.views=view.list,labels.complete=labels,labels.new=new.labels,unlabelled.ids=ids.unlabelled
                            ,weighting.threshold=weighting.threshold,no.iterations=no.iterations,expanded.output=expanded.output)

  if(expanded.output){
    final.result.list$iteration.information <- collect.iteration.lists
    final.result.list$labelling.matrix <- unlabelled.matrix
    final.result.list$labelling.matrices.views <- unlabelled.matrices.list
  }
  return(final.result.list)
}
