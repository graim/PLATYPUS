## The standard platypus function
##
## Created: Dec 2014, Kiley Graim
##    Updated: September 2015, Verena Friedl
##    Updated: March 2016, Kiley Graim
##    Last Updated: Jan 2018, Kiley Graim 

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
# -w flag for weighting the preditions by accuracy, default=FALSE
# -u flag for updating the accuracies of the single views in each iteration, default=FALSE
# -e flag for expanded output: returned result list contains a list of trained views after each iteration, default=FALSE
# -b <class_name> flag for excluding cell lines that fall into class 'class_name' for the binary drug response definition, default='intermediate'



platypus <- function(argv) {
  
  source(paste0(Sys.getenv("HOME"),'/MVL/scripts/platypus.basicFunctions.R'))
  
  ## Read in the command line arguments
  # set number of iterations
  no.iterations <- 100
  if("-i" %in% argv){
    no.iterations <- as.numeric(as.character(argv[which(argv == "-i")+1]))
    argv <- argv[-c(which(argv == "-i"),which(argv == "-i")+1)]
  }
  
  # set majority threshold
  majority.threshold.percent <- 100
  if("-m" %in% argv){
    majority.threshold.percent <- as.numeric(as.character(argv[which(argv == "-m")+1]))
    argv <- argv[-c(which(argv == "-m"),which(argv == "-m")+1)]
  }
  
  # set expanded output flag
  expanded.output = FALSE
  if("-e" %in% argv){
    expanded.output <- TRUE
    argv <- argv[-c(which(argv == "-e"))]
  }
  
  # set weighting flag
  weighting = FALSE
  if("-w" %in% argv){
    weighting <- TRUE
    argv <- argv[-c(which(argv == "-w"))]
  }
  
  # set updating flag
  updating = FALSE
  if("-u" %in% argv){
    updating <- TRUE
    argv <- argv[-c(which(argv == "-u"))]
  }
  
  # set class to ignore from training and predicting
  ignore.label <- "intermediate"
  if("-b" %in% argv){
    ignore.label <- as.character(argv[which(argv == "-b")+1])
    argv <- argv[-c(which(argv == "-b"),which(argv == "-b")+1)]
  }
  
  
  # filenames
  fn.labs <- argv[1]
  classcol.labs <- 1
  if ("-c" %in% argv) {
    col <- argv[which(argv == "-c")+1]
    if (is.na(as.numeric(as.character(col)))){
      classcol.labs <- col
    } else{
      classcol.labs <- as.numeric(as.character(col))
    }
    argv <- argv[-c(which(argv == "-c"),which(argv == "-c")+1)]
  }
  fn.views <- argv[2:length(argv)]

  ## Create view from each parameter file and store in a list of views 
  view.list <- list()
  for(view.i in 1:length(fn.views)){
    view.list[[view.i]] <- load.parameterfile(fn.views[view.i])
  }

  ## Load the data for labs and each view
  labs <- load.label.data(fn.labs,classcol.labs)
  
  # Get the the two labels
  unique.labels <- get.unique.labels(labs[,classcol.labs],ignore.label)
  
  for(view.i in 1:length(view.list)){
    view <- load.data(view.list[[view.i]])
    view.list[[view.i]] <- view
  }

  
  ## Take all IDs for each data type
  all.ids <- c()
  for(view in view.list){
    all.ids <- unique(c(all.ids,rownames(view$data.matrix))) 
  }  
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
  if(weighting){
    
    view.list <- normalize.accuracies(view.list)
    sum.acc <- 0
    for(view.i in 1:length(view.list)){
      sum.acc = sum.acc + view.list[[view.i]]$acc.norm
    }
    
    weighting.threshold <- majority.threshold.percent*sum.acc/100
    #starting values for all views agree
    weighting.threshold.upper <- sum.acc
    weighting.threshold.lower <- 0
    
  } else {
    majority <- length(view.list)
    majority.missingData <- length(view.list)
    majority.threshold <- ceiling(majority.threshold.percent*majority/100)
  }

  
  # counter
  iterations.woChange <- 0
  
  # collect information about each iteration for expanded output
  if(expanded.output){
    collect.iteration.lists <- list()
    
    unlabelled.matrix <- matrix(data=NA,nrow=length(ids.unlabelled),ncol=no.iterations,dimnames=list(ids.unlabelled))
    unlabelled.matrices.list <- list()
    for (view.i in 1:length(view.list)) {
      unlabelled.matrices.list[[view.i]] <- matrix(data=NA,nrow=length(ids.unlabelled),ncol=no.iterations,dimnames=list(ids.unlabelled))
    }
  }
  
  ## Repeat until all labels are 'known' OR we hit the iterations limit OR no new labels are added anymore
  for (i in 1:no.iterations) {
    
    ## Quit if there aren't any unlabelled IDs left
    if (length(ids.unlabelled) <= 0) { 
      break 
    }

    ## Collect starting information for the iteration
    if(expanded.output){
      if(weighting){
        iteration.list <- list(iteration=i,weighting.threshold=weighting.threshold
                               ,no.ids.labelled=dim(labels)[[1]],no.ids.unlabelled=length(ids.unlabelled))
      } else {
        iteration.list <- list(iteration=i,iterations.woChange=iterations.woChange,majority=majority,majority.missingData=majority.missingData,majority.threshold=majority.threshold
                               ,no.ids.labelled=dim(labels)[[1]],no.ids.unlabelled=length(ids.unlabelled))
      }
    }
    
    ## Train each view
    for (view.i in 1:length(view.list)) {
      view <- view.train(labels,view.list[[view.i]])
      view.list[[view.i]] <- view
    }
    
    ## Test on the unknown labels
    predictions <- matrix(data=NA, nrow=length(ids.unlabelled), ncol=length(view.list),dimnames=list(ids.unlabelled, 1:length(view.list)))
    for(view.i in 1:length(view.list)){
      ids <- intersect(ids.unlabelled,rownames(view.list[[view.i]]$data.matrix))
      predictions[ids,view.i] <- view.predict(ids.unlabelled,view.list[[view.i]]) 
    }


    ## Add unknown labels to known, where view agreement meets requirements
    if(weighting){
      
      if(updating){
          ## update accuracies
          view.list <- update.accuracies(view.list, known.labels)
          ## update weighting.threshold
          sum.acc <- 0
          for(view.i in 1:length(view.list)){
            sum.acc = sum.acc + view.list[[view.i]]$acc.norm
            print( view.list[[view.i]]$acc.norm)
          }
          weighting.threshold <- majority.threshold.percent*sum.acc/100
      }
      
      
      new.labelled.list <- get.new.labels.majorityWeighted(predictions,view.list,unique.labels)
      
      # Quit, if the maximal prediction value reached is lower than the given threshold
      if(new.labelled.list$weighting.threshold.upper < weighting.threshold){
        break
      } else{
        new.labelled <- new.labelled.list$new.labelled
        weighting.threshold.upper <- new.labelled.list$weighting.threshold.upper
        weighting.threshold.lower <- new.labelled.list$weighting.threshold.lower
      }
      
      
    } else {
      new.labelled <- get.new.labels.majorityCount(predictions,majority,majority.missingData)
    }
    
    
    colnames(new.labelled) <- colnames(labels)
    
    ## Add new labels
    labels <- rbind(labels,new.labelled)
    new.labels <- rbind(new.labels, new.labelled)
    
    ## Remove newly found labels from unlabelled IDs
    ids.unlabelled <- ids.unlabelled[!(ids.unlabelled %in% rownames(labels))]

    ## Collect information about the iterations for expanded output
    if(expanded.output){
      if(weighting){
        iteration.list$weighting.threshold.upper <- weighting.threshold.upper
        iteration.list$weighting.threshold.lower <- weighting.threshold.lower
      }
      iteration.list$no.new.labelled <- length(new.labelled)
      iteration.list$no.ids.left.unlabelled <- length(ids.unlabelled)
      iteration.list$view.list <- view.list
      collect.iteration.lists[[i]] <- iteration.list
      
      unlabelled.matrix[rownames(new.labels),i] <- new.labels
      for (view.i in 1:length(view.list)) {
        unlabelled.matrices.list[[view.i]][rownames(predictions),i] <- predictions[,view.i]
      }
    }

    

    
      
      
    if(!weighting) {
      ## Adjust majority value if we didn't learn any new labels for three iterations
      if( length(new.labelled) <= 0 ) {
        iterations.woChange <- iterations.woChange + 1
      } else{
        iterations.woChange <- 0 #reset if new labels were learned
      }
      ## Majority value for missing data is reduced first (prefer missing data over contrary predictions)
      ## Quit, if both majority values are already at the threshold 
      if(iterations.woChange >= 3){ 
        if(majority > majority.threshold){
          if(majority != majority.missingData){ # adjust overall majority value
            majority <- majority - 1
          } else{
            majority.missingData <- majority.missingData - 1  # adjust the majority value for missing data first
          }
          iterations.woChange <- 0
        } else{
          break
        }
      }
    }
    

    

    
  }

  ## Collect information for the returned result
  if(weighting){
    final.result.list <- list(final.views=view.list,labels.complete=labels,labels.new=new.labels,unlabelled.ids=ids.unlabelled
                              ,weighting.threshold=weighting.threshold,no.iterations=no.iterations,expanded.output=expanded.output)
  } else {
    final.result.list <- list(final.views=view.list,labels.complete=labels,labels.new=new.labels,unlabelled.ids=ids.unlabelled
                              ,majority.threshold=majority.threshold,no.iterations=no.iterations,expanded.output=expanded.output)
  }

  if(expanded.output){
    final.result.list$iteration.information <- collect.iteration.lists
    final.result.list$labelling.matrix <- unlabelled.matrix
    final.result.list$labelling.matrices.views <- unlabelled.matrices.list
  }
  return(final.result.list)
}


#platypus(commandArgs(TRUE))

