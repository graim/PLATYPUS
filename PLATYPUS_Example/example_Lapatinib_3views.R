## TODO: Build this up into a full example of how to take 1 data set and create interpreted views from it
## So create the interpreted data, run function to do parameter sweep & create config file, then use that in an MVL run

###############
## Initial data (provided with example)
######

## Load libraries
library(devtools)
install_github('graim/PLATYPUS')
library(PLATYPUS)

## Baseline view data: TODO 

## Labels to use: TODO

###############
## Create a new interpreted view
######

# Convert the PanCanAtlas data to the same views
# TODO: Flesh this out
#generate.feature.data.summary(gene.sample.data, feature.gene.matrix, type = 'expression', value = 'median', num.gene.threshold = 5)

###############
## Create configuration files for platypus using already created views 
## Create a config file for a set of tasks and views (e.g. data)
######

## Views to use
fns.data <- c(
  '2015-06-20_DrugTargetsExpression.txt', 
  '2015-07-02_TissueByMetabolicEnzymeExpression.txt', 
  '2015-09-08_DrugTargetNearestNeighboursPathway_HallmarkAllExpression.txt', 
  '2015-07-20_AllSummaryMetrics.txt', 
  'CCLE_expression_5000_featurematrix.tab', 
  'DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt' 
)
fns.data <- file.path('data',fns.data)

## Task file - each column is 1 task. Requires 1 or more tasks
fn.tasks <- 'CCLE_responselabel_binary_3cat_Lapatinib.tab'
fn.tasks <- file.path('data',fn.tasks)

## Generate config files - returns a list of the config filenames for using in platypus calls
config.files.en <- gen.config(fns.data[1:3], fn.tasks, model.type='en', delim='\t', config.loc='config') # elastic net configs
config.files.rf <- gen.config(fns.data[4:length(fns.data)], fn.tasks, model.type='rf', delim='\t', config.loc='config') # random forest configs
config.files    <- c(config.files.en, config.files.rf)

###############
## Run platypus using provided data
######

# rf, rf, rf, en
# TODO: Update this to have 1 baseline view, 1 interpreted view, 1 RF view, and 1 EN view. Any combo of those things
config.files <- list('DrugTargets_Lapatinib_configFile.txt','DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt','MetabolicEnzymes_Lapatinib_configFile.txt')
config.files <- file.path('config', config.files)

print('Views in use:'); flush.console()
print(config.files)

## Set parameters
n.iters  <- 10                                             # Run for this many iterations of label learning
m.thresh <- 95                                             # Threshold for learning a label (0-1 range)
of.name  <- 'platypus_output'                         # Directory where results should be stored
fn.labs  <- 'CCLE_responselabel_binary_3cat_Lapatinib.tab' # Filename of the labels being predicted

## Call multiview learning function
print('Running platypus')
platypus.res <- platypus(fn.views=config.files,fn.labs=fn.labs,e=TRUE,i=n.iters,m=m.thresh)
## TODO: Add a couple of lines highlighting how to use this

## Call label learning validation (LLV)
print('Running llv')
llv.platypus.res <- llv.platypus(fn.views=config.files,fn.labs=fn.labs,no.iterations=n.iters,majority.threshold.percent=m.thresh,output.folder=of.name)
 
## Call cross validation (CV)
print('Running cv')
cv.platypus.res <- cv.platypus(fn.views=config.files,fn.labs=fn.labs,no.iterations=n.iters,majority.threshold.percent=m.thresh,output.folder=of.name,expanded.output=TRUE)

## After running LLV and CV, generate the PLATYPUS performance plots
## TODO: R reads these as plot() functions for data types llv and cv. Would like to update my code to make that work!
plot.llv( fn.labs, 'platypus_output' ) # Plot label learning
plot.cv(fn.labs, 'platypus_output') # Plot cross validation - these plots pair with the llv plot
