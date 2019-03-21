#!/usr/bin/env Rscript
#SBATCH -p ccb --qos=ccb -N3 --exclusive -o ./logs/CCLE_PLATYPUS_Example.%j.out -e ./logs/CCLE_PLATYPUS_Example.%j.err

###############
## Run platypus using provided data
######

## Load functions
source('../scripts/platypus.basicFunctions.R')
source('../scripts/platypus.R')
source('../scripts/llv.platypus.R')
source('../scripts/cv.platypus.R')
print('Scripts loaded.'); flush.console()

## Set up an MVL call for the 3 view Lapatinib
# rf, rf, rf, en
config.files <- list('DrugTargets_Lapatinib_configFile.txt','DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt','MetabolicEnzymes_Lapatinib_configFile.txt','AllSummaryMetrics_Lapatinib_configFile.txt')
config.files <- file.path('config_CCLE', config.files)

print('Views in use:'); flush.console()
print(config.files)

## Set parameters
n.iters  <- 10                                             # Run for this many iterations of label learning
m.thresh <- 95                                             # Threshold for learning a label (0-1 range)
of.name  <- 'platypus_output_CCLE'                         # Directory where results should be stored
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

## TODO: Build this up into a full example of how to take 1 data set and create interpreted views from it
## So create the interpreted data, run function to do parameter sweep & create config file, then use that in an MVL run

## After running LLV and CV, generate the PLATYPUS performance plots
source('../scripts/plot.llv.R')
source('../scripts/plot.cv.R')
plot.llv( fn.labs, 'platypus_output_CCLE' ) # Plot label learning
plot.cv(fn.labs, 'platypus_output_CCLE') # Plot cross validation - these plots pair with the llv plot

###############
## Create configuration files for platypus using already created views 
## Create a config file for a set of tasks and views (e.g. data)
######

## Load scripts
source('../scripts/gen_config.R')
source('../scripts/single_predictor.R')

## Views to use
fns.data <- c(
  '2015-06-20_DrugTargetsExpression.txt', # en
  '2015-07-20_AllSummaryMetrics.txt', # en
  'CCLE_clinical_featurematrix.tab', # en TODO: this is the old version, not used in PSB results. still has the 'NOT.' features
  '2015-07-02_TissueByMetabolicEnzymeExpression.txt', # en
  '2015-09-08_DrugTargetNearestNeighboursPathway_HallmarkAllExpression.txt', # en
  'CCLE_expression_5000_featurematrix.tab', # en
  'DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt' # rf
)

## Task file - each column is 1 task. Requires 1 or more tasks
fn.tasks <- 'CCLE_responselabel_binary_3cat_Lapatinib.tab'

## Append the correct location for the files above
fn.tasks <- file.path('data_CCLE',fn.tasks)
fns.data <- file.path('data_CCLE',fns.data)

## Generate config files - returns a list of the config filenames for using in platypus calls
config.files <- gen.config(fns.data, fn.tasks, model.type='en', delim='\t', config.loc='config_CCLE')
config.files <- gen.config(fns.data, fn.tasks, model.type='rf', delim='\t', config.loc='config_CCLE')

## Using the configs we've generated, run platypus
platypus.res <- platypus(fn.views=config.files,fn.labs=fn.tasks)

###############
## Create a new interpreted view
######

# Convert the PanCanAtlas data to the same views
generate.feature.data.summary(gene.sample.data, feature.gene.matrix, type = 'expression', value = 'median', num.gene.threshold = 5)

