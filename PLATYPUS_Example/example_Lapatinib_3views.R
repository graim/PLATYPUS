#!/usr/bin/env Rscript
#SBATCH -p ccb --qos=ccb -N3 --exclusive -o ./logs/CCLE_PLATYPUS_Example.%j.out -e ./logs/CCLE_PLATYPUS_Example.%j.err

## Set up an MVL call for the 3 view Lapatinib
#config.files <- list('DrugTargets_Lapatinib_configFile.txt','DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt','MetabolicEnzymes_Lapatinib_configFile.txt')
config.files <- list('DrugTargets_Lapatinib_configFile.txt','DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt','MetabolicEnzymes_Lapatinib_configFile.txt','AllSummaryMetrics_Lapatinib_configFile.txt')
config.files <- paste0('./config_CCLE/', config.files)
print('Views in use:'); flush.console()
print(config.files);flush.console()
print('Loaded config files');flush.console()

## Load functions
source('../scripts/platypus.basicFunctions.R')
source('../scripts/platypus.R')
source('../scripts/llv.platypus.R')
source('../scripts/cv.platypus.R')
print('Scripts loaded.'); flush.console()


## Set parameters
n.iters  <- 100 
m.thresh <- 90
of.name  <- 'platypus_output_CCLE'
fn.labs  <- 'CCLE_responselabel_binary_3cat_Lapatinib.tab'

## Call MVL
print('Running platypus');flush.console()
platypus.res <- platypus(fn.views=config.files,fn.labs=fn.labs,w=TRUE,e=TRUE,i=n.iters,m=m.thresh)
## TODO: Add a couple of lines highlighting how to use this

## Call LLV
print('Running llv');flush.console()
llv.platypus.res <- llv.platypus(fn.views=config.files,fn.labs=fn.labs,no.iterations=n.iters,majority.threshold.percent=m.thresh,output.folder=of.name, weighting=TRUE)

## Call CV
print('Running cv');flush.console()
cv.platypus.res <- cv.platypus(fn.views=config.files,fn.labs=fn.labs,no.iterations=n.iters,majority.threshold.percent=m.thresh,weighting=TRUE,output.folder=of.name,expanded.output=TRUE)
## TODO: add a line to make the boxplot of single view vs PLATYPUS results for this
## TODO: this returns NULL. 

print('Running llv');flush.console()
#save(platypus.res, llv.platypus.res, cv.platypus.res, file='platypus_example_Lap.RData')

## TODO: Build this up into a full example of how to take 1 data set and create interpreted views from it
## So create the interpreted data, run function to do parameter sweep & create config file, then use that in an MVL run

# Convert the PanCanAtlas data to the same views
#generate.feature.data.summary(gene.sample.data, feature.gene.matrix, type = 'expression', value = 'median', num.gene.threshold = 5)


## After running LLV and CV, generate the PLATYPUS performance plots
source('../scripts/plot.llv.R')
source('../scripts/plot.cv.R')
plot.llv( fn.labs, 'platypus_output_CCLE' ) # Plot label learning
plot.cv(fn.labs, 'platypus_output_CCLE') # Plot cross validation - these plots pair with the llv plot

