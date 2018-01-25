## Set up an MVL call for the 3 view Lapatinib
config.files <- list('DrugTargets_Lapatinib_configFile.txt','DrugTargetPathway_HallmarkPathway_Lapatinib_configFile.txt','MetabolicEnzymes_Lapatinib_configFile.txt')
config.files <- paste0('./config/', config.files)

## Load functions
source('../scripts/platypus.basicFunctions.R')
source('../scripts/platypus.R')
source('../scripts/llv.platypus.R')
source('../scripts/cv.platypus.R')


## Call MVL
platypus.res <- platypus(fn.views=config.files,fn.labs='data/CCLE_responselabel_binary_3cat_Lapatinib.tab',k=5,w=TRUE,e=TRUE,i=5,m=50)
## TODO: Add a couple of lines highlighting how to use this

## Call LLV
llv.platypus.res <- llv.platypus(fn.views=config.files,fn.labs='data/CCLE_responselabel_binary_3cat_Lapatinib.tab',no.iterations=5,majority.threshold.percent=70,output.folder='platypus_output')
## TODO: add a line to plot the results

## Call CV
cv.platypus.res <- cv.platypus(fn.views=config.files,fn.labs='data/CCLE_responselabel_binary_3cat_Lapatinib.tab',no.iterations=5,majority.threshold.percent=70,weighting=TRUE)
## TODO: add a line to make the boxplot of single view vs PLATYPUS results for this
## TODO: this returns NULL. 


## TODO: Build this up into a full example of how to take 1 data set and create interpreted views from it
## So create the interpreted data, run function to do parameter sweep & create config file, then use that in an MVL run

# Convert the PanCanAtlas data to the same views
generate.feature.data.summary(gene.sample.data, feature.gene.matrix, type = 'expression', value = 'median', num.gene.threshold = 5)

