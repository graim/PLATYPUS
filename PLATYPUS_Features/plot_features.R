## Kiley Graim
## Feb 2017


## Load the summmary table
#vsums <- read.table('view_summary_differences.tab', sep='\t', header=T, row.names=3, check.names=F)
vsums <- read.table('MVL_Changes_AUC_features.tab', sep='\t', header=T, row.names=1, check.names=F)
vsums.en <- vsums[vsums$model=='en',]
vsums.rf <- vsums[vsums$model=='rf',]

## Load the PD0325901 AUCs 
aucs <- read.table('Model_AUCs_PD0225901.tab', sep='\t', header=T, row.names=1, check.names=F)

## Plot per view changes
library(ggplot2)
library(ggrepel)

ggplot(vsums, aes(x=mean, y=sum,color=model)) + 
  geom_point() +
  geom_text_repel(aes(mean,sum,label=rownames(vsums),color=model,point.padding = 0.5)) +
  theme_minimal()
ggsave('Views_mean_vs_sum.pdf')

## Load the view with the largest average change, drop Intercept because we don't care about it
dat <- read.table('en_coefficients_PathwayOncogenicGeneSets.differences.txt', row.names=1, check.names=F)
dat <- dat[-1,,drop=F]


## Plot the weight changes per feature
plot(dat[,1], type='l', lwd=3, col='navy', ylab='Feature Weight Change')


## Feature types
## Of the top 20 features with largest weight changes:
## 11 are kurtosis
## 8 are median
## 1 is variance

#MTOR_UP.V1_UP_KURTOSIS          1.10680
#CAMP_UP.V1_DN_KURTOSIS          0.90448
#MTOR_UP.N4.V1_UP_KURTOSIS       0.76098
#S100A12_MEDIAN                  0.72873
#NFE2L2.V2_KURTOSIS              0.61603
#IL2_UP.V1_UP_KURTOSIS           0.56041
#E2F3_UP.V1_UP_KURTOSIS          0.53066
#RB_P107_DN.V1_UP_KURTOSIS       0.49221
#HOXA9_DN.V1_DN_MEDIAN           0.48767
#P53_DN.V2_DN_MEDIAN             0.45707
#HOXA9_DN.V1_UP_KURTOSIS         0.45231
#ESC_V6.5_UP_EARLY.V1_UP_MEDIAN  0.45201
#PTEN_DN.V1_UP_MEDIAN            0.44405
#RAPA_EARLY_UP.V1_DN_KURTOSIS    0.41935
#PRC2_EZH2_UP.V1_UP_MEDIAN       0.41825
#DCA_UP.V1_UP_MEDIAN             0.38823
#RAF_UP.V1_DN_KURTOSIS           0.38393
#CAHOY_OLIGODENDROCUTIC_VARIANCE 0.38064
#RAPA_EARLY_UP.V1_UP_MEDIAN      0.36768
#SNF5_DN.V1_UP_KURTOSIS          0.36215


## Load the data for the most changing view
iter1 <- 'trainedViews_PD.0325901/trainedViews_Iteration1/'
iter2 <- 'trainedViews_PD.0325901/trainedViews_Iteration3/'
views <- unlist(read.table('view_names.txt', sep='\t', header=F, stringsAsFactors=F))

view <- views[5]
temp <- readRDS( paste0(iter1,view) )
view1 <- as(temp, 'matrix')
temp <- readRDS( paste0(iter2,view) )
view2 <- as(temp, 'matrix')
view.diff <- view1[,1] - view2[,1]

view.both <- cbind(view1, view2, view.diff)
view.both <- view.both[-1,,drop=F]
colnames(view.both) <- c('ENSEMBLE', 'PLATYPUS', 'difference')

## Plot ensemble vs PLATYPUS weights
plot(view.both[,1:2])

## Segments version of the first
library(ggplot2)
library(reshape2)

df <- data.frame(view.both[which(view.both[,1]!=0 | view.both[,2]!=0),])

ggplot(df, aes(x=ENSEMBLE, y=PLATYPUS)) + 
  geom_point(colour='navy',stroke=2, shape=4) +
  geom_abline(intercept = 0, slope = 1, color="grey60", linetype="dashed", size=1) + 
  theme_minimal() 
ggsave('PathwayOncogenicGeneSets.differences_scatter.svg')  
ggsave('PathwayOncogenicGeneSets.differences_scatter.pdf')  




## Grabbing the feature names, sorted
temp <- df$difference
names(temp) <- rownames(df)
write.table(sort(abs(temp), decreasing=T), file='sortedWeightChanges_PathwayOncogenicGeneSets.txt', sep='\t', col.names=T, row.names=T, quote=F)



## Version that plots model AUC vs mean feature weight changes
vsums <- read.table('AUC_vs_MeanFeatureChanges_Renamed.txt', sep='\t', header=T,check.names=F)
# rf_DrugTargets.RDS, rf_DTPathwayOncogenic.RDS, rf_Clinical.RDS, rf_Mutations.RDS

vsums[vsums$model=='rf','mean'] <-  vsums[vsums$model=='rf','mean']*10
vsums[vsums$model=='rf','mean'] <-  vsums[vsums$model=='rf','mean']/50

library(ggplot2)
library(ggrepel)

ggplot(vsums, aes(x=mean, y=AUC,color=model)) +
  geom_point() +
  geom_text_repel(aes(mean,AUC,label=rownames(vsums),color=model,point.padding = 0.5)) +
  theme_minimal()
ggsave('Views_AUC_vs_sum.pdf')



## Version that uses the AUC differences from Verena
load('vsums.RData')
library(ggplot2)
library(ggrepel)
library(RColorBrewer)


ggplot(vsums, aes(x=mean, y=AUCChange,color=model)) +
  geom_point() +
  geom_text_repel(aes(mean,AUCChange,label=rownames(vsums),color=model)) +
  theme_minimal() +
  xlab('Mean Weight Change') +
  ylab('(iter3 - iter1) / (1 - iter1)')
#ggsave('Views_AUCdiff_vs_AvgFeatChange.pdf')

cols <- brewer.pal(4, 'Set3')
ggplot(vsums.en, aes(x=mean, y=AUCChange,color=cols[4])) +
  geom_point(color=cols[4]) +
  geom_text_repel(aes(mean,AUCChange,label=rownames(vsums.en)),color=cols[4]) +
  theme_minimal() +
  xlab('Average Feature Weight Change') +
  ylab('Relative AUC Change') +
  geom_hline(aes(yintercept=0),linetype='dashed',color='grey60')
ggsave('ViewsEN_AUCdiff_vs_AvgFeatChange.pdf')

cols <- brewer.pal(4, 'Dark2')
ggplot(vsums.rf, aes(x=mean, y=AUCChange,color=cols[1])) +
  geom_point(color=cols[1]) +
  geom_text_repel(aes(mean,AUCChange,label=rownames(vsums.rf)),color=cols[1]) +
  theme_minimal() +
  xlab('Average Gini Decrease') +
  ylab('Relative AUC Change') +
  geom_hline(aes(yintercept=0),linetype='dashed',color='grey60')
ggsave('ViewsRF_AUCdiff_vs_AvgFeatChange.pdf')



## Plot the change in drug targets views - DrugTargets RF version
dat.iter1 <- read.table('rf_DrugTargets.RDS_iter1_differences.txt', row.names=1, check.names=F)
dat.iter3 <- read.table('rf_DrugTargets.RDS_iter3_differences.txt', row.names=1, check.names=F)

ids <- intersect( rownames(dat.iter1), rownames(dat.iter3) )
plot( dat.iter1[ids,1], dat.iter3[ids,1], col='navy')

view.both <- cbind(dat.iter1[ids,1], dat.iter3[ids,1])
colnames(view.both) <- c('ENSEMBLE', 'PLATYPUS')
rownames(view.both) <- ids

library(ggplot2)
library(reshape2)

df <- as.data.frame(view.both)
df$name <- rownames(df)

ggplot(df, aes(x=ENSEMBLE, y=PLATYPUS)) +
  geom_point(colour='navy',stroke=1.5, shape=4) +
  geom_abline(intercept = 0, slope = 1, color="grey60", linetype="dashed", size=1) +
  coord_cartesian(xlim=range(view.both), ylim=range(view.both)) +
  theme_minimal() +
  geom_text_repel(data = subset(df, PLATYPUS > 3 | ENSEMBLE > 2), aes(label = name),color='grey40')
ggsave('RF_DrugTargets_GiniDifferences.scatter.pdf')



## Plot the change in drug targets views - DrugTargets RF version
dat <- read.table('en_coefficients_PathwayOncogenicGeneSets.differences.txt', row.names=1, check.names=F)

dat.iter1 <- read.table('rf_DrugTargets.RDS_iter1_differences.txt', row.names=1, check.names=F)
dat.iter3 <- read.table('rf_DrugTargets.RDS_iter3_differences.txt', row.names=1, check.names=F)

ids <- intersect( rownames(dat.iter1), rownames(dat.iter3) )
plot( dat.iter1[ids,1], dat.iter3[ids,1], col='navy')

view.both <- cbind(dat.iter1[ids,1], dat.iter3[ids,1])
colnames(view.both) <- c('ENSEMBLE', 'PLATYPUS')
rownames(view.both) <- ids

library(ggplot2)
library(reshape2)

df <- as.data.frame(view.both)
df$name <- rownames(df)

ggplot(df, aes(x=ENSEMBLE, y=PLATYPUS)) +
  geom_point(colour='navy',stroke=1.5, shape=4) +
  geom_abline(intercept = 0, slope = 1, color="grey60", linetype="dashed", size=1) +
  coord_cartesian(xlim=range(view.both), ylim=range(view.both)) +
  theme_minimal() +
  geom_text_repel(data = subset(df, PLATYPUS > 3 | ENSEMBLE > 2), aes(label = name),color='grey40')
ggsave('RF_DrugTargets_GiniDifferences.scatter.pdf')

## Plot feature weight changes in View_PathwayOncogenicGeneSets_summary_temp.tab view


dat <- read.table('View_PathwayOncogenicGeneSets_summary_temp.tab', row.names=1, check.names=F, header=T)

library(ggplot2)
library(reshape2)

df <- as.data.frame(dat)
df$name <- rownames(df)

ggplot(df, aes(x=ENSEMBLE, y=PLATYPUS)) +
  geom_point(colour='navy',stroke=1.5, shape=4) +
  geom_abline(intercept = 0, slope = 1, color="grey60", linetype="dashed", size=1) +
  coord_cartesian(xlim=range(df[,1:2]), ylim=range(df[,1:2])) +
  theme_minimal() +
  geom_text_repel(data = subset(df, difference > 0.5 | difference < -0.5), aes(label = name),color='grey40')
ggsave('EN_PathwayOncogenicGeneSets_Differences.scatter.pdf')









