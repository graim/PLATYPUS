## Kiley Graim
## April 2017


## Libraries
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

## Load the data
load('vsums2.RData')
#vsums$PLATYPUS <- as.vector(vsums$PLATYPUS)
#vsums$AUCdiff <- as.vector(vsums$AUCdiff)
#aucs <- read.table('AUC_vs_MeanFeatureChanges_Renamed.txt', sep='\t', header=T, row.names=1)

## Recalculate the AUC for PLATYPUS and ensemble, per view
#vsums$Ensemble <- as.vector(vsums$AUC)
#vsums$PLATYPUS <- as.vector(vsums$AUC) + as.vector(vsums$AUCdiff)
#vsums$View <- rownames(vsums)

## Drop the extra columns
#vsums <- vsums[,2:6]

## Will plot model types separately, so store full info
vsums.stored <- vsums

## Random Forest
vsums <- vsums.stored
vsums <- vsums[vsums$model=='rf',]
vsums$mean <- vsums$mean * 10
ggplot(vsums, aes(x=mean, y=Ensemble)) +
  ylim(0.75,0.85) +
  labs(y='AUC', x='Avg Change in Feature Importance') +
  geom_segment(aes(xend=vsums$mean, yend=vsums$PLATYPUS, color=Ensemble<PLATYPUS), size=1.5, arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  scale_color_brewer(palette="Set2") +
  geom_point(aes(color=Ensemble<PLATYPUS),size=5) +
#  geom_text_repel(aes(y=PLATYPUS, label=rownames(vsums)),size=6) +
  theme_linedraw(base_size = 18) +
  theme(legend.title=element_blank(),legend.position="none")
ggsave('AUCChanges_RF_noNames.pdf')

## Elastic Net 
vsums <- vsums.stored
vsums <- vsums[vsums$model=='en',]
#dev.new()
ggplot(vsums, aes(x=mean, y=Ensemble)) +
  ylim(0.75,0.85) +
#  ylim(0.745,0.81) +
  scale_color_brewer(palette="Set2") +
  labs(y='AUC', x='Avg Change in Feature Importance') +
  geom_segment(aes(xend=vsums$mean, yend=vsums$PLATYPUS, color=Ensemble<PLATYPUS), size=1.5, arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_point(aes(color=Ensemble<PLATYPUS),size=5) + 
#  geom_text_repel(aes(y=PLATYPUS, label=rownames(vsums)),size=6) +
  theme_linedraw(base_size = 18) +
  theme(legend.title=element_blank(),legend.position="none")
ggsave('AUCChanges_EN_noNames.pdf')



### Moving on to part (c) - drugtargets

#rf_importance_DrugTargets.differences.txt

dat <- read.table('rf_importance_DrugTargets.differences.txt', row.names=1, check.names=F, header=F)

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
#ggsave('RF_DrugTargets_Differences.scatter.pdf')


## Try 2
library(plyr)
aucs$View <- rownames(aucs)
vsums$View <- rownames(vsums)

temp <- join(vsums, aucs, by='View')
colnames(temp)[10] <- 'Ensemble'
colnames(temp)[1] <- 'Ensemble'
temp$AUCdiff <- temp$AUCRawChange
temp$PLATYPUS <- temp$Ensemble + temp$AUCRawChange
vsums <- temp
