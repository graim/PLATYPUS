## Kiley Graim
## March 2017

## Load the data, extract iterations 1 and 3
dat <- read.table('perf_mvl_expanded.tab', sep='\t', header=T)

dat.iters <- dat[,c(1:2,45:55)]
dat.iters <- dat.iters[dat.iters$iteration==1 | dat.iters$iteration==3,]

dat.iter1 <- dat.iters[dat.iters$iteration==1,]
iter1 <- apply(dat.iter1, 2, mean)

dat.iter3 <- dat.iters[dat.iters$iteration==3,]
iter3 <- apply(dat.iter3, 2, mean)

iters <- rbind(iter1, iter3)

## Save the AUC changes
write.table(iters[,3:13], file='AUC_changes_iter1v3.txt', sep=',', col.names=T, row.names=T, quote=F)

## Calculate relative change in AUC
iterdelta <- (iter3-iter1)/(1-iter1)
iterdelta <- iterdelta[4:13]

## Calculate raw change in AUC
rawdelta <- iter3-iter1
rawdelta <- rawdelta[4:13]

## Make the view names match
names(iterdelta) <- c('Mutations', 'DrugTargets', 'GeneExpression', 'PathwayOncogenicGeneSets', 'PathwayMotifGeneSets', 'PathwayImmunologicalGeneSets', 'DrugTargetPathway', 'PathwayTranscriptionFactorTargets', 'ClinicalData', 'DruggableGenes')
names(rawdelta) <- c('Mutations', 'DrugTargets', 'GeneExpression', 'PathwayOncogenicGeneSets', 'PathwayMotifGeneSets', 'PathwayImmunologicalGeneSets', 'DrugTargetPathway', 'PathwayTranscriptionFactorTargets', 'ClinicalData', 'DruggableGenes')
write.table(signif(iterdelta,digits=3), file='AUC_Changes_Iter3minus1.txt', sep=',', col.names=TRUE, quote=F)


## Load the feature weight change data
feat.delta <- read.table('temp.tab', sep='\t', header=T, row.names=2)
names(iterdelta)[7] <- 'DTPathwayOncogenic'
names(iterdelta)[9] <- 'Clinical'
names(rawdelta)[7] <- 'DTPathwayOncogenic'
names(rawdelta)[9] <- 'Clinical'
iterdelta <- iterdelta[c('PathwayOncogenicGeneSets','PathwayTranscriptionFactorTargets','PathwayMotifGeneSets','DruggableGenes','PathwayImmunologicalGeneSets','GeneExpression','DrugTargets','DTPathwayOncogenic','Clinical','Mutations')]
rawdelta <-   rawdelta[c('PathwayOncogenicGeneSets','PathwayTranscriptionFactorTargets','PathwayMotifGeneSets','DruggableGenes','PathwayImmunologicalGeneSets','GeneExpression','DrugTargets','DTPathwayOncogenic','Clinical','Mutations')]

## Combine and save the changes
write.table(signif(iterdelta,digits=3), file='AUC_Changes_Iter3minus1.txt', sep=',', col.names=TRUE, quote=F)
feat.delta$AUCChange <- iterdelta[rownames(feat.delta)]
feat.delta$AUCRawChange <- rawdelta[rownames(feat.delta)]

write.table(feat.delta, 'MVL_Changes_AUC_features.tab', sep='\t', col.names=T, row.names=T, quote=F)

library( ggplot2)
library( reshape2 )

plot(feat.delta$mean, feat.delta$AUCChange, col=feat.delta$model)
plot(feat.delta$mean, feat.delta$AUCChange, col=feat.delta$model, xlim=c(0,0.01))
plot(feat.delta$mean, feat.delta$AUCChange, col=feat.delta$model, xlim=c(0,0.1))
