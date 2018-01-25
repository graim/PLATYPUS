## Aug 2016


## Read in each of the drug view scores files
## Output a matrix of all combined
load('ids.RData')

dat <- matrix(NA, nrow=24, ncol=26)
rownames(dat) <- drug.ids
colnames(dat) <- view.ids

for(dn in drug.ids) {
	fn <- paste0(dn,'.csv')
	temp <- read.table(fn, sep=',', row.names=1, stringsAsFactors=F)
	colnames(temp) <- temp[1,]
	dat.view <- temp[2,,drop=F]
	dat[dn,] <- as.numeric(unlist(dat.view[,colnames(dat)]))
}

write.table(dat, file='View_Scores.tab', sep='\t', col.names=T, row.names=T, quote=F)


## Read in the view scores file, print a boxplot highlighting MVL scores

library( ggplot2)
library( reshape2 )

dat <- as.data.frame(t(read.table("View_Scores.tab", sep='\t', header=T, row.names=1)))
dat <- dat[,names(sort(apply(dat,2,mean)))]

df <- melt(dat)
colnames(df) <- c('drug','AUC')
temp <- dat[rownames(dat)[27:30],]
df$isMVL <- 'blue'
df$isMVL[df$AUC %in% unlist(temp)] <- 'red'

ggplot(df, aes(drug,AUC)) +
	geom_boxplot(aes(fill=drug),colour=rainbow(24)) + 
        theme_minimal() + 
	theme(axis.text.x=element_text(angle=45), legend.position='none') + 
	geom_point(position = position_jitter(width = 0.2), colour=df$isMVL) +
	geom_point(aes(y = dat['mvl_10view',1]), colour='red')



## Using boxplot

## Sort by best MVL score
temp <- apply(dat[rownames(dat)[27:30],],2,max)
dat <- dat[,names(sort(temp))]

## Plot and save to pdf
pdf('Drug_scores.pdf', width=7, height=4)
cols <- c('red','grey60','lightsteelblue1')
col.axes <- 'grey40'
sym <- c(13,1,20)
stripchart(dat[1:26,], vertical = TRUE, method = "jitter", pch = sym[3], col = cols[3],ylim=c(0.4,1), xlab = "",axes=FALSE,ylab='AUC',col.lab=col.axes)
boxplot(dat[1:26,],add=TRUE, xlab = "",axes=FALSE, border='midnightblue',outline=FALSE)

axis(1, labels = FALSE, col=col.axes)
text(x =  seq_along(colnames(dat)), y = par("usr")[3]-0.05, srt = 45, adj = 1, labels = colnames(dat), xpd = TRUE, tck=1, col=col.axes)
axis(2, seq(0.4,1,by=0.1), las=2, col=col.axes, col.axis=col.axes)

legend(1,0.99, c('MVL','Single View'), pch=c(13,20), col=cols[c(1,3)], lwd=2,box.col='grey90', lty=0)
#legend(1,0.95, c('Best MVL','MVL', 'Single View'), pch=sym, col=cols, lwd=2,box.col='grey90',lty=0)

# Only include if I want ALL mvl results plotted
#points(1:24, dat[27,], col=cols[2], pch=sym[2],lwd=1)
#points(1:24, dat[28,], col=cols[2], pch=sym[2],lwd=1)
#points(1:24, dat[29,], col=cols[2], pch=sym[2],lwd=1)
#points(1:24, dat[30,], col=cols[2], pch=sym[2],lwd=1)

# plots best MVL result
points(1:24, apply(dat[rownames(dat)[27:30],],2,max), col=cols[1], pch=sym[1],lwd=2)
dev.off()




## Plot just 1 drug, color by type of view
ids.baseline <- c('AllMutation','AllMutation','AllCNV','ClinicalAggregationBinned','AllExpression')
ids.geneset.muts <- c('PositionalMutations','OncogenicMutations','ImmunogenicMutations','MotifMutations','HallmarkMutations','TranscritionMutations')
ids.genesets <- c('HallmarkAll','MetabolicEnzymes','GeneEssentialityAllCellLines','DruggableGenes','TranscriptionAll','OncogenicAll','DrugTarget',
	'ChromatinModifyingEnzymes','MultiDrugResistant','DrugTargetPathwayImmunologic','PositionalAll','DrugTargetPathwayOncogenic','MotifAll')
ids.platypus <- c('PLATYPUSBaseline','PLATYPUSAll')

view.types <- c("Baseline","Mutation","Gene Set","Summary", "PLATYPUS")

cols <- c('mediumorchid3','dodgerblue3', 'mediumseagreen','goldenrod', 'red')
dat.tmp <- dat[1:28,'PD.0325901',drop=F]
rownames(dat.tmp)[27:28] <- c('PLATYPUSBaseline','PLATYPUSAll')
dat.tmp['PLATYPUSAll',1] <- max(dat[27:30,'PD.0325901',drop=F])
dat.tmp['PLATYPUSBaseline',1] <- 0.941 # From Verena
colnames(dat.tmp)[1] <- 'AUC'
dat.tmp <- dat.tmp[with(dat.tmp, order(AUC)), ,drop=F ]
dat.tmp$Color <- NA
dat.tmp[ids.baseline,'Color'] <- cols[1] 
dat.tmp[ids.geneset.muts,'Color'] <- cols[2]
dat.tmp[ids.genesets,'Color'] <- cols[3]
dat.tmp[ids.platypus,'Color'] <- cols[5]
dat.tmp[is.na(dat.tmp$Color),'Color'] <- cols[4]

pdf('PD0325901_view_scores.pdf')
par(mar=c(14,4,4,2))
plot(1:nrow(dat.tmp), dat.tmp[,1], ylab='AUC',xlab='', pch=20,lwd=7,xaxt='n', las=2,ylim=c(0.5,1),col=dat.tmp$Color)
axis(1, at=seq(1, nrow(dat.tmp), by=1), labels = rownames(dat.tmp), las=2)
abline(v=16.5, lty=3, col='grey60',lwd=2)
abline(h=max(dat.tmp$AUC), col='red', lty=3, lwd=2)
legend('topleft',title='View Type', col=cols, view.types,pch=20,pt.cex=2, box.col='black', bg='white')
dev.off()




PositionalMutations
HallmarkMutations           
OncogenicMutations
AllCNV                      
ImmunogenicMutations
TranscritionMutations       
MotifMutations
ChromatinModifyingEnzymes   
Viper
MultiDrugResistant          
HallmarkAll
DrugTargetPathwayImmunologic
MetabolicEnzymes
PositionalAll               
GeneEssentialityAllCellLines
DrugTargetPathwayHallmark   
DruggableGenes
ClinicalAggregationBinned   
TranscriptionAll
DrugTargetPathwayOncogenic  
ImmunologicMedianVariance
MotifAll                    
OncogenicAll
AllExpression
DrugTarget
AllMutation
