# For each column, return the top n features in data frame dat
get.top <- function(dat,n) {
	dat <- apply(dat, 2, function(x) {rank(abs(x))} )
	#dat <- t(t(dat)/apply(t(dat), 1, max))
        return(apply(dat, 2, function(x) { names(rev(sort(abs(x))))[1:n] } ))
}

get.top <- function(dat,n) {
	dat <- apply(dat, 2, function(x) {abs(x)} )
	dat <- dat/max(dat)
        return(apply(dat, 2, function(x) { rev(sort(abs(x)))[1:n] } ))
}
# All files with scores
fns <- c(
	'ChromatinModifyingEnzymes_Importance_randomForest.txt','OncogenicGeneSets_Importance_randomForest.txt','Pathway_TranscriptionFactorTargets_Importance_randomForest.txt','Pathway_MotifGeneSets_Importance_randomForest.txt','Pathway_PositionalGeneSets_Importance_randomForest.txt','Pathway_HallmarkGeneSets_Importance_randomForest.txt','DrugTargets_Importance_randomForest.txt','GeneExpression_Importance_randomForest.txt','GeneExpression_LandmarkGenes_Importance_randomForest.txt','MultiDrugResistant_Importance_randomForest.txt','CTDD_AZD6224_GeneExpression_Importance_randomForest.txt','ClinicalData_Importance_randomForest.txt','DiffusedMutations_Importance_randomForest.txt'
#	'Pathway_OncogenicGeneSets_Coefficients_elasticNet.txt','Pathway_TranscriptionFactorTargets_Coefficients_elasticNet.txt','Pathway_MotifGeneSets_Coefficients_elasticNet.txt','Pathway_PositionalGeneSets_Coefficients_elasticNet.txt','Pathway_HallmarkGeneSets_HallmarkCoefficients_elasticNet.txt','GeneExpression_Coefficients_elasticNet.txt','DiffusedMutations_Coefficients_elasticNet.txt'
)


get_top_X <- function(fn, x) {
        dat <- read.table(fn, sep='\t', header=T, row.names=1, check.names=F)
        dat.scaled <- dat/apply(dat, 1, max)  # rescale 0-1
        top <- apply(dat, 2, function(y) {names(head(sort(y)))} )
}

top <- matrix(nrow=24, ncol=0)
rownames(top) <- c('AEW541','Nilotinib','X17.AAG','PHA.665752','Lapatinib','Nutlin.3','AZD0530',
	'PF2341066','L.685458','ZD.6474','Panobinostat','Sorafenib','Irinotecan','Topotecan',
	'LBW242','PD.0325901','PD.0332991','Paclitaxel','AZD6244','PLX4720','RAF265','TAE684','TKI258','Erlotinib')

dat.all <- matrix(nrow=24, ncol=0)
rownames(dat.all) <- c('AEW541','Nilotinib','X17.AAG','PHA.665752','Lapatinib','Nutlin.3','AZD0530',
	'PF2341066','L.685458','ZD.6474','Panobinostat','Sorafenib','Irinotecan','Topotecan',
	'LBW242','PD.0325901','PD.0332991','Paclitaxel','AZD6244','PLX4720','RAF265','TAE684','TKI258','Erlotinib')

for(fn in fns) {
	dat <- read.table(fn, sep='\t', header=T, row.names=1, check.names=F)
	top <- cbind(top,t(get.top(dat,10)))
	dat.all <- cbind(dat.all, t(dat))
}

temp <- dat.all[,unique(unlist(top))]
boxplot(t(temp), las=2, ylim=c(0,1.5))

for( rn in rownames(top) ) {
	write.table(sort(unique(top[rn,])), file=paste0(rn,'_top10FeaturesPerView.txt'), sep=',', col.names=F, row.names=F, quote=F)
}

# Plot the top N features for each drug across all models
## Plot and save to pdf

#pdf('Drug_scores.pdf', width=7, height=4)
cols <- c('red','grey60','lightsteelblue1')
col.axes <- 'grey40'
sym <- c(13,1,20)
stripchart(dat, vertical = TRUE, method = "jitter", pch = sym[3], col = cols[3], xlab = "",axes=FALSE,ylab='Absolute Weight',col.lab=col.axes)
boxplot(dat,add=TRUE, xlab = "",axes=FALSE, border='midnightblue',outline=FALSE)
#boxplot(dat, xlab = "",axes=FALSE, border='midnightblue',outline=FALSE)

axis(1, labels = FALSE, col=col.axes)
text(x =  seq_along(colnames(dat)), y = par("usr")[3]-1, srt = 45, adj = 1, labels = colnames(dat), xpd = TRUE, tck=1, col=col.axes)
axis(2, seq(as.integer(min(dat)),as.integer(max(dat)),by=5), las=2, col=col.axes, col.axis=col.axes)

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


