## Kiley Graim
## Code to make Fig 3

## Load the data
load('Wrangled_Scores.RData')

## Make the plot
pdf('Drug_scores.pdf', width=7, height=4)
stripchart(dat[1:26,], vertical = TRUE, method = "jitter", pch = sym[3], col = cols[3],ylim=c(0.4,1), xlab = "",axes=FALSE,ylab='AUC',col.lab=col.axes) 
boxplot(dat[1:26,],add=TRUE, xlab = "",axes=FALSE, border='midnightblue',outline=FALSE)
axis(1, labels = FALSE, col=col.axes)
text(x =  seq_along(colnames(dat)), y = par("usr")[3]-0.05, srt = 45, adj = 1, labels = colnames(dat), xpd = TRUE, tck=1, col=col.axes)
axis(2, seq(0.4,1,by=0.1), las=2, col=col.axes, col.axis=col.axes)
legend(1,0.99, c('PLATYPUS','Single View'), pch=c(13,20), col=cols[c(1,3)], lwd=2,box.col='grey90', lty=0)
points(1:24, apply(dat[rownames(dat)[27:30],],2,max), col=cols[1], pch=sym[1],lwd=2)
text(x =  seq_along(colnames(dat)), y = par("usr")[3]-0.05, srt = 45, adj = 1, labels = colnames(dat), xpd = TRUE, tck=1, col=col.axes)
dev.off()
