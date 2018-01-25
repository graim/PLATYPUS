require(gplots)

#####
# ALL Together - new runs with all compounds

TKI <- c("AZD0530","AEW541","Nilotinib","Erlotinib","TKI258","Sorafenib","Lapatinib","ZD.6474")
BRAFi <- c("PLX4720","RAF265")
MEKi <- c("PD.0325901","AZD6244")
CDKi <- c("PD.0332991")
ALKfusion <- c("TAE684","PF2341066")
TOP1i <- c("Irinotecan","Topotecan")
IAPi <- c("LBW242")
HSP90i <- c("X17.AAG")
HDACi <- c("Panobinostat")
METi <- c("PHA.665752")
tubulin <- c("Paclitaxel")
p53 <- c("Nutlin.3")
PSEN <- c("L.685458")
compounds.sorted <- c(TKI
                      ,BRAFi
                      ,MEKi
                      ,CDKi
                      ,ALKfusion
                      ,TOP1i
                      ,IAPi
                      ,HSP90i
                      ,HDACi
                      ,METi
                      ,tubulin
                      ,p53
                      ,PSEN
)
#######################
### sorted compounds
TKI <- c("AZD0530","AEW541","Nilotinib","Erlotinib","TKI258","Sorafenib","Lapatinib","ZD-6474")
BRAFi <- c("PLX4720","RAF265")
MEKi <- c("PD-0325901","AZD6244")
CDKi <- c("PD-0332991")
ALKfusion <- c("TAE684","PF2341066")
TOP1i <- c("Irinotecan","Topotecan")
IAPi <- c("LBW242")
HSP90i <- c("17-AAG")
HDACi <- c("Panobinostat")
METi <- c("PHA-665752")
tubulin <- c("Paclitaxel")
p53 <- c("Nutlin-3")
PSEN <- c("L-685458")
compounds.sorted.old <- c(TKI
                          ,BRAFi
                          ,MEKi
                          ,CDKi
                          ,ALKfusion
                          ,TOP1i
                          ,IAPi
                          ,HSP90i
                          ,HDACi
                          ,METi
                          ,tubulin
                          ,p53
                          ,PSEN
)

acc.matrix.tissue <- read.table("~/MVL/MVLruns_blood/2015-09-08_MethodAccuracyBlood.tab", sep='\t',header=T,row.names=1)
acc.matrix.tissue <- acc.matrix.tissue[compounds.sorted,]

acc.matrix <- read.table("~/MVL/MVLruns/2015-08-27_MethodAccuracy.txt", sep='\t',header=T,row.names=1)
acc.matrix <- acc.matrix[compounds.sorted,]



ds.views.tissue <- c("Clinical","AllExpression","AllMutations")
acc.matrix.ds.views.tissue <- as.matrix(acc.matrix.tissue[,ds.views.tissue])
ds.views <- c("ClinicalAggregationBinned","AllExpression","AllMutation")
acc.matrix.ds.views <- as.matrix(acc.matrix[,ds.views])
acc.matrix.ds <- rbind(acc.matrix.ds.views.tissue,acc.matrix.ds.views)

gs.views.tissue <- c("MetabolicEnzymes","MultiDrugResistant","DrugTargets","ChromatinModifyingEnzymes","GeneEssentiality")
acc.matrix.gs.view.tissues <- as.matrix(acc.matrix.tissue[,gs.views.tissue])
gs.views <- c("MetabolicEnzymes","MultiDrugResistant","DrugTarget","ChromatinModifyingEnzymes","GeneEssentialityAllCellLines")
acc.matrix.gs.views <- as.matrix(acc.matrix[,gs.views])
acc.matrix.gs <- rbind(acc.matrix.gs.view.tissues,acc.matrix.gs.views)

msig.views.tissue <- c("Hallmark","Motif","Transcription","Positional","Oncogenic","Immunologic")
acc.matrix.msig.views.tissue <- as.matrix(acc.matrix.tissue[,msig.views.tissue])
msig.views <- c("HallmarkAll","MotifAll","TranscriptionAll","PositionalAll","OncogenicAll","ImmunologicMedianVariance")
acc.matrix.msig.views <- as.matrix(acc.matrix[,msig.views])
acc.matrix.msig <- rbind(acc.matrix.msig.views.tissue,acc.matrix.msig.views)

msig.dt.views.tissue <- c("DrugTargetHallmarkPathway","DrugTargetsOncogenicPathway","DrugTargetImmunologicPathway")
acc.matrix.msig.dt.view.tissues <- as.matrix(acc.matrix.tissue[,msig.dt.views.tissue])
msig.dt.views <- c("DrugTargetPathwayHallmark","DrugTargetPathwayOncogenic","DrugTargetPathwayImmunologic")
acc.matrix.msig.dt.views <- as.matrix(acc.matrix[,msig.dt.views])
acc.matrix.msig.dt <- rbind(acc.matrix.msig.dt.view.tissues,acc.matrix.msig.dt.views)

viper.tissue <- c("Viper")
acc.matrix.viper.tissue <- as.matrix(acc.matrix.tissue[,viper.tissue])
viper <- c("Viper")
acc.matrix.viper <- as.matrix(acc.matrix[,viper])
acc.matrix.viper.bind <- rbind(acc.matrix.viper.tissue,acc.matrix.viper)

plot_matrix <- cbind(acc.matrix.ds
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.gs
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.msig
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.msig.dt
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.viper.bind)

n <- 24
alternating <- c(1,1+n,2,2+n,3,3+n,4,4+n,5,5+n,6,6+n,7,7+n,8,8+n,9,9+n,10,10+n,11,11+n,12,12+n,13,13+n,14,14+n,15,15+n,16,16+n,17,17+n,18,18+n,19,19+n,20,20+n,21,21+n,22,22+n,23,23+n,24,24+n)
plot_matrix <- plot_matrix[alternating,]

collabels <- c(paste0(compounds.sorted.old," Blood"),paste0(compounds.sorted.old," General"))
collabels <- collabels[alternating]


###### rotated 
pdf(file="~/MVL/SingleViews/single.views_tissue_rotated.pdf",width=7,height=14)

breaks.vec <- c(0,seq(0.5,1,length.out=50))
colfunc <- colorRampPalette(c("white","darkblue"))

# data specific
heatmap.2(plot_matrix
          ,dendrogram="none"
          ,Rowv=F
          ,Colv=F
          ,na.rm=T
          ,breaks=sort(breaks.vec)
          ,col=colfunc(length(breaks.vec)-1)
          
          ,na.color="grey"
          
          ,labCol=c("Clinical","Expression","Mutation",""
                    ,"Metabolic Enzymes","Multi-Drug Resistant","Drug Targets","Chromatin Modifying","Essential Genes",""
                    ,"Hallmark GS Expr","Motif GS Expr","TF Targets Expr","Positional GS Expr","Oncogenic GS Expr","Immunologic GS Expr",""
                    ,"DT GS Hallmark","DT GS Oncogenic","DT GS Immunologic",""
                    ,"Viper")
          ,labRow=collabels
          
          ,cexCol=1.1
          ,cexRow=1.1
          
          ,key=T
          ,density.info="none"
          ,key.xlab="Accuracy"
          ,key.title=NA
          ,key.par=list(mar=c(4,2,1,11))
          
          ,margins=c(9.5,11)
          ,trace="none"
          ,sepwidth=c(0.01,0.01)
          ,sepcolor="black"
          ,rowsep=c(0,(1:24)*2,nrow(plot_matrix))
          ,colsep=c(0,22)
          
          ,lmat = rbind(c(0,3),c(2,1),c(0,4))
          ,lwid = c(0.05,4)
          ,lhei = c(0.05,4,0.3)
)
dev.off()



#####
# ALL Together - new runs with 2 compounds


compounds.sorted <- c("AEW541","AZD6244")

compounds.sorted.old <- c("AEW541","AZD6244")

acc.matrix.tissue <- read.table("~/MVL/MVLruns_blood/2015-09-08_MethodAccuracyBlood.tab", sep='\t',header=T,row.names=1)
acc.matrix.tissue <- acc.matrix.tissue[compounds.sorted,]

acc.matrix <- read.table("~/MVL/MVLruns/2015-08-27_MethodAccuracy.txt", sep='\t',header=T,row.names=1)
acc.matrix <- acc.matrix[compounds.sorted,]



ds.views.tissue <- c("Clinical","AllExpression","AllMutations")
acc.matrix.ds.views.tissue <- as.matrix(acc.matrix.tissue[,ds.views.tissue])
ds.views <- c("ClinicalAggregationBinned","AllExpression","AllMutation")
acc.matrix.ds.views <- as.matrix(acc.matrix[,ds.views])
acc.matrix.ds <- rbind(acc.matrix.ds.views.tissue,acc.matrix.ds.views)

gs.views.tissue <- c("MetabolicEnzymes","MultiDrugResistant","DrugTargets","ChromatinModifyingEnzymes","GeneEssentiality")
acc.matrix.gs.view.tissues <- as.matrix(acc.matrix.tissue[,gs.views.tissue])
gs.views <- c("MetabolicEnzymes","MultiDrugResistant","DrugTarget","ChromatinModifyingEnzymes","GeneEssentialityAllCellLines")
acc.matrix.gs.views <- as.matrix(acc.matrix[,gs.views])
acc.matrix.gs <- rbind(acc.matrix.gs.view.tissues,acc.matrix.gs.views)

msig.views.tissue <- c("Hallmark","Motif","Transcription","Positional","Oncogenic","Immunologic")
acc.matrix.msig.views.tissue <- as.matrix(acc.matrix.tissue[,msig.views.tissue])
msig.views <- c("HallmarkAll","MotifAll","TranscriptionAll","PositionalAll","OncogenicAll","ImmunologicMedianVariance")
acc.matrix.msig.views <- as.matrix(acc.matrix[,msig.views])
acc.matrix.msig <- rbind(acc.matrix.msig.views.tissue,acc.matrix.msig.views)

msig.dt.views.tissue <- c("DrugTargetHallmarkPathway","DrugTargetsOncogenicPathway","DrugTargetImmunologicPathway")
acc.matrix.msig.dt.view.tissues <- as.matrix(acc.matrix.tissue[,msig.dt.views.tissue])
msig.dt.views <- c("DrugTargetPathwayHallmark","DrugTargetPathwayOncogenic","DrugTargetPathwayImmunologic")
acc.matrix.msig.dt.views <- as.matrix(acc.matrix[,msig.dt.views])
acc.matrix.msig.dt <- rbind(acc.matrix.msig.dt.view.tissues,acc.matrix.msig.dt.views)

viper.tissue <- c("Viper")
acc.matrix.viper.tissue <- as.matrix(acc.matrix.tissue[,viper.tissue])
viper <- c("Viper")
acc.matrix.viper <- as.matrix(acc.matrix[,viper])
acc.matrix.viper.bind <- rbind(acc.matrix.viper.tissue,acc.matrix.viper)

plot_matrix <- cbind(acc.matrix.ds
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.gs
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.msig
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.msig.dt
                     ,rep(0,dim(acc.matrix.ds)[1])
                     ,acc.matrix.viper.bind)


alternating <- c(1,3,2,4)
plot_matrix <- plot_matrix[alternating,]

collabels <- c(paste0(compounds.sorted.old,"\nBlood"),paste0(compounds.sorted.old,"\nGeneral"))
collabels <- collabels[alternating]


pdf(file="~/MVL/SingleViews/single.views_tissue_2vertical.pdf",width=8,height=3.5)

collabels <- c(paste0(compounds.sorted.old," Blood"),paste0(compounds.sorted.old," General"))
collabels <- collabels[alternating]

breaks.vec <- c(0,seq(0.5,1,length.out=50))
colfunc <- colorRampPalette(c("white","darkblue"))

# data specific
heatmap.2(plot_matrix
          ,dendrogram="none"
          ,Rowv=F
          ,Colv=F
          ,na.rm=T
          ,breaks=sort(breaks.vec)
          ,col=colfunc(length(breaks.vec)-1)
          
          ,na.color="grey"
          
          ,labCol=c("Clinical","Expression","Mutation",""
                    ,"Metabolic Enzymes","Multi-Drug Resistant","Drug Targets","Chromatin Modifying","Essential Genes",""
                    ,"Hallmark GS Expr","Motif GS Expr","TF Targets Expr","Positional GS Expr","Oncogenic GS Expr","Immunologic GS Expr",""
                    ,"DT GS Hallmark","DT GS Oncogenic","DT GS Immunologic",""
                    ,"Viper")
          ,labRow=collabels
          
          ,cexCol=1.1
          ,cexRow=1.1
          
          ,key=T
          ,density.info="none"
          ,key.xlab="Accuracy"
          ,key.title=NA
          ,key.par=list(mar=c(4,2,1,11))
          
          ,margins=c(10,10)
          ,trace="none"
          ,sepwidth=c(0.01,0.01)
          ,sepcolor="black"
          ,rowsep=c(0,nrow(plot_matrix))
          ,colsep=c(0,3,4,9,10,16,16,17,20,21,22)
          
          ,lmat = rbind(c(0,3),c(2,1),c(0,4))
          ,lwid = c(0.1,4)
          ,lhei = c(0.05,2.625,1)
)
dev.off()

