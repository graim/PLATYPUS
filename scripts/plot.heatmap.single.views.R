require(gplots)

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
compounds.sorted.oldID <- c(TKI
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

acc.matrix <- read.table("~/MVL/SingleViews/2015-08-27_MethodAccuracy.txt", sep='\t',header=T,row.names=1)
acc.matrix <- acc.matrix[compounds.sorted,]


#####
# ALL Together
ds.views <- c("ClinicalAggregationBinned","AllExpression","AllCNV","AllMutation")
acc.matrix.ds.views <- as.matrix(acc.matrix[,ds.views])

gs.views <- c("MetabolicEnzymes","MultiDrugResistant","DrugTarget","ChromatinModifyingEnzymes","DruggableGenes","GeneEssentialityAllCellLines")
acc.matrix.gs.views <- as.matrix(acc.matrix[,gs.views])
DT_mut <- read.table(file="~/MVL/SingleViews/coefficients/data_views/DrugTargets/Mutations/2015-08-18_SampleFeatureMutation_randomForest_BestParametersReRun.txt"
                     ,header=T, row.names=1)
dt_mut_acc <- DT_mut[compounds.sorted,"accuracy"]

acc.matrix.gs.views.mut <- as.matrix(cbind(acc.matrix.gs.views[,1:3],dt_mut_acc,acc.matrix.gs.views[,4:6]))

msig.views <- c("HallmarkAll","MotifAll","TranscriptionAll","PositionalAll","OncogenicAll","ImmunologicMedianVariance")
acc.matrix.msig.views <- as.matrix(acc.matrix[,msig.views])
msig.mut.views <- c("HallmarkMutations","MotifMutations","TranscritionMutations","PositionalMutations","OncogenicMutations","ImmunogenicMutations")
acc.matrix.msig.mut.views <- as.matrix(acc.matrix[,msig.mut.views])
msig.dt.views <- c("DrugTargetPathwayHallmark","DrugTargetPathwayOncogenic","DrugTargetPathwayImmunologic")
acc.matrix.msig.dt.views <- as.matrix(acc.matrix[,msig.dt.views])
viper <- c("Viper")
acc.matrix.viper <- as.matrix(acc.matrix[,viper])

simple.pred <- read.table(file="~/CCLE/simplePredictors/annotated_target_predictor_accuracy.tab", sep='\t',header=T,row.names=1)
simple.pred <- simple.pred[compounds.sorted.oldID,]
rownames(simple.pred) <- compounds.sorted
max_accuracy <- apply(simple.pred,1,max)
simple.pred <- cbind(simple.pred,max_accuracy)
simple.pred <- as.matrix(simple.pred[,3,drop=F])


pdf(file="~/MVL/SingleViews/single.views_simplepred.pdf",width=9,height=12)


breaks.vec <- c(0,seq(0.5,1,length.out=50))
colfunc <- colorRampPalette(c("white","darkblue"))

# data specific
heatmap.2(t(cbind(simple.pred
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.ds.views
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.gs.views.mut
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.msig.views
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.msig.mut.views
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.msig.dt.views
                  ,rep(0,dim(acc.matrix)[1])
                  ,acc.matrix.viper))
          ,dendrogram="none"
          ,Rowv=F
          ,Colv=F
          ,na.rm=T
          ,breaks=sort(breaks.vec)
          ,col=colfunc(length(breaks.vec)-1)
          
          ,na.color="grey"
          
          ,labRow=c("Annotated Target Mut",""
                    ,"Clinical","Expression","CNV","Mutation",""
                    ,"Metabolic Enzymes","Multi-Drug Resistant","Drug Targets","Drug Targets Mut","Chromatin Modifying","Druggable Genes","Essential Genes",""
                    ,"Hallmark GS Expr","Motif GS Expr","TF Targets Expr","Positional GS Expr","Oncogenic GS Expr","Immunologic GS Expr",""
                    ,"Hallmark GS Mut","Motif GS Mut","TF Targets Mut","Positional GS Mut","Oncogenic GS Mut","Immunologic GS Mut",""
                    ,"DT GS Hallmark","DT GS Oncogenic","DT GS Immunologic",""
                    ,"Viper")
          ,labCol=compounds.sorted.oldID
          ,cexCol=1.1
          ,cexRow=1.1
          
          ,key=T
          ,density.info="none"
          #,keysize = 1.5
          ,key.xlab="Accuracy"
          ,key.title=NA
          ,key.par=list(mar=c(4,2,1,11))
          
          ,margins=c(6,11)
          ,trace="none"
          ,sepwidth=c(0.01,0.01)
          ,sepcolor="black"
          ,colsep=c(0,ncol(t(acc.matrix.ds.views)))
          ,rowsep=c(0,1,c(0,4,5)+2,c(11,12,18,19,25,26,29,30,31)+3)
          
          ,lmat = rbind(c(0,3),c(2,1),c(0,4))
          ,lwid = c(0.2,4)
          ,lhei = c(0.05,4,0.4)
)
dev.off()
