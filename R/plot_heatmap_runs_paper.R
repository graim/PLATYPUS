
#######################
### sorted compounds
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

root <- c("~/MVL/MVLruns/outfiles/")

plot <- c("~/MVL/MVLruns/plots/")
no.view.vec <- c("3","5","7","10")
runs.vec <- paste0(no.view.vec,"views")



matrix.list <- list()
metric.vec <- c("balanced.accuracy.all.agree","no.correctly.predicted.balanced.all","balanced.accuracy.majority.agree","no.correctly.predicted.balanced.majority")
for(metric in metric.vec){
  platypus.first.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  platypus.final.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  platypus.best.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  
  for(compound in compounds.sorted){
    
    for(run in runs.vec){
      dir <- paste0(root,compound,"_",run,"/")
      filename <- paste0(dir,"perf_platypus_expanded.tab")
      
      if(file.exists(filename)){
        accuracy.platypus.iterations <- read.table(filename, sep='\t',header=T)
        
        cv.folds <- max(accuracy.platypus.iterations[,"cv.fold"])
        
        ## calculate no. correct predictions from accuracy and coverage for all.agree and majority.agree
        no.correctly.predicted.balanced.all <- accuracy.platypus.iterations[,"coverage.all.agree"] * accuracy.platypus.iterations[,"balanced.accuracy.all.agree"]
        accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.all)
        
        no.correctly.predicted.balanced.majority <- accuracy.platypus.iterations[,"coverage.majority.agree"] * accuracy.platypus.iterations[,"balanced.accuracy.majority.agree"]
        accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.majority)
        
        ## extend accuracy data table to maximal number of iterations
        no.iterations <- max(accuracy.platypus.iterations[,"iteration"])
        
        data.extended <- c()  
        for(cv in 1:cv.folds){
          
          iterations.in.fold <- max(accuracy.platypus.iterations[accuracy.platypus.iterations[,"cv.fold"]==cv,"iteration"])
          
          accuracy.platypus.iterations.cvfold <- accuracy.platypus.iterations[which(accuracy.platypus.iterations[,"cv.fold"] == cv),]
          # fill in the empty iterations with the result from the last iteration
          if(iterations.in.fold < no.iterations){
            for(i in (iterations.in.fold + 1):no.iterations){
              extended <- accuracy.platypus.iterations.cvfold[dim(accuracy.platypus.iterations.cvfold)[[1]],]
              extended[dim(extended)[[1]],"iteration"] <- i
              accuracy.platypus.iterations.cvfold <- rbind(accuracy.platypus.iterations.cvfold,extended)
            }
          }
          data.extended <- rbind(data.extended,accuracy.platypus.iterations.cvfold)
        }
        
        mean.vec <- c()
        for(i in 1:no.iterations){
          mean.vec <- c(mean.vec, mean(data.extended[data.extended[,"iteration"] == i,metric]))
        }
        
        platypus.first.mat[compound,run] <- mean.vec[1]
        platypus.final.mat[compound,run] <- mean.vec[no.iterations]
        platypus.best.mat[compound,run] <- max(mean.vec) 
      }
    }
  }
  

  matrix.list[[metric]] <- list()
  matrix.list[[metric]][["platypus.first"]] <- platypus.first.mat
  matrix.list[[metric]][["platypus.final"]] <- platypus.final.mat
  matrix.list[[metric]][["platypus.best"]] <- platypus.best.mat
  
  
  
  
  ##############
  # single heatmaps
  #############
  
  if(metric == "balanced.accuracy.all.agree"){
    keyname = "MVL All Agree [Balanced Accuracy]"
    breaks.vec <- c(0,seq(0.5,1,length.out=50))
  } else{
    keyname = "MVL All Agree [# correct]"
    breaks.vec <- c(0,seq(0.25,1,length.out=50))
  }
  
  if(metric == "balanced.accuracy.majority.agree"){
    keyname = "MVL Majority Agree [Balanced Accuracy]"
    breaks.vec <- c(0,seq(0.5,1,length.out=50))
  }
  if(metric == "no.correctly.predicted.balanced.majority"){
    keyname = "MVL Majority Agree [# correct]"
    breaks.vec <- c(0,seq(0.25,1,length.out=50))
  }
  
  single.mat <- cbind(platypus.first.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.best.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.final.mat)
 
  pdf(file=paste0(plot,"single_",metric,".pdf"),width=7,height=5)
  
  colfunc <- colorRampPalette(c("white","darkblue"))
  heatmap.2(t(single.mat)
            ,dendrogram="none"
            ,Rowv=F
            ,Colv=F
            ,na.rm=T
            ,breaks=sort(breaks.vec)
            ,col=colfunc(length(breaks.vec)-1)
            
            ,labCol=NA
            
            ,na.color="grey"
            
            ,ylab="# Views                                                             "
            ,labRow=c("Ensemble  3","Ensemble  5","Ensemble  7","Ensemble 10",""
                      ,"Best  3","Best  5","Best  7","Best 10",""
                      ,"Final  3","Final  5","Final  7","Final 10")
            ,cexRow=1
            
            ,key=T
            ,density.info="none"
            #,keysize = 1.5
            ,key.xlab=keyname
            ,key.title=NA
            ,key.par=list(mar=c(4,2,1,10))
            
            ,margins=c(0.1,9)
            ,trace="none"
            ,sepwidth=c(0.01,0.01)
            ,sepcolor="black"
            ,colsep=c(0,ncol(t(single.mat)))
            ,rowsep=c(0,4,5,9,10,nrow(t(single.mat)))
            
            ,lmat = rbind(c(0,4),c(2,1),c(0,3))
            ,lwid = c(0.1,4)
            ,lhei = c(1.5,5,0.1)
  )
  dev.off()
  
  #################
  # comparison heatmaps
  #################

  
  platypus.final.first.mat <- platypus.final.mat - platypus.first.mat  
  platypus.best.first.mat <- platypus.best.mat - platypus.first.mat
  
  comp_mat <- cbind(platypus.best.first.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.final.first.mat)
  
  pdf(file=paste0(plot,"comp_",metric,".pdf"),width=7,height=4.22)
  
  breaks.vec <- seq(-0.5,0.5,length.out=50)
  colfunc <- colorRampPalette(c("magenta3","white","green3"))
  heatmap.2(t(comp_mat)
            ,dendrogram="none"
            ,Rowv=F
            ,Colv=F
            ,na.rm=T
            ,breaks=sort(breaks.vec)
            ,col=colfunc(length(breaks.vec)-1)
            
            ,na.color="grey"
            
            ,ylab=" "
            ,labRow=c("Best - Ensemble  3","Best - Ensemble  5","Best - Ensemble  7","Best - Ensemble 10",""
                      ,"Final - Ensemble  3","Final - Ensemble  5","Final - Ensemble  7","Final - Ensemble 10")
            ,cexRow=1
            
            ,key=T
            ,density.info="none"
            ,key.xlab=paste0("Diff( ",keyname," )")
            ,key.title=NA
            ,key.par=list(mar=c(4,2,1,10))
            
            ,margins=c(5,9)
            ,trace="none"
            ,sepwidth=c(0.01,0.01)
            ,sepcolor="black"
            ,colsep=c(0,ncol(t(comp_mat)))
            ,rowsep=c(0,4,5,nrow(t(comp_mat)))
            
            ,lmat = rbind(c(0,3),c(2,1),c(0,4))
            ,lwid = c(0.1,4)
            ,lhei = c(0.1,4,1.5)
  )
  dev.off()
  
}





############
# DATA specific


root <- c("~/MVL/MVLruns_dataspecific/outfiles/")

plot <- c("~/MVL/MVLruns_dataspecific/plots/")
no.view.vec <- c("3")
runs.vec <- paste0(no.view.vec,"views")



metric.vec <- c("balanced.accuracy.all.agree","no.correctly.predicted.balanced.all","balanced.accuracy.majority.agree","no.correctly.predicted.balanced.majority")
for(metric in metric.vec){
  platypus.first.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  platypus.final.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  platypus.best.mat <- matrix(data=NA,nrow=length(compounds.sorted),ncol=length(runs.vec),dimnames=list(compounds.sorted,runs.vec))
  
  for(compound in compounds.sorted){
    
    for(run in runs.vec){
      dir <- paste0(root,compound,"_",run,"/")
      filename <- paste0(dir,"perf_platypus_expanded.tab")
      
      if(file.exists(filename)){
        accuracy.platypus.iterations <- read.table(filename, sep='\t',header=T)
        
        cv.folds <- max(accuracy.platypus.iterations[,"cv.fold"])
        
        ## calculate no. correct predictions from accuracy and coverage for all.agree and majority.agree
        no.correctly.predicted.balanced.all <- accuracy.platypus.iterations[,"coverage.all.agree"] * accuracy.platypus.iterations[,"balanced.accuracy.all.agree"]
        accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.all)
        
        no.correctly.predicted.balanced.majority <- accuracy.platypus.iterations[,"coverage.majority.agree"] * accuracy.platypus.iterations[,"balanced.accuracy.majority.agree"]
        accuracy.platypus.iterations <- cbind(accuracy.platypus.iterations,no.correctly.predicted.balanced.majority)
        
        ## extend accuracy data table to maximal number of iterations
        no.iterations <- max(accuracy.platypus.iterations[,"iteration"])
        
        data.extended <- c()  
        for(cv in 1:cv.folds){
          
          iterations.in.fold <- max(accuracy.platypus.iterations[accuracy.platypus.iterations[,"cv.fold"]==cv,"iteration"])
          
          accuracy.platypus.iterations.cvfold <- accuracy.platypus.iterations[which(accuracy.platypus.iterations[,"cv.fold"] == cv),]
          # fill in the empty iterations with the result from the last iteration
          if(iterations.in.fold < no.iterations){
            for(i in (iterations.in.fold + 1):no.iterations){
              extended <- accuracy.platypus.iterations.cvfold[dim(accuracy.platypus.iterations.cvfold)[[1]],]
              extended[dim(extended)[[1]],"iteration"] <- i
              accuracy.platypus.iterations.cvfold <- rbind(accuracy.platypus.iterations.cvfold,extended)
            }
          }
          data.extended <- rbind(data.extended,accuracy.platypus.iterations.cvfold)
        }
        
        mean.vec <- c()
        for(i in 1:no.iterations){
          mean.vec <- c(mean.vec, mean(data.extended[data.extended[,"iteration"] == i,metric]))
        }
        
        platypus.first.mat[compound,run] <- mean.vec[1]
        platypus.final.mat[compound,run] <- mean.vec[no.iterations]
        platypus.best.mat[compound,run] <- max(mean.vec) 
      }
    }
  }
  
  
  platypus.first.mat <- cbind(platypus.first.mat,matrix.list[[metric]][["platypus.first"]])
  platypus.final.mat <- cbind(platypus.final.mat,matrix.list[[metric]][["platypus.final"]])
  platypus.best.mat <- cbind(platypus.best.mat,matrix.list[[metric]][["platypus.best"]])
  
  ##############
  # single heatmaps
  #############
  
  if(metric == "balanced.accuracy.all.agree"){
    keyname = "MVL All Agree [Balanced Accuracy]"
    breaks.vec <- c(0,seq(0.5,1,length.out=50))
  } else{
    keyname = "MVL All Agree [# correct]"
    breaks.vec <- c(0,seq(0.25,1,length.out=50))
  }
  
  if(metric == "balanced.accuracy.majority.agree"){
    keyname = "MVL Majority Agree [Balanced Accuracy]"
    breaks.vec <- c(0,seq(0.5,1,length.out=50))
  }
  if(metric == "no.correctly.predicted.balanced.majority"){
    keyname = "MVL Majority Agree [# correct]"
    breaks.vec <- c(0,seq(0.25,1,length.out=50))
  }
  
  single.mat <- cbind(platypus.first.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.best.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.final.mat)
  
  pdf(file=paste0(plot,"together_single_",metric,".pdf"),width=7,height=5)
  
  
  
  colfunc <- colorRampPalette(c("white","darkblue"))
  heatmap.2(t(single.mat)
            ,dendrogram="none"
            ,Rowv=F
            ,Colv=F
            ,na.rm=T
            ,breaks=sort(breaks.vec)
            ,col=colfunc(length(breaks.vec)-1)
            
            ,labCol=NA
            
            ,na.color="grey"
            
            ,ylab="# Views                                                             "
            ,labRow=c("Ensemble ds","Ensemble  3","Ensemble  5","Ensemble  7","Ensemble 10",""
                      ,"Best ds","Best  3","Best  5","Best  7","Best 10",""
                      ,"Final ds","Final  3","Final  5","Final  7","Final 10")
            ,cexRow=1
            
            ,key=T
            ,density.info="none"
            ,key.xlab=keyname
            ,key.title=NA
            ,key.par=list(mar=c(4,2,1,10))
            
            ,margins=c(0.1,9)
            ,trace="none"
            ,sepwidth=c(0.01,0.01)
            ,sepcolor="black"
            ,colsep=c(0,ncol(t(single.mat)))
            ,rowsep=c(0,5,6,11,12,nrow(t(single.mat)))
            
            ,lmat = rbind(c(0,4),c(2,1),c(0,3))
            ,lwid = c(0.1,4)
            ,lhei = c(1.5,5,0.1)
  )
  dev.off()
  
  #################
  # comparison heatmaps
  #################
  
  
  platypus.final.first.mat <- platypus.final.mat - platypus.first.mat  
  platypus.best.first.mat <- platypus.best.mat - platypus.first.mat
  
  comp_mat <- cbind(platypus.best.first.mat,rep(0,dim(platypus.first.mat)[[1]]),platypus.final.first.mat)
  
  pdf(file=paste0(plot,"together_comp_",metric,".pdf"),width=7,height=4.22)
  
  
  breaks.vec <- seq(-0.5,0.5,length.out=50)
  colfunc <- colorRampPalette(c("magenta3","white","green3"))
  heatmap.2(t(comp_mat)
            ,dendrogram="none"
            ,Rowv=F
            ,Colv=F
            ,na.rm=T
            ,breaks=sort(breaks.vec)
            ,col=colfunc(length(breaks.vec)-1)
            
            ,na.color="grey"
            
            ,ylab=" "
            ,labRow=c("Best - Ensemble  ds","Best - Ensemble  3","Best - Ensemble  5","Best - Ensemble  7","Best - Ensemble 10",""
                      ,"Final - Ensemble  ds","Final - Ensemble  3","Final - Ensemble  5","Final - Ensemble  7","Final - Ensemble 10")
            ,cexRow=1
            ,labCol=compounds.sorted.old
            
            ,key=T
            ,density.info="none"
            ,key.xlab=paste0("Diff( ",keyname," )")
            ,key.title=NA
            ,key.par=list(mar=c(4,2,1,10))
            
            ,margins=c(5,9)
            ,trace="none"
            ,sepwidth=c(0.01,0.01)
            ,sepcolor="black"
            ,colsep=c(0,ncol(t(comp_mat)))
            ,rowsep=c(0,5,6,nrow(t(comp_mat)))
            
            ,lmat = rbind(c(0,3),c(2,1),c(0,4))
            ,lwid = c(0.1,4)
            ,lhei = c(0.1,4,1.5)
  )
  dev.off()
  
}
