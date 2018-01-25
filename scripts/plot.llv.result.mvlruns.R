
addX <- function(vector){
  vector[grep("^[0-9].",vector)] <- paste0("X",vector[grep("^[0-9].",vector)])
  return(vector)
}

require(MASS)

classcol.labs <- 1

plot.single.llvfolds <- FALSE
#plot.single.cvfolds <- TRUE



###################
## sorted compounds

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

plot <- c("~/MVL/MVLruns/plots/llv/")
no.view.vec <- c("3","5","7","10")
runs.vec <- paste0(no.view.vec,"views")




for(i.compound in 1:length(compounds.sorted)){
  compound <- compounds.sorted[i.compound]
  compound.o <- compounds.sorted.old[i.compound]
  
  for(run in runs.vec){
    folder <- paste0(root,compound,"_",run,"/")
    pdf.name <- paste0(compound,"_",run)
    fn.labs <- paste0("~/CCLE/DATA/DRUG_RESPONSE/CCLE_responselabel_binary_",compound.o,".tab")
    
    if(file.exists(paste0(folder,"perf_platypus_expanded.tab")) & file.exists(paste0(folder,"perf_llv_expanded.tab"))){
    
    ###
    # FROM CV OUTPUT! to get the weights
    ###
    accuracy.platypus.iterations <- read.table( paste0(folder,"perf_platypus_expanded.tab"), sep='\t',header=T)
    # only deals with starting accuracy
    acc.norm.views <- accuracy.platypus.iterations[1, grep('^weighting.norm.view.', colnames(accuracy.platypus.iterations))] 
    ####
    
    accuracy.llvfolds <- read.table(paste0(folder,"perf_llv.tab"), sep='\t',header=T)
    accuracy.platypus.iterations <- read.table(paste0(folder,"perf_llv_expanded.tab"), sep='\t',header=T)
    
    load(paste0(folder,"labelling.matrix.llvlist.Rdata"))
    load(paste0(folder,"labelling.matrices.views.llvlist.Rdata"))
    
    llv.folds <- max(accuracy.llvfolds[,"llv.fold"])
    no.views <- length(labelling.matrices.views.llvlist[[1]])
    
    
    # read in labels
    
    labs <- read.table(fn.labs, sep='\t',header=T, row.names=1, check.names=F, stringsAsFactors = FALSE)
    # take out NA values
    labs <- labs[which(!(is.na(labs[,classcol.labs]))),,drop=F]
    rownames(labs) <- addX(rownames(labs))

    
    # reduce labelling matrices to known labels
    for(llv in 1:llv.folds){
      ids <- intersect(rownames(labelling.matrix.llvlist[[llv]]),rownames(labs))
      labelling.matrix.llvlist[[llv]] <- labelling.matrix.llvlist[[llv]][ids,]
      for(view.i in 1:no.views){
        labelling.matrices.views.llvlist[[llv]][[view.i]] <- labelling.matrices.views.llvlist[[llv]][[view.i]][ids,]
      }
      
    }
    
    # extend labelling matrices like the accuracy table is extended (together with data.extended)
    ### HERE: because of the maximum approach, there is no clear threshold jumps any more. Just extend the matrices at the end
    no.iterations <- max(accuracy.platypus.iterations[,"iteration"])
    data.extended <- c()
    for(llv in 1:llv.folds){
      iterations.in.fold <- max(accuracy.platypus.iterations[accuracy.platypus.iterations[,"llv.fold"]==llv,"iteration"])
      # reduce the labelling matrices to the maximum number of iterations
      labelling.matrix.llvlist[[llv]] <- labelling.matrix.llvlist[[llv]][,1:no.iterations]
      for(view.i in 1:no.views){
        labelling.matrices.views.llvlist[[llv]][[view.i]] <- labelling.matrices.views.llvlist[[llv]][[view.i]][,1:no.iterations]
      }
      accuracy.platypus.iterations.llvfold <- accuracy.platypus.iterations[which(accuracy.platypus.iterations[,"llv.fold"] == llv),]
      # fill in the empty iterations with the result from the last iteration
      if(iterations.in.fold < no.iterations){
        for(i in (iterations.in.fold + 1):no.iterations){
          labelling.matrix.llvlist[[llv]][,i] <- labelling.matrix.llvlist[[llv]][,i-1]
          for(view.i in 1:no.views){
            labelling.matrices.views.llvlist[[llv]][[view.i]][,i] <- labelling.matrices.views.llvlist[[llv]][[view.i]][,i-1]
          }
          accuracy.platypus.iterations.llvfold <- rbind(accuracy.platypus.iterations.llvfold,accuracy.platypus.iterations.llvfold[dim(accuracy.platypus.iterations.llvfold)[[1]],])
        }
      }
      data.extended <- rbind(data.extended,accuracy.platypus.iterations.llvfold)
    }
    
    
    
    ## The matrices for single llv-folds can be put together now
    labelling.matrix.complete <- c()
    for(llv in 1:llv.folds){
      labelling.matrix.complete <- rbind(labelling.matrix.complete,labelling.matrix.llvlist[[llv]])
    }
    
    labelling.matrices.views.complete <- list()
    for(view.i in 1:no.views){
      labelling.matrix.view <- c()
      for(llv in 1:llv.folds){
        labelling.matrix.view <- rbind(labelling.matrix.view,labelling.matrices.views.llvlist[[llv]][[view.i]])
      }
      labelling.matrices.views.complete[[view.i]] <- labelling.matrix.view
    }
    
    
    
    # accuracy of MVL labelling
    ll.acc.avg <- c()
    ll.acc.sd <- c()
    for(i in 1:no.iterations){
      acc.vec <- c()
      for(llv in 1:llv.folds){
        acc.vec <- c(acc.vec,data.extended[which(data.extended[,"llv.fold"] == llv),"accuracy"][i])
      }
      ll.acc.avg <- c(ll.acc.avg, mean(acc.vec))
      ll.acc.sd <- c(ll.acc.sd, sd(acc.vec))
    }
    
    # coverage of MVL labelling
    ll.cov.avg <- c()
    ll.cov.sd <- c()
    for(i in 1:no.iterations){
      cov.vec <- c()
      for(llv in 1:llv.folds){
        cov.vec <- c(cov.vec,data.extended[which(data.extended[,"llv.fold"] == llv),"coverage"][i])
      }
      ll.cov.avg <- c(ll.cov.avg, mean(cov.vec))
      ll.cov.sd <- c(ll.cov.sd, sd(cov.vec))
    }
    
    
    # thresholds of MVL labelling
    ll.upper.avg <- c()
    ll.upper.sd <- c()
    for(i in 1:no.iterations){
      upper.vec <- c()
      for(llv in 1:llv.folds){
        upper.vec <- c(upper.vec,data.extended[which(data.extended[,"llv.fold"] == llv),"weighting.threshold.upper"][i])
      }
      ll.upper.avg <- c(ll.upper.avg, mean(upper.vec))
      ll.upper.sd <- c(ll.upper.sd, sd(upper.vec))
    }
    
    ll.lower.avg <- c()
    ll.lower.sd <- c()
    for(i in 1:no.iterations){
      lower.vec <- c()
      for(llv in 1:llv.folds){
        lower.vec <- c(lower.vec,data.extended[which(data.extended[,"llv.fold"] == llv),"weighting.threshold.lower"][i])
      }
      ll.lower.avg <- c(ll.lower.avg, mean(lower.vec))
      ll.lower.sd <- c(ll.lower.sd, sd(lower.vec))
    }
    
    
    #############
    # PLOTTING
    #############
    
    
    ##########################
    ## sample labelling plot
    require(gplots)
    
    unique.labels <- unique(as.vector(labelling.matrix.complete))
    unique.labels <- unique.labels[!is.na(unique.labels)]
    unique.labels <- sort(unique.labels,decreasing=T)
    
    numeric.labelling.matrix <- labelling.matrix.complete
    
    numeric.labelling.matrix[which(numeric.labelling.matrix==unique.labels[1])] <- 1
    numeric.labelling.matrix[which(numeric.labelling.matrix==unique.labels[2])] <- -1
    class(numeric.labelling.matrix) <- "numeric"
    
    
    ## get if the labelling is correct
    ids <- intersect(rownames(labelling.matrix.complete),rownames(labs))
    correct.labelling.vec <- labelling.matrix.complete[ids,dim(labelling.matrix.complete)[[2]]] == labs[ids,1]
    
    correct.labelling.vec.num <- correct.labelling.vec
    correct.labelling.vec.num[which(correct.labelling.vec)] <- 1
    correct.labelling.vec.num[which(!correct.labelling.vec)] <- -1
    
    
    correct.labelling.vec.col <- correct.labelling.vec
    correct.labelling.vec.col[which(correct.labelling.vec)] <- "green"
    correct.labelling.vec.col[which(!correct.labelling.vec)] <- "black" 
    
    ## create the matrix basis for the heatmap
    view.label.matrix <- numeric.labelling.matrix
    view.label.matrix[] <- NA
    not.enough.data <- c()
    for(cell in 1:length(view.label.matrix)){
      count1 <- 0
      count2 <- 0
      countNA <- 0
      for(view.i in 1:no.views){
        if(is.na(labelling.matrices.views.complete[[view.i]][cell])){
          countNA <- countNA + 1
        } else{
          if(labelling.matrices.views.complete[[view.i]][cell] == unique.labels[1]){
            count1 <- count1 + unlist(acc.norm.views[view.i])
          }
          if(labelling.matrices.views.complete[[view.i]][cell] == unique.labels[2]){
            count2 <- count2 + unlist(acc.norm.views[view.i])
          }
        }
      }
      if(!(countNA == no.views)){
        if(count1 == count2){
          view.label.matrix[cell] <- 0
        }
        if(count1 > count2){
          view.label.matrix[cell] <- count1
        }
        if(count2 > count1){
          view.label.matrix[cell] <- -count2
        }
      }
      
      # add the labelling information
      if(!(is.na(numeric.labelling.matrix[cell]))){
        if(is.na(view.label.matrix[cell])){
          view.label.matrix[cell] <- 2*numeric.labelling.matrix[cell]*sum(acc.norm.views)
        }
      }

    }
    
    
    not.enough.data <- not.enough.data[not.enough.data<= dim(view.label.matrix)[[1]]]
    exclude.rows <- rownames(view.label.matrix[not.enough.data,])
    
    
    ## color palette / breaks depend on no.views
    breaks.vec <- c()
    for(n in no.views:3){
      for(z in (n-1):2){
        b <- z/n
        if(b > 0.5 & !(b %in% breaks.vec)){
          breaks.vec <- c(breaks.vec,b)
        }
      }
    }
    breaks.vec <- sort(breaks.vec)
    breaks.vec.for.key <- as.character(fractions(breaks.vec))
    
    breaks.vec <- c(-2,-1,-(rev(breaks.vec)),0,breaks.vec,1,2)
    breaks.vec.inthemiddle <- c()
    for(b in 1:length(breaks.vec)-1){
      m <- (breaks.vec[b] + breaks.vec[b+1]) /2
      breaks.vec.inthemiddle <- c(breaks.vec.inthemiddle,m)
    }
    breaks.vec.inthemiddle <- c(-2,breaks.vec.inthemiddle,2)
    
    
    breaks.vec.inthemiddle.middle <- c()
    for(b in 1:length(breaks.vec.inthemiddle)-1){
      m <- (breaks.vec.inthemiddle[b] + breaks.vec.inthemiddle[b+1]) /2
      breaks.vec.inthemiddle.middle <- c(breaks.vec.inthemiddle.middle,m)
    }
    
    plot.matrix <- view.label.matrix
    plot.matrix <- plot.matrix[order(rowSums(cbind(plot.matrix,correct.labelling.vec.num),na.rm=T),decreasing=T),]
    
    
    breaks.vec <- unique(c(-(2*sum(acc.norm.views)),-(2*sum(acc.norm.views) - 0.5 * sum(acc.norm.views) )
                           , seq(-sum(acc.norm.views),0,length.out=20),  seq(0,sum(acc.norm.views),length.out=20),(2*sum(acc.norm.views) - 0.5 * sum(acc.norm.views) ), (2*sum(acc.norm.views)) ))
    label.vec.key <- unique(c( seq(-round(sum(acc.norm.views)),0,length.out=round(sum(acc.norm.views))+1),  seq(0,round(sum(acc.norm.views)),length.out=round(sum(acc.norm.views))+1) ))
    
    if(pdf.name == "PD.0325901_10views" | pdf.name == "Panobinostat_10views"){
      pdf(file=paste0(plot,"llv_labelling_",pdf.name,".pdf"),width=10,height=10)
      
      colfunc <- colorRampPalette(c("blue3","thistle1","red3"))
      heatmap.2(
        plot.matrix
        ,dendrogram="none"
        ,Rowv=F
        ,Colv=F
        ,na.rm=T
        ,breaks=breaks.vec
        ,col=c("blue4",colfunc(length(breaks.vec)-3),"red4")
        
        ,key=T
        ,density.info="none"
        ,key.xtickfun=function() {
          
          label.vec <- c("added","all.agree",c("","8","","6","","4","","2","","0","","2","","4","","6","","8",""),"all.agree","added")
          
          breaks.vec.mod <- c((breaks.vec[1]+breaks.vec[2]) /2 , (breaks.vec[2] + breaks.vec[3] ) /2
                              ,label.vec.key
                              , -((breaks.vec[2] + breaks.vec[3] ) /2), -((breaks.vec[1]+breaks.vec[2]) /2)
          ) + -breaks.vec[1]
          breaks.vec.mod <-  breaks.vec.mod / ( breaks.vec.mod[length(breaks.vec.mod)] + -breaks.vec[1] -  -((breaks.vec[1]+breaks.vec[2]) /2) )
          return(list(at=breaks.vec.mod,labels=label.vec,tick=T))
          
        }
        ,keysize = 1.5
        ,key.xlab=paste0(unique.labels[2],"   <-                 Predictor Agreement [SUM(votes)]                 ->  ",unique.labels[1])
        ,key.title=NA
        ,key.par=list(mar=c(4, 5, 1, 3) + 0.1)
        
        ,lmat = rbind(c(3,0,5),c(4,1,2))
        ,lwid = c(1,0.25,6)
        ,lhei = c(0.5,4)
        
        ,xlab="Iterations"
        ,ylab="Samples"
        ,labRow = NA,
        ,labCol = NA,
        ,margins=c(3,3)
        ,trace="none"
        ,RowSideColors=correct.labelling.vec.col[rownames(plot.matrix)]
      )
      legend("topleft",inset=-0.01,legend=c("Correctly Labelled","Incorrectly Labelled"),fill=c("green","black")
             ,border=NA,bty='n',cex=0.75)
      dev.off()
    } else{
      pdf(file=paste0(plot,"llv_labelling_",pdf.name,".pdf"),width=10,height=10)
      
      colfunc <- colorRampPalette(c("blue3","thistle1","red3"))
      heatmap.2(
        plot.matrix
        ,dendrogram="none"
        ,Rowv=F
        ,Colv=F
        ,na.rm=T
        ,breaks=breaks.vec
        ,col=c("blue4",colfunc(length(breaks.vec)-3),"red4")
        
        ,key=T
        ,density.info="none"
        ,key.xtickfun=function() {
          
          label.vec <- c("added","all.agree",abs(label.vec.key),"all.agree","added")
          
          breaks.vec.mod <- c((breaks.vec[1]+breaks.vec[2]) /2 , (breaks.vec[2] + breaks.vec[3] ) /2
                              ,label.vec.key
                              , -((breaks.vec[2] + breaks.vec[3] ) /2), -((breaks.vec[1]+breaks.vec[2]) /2)
          ) + -breaks.vec[1]
          breaks.vec.mod <-  breaks.vec.mod / ( breaks.vec.mod[length(breaks.vec.mod)] + -breaks.vec[1] -  -((breaks.vec[1]+breaks.vec[2]) /2) )
          return(list(at=breaks.vec.mod,labels=label.vec,tick=T))
          
        }
        ,keysize = 1.5
        ,key.xlab=paste0(unique.labels[2],"   <-                 Predictor Agreement [SUM(votes)]                 ->  ",unique.labels[1])
        ,key.title=NA
        ,key.par=list(mar=c(4, 5, 1, 3) + 0.1)
        
        ,lmat = rbind(c(3,0,5),c(4,1,2))
        ,lwid = c(1,0.25,6)
        ,lhei = c(0.5,4)
        
        ,xlab="Iterations"
        ,ylab="Samples"
        ,labRow = NA,
        ,labCol = NA,
        ,margins=c(3,3)
        ,trace="none"
        ,RowSideColors=correct.labelling.vec.col[rownames(plot.matrix)]
      )
      legend("topleft",inset=-0.01,legend=c("Correctly Labelled","Incorrectly Labelled"),fill=c("green","black")
             ,border=NA,bty='n',cex=0.75)
      dev.off()
    }
    
    
    
    }
    
  }
}


