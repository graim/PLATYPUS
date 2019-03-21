
plot.single.cvfolds <- FALSE
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


root <- c("~/MVL/MVLruns/outfiles/")

plot <- c("~/MVL/MVLruns/plots/cv/")
no.view.vec <- c("3","5","7","10")
runs.vec <- paste0(no.view.vec,"views")




  
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
              accuracy.platypus.iterations.cvfold <- rbind(accuracy.platypus.iterations.cvfold,accuracy.platypus.iterations.cvfold[dim(accuracy.platypus.iterations.cvfold)[[1]],])
            }
          }
          data.extended <- rbind(data.extended,accuracy.platypus.iterations.cvfold)
        }
        
        
        # accuracy of MVL
        platypus.acc.avg.all <- c()
        platypus.acc.avg.majority <- c()
        platypus.acc.sd.all <- c()
        platypus.acc.sd.majority <- c()
        for(i in 1:no.iterations){
          acc.vec.all <- c()
          acc.vec.majority <- c()
          for(cv in 1:cv.folds){
            acc.vec.all <- c(acc.vec.all,data.extended[which(data.extended[,"cv.fold"] == cv),"balanced.accuracy.all.agree"][i])
            acc.vec.majority <- c(acc.vec.majority,data.extended[which(data.extended[,"cv.fold"] == cv),"balanced.accuracy.majority.agree"][i])
          }
          platypus.acc.avg.all <- c(platypus.acc.avg.all, mean(acc.vec.all))
          platypus.acc.sd.all <- c(platypus.acc.sd.all, sd(acc.vec.all))
          platypus.acc.avg.majority <- c(platypus.acc.avg.majority, mean(acc.vec.majority))
          platypus.acc.sd.majority <- c(platypus.acc.sd.majority, sd(acc.vec.majority))
        }
        
        ## accuracy of views
        views.acc.avg <- c()
        views.acc.sd <- c()
        for(view.i in 1:length(grep("balanced.accuracy.view.", colnames(data.extended)))){ #TODO
          view.acc.avg <- c()
          view.acc.sd <- c()
          for(i in 1:no.iterations){
            acc.vec <- c()
            for(cv in 1:cv.folds){
              acc.vec <- c(acc.vec,data.extended[which(data.extended[,"cv.fold"] == cv),paste0("balanced.accuracy.view.",view.i)][i])
            }
            view.acc.avg <- c(view.acc.avg, mean(acc.vec))
            view.acc.sd <- c(view.acc.sd, sd(acc.vec))
          }
          views.acc.avg <- cbind(views.acc.avg,view.acc.avg)
          views.acc.sd <- cbind(views.acc.sd,view.acc.sd)
        }
        
        ## prediction coverage
        platypus.cov.avg.all <- c()
        platypus.cov.avg.majority <- c()
        platypus.cov.sd.all <- c()
        platypus.cov.sd.majority <- c()
        for(i in 1:no.iterations){
          cov.vec.all <- c()
          cov.vec.majority <- c()
          for(cv in 1:cv.folds){
            cov.vec.all <- c(cov.vec.all,data.extended[which(data.extended[,"cv.fold"] == cv),"coverage.all.agree"][i])
            cov.vec.majority <- c(cov.vec.majority,data.extended[which(data.extended[,"cv.fold"] == cv),"coverage.majority.agree"][i])
          }
          platypus.cov.avg.all <- c(platypus.cov.avg.all, mean(cov.vec.all))
          platypus.cov.sd.all <- c(platypus.cov.sd.all, sd(cov.vec.all))
          platypus.cov.avg.majority <- c(platypus.cov.avg.majority, mean(cov.vec.majority))
          platypus.cov.sd.majority <- c(platypus.cov.sd.majority, sd(cov.vec.majority))
        }
        
        ## no.correctly.predicted = acc*cov
        platypus.corr.avg.all <- c()
        platypus.corr.avg.majority <- c()
        platypus.corr.sd.all <- c()
        platypus.corr.sd.majority <- c()
        for(i in 1:no.iterations){
          corr.vec.all <- c()
          corr.vec.majority <- c()
          for(cv in 1:cv.folds){
            corr.vec.all <- c(corr.vec.all,data.extended[which(data.extended[,"cv.fold"] == cv),"no.correctly.predicted.balanced.all"][i])
            corr.vec.majority <- c(corr.vec.majority,data.extended[which(data.extended[,"cv.fold"] == cv),"no.correctly.predicted.balanced.majority"][i])
          }
          platypus.corr.avg.all <- c(platypus.corr.avg.all, mean(corr.vec.all))
          platypus.corr.sd.all <- c(platypus.corr.sd.all, sd(corr.vec.all))
          platypus.corr.avg.majority <- c(platypus.corr.avg.majority, mean(corr.vec.majority))
          platypus.corr.sd.majority <- c(platypus.corr.sd.majority, sd(corr.vec.majority))
        }
        
        
        ## inclusion of unlabelled data in training
        platypus.inclusion.avg <- c()
        platypus.inclusion.sd <- c()
        for(i in 1:no.iterations){
          inclusion.vec <- c()
          for(cv in 1:cv.folds){
            inclusion.vec <- c(inclusion.vec,data.extended[which(data.extended[,"cv.fold"] == cv),"no.ids.labelled"][i])
          }
          platypus.inclusion.avg <- c(platypus.inclusion.avg, mean(inclusion.vec))
          platypus.inclusion.sd <- c(platypus.inclusion.sd, sd(inclusion.vec))
        }
        
        
        # thresholds of MVL labelling
        ll.upper.avg <- c()
        ll.upper.sd <- c()
        for(i in 1:no.iterations){
          upper.vec <- c()
          for(cv in 1:cv.folds){
            upper.vec <- c(upper.vec,data.extended[which(data.extended[,"cv.fold"] == cv),"weighting.threshold.upper"][i])
          }
          ll.upper.avg <- c(ll.upper.avg, mean(upper.vec))
          ll.upper.sd <- c(ll.upper.sd, sd(upper.vec))
        }
        
        ll.lower.avg <- c()
        ll.lower.sd <- c()
        for(i in 1:no.iterations){
          lower.vec <- c()
          for(cv in 1:cv.folds){
            lower.vec <- c(lower.vec,data.extended[which(data.extended[,"cv.fold"] == cv),"weighting.threshold.lower"][i])
          }
          ll.lower.avg <- c(ll.lower.avg, mean(lower.vec))
          ll.lower.sd <- c(ll.lower.sd, sd(lower.vec))
        }
        
        
        #############
        # PLOTTING
        #############
        
        pdf(file=paste0(plot,"cv_performance_",compound,"_",run,".pdf"),width=8,height=7)
        
        green.transparent <- rgb(col2rgb("green")[1,1],col2rgb("green")[2,1],col2rgb("green")[3,1],100,maxColorValue=255)
        
        par(mar=c(4, 4, 0, 0) + 0.1,par(xaxs='i',yaxs='i')) 
        
        layout(as.matrix(rbind(1,2,3,4,5)), 
               widths=c(1,1,1,1,1), heights=c(1,1,1,1,1.25))
        
         par(mar=c(1, 4, 0.6, 0.1) + 0.1)
        ## prediction accuracy plot
        plot(1:no.iterations,platypus.acc.avg.all
             ,ylim=c(0.5,1)
             ,xlim=c(1,no.iterations)
             ,type='l'
             #,pch=20
             ,lwd=2
             ,col="red3"
             ,xlab=""
             ,ylab="CV Accuracy"
             ,xaxt='n'
        )
        grid()

        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"balanced.accuracy.all.agree"],type='l',col="tomato")
          }
        }
        
        
        #single views
        for(view.i in 1:dim(views.acc.avg)[[2]]){
          lines(1:no.iterations,views.acc.avg[,view.i]
                ,pch=20
                ,col=green.transparent
                ,lwd=2
          )
        }
        
        
        upper <- platypus.acc.avg.majority + platypus.acc.sd.majority
        lower <- platypus.acc.avg.majority - platypus.acc.sd.majority
        lightblue.transparent <- rgb(col2rgb("lightblue")[1,1],col2rgb("lightblue")[2,1],col2rgb("lightblue")[3,1],100,maxColorValue=255)
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = lightblue.transparent, border = NA)
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"balanced.accuracy.majority.agree"],type='l',col="lightblue")
          }
        }
        
       
        #red polygon again
        upper <- platypus.acc.avg.all + platypus.acc.sd.all
        lower <- platypus.acc.avg.all - platypus.acc.sd.all
        tomato.transparent <- rgb(col2rgb("tomato")[1,1],col2rgb("tomato")[2,1],col2rgb("tomato")[3,1],100,maxColorValue=255)
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = tomato.transparent, border = NA)
        
        # lines again
        lines(1:no.iterations,platypus.acc.avg.majority,col="blue",lwd=2)
        lines(1:no.iterations, platypus.acc.avg.all, col = 'red3',lwd=2)
        
        
        legend("bottomright",legend=c("all agree","majority agrees","single views"), lty=c(1,1,1),lwd=c(2,2,2), col=c("red3","blue","green"), bty='n')
        
        
        ## prediction coverage plot
         par(mar=c(1, 4, 0.6, 0.1) + 0.1)
        plot(1:no.iterations,platypus.cov.avg.all
             ,ylim=c(0,1)
             ,type='l'
             #,pch=20
             ,lwd=2
             ,col="red3"
             ,xlab=""
             ,ylab="CV Coverage"# \npredictions made on hold-out samples"
             ,xaxt='n'
        )
        grid()
        

        
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"coverage.all.agree"],type='l',col="tomato")
          }
        }
        
        
        upper <- platypus.cov.avg.majority + platypus.cov.sd.majority
        lower <- platypus.cov.avg.majority - platypus.cov.sd.majority
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = lightblue.transparent, border = NA)
        
        upper <- platypus.cov.avg.all + platypus.cov.sd.all
        lower <- platypus.cov.avg.all - platypus.cov.sd.all
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = tomato.transparent, border = NA)
        
        
        lines(1:no.iterations,platypus.cov.avg.majority
              ,type='l'
              ,lwd=2
              ,col="blue"
        )
        lines(1:no.iterations,platypus.cov.avg.all,col="red3",lwd=2)
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"coverage.majority.agree"],type='l',col="lightblue")
          }
        }
        
        legend("bottomright",legend=c("all agree","majority agrees"), lty=c(1,1),lwd=c(2,2), col=c("red3","blue"), bty='n')
        
        
        
        ## no.correctly.predicted plot
         par(mar=c(1, 4, 0.6, 0.1) + 0.1)
        plot(1:no.iterations,platypus.corr.avg.all
             ,ylim=c(0,1)
             ,type='l'
             ,lwd=2
             ,col="red3"
             ,xlab=""
             ,ylab="# Correct"#"Correctly Predicted Samples \nAccuracy * Coverage"
             ,xaxt='n'
        )
        grid()
        

        

        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"coverage.all.agree"],type='l',col="tomato")
          }
        }
        
        
        upper <- platypus.corr.avg.majority + platypus.corr.sd.majority
        lower <- platypus.corr.avg.majority - platypus.corr.sd.majority
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = lightblue.transparent, border = NA)
        
        upper <- platypus.corr.avg.all + platypus.corr.sd.all
        lower <- platypus.corr.avg.all - platypus.corr.sd.all
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = tomato.transparent, border = NA)
        
        
        lines(1:no.iterations,platypus.corr.avg.majority
              ,type='l'
              ,lwd=2
              ,col="blue"
        )
        lines(1:no.iterations,platypus.corr.avg.all,col="red3",lwd=2)
        
        legend("bottomright",legend=c("all agree","majority agrees"), lty=c(1,1),lwd=c(2,2), col=c("red3","blue"), bty='n')
        
        
        
        ## inclusion of unlabelled data plot
         par(mar=c(1, 4, 0.6, 0.1) + 0.1)
        plot(1:no.iterations,platypus.inclusion.avg
             ,ylim=c(0,max(platypus.inclusion.avg))
             ,type='l'
             #,pch=20
             ,lwd=2
             ,col="black"
             ,xlab="Iterations"
             ,ylab="# Labelled Samples"
             ,xaxt='n'
        )
        grid()
        
        upper <- platypus.inclusion.avg + platypus.inclusion.sd
        lower <- platypus.inclusion.avg - platypus.inclusion.sd
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = "grey", border = NA)
        lines(1:no.iterations,platypus.inclusion.avg,col="black",lwd=2)
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"no.ids.labelled"],type='l',col="grey")
          }
        }
        
        
        
        ### Threshold plot
        par(mar=c(4, 4, 0.6, 0.1) + 0.1)
        
        magenta.transparent <- rgb(col2rgb("magenta")[1,1],col2rgb("magenta")[2,1],col2rgb("magenta")[3,1],100,maxColorValue=255)
        
        plot(1:no.iterations,ll.upper.avg
             ,ylim=c(0,max(ll.upper.avg+ll.upper.sd))
             ,type='l'
             ,lwd=2
             ,col="magenta3"
             ,xlab="Iteration"
             ,ylab="Voting Sum"
        )
        grid()
        
        upper <- ll.upper.avg + ll.upper.sd
        lower <- ll.upper.avg - ll.upper.sd
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = magenta.transparent, border = NA)
        lines(1:no.iterations,ll.upper.avg,col="red3",lwd=2)
        
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"weighting.threshold.upper"],type='l',col="tomato")
          }
          lines(1:no.iterations, ll.upper.avg, col = 'red3',lwd=2)
        }
        
        
        
        
        #lower
        upper <- ll.lower.avg + ll.lower.sd
        lower <- ll.lower.avg - ll.lower.sd
        lightblue.transparent <- rgb(col2rgb("lightblue")[1,1],col2rgb("lightblue")[2,1],col2rgb("lightblue")[3,1],100,maxColorValue=255)
       
       
        
        polygon(c(1:no.iterations, rev(1:no.iterations)), c(upper, rev(lower)),
                col = green.transparent, border = NA)
        
        if(plot.single.cvfolds){
          for(cv in 1:cv.folds){
            lines(1:no.iterations,data.extended[which(data.extended[,"cv.fold"] == cv),"weighting.threshold.lower"],type='l',col="lightblue")
          }
        }
        
        lines(1:no.iterations,ll.lower.avg
              ,type='l'
              ,lwd=2
              ,col="green3")
        
        
        legend("left",legend=c("min. voting for favored label","max. voting for contrary label"), lty=c(1,1),lwd=c(2,2), col=c("magenta3","green3"), bty='n')
        
        
        dev.off()
        
        
      }
    }
  }
  
  
