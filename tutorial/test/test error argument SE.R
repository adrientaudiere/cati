

\section{Handling measurement Error}

Cati propose a simple argument to take into account for measurement error on "individual" values.
Indeed, most of times, one measure of a trait is assume to be representative of the real 

<<Handling_measurement_Error, eval = FALSE>>= 
  nperm <- 10
SE_vector <- c(0.1, 1, 10, 100, 10000)

res.simu5 <- list()
res.simu5.pval <- list()

for (i in 1: length(SE_vector)){
  res.simu5[[i]] <- list()
  res.simu5.pval[[i]] <- list()
  for(n in 1:nperm){
    res.simu5[[i]][[n]] <- Tstats(res.simu.traits2[[n]], ex.indplot2, ex.sp2,
                                  SE = SE_vector[i], printprogress = FALSE)
    
    res.simu5.pval[[i]][[n]] <- sum_Tstats(res.simu5[[i]][[n]], type = "p.value")
    
    print(paste("---", round(n*length(SE_vector)/(nperm*length(SE_vector)), 2) * 100, "%", sep = " "))
  }
}#end of calcul

@
  
  <<Handling_measurement_Error_results_SES>>= 
  #T_IP.IC
  meanSES.5.T_IP.IC.distriNorm <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_1$ses[, 1])))

meanSES.5.T_IP.IC.distriUni <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_1$ses[, 2])))

#T_IC.IR
meanSES.5.T_IC.IR.distriNorm <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_2$ses[, 1])))

meanSES.5.T_IC.IR.distriUni <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_2$ses[, 2])))

#T_PC.PR
meanSES.5.T_PC.PR.distriNorm <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_3$ses[, 1])))

meanSES.5.T_PC.PR.distriUni <- 
  lapply(res.simu5, function(y) 
    lapply(y, function(x) mean(ses.listofindex(as.listofindex(x))
                               $index_1_3$ses[, 2])))
@
  
  <<Handling_measurement_Error_plots_SES>>= 
  par(mfrow = c(2, 3))
boxplot(unlist(meanSES.5.T_IP.IC.distriNorm[[1]]), 
        unlist(meanSES.5.T_IP.IC.distriNorm[[2]]),
        unlist(meanSES.5.T_IP.IC.distriNorm[[3]]),
        unlist(meanSES.5.T_IP.IC.distriNorm[[4]]),
        unlist(meanSES.5.T_IP.IC.distriNorm[[5]]))

boxplot(unlist(meanSES.5.T_IP.IC.distriUni[[1]]), 
        unlist(meanSES.5.T_IP.IC.distriUni[[2]]),
        unlist(meanSES.5.T_IP.IC.distriUni[[3]]),
        unlist(meanSES.5.T_IP.IC.distriUni[[4]]))

boxplot(unlist(meanSES.5.T_IC.IR.distriNorm[[1]]), 
        unlist(meanSES.5.T_IC.IR.distriNorm[[2]]),
        unlist(meanSES.5.T_IC.IR.distriNorm[[3]]),
        unlist(meanSES.5.T_IC.IR.distriNorm[[4]]))

boxplot(unlist(meanSES.5.T_IC.IR.distriUni[[1]]), 
        unlist(meanSES.5.T_IC.IR.distriUni[[2]]),
        unlist(meanSES.5.T_IC.IR.distriUni[[3]]),
        unlist(meanSES.5.T_IC.IR.distriUni[[4]]))

boxplot(unlist(meanSES.5.T_PC.PR.distriNorm[[1]]), 
        unlist(meanSES.5.T_PC.PR.distriNorm[[2]]),
        unlist(meanSES.5.T_PC.PR.distriNorm[[3]]),
        unlist(meanSES.5.T_PC.PR.distriNorm[[4]]))

boxplot(unlist(meanSES.5.T_PC.PR.distriUni[[1]]), 
        unlist(meanSES.5.T_PC.PR.distriUni[[2]]),
        unlist(meanSES.5.T_PC.PR.distriUni[[3]]),
        unlist(meanSES.5.T_PC.PR.distriUni[[4]]))

par(mfrow = c(1, 1))
@
  
  <<>>= 
  xx.1.1 <- log10(sort(unlist(lapply(res.simu5.pval[[1]], function(x) x[1:10, 1]))))
xx.1.2 <- log10(sort(unlist(lapply(res.simu5.pval[[2]], function(x) x[1:10, 1]))))
xx.1.3 <- log10(sort(unlist(lapply(res.simu5.pval[[3]], function(x) x[1:10, 1]))))
xx.1.4 <- log10(sort(unlist(lapply(res.simu5.pval[[4]], function(x) x[1:10, 1]))))

boxplot(xx.1.1, xx.1.2, xx.1.3, xx.1.4)

xx.2.1 <- log10(sort(unlist(lapply(res.simu5.pval[[1]], function(x) x[1:10, 2]))))
xx.2.2 <- log10(sort(unlist(lapply(res.simu5.pval[[2]], function(x) x[1:10, 2]))))
xx.2.3 <- log10(sort(unlist(lapply(res.simu5.pval[[3]], function(x) x[1:10, 2]))))
xx.2.4 <- log10(sort(unlist(lapply(res.simu5.pval[[4]], function(x) x[1:10, 2]))))

boxplot(xx.2.1, xx.2.2, xx.2.3, xx.2.4)
@
  
  
  funct <- c("mean(x, na.rm = TRUE)", "sd(x, na.rm = TRUE)", "max(x, na.rm = TRUE) - min(x, na.rm = TRUE)" )

res0 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, nullmodels = "regional.ind", ind.plot = ind.plot.finch, nperm = 99)
res1 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, nullmodels = "regional.ind", ind.plot = ind.plot.finch, nperm = 99, SE = 1)
res100 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, nullmodels = "regional.ind", ind.plot = ind.plot.finch, nperm = 99, SE = 100)
res <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, nullmodels = "regional.ind", ind.plot = ind.plot.finch, nperm = 99, SE = apply(traits.finch, 2, function(x) sd(x, na.rm = T)))

plot(res)
plot(res0)
plot(res1)
plot(res100)
plot(as.listofindex(list(res0, res1, res, res100)))


