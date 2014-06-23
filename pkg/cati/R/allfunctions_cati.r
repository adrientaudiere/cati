#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__partvar

# traits is the Matrix of traits and factors are the nested factors to take into account in the partition of variance
partvar<-function(traits, factors, printprogress=TRUE){
	
	traits<-as.matrix(traits)
	factors<-as.matrix(factors)
	nfactors <- ncol(factors)
	ntraits  <- ncol(traits)
	res<-matrix(0, nrow=nfactors+1, ncol=ntraits)
	colnames(res)<-colnames(traits)
	
	if(!is.null(colnames(factors)))
		{rownames(res)<-c(colnames(factors), "within")
	} 
	
	else {
		rownames(res)<-c(paste("factor",1:(nfactors),sep=""),  "within") ; colnames(factors)<-c(paste("factor", 1:(nfactors),sep="")) 
	} 
	
	factors<-as.data.frame(factors)
	
	attach(factors)
	on.exit(expr = detach(factors))
	
	for (t in 1 : ntraits) {
		trait<-traits[,t]
		functionlme= paste('varcomp(lme(trait~1, random=~1|', paste(colnames(factors), collapse='/'), ",na.action=na.omit),1)", sep="")
		res[,t]<-as.vector(eval(parse(text=functionlme)))
		
		if(printprogress==TRUE) 
			{print(paste(round(t/ntraits*100, 2), "%", sep=" ")) }
		else{}
	}
	
	class(res)<-"partvar"
	
	print(res)
}

#Pie of variance partitioning
pie_partvar<-function(partvar, col.pie=NA){
	nfactors <- nrow(partvar)
	ntraits  <- ncol(partvar)
	for (t in 1 : ntraits) {
		pie(partvar[,t], main= colnames(partvar)[t], col=col.pie , labels=rownames(partvar))
	}
}

#barplot of variance partitioning
bar_partvar<-function(partvar,  col.bar=NA, leg=FALSE){
	
	oldpar<-par(no.readonly = TRUE)
	par(mar=c(5,6.5,4,2), cex=0.7)
	
	if(leg==FALSE)	{
		barplot(partvar, col=col.bar, las=1, horiz=T, xlab="% of variance")
	}
	
	else if(leg==TRUE) {
		barplot(partvar, col=col.bar, las=1, horiz=T, legend.=rownames(partvar), xlab="% of variance")
	}
	
	else{stop("leg must be TRUE or FALSE")}
	 par(oldpar)
}




#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__Tstats

### Function to calcul Tstats 
Tstats<-function(Ttraits, ind_plot, sp, reg.pool=NULL, nperm=NULL, printprogress=TRUE, p.value=TRUE){
	#6 variances: I: individual, P: population, C: community, R: region
	#IP; IC; IR; PC; PR; CR
	
	#Ttraits is the matrix of individual traits, ind_plot is the name of the plot in which the individual is (factor type), and sp is the species name of each individual
	
	names_sp_ind_plot<-as.factor(paste(sp, ind_plot, sep="@")) 
	Tplosp=unlist(strsplit(levels(names_sp_ind_plot),split="@"))[2*(1:nlevels(names_sp_ind_plot))]; names(Tplosp)=levels(names_sp_ind_plot);
	#Tplosp is the plot in wich the population is
	
	
	
  
  
	######################################## 
	####	Calcul of observed values	####
	######################################## 

	#________________________________________
	#Objects creation
	mean_IP<-matrix(nrow=nlevels(names_sp_ind_plot), ncol=ncol(Ttraits)); rownames(mean_IP)=levels(names_sp_ind_plot);
	mean_PC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_IP<-matrix(nrow=nlevels(names_sp_ind_plot), ncol=ncol(Ttraits))
	var_PC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_CR<-vector()
	var_IC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_PR<-vector()
	var_IR<-vector()
	T_IP.IC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	T_IC.IR<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	T_PC.PR<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
  
	for (t in 1: ncol(Ttraits)){
		mean_IP[,t]<-tapply(Ttraits[,t], names_sp_ind_plot  ,mean, na.rm=T)
		mean_PC[,t]<-tapply(mean_IP[,t], Tplosp , mean, na.rm=T)
		
		var_IP[,t]<-tapply(Ttraits[,t], names_sp_ind_plot, var, na.rm=T)
		var_PC[,t]<-tapply(mean_IP[,t], Tplosp  ,var, na.rm=T)
		var_CR[t]<-var(mean_PC[,t], na.rm=T)
		var_IC[,t]<-tapply(Ttraits[,t], ind_plot  ,var, na.rm=T)
		var_PR[t]<-var(as.vector(mean_IP[,t]), na.rm=T)
		var_IR[t]<-var(Ttraits[,t], na.rm=T)
		  
		for(s in 1 : nlevels(ind_plot)){
			T_IP.IC[s,t]<-mean(var_IP[grepl(levels(ind_plot)[s],Tplosp),t], na.rm=T)/var_IC[s,t]
			T_IC.IR[s,t]<-var_IC[s,t]/var_IR[t]
			T_PC.PR[s,t]<-var_PC[s,t]/var_PR[t]
		}
	}
	
	#________________________________________
	
	######################################### 
	#### 	  Creating null models  	 ####
	######################################### 
	
	if(is.numeric(nperm)){
		
		var_IP_nm1<-array(dim=c(nperm,ncol(Ttraits),nrow=length(Tplosp)))
		var_PC_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_IC_nm1<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_IC_nm2<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_PR_nm3<-array(dim=c(nperm,ncol(Ttraits)))
		var_IR_nm2<-array(dim=c(nperm,ncol(Ttraits)))
       
		mean_IP_nm3<-array(dim=c(nperm,ncol(Ttraits),length(Tplosp)))
		mean_PC_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
       
		Ttraits.nm1<-list()
		Ttraits.nm2<-list()
		Ttraits.nm3<-list()
              
		T_IP.IC_nm1<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		T_IC.IR_nm2<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		T_PC.PR_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		
		#Creation of the regional pool if not inform
		if(is.null(reg.pool)) {
			reg.pool<-Ttraits
		}
		  
		#Creation of three null models 
		if(printprogress==T){print("creating null models")}
		
		#________________________________________
		#modèle nul 1: permutations des valeurs de traits des individus dans la communauté   
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm1[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)) {
				Ttraits.nm1[[t]][[s]]<-list()
				for(i in 1:nperm){
					if (length(Ttraits[ind_plot==levels(ind_plot)[s], t]) != 1) {
						perm_ind_plot1<-sample(Ttraits[ind_plot==levels(ind_plot)[s], t], table(ind_plot)[s])
						Ttraits.nm1[[t]][[s]][[i]]<-perm_ind_plot1
					}
					else {Ttraits.nm1[[t]][[s]][[i]]<-"NA"}
				}
			} 
			if(printprogress==T){print(paste(round(t/ncol(Ttraits)/3*100,2),"%")) } else {}
		}
		
		#________________________________________
		#modèle nul 2: permutations des valeurs de traits des individus de la région    
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm2[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)) {
				Ttraits.nm2[[t]][[s]]<-list()
				for(i in 1:nperm){
					perm_ind_plot2<-sample(reg.pool[, t], table(ind_plot)[s])
					Ttraits.nm2[[t]][[s]][[i]]<-perm_ind_plot2
				}
			}
			if(printprogress==T){print(paste(round(33.3+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		
		#________________________________________  
		#modèle nul 3: permutations des espèces au niveau de la région   
		Ttraits_by_sp<-apply(Ttraits,2,function(x) tapply(x,names_sp_ind_plot,mean))  
		Ttraits_by_pop<-Ttraits_by_sp[match(names_sp_ind_plot,rownames(Ttraits_by_sp)),]
		#Ttraits_by_sp<-aggregate(Ttraits, by = list(names_sp_ind_plot), mean, na.rm = T)[,-1] 
				
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm3[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)){
				Ttraits.nm3[[t]][[s]]<-list()
				for(i in 1:nperm){
					perm_ind_plot3<-sample(Ttraits_by_pop, table(ind_plot)[s])
					Ttraits.nm3[[t]][[s]][[i]]<-perm_ind_plot3
				}
			} 
			if(printprogress==T){print(paste(round(66.6+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		
		#________________________________________
	
		######################################### 
		#### calcul of Tstats on null models ####
		######################################### 

		if(printprogress==T){print("calcul of Tstats using null models")}
		
		yy<-length(names_sp_ind_plot)
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){ 
				mean_IP_nm3[i,t,]<-tapply(unlist(Ttraits.nm3[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind_plot  ,function(x) mean(x, na.rm=T))
				mean_PC_nm3[i,t,]<-tapply(mean_IP_nm3[i,t,], Tplosp, mean, na.rm=T)
			}
			if(printprogress==T){print(paste(round(t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		} 
		   
		   
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){
				var_IP_nm1[i,t,]<-tapply(unlist(Ttraits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind_plot  ,function(x) var(x, na.rm=T))
				var_PC_nm3[i,t,]<-tapply(mean_IP_nm3[i,t,], Tplosp  ,var, na.rm=T)
				var_IC_nm1[i,t,]<-tapply(unlist(Ttraits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], ind_plot  ,function(x) var(x, na.rm=T))
				var_IC_nm2[i,t,]<-tapply(unlist(Ttraits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], ind_plot  ,function(x) var(x, na.rm=T))
				var_PR_nm3[i,t]<-var(as.vector(mean_IP_nm3[i,t,]), na.rm=T)
				var_IR_nm2[i,t]<-var(unlist(Ttraits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], na.rm=T)
			}
			if(printprogress==T){print(paste(round(33.3+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		 
		   
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){
				for(s in 1 : nlevels(ind_plot)){
					T_IP.IC_nm1[i,t,s]<-mean(var_IP_nm1[i,t,grepl(levels(ind_plot)[s],Tplosp)], na.rm=T)/var_IC_nm1[i,t,s] 
					T_IC.IR_nm2[i,t,s]<-var_IC_nm2[i,t,s]/var_IR_nm2[i,t]
					T_PC.PR_nm3[i,t,s]<-var_PC_nm3[i,t,s]/var_PR_nm3[i,t]
				}
			} 
			if(printprogress==T){print(paste(round(66.6+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}       
		      
	}#end of calcul of Tstats using null models
         
	colnames(T_IP.IC)<-colnames(Ttraits)
    colnames(T_IC.IR)<-colnames(Ttraits)
    colnames(T_PC.PR)<-colnames(Ttraits)
  
	if(is.numeric(nperm)){
		colnames(T_IP.IC_nm1)<-colnames(Ttraits)
		colnames(T_IC.IR_nm2)<-colnames(Ttraits)
		colnames(T_PC.PR_nm3)<-colnames(Ttraits)
	}
	
	rownames(T_IP.IC)<-levels(as.factor(Tplosp))
  	rownames(T_IC.IR)<-levels(as.factor(Tplosp))
 	rownames(T_PC.PR)<-levels(as.factor(Tplosp))
 
	
	#________________________________________
    res<-list()
    res$T_IP.IC<-T_IP.IC
    res$T_IC.IR<-T_IC.IR
    res$T_PC.PR<-T_PC.PR
    
    res$variances<-list()
    
    res$variances$var_IP<-var_IP
    res$variances$var_PC<-var_PC
    res$variances$var_CR<-var_CR
	res$variances$var_IC<-var_IC
    res$variances$var_PR<-var_PR
    res$variances$var_IR<-var_IR
    
    res$variances$var_IP_nm1<-var_IP_nm1
    res$variances$var_PC_nm3<-var_PC_nm3
    res$variances$var_IC_nm1<-var_IC_nm1
	res$variances$var_IC_nm2<-var_IC_nm2
    res$variances$var_PR_nm3<-var_PR_nm3
    res$variances$var_IR_nm2<-var_IR_nm2
    	
	if(is.numeric(nperm)){	 
		res$T_IP.IC_nm<-T_IP.IC_nm1
       	res$T_IC.IR_nm<-T_IC.IR_nm2
        res$T_PC.PR_nm<-T_PC.PR_nm3
    }   
    else{}
 	
 	#________________________________________
 	
 	######################################### 
	####		 calcul of p.value		 ####
	######################################### 
 
 	if(printprogress==T){print("calcul of p.value")}
 	
 	if(p.value==T){
		p.valueT_IP.IC.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_IC.IR.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_PC.PR.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		
		p.valueT_IP.IC.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_IC.IR.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_PC.PR.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		
		for (t in 1: ncol(Ttraits)){
			for(s in 1:  nlevels(ind_plot)){
 				p.valueT_IP.IC.sup[s,t]<-(sum(res$T_IP.IC[s,t]<res$T_IP.IC_nm[,t,s], na.rm=T)+1)/(1+length(res$T_IP.IC_nm[,t,s]))
 				p.valueT_IC.IR.sup[s,t]<-(sum(res$T_IC.IR[s,t]<res$T_IC.IR_nm[,t,s], na.rm=T)+1)/(1+length(res$T_IC.IR_nm[,t,s]))
 				p.valueT_PC.PR.sup[s,t]<-(sum(res$T_PC.PR[s,t]<res$T_PC.PR_nm[,t,s], na.rm=T)+1)/(1+length(res$T_PC.PR_nm[,t,s]))
		
				p.valueT_IP.IC.inf[s,t]<-(sum(res$T_IP.IC[s,t]>res$T_IP.IC_nm[,t,s], na.rm=T)+1)/(1+length(res$T_IP.IC_nm[,t,s]))
				p.valueT_IC.IR.inf[s,t]<-(sum(res$T_IC.IR[s,t]>res$T_IC.IR_nm[,t,s], na.rm=T)+1)/(1+length(res$T_IC.IR_nm[,t,s]))
				p.valueT_PC.PR.inf[s,t]<-(sum(res$T_PC.PR[s,t]>res$T_PC.PR_nm[,t,s], na.rm=T)+1)/(1+length(res$T_PC.PR_nm[,t,s]))
			}
		}	
	    
		colnames(p.valueT_IP.IC.sup)<-colnames(Ttraits)
		colnames(p.valueT_IC.IR.sup)<-colnames(Ttraits)
		colnames(p.valueT_PC.PR.sup)<-colnames(Ttraits)
		
		rownames(p.valueT_IP.IC.sup)<-levels(Tplosp)
		rownames(p.valueT_IC.IR.sup)<-levels(Tplosp)
		rownames(p.valueT_PC.PR.sup)<-levels(Tplosp)
  	
		colnames(p.valueT_IP.IC.inf)<-colnames(Ttraits)
		colnames(p.valueT_IC.IR.inf)<-colnames(Ttraits)
		colnames(p.valueT_PC.PR.inf)<-colnames(Ttraits)
		
		rownames(p.valueT_IP.IC.inf)<-levels(Tplosp)
		rownames(p.valueT_IC.IR.inf)<-levels(Tplosp)
		rownames(p.valueT_PC.PR.inf)<-levels(Tplosp)
		
		res$pval<-list()
		
		res$pval$T_IP.IC.inf<-p.valueT_IP.IC.inf
		res$pval$T_IC.IR.inf<-p.valueT_IC.IR.inf
		res$pval$T_PC.PR.inf<-p.valueT_PC.PR.inf
		
		res$pval$T_IP.IC.sup<-p.valueT_IP.IC.sup
		res$pval$T_IC.IR.sup<-p.valueT_IC.IR.sup
		res$pval$T_PC.PR.sup<-p.valueT_PC.PR.sup
	}
	else{}
	
    class(res)<-"Tstats"
    
    return(res)
}

### Function to represent standardised effect size of Tstats using null models
plot.Tstats<-function(x, val.quant=c(0.025,0.975), col.Tstats=c("red","purple","green"), type="normal", add.conf=TRUE, ylim=NULL, xlim=NULL, ...){
	#possible type = "color_cond", "simple", "simple_sd", "normal" and "barplot"	
	
	tstats<-x
	
	#________________________________________
	#Calcul of standardised effect size
	ses.T_IP.IC<-(tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR<-(tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR<-(tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.inf<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.inf<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.inf<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.sup<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.sup<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.sup<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	
	#________________________________________
	#Condition to be significantly different from null models with respect to values of quantile choose
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
		
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	all=c(ses.T_IP.IC,ses.T_IC.IR,ses.T_PC.PR)
	if(is.null(ylim)) {ylim=c(5*dim(tstats$T_IP.IC)[2]+3,3)}
	if(is.null(xlim)) {xlim=c(min(all, na.rm=T),max(all, na.rm=T))}
	
	par(mar=c(5, 7, 4, 2))
	plot(0,0, ylab="Traits",yaxt= "n", xlab="Tstats Standardized Effect Size",  col="black", type="l", xlim=xlim, ylim=ylim, ...)
	axis(side=2, seq(from=5.5, to=4*dim(tstats$T_IP.IC)[2]+1.5, by=4), labels=colnames(tstats$T_IP.IC), las=1, cex.axis=0.7 ) 
	legend("bottom", inset=.005, title="Tstats", c("T_IP.IC","T_IC.IR","T_PC.PR"), fill=col.Tstats, horiz=TRUE, cex=0.7, bty="n")
	
	#________________________________________
	#plot : possible type = "color_cond", "simple", "simple_range", "normal" and "barplot"	
	
	#__________
	if(type=="color_cond"){
		
		if(length(col.Tstats)==3) {col.Tstats[4:6]<-"grey"} 
		if(length(col.Tstats)!=6) {print("Warnings: plot type color_cond need 3 or 6 colors in the argument col.Tstats")}
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			points(ses.T_IP.IC[,t], rep(t*4, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[4])
			points(ses.T_IC.IR[,t], rep(t*4+1, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[5])
			points(ses.T_PC.PR[,t], rep(t*4+2, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[6])
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			points(ses.T_IP.IC[,t][cond.T_IP.IC.inf[,t]], rep(t*4,times=length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])),pch=16, col=col.Tstats[1])
			points(ses.T_IC.IR[,t][cond.T_IC.IR.inf[,t]], rep(t*4+1,times=length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])), pch=16, col=col.Tstats[2])
			points(ses.T_PC.PR[,t][cond.T_PC.PR.inf[,t]], rep(t*4+2, times=length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])), pch=16, col=col.Tstats[3])
			
			points(ses.T_IP.IC[,t][cond.T_IP.IC.sup[,t]], rep(t*4,times=length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])), pch=16, col=col.Tstats[1])
			points(ses.T_IC.IR[,t][cond.T_IC.IR.sup[,t]], rep(t*4+1,times=length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])), pch=16, col=col.Tstats[2])
			points(ses.T_PC.PR[,t][cond.T_PC.PR.sup[,t]], rep(t*4+2, times=length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])), pch=16, col=col.Tstats[3])
			
			if (length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])>0) 	{text(ses.T_IP.IC[,t][cond.T_IP.IC.inf[,t]], rep(t*4-2,times=length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IP.IC.inf[,t]], cex=0.6, srt=45, col=col.Tstats[1], pos=3)}
			if (length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])>0) 	{text(ses.T_IC.IR[,t][cond.T_IC.IR.inf[,t]], rep(t*4,times=length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IC.IR.inf[,t]], cex=0.6, srt=45, col=col.Tstats[2], pos=2)}
			if (length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])>0) 	{text(ses.T_PC.PR[,t][cond.T_PC.PR.inf[,t]], rep(t*4+4, times=length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_PC.PR.inf[,t]], cex=0.6, srt=45, col=col.Tstats[3], pos=1)}
			
			if (length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])>0) 	{text(ses.T_IP.IC[,t][cond.T_IP.IC.sup[,t]], rep(t*4-2,times=length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IP.IC.sup[,t]], cex=0.6, srt=45, col=col.Tstats[1], pos=3)}
			if (length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])>0) 	{text(ses.T_IC.IR[,t][cond.T_IC.IR.sup[,t]], rep(t*4,times=length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IC.IR.sup[,t]], cex=0.6, srt=45, col=col.Tstats[2], pos=2)}
			if (length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])>0) 	{text(ses.T_PC.PR[,t][cond.T_PC.PR.sup[,t]], rep(t*4+4, times=length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_PC.PR.sup[,t]], cex=0.6, srt=45, col=col.Tstats[3], pos=1)}
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
	}

	#__________
	else if(type=="simple"){
	
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			
			points(ses.T_IP.IC[,t], rep(t*4, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[1])
			points(ses.T_IC.IR[,t], rep(t*4+1, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[2])
			points(ses.T_PC.PR[,t], rep(t*4+2, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.sup[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.sup[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.sup[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.inf[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.inf[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.inf[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])	
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
	}
	
	#__________
	else if(type=="simple_sd"){
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
					
			points(mean(ses.T_IP.IC.sup[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.sup[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.sup[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.inf[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.inf[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.inf[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])	
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}	
	}
	
	#__________
	else if(type=="normal"){
		
		for(t in 1:dim(tstats$T_IP.IC)[2]){
		
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
	   
			points(min(ses.T_IP.IC[,t], na.rm=T), t*4, pch="*", col=col.Tstats[1])
			points(min(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch="*", col=col.Tstats[2])
			points(min(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch="*", col=col.Tstats[3])
			
			points(max(ses.T_IP.IC[,t], na.rm=T), t*4, pch="*", col=col.Tstats[1])
			points(max(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch="*", col=col.Tstats[2])
			points(max(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch="*", col=col.Tstats[3])
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
		
		if(add.conf==T){
			points(colMeans(ses.T_IP.IC.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])
			
			points(colMeans(ses.T_IP.IC.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])  
		}	
		else {}
	}

	#__________
	else if(type=="barplot"){
	
		for(t in 1:dim(tstats$T_IP.IC)[2]){
		
			segments(mean(ses.T_IP.IC[,t], na.rm=T), t*4  , 0, t*4, pch=17, col=col.Tstats[1], lwd=8)
			segments(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1, 0, t*4+1, pch=17, col=col.Tstats[2], lwd=8)
			segments(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, 0, t*4+2, pch=17, col=col.Tstats[3], lwd=8)
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])

			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
		
		if(add.conf==T){
			points(colMeans(ses.T_IP.IC.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])
			
			points(colMeans(ses.T_IP.IC.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])  
		}
		else {}
	}
	else{print(paste("Error:",type,"is not a valid type of plot"))}
	
	par(mar=c(5, 4, 4, 2) + 0.1) #return to default parameter
}

### Function to summarize traits and community which show a significant difference between observed and simulated value
summary_Tstats<-function(x, val.quant=c(0.025,0.975), type="all") {
	
	tstats<-x
	#________________________________________
	ses.T_IP.IC<-(tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR<-(tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR<-(tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.inf<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.inf<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.inf<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.sup<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.sup<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.sup<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	
	ses.T_IP.IC.mean<-t(colMeans((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_IC.IR.mean<-t(colMeans((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_PC.PR.mean<-t(colMeans((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	  
	ses.T_IP.IC.inf.mean<-apply(ses.T_IP.IC.inf,2, mean)
	ses.T_IC.IR.inf.mean<-apply(ses.T_IC.IR.inf,2, mean)
	ses.T_PC.PR.inf.mean<-apply(ses.T_PC.PR.inf,2, mean)	  
	
	ses.T_IP.IC.sup.mean<-apply(ses.T_IP.IC.sup,2, mean)
	ses.T_IC.IR.sup.mean<-apply(ses.T_IC.IR.sup,2, mean)
	ses.T_PC.PR.sup.mean<-apply(ses.T_PC.PR.sup,2, mean)
	
	#________________________________________
	#Condition to be significantly different from null models with respect to values of quantile choosen
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
		
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	cond.T_IP.IC.inf.mean<-ses.T_IP.IC.mean<ses.T_IP.IC.inf.mean
	cond.T_IC.IR.inf.mean<-ses.T_IC.IR.mean<ses.T_IC.IR.inf.mean
	cond.T_PC.PR.inf.mean<-ses.T_PC.PR.mean<ses.T_PC.PR.inf.mean
		
	cond.T_IP.IC.sup.mean<-ses.T_IP.IC.mean>ses.T_IP.IC.sup.mean
	cond.T_IC.IR.sup.mean<-ses.T_IC.IR.mean>ses.T_IC.IR.sup.mean
	cond.T_PC.PR.sup.mean<-ses.T_PC.PR.mean>ses.T_PC.PR.sup.mean

	
	#________________________________________
	if(type=="binary"){
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)
	}
	
	#________________________________________
	else if(type=="percent"){
		
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats[1,]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep="")
		summ.Tstats[2,]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep="")
		summ.Tstats[3,]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep="")
		summ.Tstats[4,]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep="")
		summ.Tstats[5,]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep="")
		summ.Tstats[6,]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep="")
				
		summ.Tstats[1,][cond.T_IP.IC.inf.mean]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep="")
		summ.Tstats[2,][cond.T_IP.IC.sup.mean]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep="")
		summ.Tstats[3,][cond.T_IC.IR.inf.mean]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep="")
		summ.Tstats[4,][cond.T_IC.IR.sup.mean]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep="")
		summ.Tstats[5,][cond.T_PC.PR.inf.mean]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep="")
		summ.Tstats[6,][cond.T_PC.PR.sup.mean]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep="")	
	
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)

	}
	
	#________________________________________
	else if(type=="site"){
	
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){
			
			if(sum(cond.T_IP.IC.inf[,t], na.rm=T)>0) 
				{summ.Tstats[1,t]<-paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ") }
			else{summ.Tstats[1,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IP.IC.sup[,t], na.rm=T)>0) 
				{summ.Tstats[2,t]<-paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse=" ")  }
			else{summ.Tstats[2,t]<-"H0 not rejected"}
								
			if(sum(cond.T_IC.IR.inf[,t], na.rm=T)>0)
				{summ.Tstats[3,t]<-paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ")  }
			else{summ.Tstats[3,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IC.IR.sup[,t], na.rm=T)>0)
				{summ.Tstats[4,t]<-paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse=" ")	}
			else{summ.Tstats[4,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.inf[,t], na.rm=T)>0)
				{summ.Tstats[5,t]<-paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse=" ") 	}
			else{summ.Tstats[5,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.sup[,t], na.rm=T)>0)
				{summ.Tstats[6,t]<-paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse=" ") 	}
			else{summ.Tstats[6,t]<-"H0 not rejected"}		
		}
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)

	}
	
	
	#________________________________________
	else if(type=="p.value"){
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(tstats$pval$T_IP.IC.inf, tstats$pval$T_IP.IC.sup , tstats$pval$T_IC.IR.inf, tstats$pval$T_IC.IR.sup , tstats$pval$T_PC.PR.inf, tstats$pval$T_PC.PR.sup)
		rownames(summ.Tstats)<-c(paste(rep("T_IP.IC.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IP.IC.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IC.IR.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IC.IR.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_PC.PR.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_PC.PR.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)))
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)
	}
	

	#________________________________________
	else if(type=="all"){
		summ.Tstats<-list()
		
		#__________
		##p.value
		summ.Tstats$p.value <-matrix("H0 not rejected", nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$p.value <- rbind(tstats$pval$T_IP.IC.inf, tstats$pval$T_IP.IC.sup , tstats$pval$T_IC.IR.inf, tstats$pval$T_IC.IR.sup , tstats$pval$T_PC.PR.inf, tstats$pval$T_PC.PR.sup)
	
		#__________
		##percent
		summ.Tstats$percent <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$percent[1,]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[2,]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep="")
		summ.Tstats$percent[3,]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[4,]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep="")
		summ.Tstats$percent[5,]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[6,]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep="")
				
		summ.Tstats$percent[1,][cond.T_IP.IC.inf.mean]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep="")
		summ.Tstats$percent[2,][cond.T_IP.IC.sup.mean]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep="")
		summ.Tstats$percent[3,][cond.T_IC.IR.inf.mean]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep="")
		summ.Tstats$percent[4,][cond.T_IC.IR.sup.mean]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep="")
		summ.Tstats$percent[5,][cond.T_PC.PR.inf.mean]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep="")
		summ.Tstats$percent[6,][cond.T_PC.PR.sup.mean]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep="")	
		
		#__________
		##sites
		summ.Tstats$sites <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){
			
			if(sum(cond.T_IP.IC.inf[,t], na.rm=T)>0) 
				{summ.Tstats$sites[1,t]<-paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ") }
			else{summ.Tstats$sites[1,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IP.IC.sup[,t], na.rm=T)>0) 
				{summ.Tstats$sites[2,t]<-paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse=" ")  }
			else{summ.Tstats$sites[2,t]<-"H0 not rejected"}
								
			if(sum(cond.T_IC.IR.inf[,t], na.rm=T)>0)
				{summ.Tstats$sites[3,t]<-paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ")  }
			else{summ.Tstats$sites[3,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IC.IR.sup[,t], na.rm=T)>0)
				{summ.Tstats$sites[4,t]<-paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse=" ")	}
			else{summ.Tstats$sites[4,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.inf[,t], na.rm=T)>0)
				{summ.Tstats$sites[5,t]<-paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse=" ") 	}
			else{summ.Tstats$sites[5,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.sup[,t], na.rm=T)>0)
				{summ.Tstats$sites[6,t]<-paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse=" ") 	}
			else{summ.Tstats$sites[6,t]<-"H0 not rejected"}		
		}
		
		#__________
		##binary
		summ.Tstats$binary <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$binary <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)
		
		#__________
		rownames(summ.Tstats$p.value)<-c(rep("T_IP.IC.inf",dim(tstats$T_IP.IC)[1]), rep("T_IP.IC.sup",dim(tstats$T_IP.IC)[1]), rep("T_IC.IR.inf",dim(tstats$T_IP.IC)[1]), rep("T_IC.IR.sup",dim(tstats$T_IP.IC)[1]), rep("T_PC.PR.inf",dim(tstats$T_IP.IC)[1]), rep("T_PC.PR.sup",dim(tstats$T_IP.IC)[1]))
		rownames(summ.Tstats$binary)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$percent)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$sites)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats$p.value) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$binary) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$sites) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$percent) <- colnames(tstats$T_IP.IC)
	}
	
	else{stop("Error: type must be 'binary', 'percent', 'p.value', 'site' or 'all'.")}
	
	return(summ.Tstats)
}

### Function to represent summarize Tstats
barplot.Tstats<-function(height, val.quant=c(0.025,0.975), col.Tstats=c("red","purple","green","white"), ylim=NULL, ...){
   
  tstats<-height
  
  T_IP.IC.inf<-apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  T_IC.IR.inf<-apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  T_PC.PR.inf<-apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  
  T_IP.IC.sup<-apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
  T_IC.IR.sup<-apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
  T_PC.PR.sup<-apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
  
  if(is.null(ylim)){
	ylim=c(min(c(T_IP.IC.inf,T_IC.IR.inf,T_PC.PR.inf,colMeans(na.omit(tstats$T_IP.IC))-apply(na.omit(tstats$T_IP.IC), 2,sd),colMeans(na.omit(tstats$T_IC.IR))-apply(na.omit(tstats$T_IC.IR), 2,sd),colMeans(na.omit(tstats$T_PC.PR))-apply(na.omit(tstats$T_PC.PR), 2,sd)), na.rm=T) , max(c(colMeans(na.omit(tstats$T_IP.IC))+apply(na.omit(tstats$T_IP.IC), 2,sd),colMeans(na.omit(tstats$T_IC.IR))+apply(na.omit(tstats$T_IC.IR), 2,sd),colMeans(na.omit(tstats$T_PC.PR))+apply(na.omit(tstats$T_PC.PR), 2,sd)), na.rm=T)) 
  } 
  
  df.bar<-barplot(rbind(colMeans(na.omit(tstats$T_IP.IC)), colMeans(na.omit(tstats$T_IC.IR)),colMeans(na.omit(tstats$T_PC.PR)),0), beside=T, plot=F)
  barplot(rbind(colMeans(na.omit(tstats$T_IP.IC)), colMeans(na.omit(tstats$T_IC.IR)),colMeans(na.omit(tstats$T_PC.PR)),0), col=col.Tstats, beside=T, ylim=ylim, ...)
  segments( df.bar[1,], colMeans(na.omit(tstats$T_IP.IC))+apply(na.omit(tstats$T_IP.IC), 2,sd),df.bar[1,],colMeans(na.omit(tstats$T_IP.IC))-apply(na.omit(tstats$T_IP.IC), 2,sd))
  segments( df.bar[2,], colMeans(na.omit(tstats$T_IC.IR))+apply(na.omit(tstats$T_IC.IR), 2,sd),df.bar[2,],colMeans(na.omit(tstats$T_IC.IR))-apply(na.omit(tstats$T_IC.IR), 2,sd))
  segments( df.bar[3,], colMeans(na.omit(tstats$T_PC.PR))+apply(na.omit(tstats$T_PC.PR), 2,sd),df.bar[3,],colMeans(na.omit(tstats$T_PC.PR))-apply(na.omit(tstats$T_PC.PR), 2,sd))
  
  points(type="l", df.bar[1,], colMeans(T_IP.IC.sup, na.rm=T), col=col.Tstats[1])
  points(type="l", df.bar[2,], colMeans(T_IC.IR.sup, na.rm=T), col=col.Tstats[2])
  points(type="l", df.bar[3,], colMeans(T_PC.PR.sup, na.rm=T), col=col.Tstats[3])
  
  points(type="l", df.bar[1,], colMeans(T_IP.IC.inf, na.rm=T), col=col.Tstats[1])
  points(type="l", df.bar[2,], colMeans(T_IC.IR.inf, na.rm=T), col=col.Tstats[2])
  points(type="l", df.bar[3,], colMeans(T_PC.PR.inf, na.rm=T), col=col.Tstats[3])
  
}




#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ com.index

#Calcul of statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models
#For each statistic this function return observed value and correspondant Null distribution
#This function implement three null models which keep unchanged the number of individual per community
#Models 1 correspond to randomization of individual values within community
#Models 2 correspond to randomization of individual values within region
#Models 3 correspond to randomization of population values within region

#In most case, model 1 and 2 correspond to index at the individual level and the model 3 to index at the species (or any other aggregate variable like genus or family) level

com.index<-function(traits=NULL, index=NULL, namesindex=NULL, nullmodels=NULL, ind.plot=NULL, sp=NULL, com=NULL, reg.pool=NULL, nperm=99, printprogress=TRUE, ind.value=TRUE, type="count"){
	
	#If data are from species or population traits, this function transform this data in a suitable format for cati
	if(!ind.value){
		if(is.null(com)) {stop("if ind.value=FALSE, you need to replace arguments ind_plot by a community matrix 'com' ")}
		
		rownames(traits)<-sp
		res.interm<-ab_to_ind(traits, com, type=type)
	
		traits<-res.interm$traits
		sp<-res.interm$sp
		ind.plot<-res.interm$ind.plot
	}
	
	if(!is.null(ind.plot) & !is.null(com)){
		warnings("If ind.plot and com are provide and ind.value=F, the function use only the argument com")
	}
	
	nindex<-length(index)
	
	if(length(nullmodels)==1){
		nullmodels<-rep(nullmodels,times=nindex)
	}	
	
	if(is.null(namesindex)) {  namesindex<-index }
	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)
	
	traits<-traits[order(ind.plot),]
	ind.plot<-ind.plot[order(ind.plot)]
	sp<-sp[order(ind.plot)]
	
	name_sp_sites=paste(sp, ind.plot, sep="_")
	comm=NULL
	comm<-t(table(ind.plot,1:length(ind.plot)))
	
	S = colSums(comm>0)
	ncom=length(S)
	
	
	#Creation of the regional pool if not inform
	if(is.null(reg.pool)) {
		reg.pool<-traits
	}
	
	if(!is.null(reg.pool) & sum(nullmodels==2)==0) {
		warnings("Custom regional pool'reg.pool' is only used in the case of null model 2")
	}
	
	if(is.numeric(nperm)){
		######################################### 
		#### 	  Calcul of null models  	 ####
		######################################### 
		#Creation of three null models 
		if(printprogress==T){ print("creating null models")}
		
		if(sum(nullmodels==1)>0){
			#________________________________________
			#modèle nul 1: permutations des valeurs de traits des individus dans la communauté   
			
			traits.nm1<-list()
			
			for (t in 1: ntr){	
				traits.nm1[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits[ind.plot==levels(ind.plot)[s], t], table(ind.plot)[s])
					}
					
					traits.nm1[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.1",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==2)>0){	
			#________________________________________
			#modèle nul 2: permutations des valeurs de traits des individus de la région    
			traits.nm2<-list()
			
			for (t in 1: ntr){	
				traits.nm2[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(reg.pool[, t], table(ind.plot)[s])
					}
					
					traits.nm2[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.2",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==3)>0){
			#________________________________________  
			#modèle nul 3: permutations des espèces au niveau de la région   
			traits.nm3<-list()
			traits_by_sp<-apply(traits,2,function(x) tapply(x,name_sp_sites,mean, na.rm=T))  
			traits_by_pop<-traits_by_sp[match(name_sp_sites,rownames(traits_by_sp)),]
			
			for (t in 1: ntr){	
				traits.nm3[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits_by_pop, table(ind.plot)[s])
					}
					
					traits.nm3[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 		
				if(printprogress==T){
					print(paste("nm.3",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		######################################## 
		####	 Calcul of random values   	####
		######################################## 
		Null<-list()
		nm_bypop<-list()
		nm_bypop.bis<-list()
		
		if(printprogress==T){print("calcul of null values using null models")}
		
		for(i in 1:nindex){
			if(nullmodels[i]==1){nm.bis<-traits.nm1[[1]]}
			else if(nullmodels[i]==2){nm.bis<-traits.nm2[[1]]}
			else if(nullmodels[i]==3){nm.bis<-traits.nm3[[1]]}
			else{print("nullmodels need 1, 2 or 3")}
			
			functionindex= eval(index[i])
			
			if(nullmodels[i]==3){
				nm_bypop.bis[[eval(namesindex[i])]]<- apply(nm.bis, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
					
				dim2<-dim(apply(nm_bypop.bis[[eval(namesindex[i])]], 2, function (x) eval(parse(text=functionindex))))[1]
				Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, dim2, nperm) )
				if(is.null(dim2)) {
					Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, 1, nperm) )
				}	
			}
				
			else{
				dim2<-dim(apply(nm.bis, 2, function (x) eval(parse(text=functionindex))))[1]
				Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, dim2, nperm) )
				if(is.null(dim2)) {
					Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, 1, nperm) )
				}		
			}
			
			for (t in 1: ntr){
			
				if(nullmodels[i]==1){nm<-traits.nm1[[t]]}
				else if(nullmodels[i]==2){nm<-traits.nm2[[t]]}
				else if(nullmodels[i]==3){nm<-traits.nm3[[t]]}
				else{print("nullmodels need to be either 1, 2 or 3")}
				
				if(nullmodels[i]==3){
					nm_bypop[[eval(namesindex[i])]]<- apply(nm, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
					Null[[eval(namesindex[i])]] [t,,] <- apply(nm_bypop[[eval(namesindex[i])]], 2, function (x) eval(parse(text=functionindex)))			
				}
				
				else{				
					Null[[eval(namesindex[i])]] [t,,] <- apply(nm, 2, function (x) eval(parse(text=functionindex)))				
				}
				
				if(printprogress==T){
					print(paste(eval(namesindex[i]), round(t/ntr*100,2),"%")) 
				} 
			}
		}
	}
		  
	######################################## 
	####	Calcul of observed values	####
	######################################## 
	obs<-list()
	
	if(printprogress==T){print("calcul of observed values")}
	
	for(i in 1:nindex){
		functionindex= eval(index[i])
		
		if(nullmodels[i]==3) {
			traits.pop<-apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
			obs[[eval(namesindex[i])]] <- array(dim=c(ntr, dim(apply(traits.pop, 2, function (x) eval(parse(text=functionindex))))[1]))
			obs[[eval(namesindex[i])]] <-  apply(traits.pop, 2, function (x) eval(parse(text=functionindex)))
		}
		
		else if(nullmodels[i]==1  |  nullmodels[i]==2) {
			obs[[eval(namesindex[i])]] <- array(dim=c(ntr, dim(apply(traits, 2, function (x) eval(parse(text=functionindex))))[1]))
			obs[[eval(namesindex[i])]] <- apply(traits, 2, function (x) eval(parse(text=functionindex)))
			#obs[[eval(namesindex[i])]] [ !is.finite(obs[[eval(namesindex[i])]] )]<-NA
		}
		if(printprogress==T){
			print(paste(round(i/nindex*100,2),"%")) 
		} 
	}
			
		
	######################################## 
	####		Create results list		####
	######################################## 
	
	com.index<-list()  
	com.index$obs<-obs
		
	if(is.numeric(nperm)){
		com.index$Null<-Null
	}
	
	com.index$list.index<-list()
	com.index$list.index.t<-list()
	name.com.index_list.index<-vector()
	
	for(i in 1:nindex){
		com.index$list.index.t[[seq(1,nindex*2,by=2)[i]]]<-t(obs[[i]])
		com.index$list.index[[seq(1,nindex*2,by=2)[i]]]<-obs[[i]]
		name.com.index_list.index[seq(1,nindex*2,by=2)[i]]<-names(obs)[i]
		
		if(is.numeric(nperm)){
			com.index$list.index[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			com.index$list.index.t[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			name.com.index_list.index[seq(1,nindex*2,by=2)[i]+1]<-paste(names(Null)[i], "nm", sep="_")
		}
	}
	
	names(com.index$list.index.t)<-name.com.index_list.index
	names(com.index$list.index)<-name.com.index_list.index
		
	com.index$sites_richness<-S
	com.index$namestraits<-namestraits
	
	class(com.index)<-"com.index"
	
	return(com.index)
}

com.index.multi<-function(traits=NULL, index=NULL, by.factor=NULL, namesindex=NULL, nullmodels=NULL, ind.plot=NULL, sp=NULL, com=NULL, reg.pool=NULL, nperm=99, printprogress=TRUE, ind.value=TRUE, type="count"){
	
	names_sp_ind_plot<-as.factor(paste(sp, ind.plot, sep="@")) 
	
	#If data are from species or population traits, this function transform this data in a suitable format for cati
	if(!ind.value){
		if(is.null(com)) {stop("if ind.value=FALSE, you need to replace arguments ind_plot by a community matrix 'com' ")}
		
		rownames(traits)<-sp
		res.interm<-ab_to_ind(traits, com, type=type)
	
		traits<-res.interm$traits
		sp<-res.interm$sp
		ind.plot<-res.interm$ind.plot
	}
	
	if(!is.null(ind.plot) & !is.null(com)){
		warnings("If ind.plot and com are provide and ind.value=F, the function use only the argument com")	
	}
	
	####
	nindex<-length(index)
	
	if(length(nullmodels)==1){
		nullmodels<-rep(nullmodels,times=nindex)
	}	
	
	if(is.null(namesindex)) {  namesindex<-index }
	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)
	
	if(is.null(by.factor)) {  by.factor=rep(1,length(names_sp_ind_plot)) }

	traits<-traits[order(ind.plot),]
	ind.plot<-ind.plot[order(ind.plot)]
	sp<-sp[order(ind.plot)]
	
	name_sp_sites=paste(sp, ind.plot, sep="_")
	comm=NULL
	comm<-t(table(ind.plot,1:length(ind.plot)))
	
	S = colSums(comm>0)
	ncom=length(S)
	
	#Creation of the regional pool if not inform
	if(is.null(reg.pool)) {
		reg.pool<-traits
	}
	
	if(!is.null(reg.pool) & sum(nullmodels==2)==0) {
		warnings("Custom regional pool'reg.pool' is only used in the case of null model 2")
	}
	
	if(is.numeric(nperm)){
		######################################### 
		#### 	  Calcul of null models  	 ####
		######################################### 
		#Creation of three null models 
		if(printprogress==T){ print("creating null models")}
		
		if(sum(nullmodels==1)>0){
			#________________________________________
			#modèle nul 1: permutations des valeurs de traits des individus dans la communauté   
			
			traits.nm1<-list()
			
			for (t in 1: ntr){	
				traits.nm1[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits[ind.plot==levels(ind.plot)[s], t], table(ind.plot)[s])
					}
					
					traits.nm1[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.1",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==2)>0){	
			#________________________________________
			#modèle nul 2: permutations des valeurs de traits des individus de la région    
			traits.nm2<-list()
			
			for (t in 1: ntr){	
				traits.nm2[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(reg.pool[, t], table(ind.plot)[s])
					}
					
					traits.nm2[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.2",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==3)>0){
			#________________________________________  
			#modèle nul 3: permutations des espèces au niveau de la région   
			traits.nm3<-list()
			traits_by_sp<-apply(traits,2,function(x) tapply(x,name_sp_sites,mean, na.rm=T))  
			traits_by_pop<-traits_by_sp[match(name_sp_sites,rownames(traits_by_sp)),]
			
			for (t in 1: ntr){	
				traits.nm3[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits_by_pop, table(ind.plot)[s])
					}
					
					traits.nm3[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 		
				if(printprogress==T){
					print(paste("nm.3",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		######################################## 
		####	 Calcul of random values   	####
		######################################## 
		Null<-list()
		
		if(printprogress==T){print("calcul of null values using null models")}
		
		for(i in 1:nindex){
		
			if(nullmodels[i]==1){nm<-array(unlist(traits.nm1),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else if(nullmodels[i]==2){nm<-array(unlist(traits.nm2),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else if(nullmodels[i]==3){nm<-array(unlist(traits.nm3),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else{print("nullmodels need 1, 2 or 3")}
			
			nm_n<-nm[,,n]
			colnames(nm_n)<-rownames(comm)
			rownames(nm_n)<-colnames(traits)
			
			functionindex= eval(index[i])
			
			dim2<-dim(by(t(nm_n), by.factor, function (x) eval(parse(text=functionindex))))[1]
			Null[[eval(namesindex[i])]] <- array(NA, dim=c(dim2, nperm) )
				
			if(is.null(dim2)) {
				Null[[eval(namesindex[i])]] <- array(NA, dim=c(1, nperm) )
			}
				
			for(n in 1:nperm){
				Null[[eval(namesindex[i])]][,n]  <- as.vector(by(t(nm[,,n]), by.factor, function (x) eval(parse(text=functionindex))))
			}
			
			if(printprogress==T){
				print(paste(eval(namesindex[i]), round(i/nindex*100,2),"%")) 
			} 
		}
		
	}
		  
	######################################## 
	####	Calcul of observed values	####
	######################################## 
	obs<-list()
	
	if(printprogress==T){print("calcul of observed values")}
	
	for(i in 1:nindex){
		if(nullmodels[i]==1){nm<-array(unlist(traits.nm1),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else if(nullmodels[i]==2){nm<-array(unlist(traits.nm2),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else if(nullmodels[i]==3){nm<-array(unlist(traits.nm3),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else{print("nullmodels need 1, 2 or 3")}
		
		nm_n<-nm[,,n]
		colnames(nm_n)<-rownames(comm)
		rownames(nm_n)<-colnames(traits)
		
		functionindex= eval(index[i])
		
		dim2<-dim(by(t(nm_n), by.factor, function (x) eval(parse(text=functionindex))))[1]
		obs[[eval(namesindex[i])]] <-rep(NA, times=dim2)
	
		if(nullmodels[i]==3) {
			traits.pop<-apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
			obs[[eval(namesindex[i])]] <- as.vector(by(t(traits.pop), by.factor, function (x) eval(parse(text=functionindex))))
		}
		
		else if(nullmodels[i]==1  |  nullmodels[i]==2) {
			obs[[eval(namesindex[i])]] <- as.vector(by(traits, by.factor, function (x) eval(parse(text=functionindex))))
		}
		
		obs[[eval(namesindex[i])]]<-as.vector(obs[[eval(namesindex[i])]])
		
		if(printprogress==T){
			print(paste(round(i/nindex*100,2),"%")) 
		} 
	}
			
		
	######################################## 
	####		Create results list		####
	######################################## 
	
	com.index<-list()  
	com.index$obs<-obs
		
	if(is.numeric(nperm)){
		com.index$Null<-Null
	}
	
	com.index$list.index<-list()
	com.index$list.index.t<-list()
	name.com.index_list.index<-vector()
	
	for(i in 1:nindex){
		com.index$list.index.t[[seq(1,nindex*2,by=2)[i]]]<-t(obs[[i]])
		com.index$list.index[[seq(1,nindex*2,by=2)[i]]]<-obs[[i]]
		name.com.index_list.index[seq(1,nindex*2,by=2)[i]]<-names(obs)[i]
		
		if(is.numeric(nperm)){
			com.index$list.index[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			com.index$list.index.t[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			name.com.index_list.index[seq(1,nindex*2,by=2)[i]+1]<-paste(names(Null)[i], "nm", sep="_")
		}
	}
	
	names(com.index$list.index.t)<-name.com.index_list.index
	names(com.index$list.index)<-name.com.index_list.index
		
	com.index$sites_richness<-S
	com.index$namestraits<-namestraits
		
	class(com.index)<-"com.index.multi"
	
	return(com.index)
}

as.listofindex<-function(x, namesindex=NULL) {
	
	if(class(x)!="list"){x<-list(x)}
	
	nlist<-length(x)
	nindex<-vector()
	for(i in 1: nlist){
		if(inherits(x[[i]], "Tstats")) {
			nindex[i]<-3
		}
		else if (inherits(x[[i]], "com.index") | inherits(x[[i]], "com.index.multi")) {
			nindex[i]<-length(x[[i]]$obs)
		}
		else{stop("x must be a list of objects of class Tstats, com.index or com.index.multi")}
	}
	
	res<-list()
	
	for(l in 1: nlist){
		if(inherits(x[[l]], "Tstats")) {
			for(i in c(1,5,2,6,3,7) ){
				res<-c(res, list(x[[l]][[i]]))
			}
		}
		
		else{
			for(i in 1: nindex[l]){
				res<-c(res, list(x[[l]]$obs[[i]]), list(x[[l]]$Null[[i]]) )
			}		
		}
	}
	
	if(is.null(namesindex)) {
		for(l in 1: nlist){
			namesindex<-c(namesindex, paste( rep("index", nindex[l]), l, 1:nindex[l], sep="_") )
		}
	}
	
	if(length(namesindex)==sum(nindex)) {		
		interm<-c()
		
		for(i in seq(1, sum(nindex))){
			interm<-c(interm, i, i+sum(nindex))
		}
		
		namesindex<-c(namesindex, paste(namesindex, "nm"))[interm]
	}
	
	names(res)<-namesindex
	
	class(res)<-"listofindex"
	return(res)
}

### Function to represent standardised effect size of all indices using null models
# index.list is a list of index associate with a list of corresponding null models in this order: [1] index 1 obs - [2] index 1 null model - [3] index 2 obs - [4] index 2 null model ...
#e.g. index.list<-list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm)
#observed matrix of values need to be of the same dimension
#You can transpose the observed matrix to represent either the ses by traits or by plots
plot.listofindex<-function(x, type="normal", col.index=c("red","purple","green"), add.conf=TRUE, color.cond=TRUE, val.quant=c(0.025,0.975), grid.v=TRUE, grid.h=TRUE, xlim=NULL, ylim=NULL, cex.text =0.8, plot.ask=FALSE, srt.text=90, bysite=FALSE,...){
	#possible type = "simple",  "simple_range", "normal" , "barplot" and "bytraits"
	
	if(!inherits(x, "listofindex")) {
		if(inherits(x[[1]], "Tstats") | inherits(x[[2]], "com.index")  | inherits(x[[3]], "com.index.multi")) {
			x<-as.listofindex(x)
		}	
		else{stop("x must be a list of objects of class Tstats, com.index or com.index.multi")}
	}
	
	index.list<-x
	
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
		
	ncom<-c()
	ntr<-c()
	for(i in seq(1, 2*nindex, by=2)){
		ncom<-c(ncom,dim(as.matrix(index.list[[i]]))[1])
		ntr<-c(ntr,dim(as.matrix(index.list[[i]]))[2])
	}
	
	if(is.null(namescommunity)) {warnings("rownames of index.list[[1]] is empty so names of plots cannot be plot")}
	if(is.null(namestraits)) {warnings("colnames of index.list[[1]] is empty so names of traits cannot be plot")}
	
	if(is.null(ncom)) {ncom<-dim(as.matrix(index.list[[1]]))[1]}
	if(is.null(ntr)) {ntr<-dim(as.matrix(index.list[[1]]))[2]}
	
	if(is.null(ncom)) {ncom=1}
	if(is.null(ntr)) {ntr=1}
	
	if(length(col.index)<nindex){
		col.index<-palette(rainbow(nindex))
	}
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}
	
	res<-lapply(res, function(x) lapply(x, as.matrix)   )
	
	nfactor<-c()
	for(i in 1:nindex){
		if(ntr[i]>1) nfactor<-c(nfactor, dim(as.matrix(res[[eval(namesindex[i])]]$ses))[1])
		if(ntr[i]==1) nfactor<-c(nfactor, length(res[[eval(namesindex[i])]]$ses))
	}
	
	
	#________________________________________
	#plot : possible type = "simple", "simple_range", "normal", "barplot" and "bytraits"
	
	par(mar=c(5, 7, 4, 2))
	
	#__________
	if(type=="bytraits"){
		par(mar=c(5, 4, 4, 2) + 0.1)
		
		if(is.null(ylim)) {  ylim=c(0,nindex+1) }
		res.total<-unlist(res)
		res.total<-res.total[is.finite(res.total)]
		if(is.null(xlim)) {xlim=c(min(c(res.total,-2), na.rm=T), max(c(res.total,2), na.rm=T))}
		
		if(is.logical(color.cond)) {color.cond=c("blue","orange")}
			
		if(bysite==T){
			
			if(length(unique(ncom))!=1){stop("if type = bytrait and bysite = T, all index need the same number of community")}
			
			for (s in 1: ncom[1]){
				plot(mean(res[[eval(namesindex.all[1])]]$ses[s,], na.rm=T),(1:nindex)[1] ,bty="n", cex.lab=0.8, yaxt="n", xlab=paste("SES", namescommunity[s]), ylim=ylim, xlim=xlim, pch=16, type="n", ...)
				abline(v=0)	
						
				for(i in 1:nindex){
					abline(h=(1:nindex)[i], lty=2, col="lightgray")	
						
					segments(mean(res[[eval(namesindex[i])]]$ses.sup[s,], na.rm=T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[s,], na.rm=T), (1:nindex)[i])
					
					points(mean(res[[eval(namesindex[i])]]$ses.sup[s,], na.rm=T), (1:nindex)[i], pch="|")
					points(mean(res[[eval(namesindex[i])]]$ses.inf[s,], na.rm=T), (1:nindex)[i], pch="|")
					
					points(res[[eval(namesindex[i])]]$ses[s,], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,]) ), pch="*")
									
					cond.sup<-res[[eval(namesindex[i])]]$ses[s,]>res[[eval(namesindex[i])]]$ses.sup[s,]
					points(res[[eval(namesindex[i])]]$ses[s,][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.sup]) ), pch="*", cex=3, col=color.cond[2])
					
					cond.inf<-res[[eval(namesindex[i])]]$ses[s,]<res[[eval(namesindex[i])]]$ses.inf[s,]
					points(res[[eval(namesindex[i])]]$ses[s,][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.inf]) ), pch="*", cex=3, col=color.cond[1])
					
					
					points(mean(res[[eval(namesindex[i])]]$ses[s,], na.rm=T), (1:nindex)[i], col="red", pch=16)
					
					text(1, (1:nindex)[i]+0.3, namesindex[i], cex=cex.text,  pos=4, font=2)
									
					chh <- par()$cxy[ 2 ]  ##  character height
					text(res[[eval(namesindex[i])]]$ses[s,], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,]) ), namestraits, cex=cex.text, srt=srt.text,)
					
				}
			}	
		}
		
		
		else if(bysite==F){
			if(length(unique(ntr))!=1){stop("if type = bytrait and bysite = F, all index need the same number of traits")}
			
			for (t in 1: ntr[1]){
								
				plot(mean(res[[eval(namesindex.all[i])]]$ses[,t], na.rm=T), (1:nindex)[i] ,bty="n", cex.lab=0.8, yaxt="n", xlab=paste("SES", namestraits[t]), ylim=ylim, xlim=xlim, pch=16, type="n", ...)
				abline(v=0)	
						
				for(i in 1:nindex){
					
					abline(h=(1:nindex)[i], lty=2, col="lightgray")		
					
					segments(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm=T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm=T), (1:nindex)[i])
					
					points(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm=T), (1:nindex)[i], pch="|")
					points(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm=T), (1:nindex)[i], pch="|")
					
					points(res[[eval(namesindex[i])]]$ses[,t], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), pch="*")
									
					cond.sup<-res[[eval(namesindex[i])]]$ses[,t]>res[[eval(namesindex[i])]]$ses.sup[,t]
					points(res[[eval(namesindex[i])]]$ses[,t][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.sup]) ), pch="*", cex=3, col=color.cond[2])
					
					cond.inf<-res[[eval(namesindex[i])]]$ses[,t]<res[[eval(namesindex[i])]]$ses.inf[,t]
					points(res[[eval(namesindex[i])]]$ses[,t][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.inf]) ), pch="*", cex=3, col=color.cond[1])
					
					
					points(mean(res[[eval(namesindex[i])]]$ses[,t], na.rm=T), (1:nindex)[i], col="red", pch=16)
					
					text(1, (1:nindex)[i]+0.3, namesindex[i], cex=0.8,  pos=4, font=2)
									
					chh <- par()$cxy[ 2 ]  ##  character height
					text(res[[eval(namesindex[i])]]$ses[,t], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), namescommunity, cex=cex.text, srt=srt.text,)
				}
			}	
		}
	}
	
	else if(type=="simple" | type=="simple_range" | type=="normal" | type=="barplot"){
		
		if(is.null(ylim)) { ylim=c(0,5.5+(nindex+1)*max(ntr)) }
		res.total<-unlist(res)
		res.total<-res.total[is.finite(res.total)]
		if(is.null(xlim)) {xlim=c(min(c(res.total,-2), na.rm=T), max(c(res.total,2), na.rm=T))}
		
		
		plot(0, ylab="Traits", yaxt= "n", xlab="Standardized Effect Size", ylim=ylim, xlim=xlim, col="black", type="l", ...)
		try(axis(side=2, seq(from=5.5, to=4.5+(nindex+1)*max(ntr), by=nindex+1)+(nindex+1)/2, labels=namestraits, las=1, cex.axis=0.7 ) )
		abline(v=0)
		
		if(grid.v==T) {
			range.<-max(c(res.total,2), na.rm=T)-min(c(res.total,-2), na.rm=T)
			
			vect.grid<-seq(min(res.total, na.rm=T),max(res.total, na.rm=T), by=round(range.,2)/9)
			for(j in vect.grid){
				abline(v=j, lty=2, col="lightgray")	
			}
		}
		
		if(grid.h==T) {
			for(j in seq(5.5,5.5+(nindex+1)*max(ntr))  ){
				abline(h=j, lty=2, col="lightgray")
			}
		}
		
		#__________
		if(type=="simple"){
				
			for(i in 1:nindex){
				for(t in 1:ntr[i]){
					
					if(color.cond==F){
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=nfactor[i]), pch=20, col=col.index[i])
					}
					
					if(color.cond==T){
						if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=nfactor[i]), pch=20, col=col.index[nindex+i])
						condition<-res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t]  |   res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
						condition[is.na(condition)]<-FALSE					
						
						if(sum(condition)>0){
							points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times=sum(condition)), pch=20, col=col.index[i])
						}
					}
									
					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
													
					if(add.conf==T){
						points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
						points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					}
			
					abline(seq(from=5.5, to=4.5+(nindex+1)*ntr[i], by=nindex+1)[t],b=0, lty=4, lwd=0.2)
					
				}		
			}
		}
		
		#__________
		else if(type=="simple_range"){
			
			for(i in 1:nindex){
				for(t in 1:ntr[i]){
					
					if(color.cond==F){
						points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
					}
					
					if(color.cond==T){
						if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
						points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[nindex+i])
						
						condition<-mean(res[[eval(namesindex[i])]]$ses [,t], na.rm=T) > mean(res[[eval(namesindex[i])]]$ses.sup [,t], na.rm=T)  |  mean( res[[eval(namesindex[i])]]$ses [,t], na.rm=T) < mean(res[[eval(namesindex[i])]]$ses.inf [,t], na.rm=T)
						if(sum(condition)>0){
							points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
						}
					}
					
					if(add.conf==T){
						points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
						points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					}
					
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
					
					points(min(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="*", col=col.index[i])
					points(max(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="*", col=col.index[i])
					
					abline(seq(from=5.5, to=4.5+(nindex+1)*ntr[i], by=nindex+1)[t],b=0, lty=4, lwd=0.2)
				}		
			}
		}
		
		#__________
		else if(type=="normal"){
			
			for(i in 1:nindex){
				for(t in 1:ntr[i]){
					
					if(color.cond==F){
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=nfactor[i]), pch=20, col=col.index[i])
					}
					
					if(color.cond==T){
						if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=nfactor[i]), pch=20, col=col.index[nindex+i])
						condition<-res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t]  |   res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
						condition[is.na(condition)]<-FALSE					
						
						if(sum(condition)>0){
							points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times=sum(condition)), pch=20, col=col.index[i])
						}
					}
			
					if(add.conf==T){
						points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
						points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					}
					
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=10, col=col.index[i])
					
					abline(seq(from=5.5, to=4.5+(nindex+1)*ntr[i], by=nindex+1)[t],b=0, lty=4, lwd=0.2)
				}		
			}
		}
	
		#__________
		else if(type=="barplot"){
			for(i in 1:nindex){
				for(t in 1:ntr[i]){
					
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, 0, 5.5+(nindex+1)*t-i, pch=17, col=col.index[i], lwd=8)
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
									
					if(add.conf==T){
						points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i], cex=2)
						points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i], cex=2)
					}
					
					abline(seq(from=5.5, to=4.5+(nindex+1)*ntr[i], by=nindex+1)[t],b=0, lty=4, lwd=0.2)
			
				}		
			}
			
		}
		legend("bottom", inset=.005, namesindex, fill=col.index, ncol=round(nindex/3) ,cex=0.6, bty="0")
	}
	
	else{print(paste("Error:",type,"is not a valid type of plot"))}
	
	par(mar=c(5, 4, 4, 2) + 0.1) #return to default parameter
	par(oldpar) #return to default parameter
}





#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Decomposition of variance

decomp_within<-function(traits=NULL , formula=~ 1, ind.plot=NULL, sp=NULL, printprogress=TRUE, ...) {
	
	form.string<-deparse(substitute(formula))
	
 	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)

	comm<-t(table(ind.plot,1:length(ind.plot)))
	S = colSums(comm>0)
	ncom=length(S)
	
	decomp<-list()
	
	moy_pop<-apply(traits, 2, function(x) tapply(x, paste(ind.plot, sp, sep="."), mean, na.rm=T)) 
	moy_sp<-apply(traits, 2, function(x) tapply(x, sp, mean, na.rm=T))
	

	
	for(t in 1:ntr){
		interm<-apply(comm, 2, function(x) tapply(x, paste(ind.plot, sp, sep=".") , sum, na.rm=T ) )[!is.na(moy_pop[,t]), 1:ncom]
		interm2<- apply(comm, 2, function(x) tapply(x, sp, sum, na.rm=T ) )[!is.na(moy_sp[,t]), 1:ncom]
		
		specif.avg<-t(moy_pop[,t][!is.na(moy_pop[,t])] %*% interm )/colSums(interm)
		const.avg<-t(moy_sp[,t][!is.na(moy_sp[,t])] %*%  interm2 )/colSums(interm2)
		
		flex= paste("trait.flex.anova(", form.string , ", specif.avg, const.avg, ...)", sep="")
		
 		decomp[[eval(namestraits[t])]]<-as.vector(eval(parse(text=flex)))
		
		if(printprogress==T){print(paste(round(t/ntr*100, 2),"%"))} else {}
	}
	
	decomp$traits<-namestraits
	class(decomp) <-"decomp_within"
	return(decomp)
}

barplot.decomp_within<-function(height, resume=TRUE, ...) { 
	
	x<-height
	
	oldpar<-par(no.readonly = TRUE)
	
	if(resume==TRUE){
		res<-matrix(ncol= dim(as.matrix(x))[1]-1, nrow=dim(x[[1]]$SumSq)[2]-1)
		for(i in 1: (dim(as.matrix(x))[1]-1)) {
			res[,i]<-plot.trait.flex(x[[i]], plot.res=FALSE, cumul=TRUE, beside=T, ...)
		}  
		colnames(res)<-x$traits
		barplot(res, beside=T)
		legend("topright", c("Turnover", "Intraspec.", "Covariation"), fill=c(gray.colors(3)) )
	}
	
	else if(resume==FALSE){
		par(cex=2/length(x$traits))
		
		for(i in 1:dim(x[[1]]$SumSq)[2]){
			plot.trait.flex(x[[i]], main=x$traits[i], ...)
			abline(v=0, lty=2)
		}  
	}
	
	par(oldpar)

}

# Leps et al. function
trait.flex.anova <-  function(formula, specif.avg, const.avg, ...)  {
    # Formula into string form
    form.string<-deparse(substitute(formula))
    
    # Check formula parameter
    form.parts<-unlist(strsplit(form.string,"~"))
    if(length(form.parts) != 2)
      stop("First parameter must be valid one-sided formula, like ~A*B");
    if(nchar(form.parts[1])>0)
      warning("Left side of the formula was ignored!");
    
    # The two average variables into string form
    spec.av <- deparse(substitute(specif.avg))
    cons.av <- deparse(substitute(const.avg))
    
    test.has.onelevel<-function(aov.summ)
    {
      (length(aov.summ) == 1);
    }
    
    test.has.resid<-function(aov.one)
    {
      if(class(aov.one)[1] != "anova")
        warning("specified object is not aov result!");
      nrows <- dim(aov.one)[1]
      if(nrows == 0)
        return( FALSE);
      last.row.lbl <- dimnames(aov.one)[[1]][nrows]
      if(unlist(strsplit(last.row.lbl, " "))[1] != "Residuals")
        return( FALSE);
      if(dim(aov.one)[2] < 5)  # P and F are missing
        return( FALSE);
      TRUE;
    }
    
    # Specific averages ANOVA
    form.1 <- as.formula(
      paste(spec.av,form.parts[2],sep="~"))
    res.1 <- summary( aov(form.1, ...))
    if(test.has.onelevel( res.1) == FALSE)
      stop("Cannot evaluate ANOVAs with multiple error levels!")
    res.1 <- res.1[[1]]
    if(test.has.resid( res.1) == FALSE)
      stop("No residual DFs left, cannot continue");
    
    # Constant averages ANOVA
    form.2 <- as.formula(
      paste(cons.av,form.parts[2],sep="~"))
    # no need to test for multilevels by now
    res.2 <- summary( aov(form.2, ...))[[1]]
    if(test.has.resid( res.2) == FALSE)
      stop("No residual DFs left in constant averages ANOVA, cannot continue");
    
    
    # Now the differences:
    spec.const.diff <- paste("I(", spec.av, "-", cons.av, ")", sep="")
    
    form.3 <- as.formula(
      paste(spec.const.diff,form.parts[2],sep="~"))
    # no need to test for multilevels by now
    res.3 <- summary( aov(form.3, ...))[[1]]
    if(test.has.resid( res.3) == FALSE)
      stop("No residual DFs left in (specif-const) ANOVA, cannot continue");
    
    
    if((dim(res.1) != dim(res.2)) || (dim(res.1) != dim(res.3)))
      stop("Tables from the three ANOVAs have incompatible sizes");
    
    # Create sum of squares table: add SS(Tot) except for null models    
    nrows <- dim(res.1)[1]
    ss.turn <- res.2[,2]
    ss.var  <- res.3[,2]
    ss.tot  <- res.1[,2]
    ss.covar<- ss.tot - ss.turn - ss.var
    ss.row.names <- dimnames(res.1)[[1]]
    if(nrows > 1)
    {
      ss.turn <- c(ss.turn, sum(ss.turn))
      ss.var  <- c(ss.var,  sum(ss.var))
      ss.tot  <- c(ss.tot,  sum(ss.tot))
      ss.covar<- c(ss.covar,sum(ss.covar))
      ss.row.names <- c(ss.row.names, "Total")
      nrows   <- nrows + 1
    }
    else
    {
      # replace row title
      ss.row.names[1] <- "Total"
    }
    SS.tab <- data.frame( Turnover = ss.turn, Intraspec. = ss.var,
                          Covariation = ss.covar, Total = ss.tot,
                          row.names = ss.row.names)
    # Calculate relative fractions
    TotalSS <- SS.tab[nrows, 4] # lower right corner
    SS.tab.rel <- SS.tab / TotalSS
    
    # Collect significances
    if(nrows > 1)  # get rid of the "Total" label again
      ss.row.names <- ss.row.names[-nrows]
    P.tab <- data.frame( Turnover = res.2[,5], Intraspec. = res.3[,5],
                         Total = res.1[,5], row.names = ss.row.names)
    
    res <- list( SumSq=SS.tab, RelSumSq=SS.tab.rel, Pvals=P.tab, 
                 anova.turnover=res.2, anova.total=res.1, anova.diff=res.3)
    class(res)<- "trait.flex"
    res
  }

print.trait.flex <-  function(x, ...)   {
    op <- options();
    cat("\n Decomposing trait sum of squares into composition turnover\n");
    cat(  " effect, intraspecific trait variability, and their covariation:\n");
    options(digits=5)
    print(x$SumSq, ...);
    cat("\n Relative contributions:\n");
    options(digits=4)
    print(x$RelSumSq, ...)
    nPvals <- dim(x$Pvals)[1]
    if(nPvals > 1)
    {
      cat("\n Significance of testable effects:\n");
      options(digits=5)
      print(x$Pvals[-nPvals,], ...);
    }
    options(op)
    invisible(x)
  }

plot.trait.flex <-  function(x, plot.total = FALSE, use.percentage = TRUE, plot.covar = FALSE, cumul = FALSE, legend.pos = if(plot.total) "topleft" else "topright", plot.res=TRUE, ...) { 
    if(use.percentage)
      SumSq <- 100 * x$RelSumSq
    else
      SumSq <- x$SumSq
    if(legend.pos == "none")
      legend.txt <- NULL
    else
      legend.txt <- colnames(SumSq)[-4]
    
    nrows <- dim(SumSq)[1]
    plot.tab <- as.matrix(SumSq)
    
    if(nrows > 1)
    {
      if(plot.covar)
      {
        if(plot.total)
          plot.tab <- plot.tab[,-4]
        else
          plot.tab <- plot.tab[-nrows,-4]
      }
      else
      {
        if(plot.total)
          plot.tab <- plot.tab[,1:2]
        else
          plot.tab <- plot.tab[-nrows, 1:2]   
        if(legend.pos != "none")
          legend.txt <- legend.txt[1:2]
      }
    }
    
    else
    {
		if(!plot.total){
			plot.tab <- t(plot.tab[,-4])
        }
        
		legend.pos <- "none"
		legend.txt <- NULL
    }
    
    if(plot.res==TRUE){
	    if(!cumul){
	    
		    if(use.percentage)
		    {
		      xpos <- barplot( plot.tab, horiz=T,
		                       ylab = "Explained variation (%)", ...) 
		    }
		    else
		      xpos <- barplot( plot.tab=T, horiz=T, 
		                       ylab = "Sum of squares of analysed trait", ...) 
		    
		    if(plot.total==FALSE){
				abline(v=as.matrix(SumSq)[,4])
		    }
		    
	    }
			
			
		if(cumul){
	    
		    if(use.percentage)
		    {
		      xpos <- barplot( t(plot.tab), horiz=T,  
		                       ylab = "Explained variation (%)", ...) 
		    }
		    else
		      xpos <- barplot( t(plot.tab), horiz=T, 
		                       ylab = "Sum of squares of analysed trait", ...) 
		    
		    if(plot.total==FALSE){
				abline(v=as.matrix(SumSq)[,4], lty=2)
		    }
		    
		    if(min(plot.tab)<0){
				abline(v=0)
			}
	    }
	}
	 
	 
    if(nrows > 1)
    {
      if(!plot.covar)
      { if(length(xpos) > 1)
        line.half <- 0.4*(xpos[2] - xpos[1])
        else
          line.half <- 0.4*(xpos[1])
        segments( xpos-line.half, SumSq[,4],
                  xpos+line.half, SumSq[,4], lwd=3)
        if(legend.pos != "none")
        { NL <- length(legend.txt)       
          legend( legend.pos, legend=c(legend.txt,"Total"),
                  fill=c(gray.colors(NL), NA), 
                  lty=c(rep( 0,NL),1),
                  lwd=c(rep( 0,NL),3))
        }
      }
      else
      {
        if(legend.pos != "none")
          legend( legend.pos, legend=legend.txt, 
                  fill=gray.colors(length(legend.txt)))
      }
    }
    
    if(!plot.res){
		return(plot.tab)
    }
}





#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Other plotting functions

### Function to represent correlations between Tstats
plot_cor.Tstats<-function(tstats=NULL, val.quant=c(0.025,0.975), add.text=FALSE, bysite=FALSE, col.obj=NULL, plot.ask=TRUE, multipanel=TRUE, ...) {
	
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	#________________________________________
	ses.T_IP.IC.moy<-t(colMeans((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_IC.IR.moy<-t(colMeans((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_PC.PR.moy<-t(colMeans((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	
	ses.T_IP.IC<-t((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR<-t((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR<-t((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	ses.T_IP.IC.inf<-t((apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR.inf<-t((apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR.inf<-t((apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	ses.T_IP.IC.sup<-t((apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR.sup<-t((apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR.sup<-t((apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
	
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	#________________________________________
	if(bysite==F){
		
		if(is.null(col.obj)) {col.obj<-rainbow(dim(tstats$T_IP.IC)[2])}
		else{}
		if(multipanel) {par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[2])+1,sqrt(dim(tstats$T_IP.IC)[2])+1)) }
		
		#__________
		#First panel of figures

		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_IC.IR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")		
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col="grey", pch=20, main=rownames(ses.T_IC.IR)[t], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			segments(rep(rowMeans(ses.T_IC.IR, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]) ,ses.T_IC.IR[t,],  ses.T_IP.IC[t,], col=col.obj[t])
		}
		
		plot.new()
		
		#__________
		#Second panel of figures

		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")		
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=rownames(ses.T_PC.PR)[t], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]) ,ses.T_PC.PR[t,],  ses.T_IP.IC[t,], col=col.obj[t])
		}
		
		plot.new()
		
		#__________
		#Third panel of figures
		
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")	
				
		for(t in 1:dim(tstats$T_IC.IR)[2]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=rownames(ses.T_PC.PR)[t], ...)
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm=T)[t], times=dim(tstats$T_IC.IR)[1]), rep(rowMeans(ses.T_IC.IR, na.rm=T)[t], times=dim(tstats$T_IC.IR)[1]) ,ses.T_PC.PR[t,],  ses.T_IC.IR[t,], col=col.obj[t])
		}
	}

	#________________________________________
	else if(bysite==T){
	
		if(is.null(col.obj)) {col.obj<-rainbow(dim(tstats$T_IP.IC)[1])}
		else{}
		
		if(multipanel) {par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[1])+1,sqrt(dim(tstats$T_IP.IC)[1])+1))}

		#__________
		#First panel of figures

		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_IC.IR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")
		
		for(s in 1:dim(tstats$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			segments(rep(colMeans(ses.T_IC.IR, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]) ,ses.T_IC.IR[,s],  ses.T_IP.IC[,s], col=col.obj[s])
		}
		
		#__________
		#Second panel of figures
		
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		
		text(0,0,"null \r\n model \r\n zone")
				
		for(s in 1:dim(tstats$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s], ... )
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]) ,ses.T_PC.PR[,s],  ses.T_IP.IC[,s], col=col.obj[s])
		}
		
		#__________
		#Third panel of figures

		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		
		text(0,0,"null \r\n model \r\n zone")
				
		for(s in 1:dim(tstats$T_IC.IR)[1]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s], ... )
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm=T)[s], times=dim(tstats$T_IC.IR)[2]), rep(colMeans(ses.T_IC.IR, na.rm=T)[s], times=dim(tstats$T_IC.IR)[2]) ,ses.T_PC.PR[,s],  ses.T_IC.IR[,s], col=col.obj[s])
		}	
	}
	
	else{print("Error: obj need to be either traits or sites")}
	
	par(oldpar)
}

# plot ses of an index against an other variable wich correspond to plot. For example species richness or a gradient variable

plot_ses.var<-function(index.list, variable=NULL, ylab="variable" ,color.traits=NULL, val.quant=c(0.025,0.975), resume=FALSE, multipanel=TRUE){

	y<-variable
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
	
	ncom<-dim(index.list[[1]])[1]
	ntr<-dim(index.list[[1]])[2]
	
	if(is.null(color.traits)){
		color.traits<-palette(terrain.colors(ntr))
	}
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}

	oldpar<-par(no.readonly = TRUE)
	if(multipanel) {par(mfrow=c(ceiling(sqrt(nindex))-1, ceiling(sqrt(nindex))))}
	
	ylim=c(min(y, na.rm=T), max(y, na.rm=T))
		
		
	for(i in seq(1,nindex*2, by=2)){
		if(resume==FALSE){xlim=c(min(c(-4,res[[eval(namesindex.all[i])]]$ses), na.rm=T), max(c(4,res[[eval(namesindex.all[i])]]$ses), na.rm=T))}
		else{xlim=c(min(c(-4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T)), na.rm=T), max(c(4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T)), na.rm=T))  }
		
		plot(0, 0 ,bty="n", cex.lab=0.8, xlab=paste("SES", namesindex.all[i]), ylab=ylab, ylim=ylim, xlim=xlim, pch=16, type="n")		
		abline(v=0, lty=1, col="black")	
		
		
		if(resume==TRUE){
			points(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T), y, pch=16)
			segments(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), y, rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), y, pch=16)
		}
		
		else{
			for(nco in 1:ncom){
				abline(h=y[nco], lty=2, pch=0.5, col="lightgray")
			}
			
			for (t in 1: ntr){
				points(res[[eval(namesindex.all[i])]]$ses[,t] , y, pch=16, col=color.traits[t])		
			}
		}	
	}
	
	if(resume!=TRUE){
		plot(0, 0 ,bty="n", cex.lab=0.8, xlab=paste("SES", namesindex.all[i]), ylim=ylim, xlim=xlim, pch=16, type="n")
		legend("center", legend=namestraits, fill=color.traits, bty="n", ncol= round(sqrt(nlevels(as.factor(namestraits)))-1 ) )
	}
	
	par(oldpar)
	
}


plot_dens<-function(traits=NULL, var.1=NULL, var.2=NULL, col.dens=NULL, plot.ask=TRUE, ylim.cex=1, cex.leg=0.8, polyg=TRUE, multipanel=TRUE, leg=TRUE)  {

	var.1<-as.factor(as.vector(var.1))
	var.2<-as.factor(as.vector(var.2))
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	namestraits<-colnames(traits)
	namescommunity<-unique(var.1)
	ncom<-length(namescommunity)
	
	ntr<-dim(traits)[2]
	
	if(is.null(col.dens)) { col.dens<-rainbow(nlevels(as.factor(var.2))) }
	
	if(multipanel) {
		par(mfrow=c(2,2))
	}
	
	
	for(co in 1:ncom){
		for(t in 1:ntr){
			if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t]))>1){
			
				plot(main= paste(namestraits[t],levels(as.factor(var.1))[co], " "), density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T), ylim=c(0,max(density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$y)*ylim.cex), xlim=c(min(traits[,t], na.rm=T),max(traits[,t], na.rm=T)), col="black")
				lines(density(traits[,t], na.rm=T), lty=2, col="grey")
				
				if(leg){
					if(mean(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)  <   mean(traits[,t], na.rm=T)  )  {pos= "topright"}
					else{pos= "topleft"}
					legend(pos, inset=0.05,levels(as.factor(var.2)), fill=col.dens, cex=cex.leg, bty="n", ncol= round(sqrt(nlevels(as.factor((var.2))))-1 ) )
				}
				
				if (polyg==T) {
					x<-density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$x
					y<-density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$y
					polygon(c(x,rev(x)), c(rep(0,length(x)),rev(y)), border=NA, col=rgb(0.5,0.5,0.5,0.7))
				}
				
				for(s in 1:nlevels(as.factor(var.2))) {
					if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]))>1) 
					{lines(density(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t], na.rm=T), col=col.dens[s])}
				
					else if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]))==1) 
					{points(0,na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]), col=col.dens[s])}
				}
			}
			
		}
	}
	par(oldpar)
}

# Ackerly & Cornwell 2007
# Plot populations values against species values
plot_sp_pop<-function(traits=NULL, ind.plot=NULL, sp=NULL, col.ind = rgb(0.5,0.5,0.5,0.5), col.pop=NULL, col.sp=NULL, col.site=NULL, resume=FALSE, p.val=0.05, min.ind.signif=10 , multipanel=TRUE, col.nonsignif.lm=rgb(0,0,0,0.5), col.signif.lm=rgb(1,0.1,0.1,0.8), silent=FALSE) {

	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)
	
	traits<-traits[order(ind.plot),]
	ind.plot<-ind.plot[order(ind.plot)]
	sp<-sp[order(ind.plot)]
	
	name_sp_sites=paste(sp, ind.plot, sep="_")
	comm<-t(table(ind.plot,1:length(ind.plot)))
	
	S = colSums(comm>0)
	ncom=length(S)

	
	plotsp=unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(3,3*nlevels(as.factor(name_sp_sites)), by=3)]
	#plosp is the plot in wich the population is
	plotind=unlist(strsplit(name_sp_sites,split="_"))[seq(3,3*length(name_sp_sites), by=3)]
	spplot=paste(unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(1,3*nlevels(as.factor(name_sp_sites)), by=3)], unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(2,3*nlevels(as.factor(name_sp_sites)), by=3)], sep="_")
	
	traits_by_pop<-apply(traits,2,function(x) tapply(x, name_sp_sites,mean, na.rm=T))  
	
	traits_by_sites<-apply(traits,2,function(x) tapply(x, ind.plot,mean, na.rm=T))  
	
	
	if(multipanel){
		par(mfrow = c(sqrt(ntr), sqrt(ntr)) )
	}
  
	
	if(is.null(col.sp)){
		col.sp<-rainbow(nlevels(sp))
	}
	
	if(length(col.sp)<nlevels(sp)){
		col.sp<-rep(col.sp ,length.out=nlevels(sp))
	}
	
	if(is.null(col.pop)){
		col.pop<-col.sp[match(spplot, levels(sp))]
	}
	
	if(is.null(col.site)){
		col.site<-rep(rgb(0.1, 0.1, 0.1 , 0.8),ncom)
	}
	
	if(!resume){
		for(t in 1:ntr){
			x.ind<-traits_by_sites[match(plotind,rownames(traits_by_sites)),t]
			y.ind<-traits[,t]
			plot(x.ind, y.ind, pch=16, col=col.ind, cex=0.5)
			
			x.pop<-traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop<-traits_by_pop[,t]
			points(x.pop, y.pop, pch=16,  col=col.pop)
					
			for(s in 1:nlevels(sp)){
				
				try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)  
				
				options(warn=-1)
				if(class(try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)) =="lm"){
					if(!is.na(summary(interm)$coefficient[,4])){
						if(summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm=1
							lwd.lm=3
							color.lm<-col.signif.lm
						}
						else{		
							lty.lm=0
							lwd.lm=0
						}
						options(warn=0)
					}
				
				   
					else{
						lty.lm=3
						lwd.lm=1
						color.lm<-col.nonsignif.lm
					}
			
					if(!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty=lty.lm, lwd=lwd.lm, col=color.lm )
					}
				}
			}
			
			options(warn=-1)
			
			try( interm2<-lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)
			
			color.lm2<-col.nonsignif.lm
			
			if(class(try(lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)) =="lm"){
				if(!is.na(summary(interm2)$coefficient[,4])){
					if(summary(interm2)$coefficients[2,4]<p.val & length(interm2$fitted.values)>min.ind.signif){
						color.lm2<-col.signif.lm
					}
					else{}
					options(warn=0)
				}
			}
			
			
			points(tapply( traits_by_pop[,t], plotsp ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], plotsp ,mean, na.rm=T) , col=col.site, pch="*", cex=3)
		
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm=T) , col=col.sp, pch=16, cex=1.5)
			abline(try(lm( traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent) , lty=2, lwd=2, col=color.lm2)
		
		}
	}
	
	if(resume){
		for(t in 1:ntr){
			x.pop<-traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop<-traits_by_pop[,t]
			plot(x.pop, y.pop, pch=16,  col=col.pop)
			abline(a=0, b=1, lty=3, lwd=2)
					
			for(s in 1:nlevels(sp)){
				
				try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)  
				
				options(warn=-1)
				if(class(try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)) =="lm"){
					if(!is.na(summary(interm)$coefficient[,4])){
						if(summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm=1
							lwd.lm=3
						}
						else{
							lty.lm=0
							lwd.lm=0
						}
						options(warn=0)
					}
				
					else{
						lty.lm=3
						lwd.lm=1
					}
			
					if(!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty=lty.lm, lwd=lwd.lm, col=col.sp[s] )
					}
				}
			}		
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm=T) , col=col.sp, pch=16, cex=1.5)
		}
	}
	
	
	par(mfrow = c(1, 1))

}



# Plot populations values against environmental variable(s)
plot_sp_var<-function(traits=NULL, ind.plot=NULL, sp=NULL, variable = NULL, col.ind = rgb(0.5,0.5,0.5,0.5), col.pop=NULL, col.sp=NULL, col.site=NULL, resume=FALSE, p.val=0.05, min.ind.signif=10 , multipanel=TRUE, col.nonsignif.lm=rgb(0,0,0,0.5), col.signif.lm=rgb(1,0.1,0.1,0.8), silent=FALSE) {

	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)
	
	traits<-traits[order(ind.plot),]
	ind.plot<-ind.plot[order(ind.plot)]
	sp<-sp[order(ind.plot)]
	
	name_sp_sites=paste(sp, ind.plot, sep="_")
	comm<-t(table(ind.plot,1:length(ind.plot)))
	
	S = colSums(comm>0)
	ncom=length(S)

	
	plotsp=unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(3,3*nlevels(as.factor(name_sp_sites)), by=3)]
	#plosp is the plot in wich the population is
	plotind=unlist(strsplit(name_sp_sites,split="_"))[seq(3,3*length(name_sp_sites), by=3)]
	spplot=paste(unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(1,3*nlevels(as.factor(name_sp_sites)), by=3)], unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[seq(2,3*nlevels(as.factor(name_sp_sites)), by=3)], sep="_")
	
	traits_by_pop<-apply(traits,2,function(x) tapply(x, name_sp_sites,mean, na.rm=T))  
	
	traits_by_sites<-variable
	
	if(is.vector(traits_by_sites)){
		traits_by_sites<-matrix(rep(traits_by_sites, times=ntr), ncol=ntr)
	}
	
	interm.for.names<-apply(traits,2,function(x) tapply(x, ind.plot,mean, na.rm=T)) 
	colnames(traits_by_sites)<-colnames(interm.for.names)
	
	rownames(traits_by_sites)<-rownames(interm.for.names)
	
	
	
	if(multipanel){
		par(mfrow = c(sqrt(ntr), sqrt(ntr)) )
	}
  
	
	if(is.null(col.sp)){
		col.sp<-rainbow(nlevels(sp))
	}
	
	if(length(col.sp)<nlevels(sp)){
		col.sp<-rep(col.sp ,length.out=nlevels(sp))
	}
	
	if(is.null(col.pop)){
		col.pop<-col.sp[match(spplot, levels(sp))]
	}
	
	if(is.null(col.site)){
		col.site<-rep(rgb(0.1, 0.1, 0.1 , 0.8),ncom)
	}
	
	if(!resume){
		for(t in 1:ntr){
			x.ind<-traits_by_sites[match(plotind,rownames(traits_by_sites)),t]
			y.ind<-traits[,t]
			plot(x.ind, y.ind, pch=16, col=col.ind, cex=0.5)
			
			x.pop<-traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop<-traits_by_pop[,t]
			points(x.pop, y.pop, pch=16,  col=col.pop)
					
			for(s in 1:nlevels(sp)){
				
				try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)  
				
				options(warn=-1)
				if(class(try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)) =="lm"){
					if(!is.na(summary(interm)$coefficient[,4])){
						if(summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm=1
							lwd.lm=3
							color.lm<-col.signif.lm
						}
						else{		
							lty.lm=0
							lwd.lm=0
						}
						options(warn=0)
					}
				
				   
					else{
						lty.lm=3
						lwd.lm=1
						color.lm<-col.nonsignif.lm
					}
			
					if(!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty=lty.lm, lwd=lwd.lm, col=color.lm )
					}
				}
			}
			
			options(warn=-1)
			
			try( interm2<-lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)
			
			color.lm2<-col.nonsignif.lm
			
			if(class(try(lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)) =="lm"){
				if(!is.na(summary(interm2)$coefficient[,4])){
					if(summary(interm2)$coefficients[2,4]<p.val & length(interm2$fitted.values)>min.ind.signif){
						color.lm2<-col.signif.lm
					}
					else{}
					options(warn=0)
				}
			}
			
			
			points(tapply( traits_by_pop[,t], plotsp ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], plotsp ,mean, na.rm=T) , col=col.site, pch="*", cex=3)
		
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm=T) , col=col.sp, pch=16, cex=1.5)
			abline(try(lm( traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent) , lty=2, lwd=2, col=color.lm2)
		
		}
	}
	
	if(resume){
		for(t in 1:ntr){
			x.pop<-traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop<-traits_by_pop[,t]
			plot(x.pop, y.pop, pch=16,  col=col.pop)
			abline(a=0, b=1, lty=3, lwd=2)
					
			for(s in 1:nlevels(sp)){
				
				try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)  
				
				options(warn=-1)
				if(class(try( interm<-lm(y.pop[spplot==levels(sp)[s]] ~ x.pop[spplot==levels(sp)[s]]), silent = silent)) =="lm"){
					if(!is.na(summary(interm)$coefficient[,4])){
						if(summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm=1
							lwd.lm=3
						}
						else{
							lty.lm=0
							lwd.lm=0
						}
						options(warn=0)
					}
				
					else{
						lty.lm=3
						lwd.lm=1
					}
			
					if(!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty=lty.lm, lwd=lwd.lm, col=col.sp[s] )
					}
				}
			}		
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm=T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm=T) , col=col.sp, pch=16, cex=1.5)
		}
	}
	
	
	par(mfrow = c(1, 1))

}



#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Other functions

#Calcul of CVNND for one trait with our without division by the range of the trait
CVNND<-function(trait, div_range=FALSE){
	
	r=sort(trait)
	if(length(r)<1){
	nnd=NA}
	
	else{nnd=vector(length=length(r)-1)
		for(j in 2:length(r)){
			nnd[j-1]=r[j]-r[j-1]
		}
	}
	
	CVNND<-sd(nnd,na.rm=T)/mean(nnd, na.rm=T)
	
	if (div_range==T) {CVNND<-CVNND/range(trait)} 
	else {}
	
	return(CVNND)
}

### Function to calcul SES on list of index
ses.listofindex<-function(index.list=NULL, val.quant=c(0.025,0.975) ){
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
	
	ncom<-dim(index.list[[1]])[1]
	ntr<-dim(index.list[[1]])[2]
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}
	
	class(res)<-"ses.list"
	return(res)
}

### Function to calcul SES
ses<-function(obs=NULL, nullmodel=NULL, val.quant=c(0.025,0.975) ){
		
	if(is.vector(obs)){
		obs<-as.matrix(obs)
	}
	
	if(length(dim(obs))!=2 ) {
		obs<-as.matrix(obs)
	}
	
	if(dim(obs)[1]==dim(obs)[2]) {
		warnings("Observed matrix have the same number of rows and columns. The function is not able to detect automatically the correspondance between dimension of observed matrix and null model. You need to be sure that the null model is in the form of an array within the first and second dimension corresespond respectively to the first and second dimension of the observed matrix and the third dimension correspond to permutations")
		cond=c(1,2)
	}
	
	else{
	
		if (class(nullmodel)=="list"){
			if (class(nullmodel[[1]])=="list"){
				nullmodel<-array(unlist(nullmodel), dim=c(nrow(nullmodel[[1]][[1]]),ncol(nullmodel[[1]][[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]][[1]])/ncol(nullmodel[[1]][[1]])))
			}
			
			else {nullmodel<-array(unlist(nullmodel), dim=c(nrow(nullmodel[[1]]),ncol(nullmodel[[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]])/ncol(nullmodel[[1]])))}
		}
		
		if (class(obs)=="list"){
			obs<-matrix(obs[[1]], nrow=nrow(obs[[1]]), ncol=ncol(obs[[1]]))
		}
		
		if(!is.null(dim(obs))) {	
			cond<-c(NA,NA)
			
			if(dim(obs)[1]==dim(nullmodel)[1]){
			cond[1]<-1
			}
			
			if(dim(obs)[1]==dim(nullmodel)[2]){
			cond[1]<-2
			}
			
			if(length(dim(nullmodel))==3){
				if(dim(obs)[1]==dim(nullmodel)[3]){
				cond[1]<-3
				}
			}
			
			if(dim(obs)[2]==dim(nullmodel)[1]){
			cond[2]<-1
			}
			
			if(dim(obs)[2]==dim(nullmodel)[2]){
			cond[2]<-2
			}
			
			if(length(dim(nullmodel))==3){
				if(dim(obs)[2]==dim(nullmodel)[3]){
				cond[2]<-3
				}	
			}
		}
	}

	cond<-na.omit(cond)
	
	res<-list()
	res$ses<-(obs-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))
	res$ses.inf<-(apply(nullmodel, cond, function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))
	res$ses.sup<-(apply(nullmodel, cond, function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))
	
	class(res)<-"ses"
	return(res)
}

RaoRel<-function(sample, dfunc, dphyl, weight=FALSE, Jost=FALSE, structure=NULL)   {
	
	
	####function Qdecomp by by VillÈger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller #####
	
	Qdecomp = function(functdist,abundances, w=TRUE) {
	
		# number and names of local communities
		c<-dim(abundances)[1] ; namescomm<-row.names(abundances)
		abundances<-as.matrix(abundances)
		
		# if necessary, transformation of functdist into matrix object
		if (is.matrix(functdist)==F) functdist<-as.matrix(functdist)
		
		# checking 'abundances' and 'functdist' dimensions
		if (dim(functdist)[1]!=dim(functdist)[2])  stop("error : 'functdist' has different number of rows and columns")
		if (dim(abundances)[2]!=dim(functdist)[1]) stop("error : different number of species in 'functdist' and 'abundances' ")
		
		# checking NA absence in 'functdist'
		if (length(which(is.na(functdist)==T))!=0)  stop("error : NA in 'functdist'")
		
		# replacement of NA by 0 in abundances
		if (is.na(sum(abundances))==T)  {
		for (i in 1:dim(abundances)[1])
		for (j in 1:dim(abundances)[2] )
		{ if(is.na(abundances[i,j])==T) abundances[i,j]<- 0 } # end of i j
		} # end of if
		
		#  species richness and total abundances in local communities
		abloc<-apply(abundances,1,sum)
		nbsploc<-apply(abundances,1,function(x) {length(which(x>0))} )
		
		# relative abundances inside each local community
		locabrel<-abundances/abloc
		
		# alpha diversity
		Qalpha=apply(locabrel, 1, function(x) t(x) %*%  functdist %*% x)
		
		#Wc
		#Wc = abloc/sum(abloc)
		relsp<-apply(abundances,1,max)
		Wc = relsp/sum(relsp)
		
		# abundance-weighted mean alpha
		mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
		#mQalpha<-as.numeric(Qalpha%*%relsp/sum(relsp) )
		
		#Villeger's correction
		if(w==T) {
			# abundance-weighted mean alpha
			mQalpha<-as.numeric(Qalpha%*%relsp/sum(relsp) )
			#mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
			#totabrel<-apply(abundances,2,sum)/sum(abundances) 
			totabrel<-apply(locabrel*Wc, 2, sum)
			Qalpha = Qalpha*Wc
		}	
		
		# Rao's original definition: mean of Pi
		else {
			mQalpha<-mean(Qalpha)
			totabrel<-apply(locabrel,2,mean)  
			}
		
		
		# gamma diversity
		Qgamma<-( totabrel %*% functdist %*% totabrel ) [1]
		
		# beta diversity
		Qbeta<-as.numeric( Qgamma-mQalpha )
		
		# standardized beta diversity
		Qbetastd<-as.numeric(Qbeta/Qgamma )
		
		# list of results
		resQ<-list(Richness_per_plot=nbsploc, Relative_abundance= locabrel, Pi=totabrel, Wc=Wc, Species_abundance_per_plot=abloc, Alpha=Qalpha, Mean_alpha=mQalpha, Gamma=Qgamma, Beta=Qbeta, Standardize_Beta =Qbetastd )
		
		return(resQ)
		
	} 


	###########function disc originally from S. Pavoine####
	 
	disc = function (samples, dis = NULL, structures = NULL, Jost = F){
	    if (!inherits(samples, "data.frame"))
	        stop("Non convenient samples")
	    if (any(samples < 0))
	        stop("Negative value in samples")
	    if (any(apply(samples, 2, sum) < 1e-16))
	        stop("Empty samples")
	    if (!is.null(dis)) {
	        if (!inherits(dis, "dist"))
	            stop("Object of class 'dist' expected for distance")
	       # if (!is.euclid(dis))
	            #stop("Euclidean property is expected for distance")
	        dis <- as.matrix(dis)
	        if (nrow(samples) != nrow(dis))
	            stop("Non convenient samples")
	    }
		if (!is.null(structures)) {
	        if (!inherits(structures, "data.frame"))
	            stop("Non convenient structures")
	        m <- match(apply(structures, 2, function(x) length(x)),
	            ncol(samples), 0)
	        if (length(m[m == 1]) != ncol(structures))
	            stop("Non convenient structures")
	        m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)),
	            function(x) is.factor(structures[, x])), TRUE, 0)
	        if (length(m[m == 1]) != ncol(structures))
	            stop("Non convenient structures")
	    }
	    Structutil <- function(dp2, Np, unit, Jost) {
	        if (!is.null(unit)) {
	            modunit <- model.matrix(~-1 + unit)
	            sumcol <- apply(Np, 2, sum)
	            Ng <- modunit * sumcol
	            lesnoms <- levels(unit)
	        }
	        else {
	            Ng <- as.matrix(Np)
	            lesnoms <- colnames(Np)
	        }
	        sumcol <- apply(Ng, 2, sum)
	        Lg <- t(t(Ng)/sumcol)
	        colnames(Lg) <- lesnoms
	        Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
	        rownames(Pg) <- lesnoms
	        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*%
	            dp2 %*% x))
	        ug <- matrix(1, ncol(Lg), 1)
	        if(Jost) {
	            #dp2 <- as.matrix(as.dist(dfunct01))
	            deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
	            X=t(Lg) %*% dp2 %*% Lg
	            alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
	            Gam = (X + alpha)/2
	            alpha = 1/(1-alpha) #Jost correction
	            Gam = 1/(1-Gam)  #Jost correction
	            Beta_add = Gam - alpha
	            Beta_mult = 100*(Gam - alpha)/Gam
	        }
	        else {
	          deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
	          X=t(Lg) %*% dp2 %*% Lg
	          alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
	          Gam = (X + alpha)/2
	          Beta_add = Gam - alpha
	          Beta_mult = 100*(Gam - alpha)/Gam
	        }
	        colnames(Beta_add) <- lesnoms
	        rownames(Beta_add) <- lesnoms
	        return(list(Beta_add = as.dist(Beta_add), Beta_mult = as.dist(Beta_mult),
	          Gamma=as.dist(Gam), Alpha=as.dist(alpha), Ng = Ng, Pg = Pg))
	    }
	    Diss <- function(dis, nbhaplotypes, samples, structures, Jost) {
	        structutil <- list(0)
	        structutil[[1]] <- Structutil(dp2 = dis, Np = samples, NULL, Jost)
	        diss <- list(structutil[[1]]$Alpha, structutil[[1]]$Gamma, structutil[[1]]$Beta_add, structutil[[1]]$Beta_mult)
	         if (!is.null(structures)) {
	            for (i in 1:length(structures)) {
	                structutil[[i + 1]] <- Structutil(as.matrix(structutil[[1]]$Beta_add), 
	                  structutil[[1]]$Ng, structures[, i], Jost)
	            }
	            diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)), 
	                function(x) as.dist(structutil[[x + 1]]$Beta_add)))
	        }    
	        return(diss)
	    }
	    nbhaplotypes <- sum(samples)
	    diss <- Diss(dis, nbhaplotypes, samples, structures, Jost)
	    if (!is.null(structures)) {
	        names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop", "Beta_region")
	        return(diss)
	    }
	    names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop")
	    return(diss)
	}


	TD<-FD<-PD<-NULL

	#Taxonomic diversity
	dS <- matrix(1, nrow(sample), nrow(sample)) - diag(rep(1, nrow(sample)))
	temp_qdec<- Qdecomp(dS,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
	TD$Richness_per_plot = temp_qdec$Richness_per_plot
	TD$Relative_abundance = temp_qdec$Relative_abundance
	TD$Pi = temp_qdec$Pi
	TD$Wc = temp_qdec$Wc
	if(Jost){
		TD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
		TD$Alpha = 1/(1-temp_qdec$Alpha)
		TD$Gamma = 1/(1-temp_qdec$Gamma)
		TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
		TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
		#Call the disc function for alpha, gamma and beta estimations for each pair of samples
		TD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
		}
	else {
	    TD$Mean_Alpha = temp_qdec$Mean_alpha
		TD$Alpha = temp_qdec$Alpha
		TD$Gamma = temp_qdec$Gamma
		TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
		TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
		#Call the disc function for alpha, gamma and beta estimations for each pair of samples
		TD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
	}
  
	#Functional diversity estimation
	if(!is.null(dfunc)){
		FD<-list()
		if(Jost){
			if(max(dfunc)>1) dfunc <- dfunc/max(dfunc)   #Make sure the distance are between 0 and 1 for the Jost correction
			temp_qdec<- Qdecomp(dfunc,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
			#  FD$Alpha = 1/(1-temp_qdec$Alpha)
			#  FD$Mean_Alpha = mean(FD$Alpha)
			FD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
			FD$Alpha = 1/(1-temp_qdec$Alpha)
			FD$Gamma = 1/(1-temp_qdec$Gamma)
			#FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
			FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			FD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
		}
		else {
			temp_qdec<- Qdecomp(dfunc,t(sample), w=weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
			FD$Mean_Alpha = temp_qdec$Mean_alpha
			FD$Alpha = temp_qdec$Alpha
			FD$Gamma = temp_qdec$Gamma
			FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
			FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
		    #FD$Beta =  temp_qdec$Beta#
		    #Call the disc function for alpha, gamma and beta estimations for each pair of samples
		    FD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
		}
	}
	
	#Phylogenetic diversity estimation
	if(!is.null(dphyl)){
	    PD<-list()
	    if(Jost){
			if(max(dphyl)>1) dphyl <- dphyl/max(dphyl)   #Make sure the distance are between 0 and 1 for the Jost correction
			temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
			PD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
			PD$Alpha = 1/(1-temp_qdec$Alpha)
			PD$Gamma = 1/(1-temp_qdec$Gamma)
			PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
			PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			PD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
		}
	    else {
			temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
			PD$Mean_Alpha = temp_qdec$Mean_alpha
			PD$Alpha = temp_qdec$Alpha
			PD$Gamma = temp_qdec$Gamma
			PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
			PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
			#PD$Beta =  temp_qdec$Beta
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			PD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
	    }
	}
	out <- list(TD, FD, PD)
	names(out) <- c("TD", "FD", "PD")
	return(out)
}

###################################################################################################################################
# 	The Rao function computes alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity with the Rao index          
# 	The script integrates two functions: "Qdecomp", by Villeger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller, and "disc", by S. Pavoine, in the package ade4.
# 	For a regional assemblage of C local communities gamma = mean(alpha) + beta, where:
#  	gamma is the diversity of the regional pool
#  	alpha are the diversities of the local communities
#  	beta is the turn over between local communities
#  	diversity is estimated with the Rao quadratic entropy index (Rao 1982)
#                                      
# INPUTS:                                                                                 
#	- "abundances": matrix of abundances (c x s) of the s species for the c local communities (or samples)           
#	- "dfunct": matrix (s x s) or dist object with pairwise functional trait distances between the s species
#	- "dphyl": as dfunct but for phylogenetic distances
#	- "weight": defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
#	- "Jost": defining if the Jost correction is applied (this paper and Jost 2007)
#	- "structure": a data frame containing the name of the group to which samples belong see                                
#      NA are not allowed in 'locabrel<-abundances/ablocist'
#      NA are automatically replaced by 0 in 'abundances'
#                                                                                         
# OUTPUTS:                                                                           
#	- The results are organized for Taxonomic diversity ($TD), Functional diversity ($FD) and phylogenetical diversity ($PD). Beta and gamma diversities are calculated for the whole data set and for each pair of samples ("Pairwise_samples") 
#	- "$Richness_per_plot"(number of species per sample)
#	- "$Relative_abundance" (species relative abundances per plot)
#	- "$Pi" (species regional relative abundance)
#	- "$Wc" (weigthing factor),                               
#	- "$Mean_Alpha" (mean aplpha diversity; for taxonomic diversity the Simpson index is calculated)                                   
#	- "$Alpha" (alpha diversity for each sample; for taxonomic diversity the Simpson index is calculated)                                       
#	- "$Gamma" (gamma diversity; for taxonomic diversity the Simpson index is calculated)                                      
#	- "$Beta_add" (Gamma-Mean_Alpha)                                          
#	- "$Beta_prop" (Beta_add*100/Gamma)                                                         
#	- "$Pairwise_samples$Alpha" (mean alpha for each pair of samples)                        
#	- "$Pairwise_samples$Gamma" (gamma for each pair of samples)
#	- "$Pairwise_samples$Beta_add" (beta for each pair of samples as Gamma-Mean_Alpha)  
#	- "$Pairwise_samples$Beta_prop" (beta for each pair of samples as Beta_add*100/Gamma)  

#####################################################################################################################################





#Function to plot result of observed indices values against null distribution
#Use the function plot(as.randtest(x)) from ade4
plot_randtest<-function(x, alter=c("greater", "less", "two-sided"), ...){
	
	if(!inherits(x, "listofindex")) {
		if(inherits(x, "Tstats") | inherits(x, "com.index")  | inherits(x, "com.index.multi")) {
			x<-as.listofindex(x)
		}	
		else{stop("x must be a list of objects of class Tstats, com.index or com.index.multi")}
	}

	index.list<-x
	
	oldpar<-par(no.readonly = TRUE)
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
		
	ncom<-c()
	ntr<-c()
	for(i in seq(1, 2*nindex, by=2)){
		ncom<-c(ncom,dim(as.matrix(index.list[[i]]))[1])
		ntr<-c(ntr,dim(as.matrix(index.list[[i]]))[2])
	}
	
	if(is.null(ncom)) {ncom<-dim(as.matrix(index.list[[1]]))[1]}
	if(is.null(ntr)) {ntr<-dim(as.matrix(index.list[[1]]))[2]}
	
	if(is.null(ncom)) {ncom=1}
	if(is.null(ntr)) {ntr=1}
	
	for (i in seq(1,nindex*2,by=2)){
		for (t in 1:ntr[1]){
			rt<-as.randtest(sim=index.list[[i+1]][t,,], obs=index.list[[i]][t], alter=alter)
			plot(rt, main=paste(namesindex.all[i], namestraits[t], "p.value = ", round(rt$pvalue, digits = 5)), ...)
		}
	}

}



#Replace a matrix with abundance of species and mean traits by pop in a pseudo-individual matrix
#Each individual take therefore the value of the population 

ab_to_ind<-function(traits, com, type="count"){
	
	if(nrow(traits) != nrow(com)){
		stop("number of rows of traits and com need to be equal")
	}
	
	#type is either count data or abundance
	#transform abundance data to number of individual
	
	if(type=="abundance"){
		if(min(com)<1){
			com<-apply(com,2, function(x) x/min(x[x>0], na.rm=T))
		}
	}
	
	ntr<-ncol(traits)
	
	#number of individual to individual traits data
	x2<-c()
	x4<-c()
	x6<-c()
	for(i in 1 : nrow(com)){
		x1<-matrix(rep(traits[i,], rowSums(com)[i]), nrow=ntr, ncol=rowSums(com)[i])
		x2<-cbind(x2,x1)	
		
		x3<-rep(rownames(com)[i], rowSums(com)[i])
		x4<-c(x4,x3)
		
		for(co in 1:ncol(com)){
			x5<-rep(colnames(com)[co], com[i,co])
			x6<-c(x5,x6)
		}
	}
	
	res<-list()
	
	res$traits<-data.frame(t(x2))
	res$sp<-as.factor(x4)
	res$ind.plot<-as.factor(x6)
	
	return(res)

}
