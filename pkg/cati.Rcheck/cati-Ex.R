pkgname <- "cati"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "cati-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('cati')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CVNND")
### * CVNND

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Neigbourhood distance metrics
### Title: Coefficient of variation, mean, minimum and standard deviation
###   of the nearest neigbourhood distance.
### Aliases: CVNND MNND MinNND SDNND SDND MND

### ** Examples

data(finch.ind)
CVNND(traits.finch[,1], na.rm = TRUE)
CVNND(traits.finch[,1], div_range =  TRUE, na.rm = TRUE)
CVNND(traits.finch, na.rm = TRUE)
CVNND(traits.finch, scale.tr = FALSE, na.rm = TRUE)
SDND(traits.finch[,1], na.rm = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CVNND", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ComIndex")
### * ComIndex

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ComIndex
### Title: Computing metrics to test and quantify the non-random assembly
###   of communities
### Aliases: ComIndex plot.ComIndex summary.ComIndex print.ComIndex

### ** Examples
	
data(finch.ind)

## Not run: 
##D #Define the functions that will be calculating
##D funct<-c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)",
##D      "max(x, na.rm = TRUE) - min(x, na.rm = TRUE)" )
##D 
##D #Test against the null model regional.ind
##D res.finch.sp_regional.ind<-ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
##D                            nullmodels = "regional.ind", ind.plot = ind.plot.finch,
##D                             nperm = 9, print = FALSE)
##D  
##D #Test against the null model regional.pop
##D #Individuals values are transformed in populational values
##D res.finch.sp_regional.pop<-ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
##D                nullmodels = "regional.pop", ind.plot = ind.plot.finch, 
##D                nperm = 9, print = FALSE)
##D 
##D 
##D #We can calculate index with or without intraspecific variance.
##D 
##D #calculate  of means by population (name_sp_site is a name of a population) 
##D #determine the site for each population (sites_bypop)
##D  
##D name_sp_sites = paste(sp.finch, ind.plot.finch,sep = "_")
##D traits.by.pop<-apply(traits.finch, 2 , 
##D            function (x) tapply(x, name_sp_sites, mean , na.rm = TRUE))
##D 
##D sites_bypop<-lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), split = "_"), 
##D           function(x) x[3])
##D 
##D #New list of function "funct"
##D 
##D funct.1<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))",
##D      "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))",
##D      "tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE)-min(x, na.rm = TRUE))", 
##D      "tapply(x, ind.plot.finch, function(x) CVNND(x, na.rm = TRUE))" )
##D 
##D fact<-unlist(sites_bypop)  
##D funct.2<-c("tapply(x, fact, function(x) mean(x, na.rm = TRUE))",
##D           "tapply(x, fact, function(x) kurtosis(x, na.rm = TRUE))",
##D           "tapply(x, fact, function(x) max(x, na.rm = TRUE)-min(x, na.rm = TRUE))", 
##D           "tapply(x, fact, function(x) CVNND(x, na.rm = TRUE))")
##D 
##D 
##D res.finch.withIV<-ComIndex(traits = traits.finch, index = funct.1,
##D                sp = sp.finch, nullmodels = "regional.ind",
##D                ind.plot = ind.plot.finch, nperm = 9, print = FALSE)
##D 
##D res.finch.withoutIV<-ComIndex(traits = traits.finch, index = funct.2, 
##D                sp = sp.finch, nullmodels = "regional.pop", 
##D                ind.plot = ind.plot.finch, nperm = 9, print = FALSE)
##D 
##D #ComIndex class are associated to S3 methods plot, print and summary.
##D 
##D res.finch.withIV
##D summary(res.finch.withIV)
##D plot(res.finch.withIV)
##D plot(res.finch.withoutIV)
##D 
##D plot(as.listofindex(list(res.finch.withIV, res.finch.withoutIV)))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ComIndex", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ComIndexMulti")
### * ComIndexMulti

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ComIndexMulti
### Title: Computing multitraits metrics to test and quantify the
###   non-random assembly of communities
### Aliases: ComIndexMulti plot.ComIndexMulti summary.ComIndexMulti
###   print.ComIndexMulti

### ** Examples

data(finch.ind)

## Not run: 
##D #For most multivariate functions we need to replace (or exclude) NA values. 
##D #For this example, we use the package mice to complete the data.
##D 
##D comm<-t(table(ind.plot.finch,1:length(ind.plot.finch)))
##D 
##D library(mice)
##D traits = traits.finch
##D mice<-mice(traits.finch)
##D traits.finch.mice<-complete(mice)
##D 
##D 
##D #A simple example to illustrate the concept of the function ComIndexMulti
##D 
##D n_sp_plot<-as.factor(paste(sp.finch, ind.plot.finch, sep = "_")) 
##D res.sum.1<-ComIndexMulti(traits.finch, 
##D               index = c("sum(scale(x), na.rm = T)", "sum(x, na.rm = T)"), 
##D               by.factor = n_sp_plot, nullmodels = "regional.ind", 
##D               ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
##D res.sum.1
##D 
##D 
##D 
##D #A more interesting example using the function hypervolume
##D library(hypervolume)
##D 
##D hv<-hypervolume(traits.finch.mice, 
##D         reps = 100,bandwidth = 0.2, 
##D         verbose = F, warnings = F)
##D plot(hv)
##D 
##D hv.1<-ComIndexMulti(traits.finch.mice, 
##D              index = c("as.numeric(try(hypervolume(na.omit(x), reps = 100, 
##D                  bandwidth = 0.2, verbose = F, warnings = F)@Volume))"),
##D              by.factor = rep(1,length(n_sp_plot)), nullmodels = "regional.ind",
##D              ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch) 
##D 
##D hv.1
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ComIndexMulti", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Fred")
### * Fred

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Fred
### Title: Functional richness, evenness and divergence following Villeger
###   et al. 2008
### Aliases: Fred

### ** Examples

	data(finch.ind)
	## Not run: 
##D 		fred<-Fred(traits.finch.mice, ind.plot.finch)
##D 	
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Fred", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MinMaxMST")
### * MinMaxMST

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MinMaxMST
### Title: Ratio of the shortest distance to the longest distance in a
###   minimum spanning tree
### Aliases: MinMaxMST

### ** Examples

## Not run: 
##D 
##D 	data(finch.ind)
##D 	
##D 	MinMaxMST(traits.finch[1:10,])
##D 	MinMaxMST(traits.finch[1:10,], gower.dist = FALSE)
##D 	MinMaxMST(traits.finch[1:10,], gower.dist = FALSE, scale.tr = FALSE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MinMaxMST", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Pval")
### * Pval

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Pval
### Title: Calcul of p-value for object of class Tstats, ComIndex,
###   ComIndexMulti and listofindex
### Aliases: Pval

### ** Examples


data(finch.ind)
res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)

Pval(res.finch)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Pval", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RandCom")
### * RandCom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RandCom
### Title: Toy model to simulate internal and/or external filtering)
### Aliases: RandCom

### ** Examples
	
	res <- RandCom()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RandCom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RaoRel")
### * RaoRel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RaoRel
### Title: Alpha, gamma and beta-components for taxonomic, functional and
###   phylogenetic diversity
### Aliases: RaoRel

### ** Examples

	data(finch.ind)
	
	## Not run: 
##D 		comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))
##D 		comm.sp <- table(sp.finch, ind.plot.finch)
##D 		class(comm.sp) <- "matrix"
##D 		
##D 		traits.finch.sp <- apply( apply(traits.finch, 2, scale ), 2, 
##D 
##D 		function(x) tapply(x, sp.finch, mean, na.rm = TRUE))
##D 		
##D 		mat.dist <- as.matrix(dist(traits.finch.sp))^2
##D 		
##D 		res.rao <- RaoRel(sample = as.matrix(comm.sp), dfunc = mat.dist, dphyl = NULL, 
##D 		weight = FALSE, Jost = FALSE, structure = NULL)
##D 
##D 		function(x) tapply(x, sp.finch, mean, na.rm=TRUE))
##D 		
##D 		mat.dist <- as.matrix(dist(traits.finch.sp))^2
##D 		
##D 		res.rao <- RaoRel(sample=as.matrix(comm.sp), dfunc=mat.dist, dphyl=NULL, 
##D 		weight=FALSE, Jost=FALSE, structure=NULL)
##D 
##D 		
##D 		witRao <- res.rao$FD$Mean_Alpha  #overall within species variance
##D 		betRao <- res.rao$FD$Beta_add    #between species variance
##D 		totRao <- res.rao$FD$Gamma       #the total variance
##D 		
##D 		witRao+betRao
##D 		totRao
##D 		
##D 		
##D 		#Now let"s take the abundance to calculate Rao diversity.
##D 		
##D 
##D 		res.rao.w <- RaoRel(sample = as.matrix(comm.sp), dfunc = mat.dist, dphyl = NULL, 
##D 		weight = TRUE, Jost = FALSE, structure = NULL)
##D 		
##D 		witRao.w <- res.rao.w$FD$Mean_Alpha  #overall within species variance
##D 		betRao.w <- res.rao.w$FD$Beta_add    #between species variance
##D 		totRao.w <- res.rao.w$FD$Gamma       #the total variance
##D 		
##D 		witRao.w
##D 		betRao.w
##D 		
##D 		
##D 		#Plot the results
##D 		
##D 		barplot(cbind(c(witRao.w, betRao.w), c(witRao, betRao)), 
##D 		names.arg  = c("abundance" ,"presence"), 
##D 		legend.text = c("within species", "between species"), 
##D 		ylab = "Rao", ylim = c(0,10))
##D 		
##D 		
##D 		#We can do this analysis for each trait separately. 
##D 		#First we need to replace (or exclude) NA values. 
##D 		#For this example, we use the package mice to complete the data.
##D 		
##D 		comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))
##D 		
##D 		library(mice)
##D 
##D 		traits = traits.finch
##D 
##D 		mice <- mice(traits.finch)
##D 		traits.finch.mice <- complete(mice)
##D 		
##D 		traits.finch.mice.sp <- apply(apply(traits.finch.mice, 2, scale ), 2, 
##D 		function(x) tapply(x, sp.finch, mean, na.rm = TRUE))
##D 
##D 
##D 		trait.rao.w <- list()
##D 		witRao.w.bytrait <- c()
##D 		betRao.w.bytrait <- c()
##D 
##D 		for (t in 1 : 4){
##D 		  trait.rao.w[[t]] <- RaoRel(sample = as.matrix(comm.sp), 
##D 		  dfunc = dist(traits.finch.mice.sp[,t]), dphyl = NULL, weight = TRUE, 
##D 		  Jost = FALSE, structure = NULL)
##D 				  
##D 		  witRao.w.bytrait <- c(witRao.w.bytrait, trait.rao.w[[t]]$FD$Mean_Alpha)
##D 		  betRao.w.bytrait <- c(betRao.w.bytrait, trait.rao.w[[t]]$FD$Beta_add)
##D 		}
##D 		
##D 		#Plot the results by traits.
##D 		
##D 		barplot(t(cbind( witRao.w.bytrait, betRao.w.bytrait)), 
##D 		names.arg = colnames(traits.finch),
##D 		legend.text = c("within species", "between species"), 
##D 		ylab = "Rao", ylim = c(0,1.5))	
##D 	
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RaoRel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SumBL")
### * SumBL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SumBL
### Title: Sum of branch length of a classification dendrogram (Petchey and
###   Gaston, 2002)
### Aliases: SumBL

### ** Examples


## Not run: 
##D 
##D data(finch.ind)
##D SumBL(traits.finch)
##D SumBL(traits.finch, gower.dist = FALSE)
##D 
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SumBL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Tstats")
### * Tstats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Tstats
### Title: Computing observed T-statistics (T for Traits) and null
###   expectations.
### Aliases: Tstats barplot.Tstats plot.Tstats summary.Tstats print.Tstats
###   sum_Tstats

### ** Examples

	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)
	
	res.finch

	#Tstats class is associated to S3 methods plot, barplot and summary
	
	plot(res.finch)
	
	## Not run: 
##D 	plot(res.finch, type = "simple")
##D 	plot(res.finch, type = "simple_range")
##D 	plot(res.finch, type = "barplot")
##D 	plot(res.finch, type = "bysites")
##D 	plot(res.finch, type = "bytraits")
##D 	
## End(Not run)
	
	attributes(sum_Tstats(res.finch))
	head(sum_Tstats(res.finch)$p.value, 10)
	
	sum_Tstats(res.finch, type = "binary")
	sum_Tstats(res.finch, type = "percent")
	sum_Tstats(res.finch, type = "site")
	sum_Tstats(res.finch, type = "p.value")
	sum_Tstats(res.finch, type = "all")
	
	barplot(res.finch)
	
	attributes(sum_Tstats(res.finch))
	head(sum_Tstats(res.finch)$p.value, 10)
	
	#### An other way to see "ses values" of T-statistics
	
	# Custom theme (from rasterVis package)
	require(rasterVis)
	
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col = my.theme$regions$col)
	
	levelplot(t(ses(res.finch$Tstats$T_IP.IC,res.finch$Tstats$T_IP.IC_nm)$ses), 
	colorkey = my.ckey, par.settings = my.theme,border = "black")
	
	
	## Not run: 
##D 		#### Use a different regional pool than the binding of studied communities
##D 		#create a random regional pool for the example
##D 	
##D 		reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
##D 	
##D 		res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch,res.finch2)))
##D     
##D     
##D 	    #### Use a different regional pool for each communities
##D 		#create a random regional pool for each communities for the example
##D 		list.reg.p <- list(
##D 		traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
##D 		traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
##D 		traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )
##D 
##D 		# Warning: the regional pool need to be larger than the observed communities
##D 		table(ind.plot.finch)
##D 		# For exemple, the third community need a regional pool of more than 981 individuals
##D 		
##D 		res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))	
##D 	
## End(Not run)
	



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Tstats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as.listofindex")
### * as.listofindex

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as.listofindex
### Title: Transform index results in a list of index
### Aliases: as.listofindex

### ** Examples

	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)

	## Not run: 
##D 		#### Use a different regional pool than the binding of studied communities
##D 		#create a random regional pool for the example
##D 	
##D 		reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
##D 	
##D 		res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch,res.finch2)))
##D     
##D     
##D 	    #### Use a different regional pool for each communities
##D 		#create a random regional pool for each communities for the example
##D 		list.reg.p <- list(
##D 		traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
##D 		traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
##D 		traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )
##D 
##D 		# Warning: the regional pool need to be larger than the observed communities
##D 		table(ind.plot.finch)
##D 		# For exemple, the third community need a regional pool of more than 981 individuals
##D 		
##D 		res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))	
##D 	
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as.listofindex", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("decompCTRE")
### * decompCTRE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: decompCTRE
### Title: Variance partitioning for multiple traits
### Aliases: decompCTRE barplot.decompCTRE

### ** Examples

data(finch.ind)

res.decomp <- decompCTRE(traits = traits.finch, sp = sp.finch, 
ind.plot = ind.plot.finch, print = FALSE)

barplot.decompCTRE(res.decomp)

par(mfrow = c(2,2))
barplot.decompCTRE(res.decomp, resume = FALSE)
par(mfrow = c(1,1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("decompCTRE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("finch.ind")
### * finch.ind

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: finch.ind
### Title: Finch morphological data
### Aliases: finch.ind ind.plot.finch sp.finch traits.finch .Random.seed

### ** Examples

	data(finch.ind)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("finch.ind", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("partvar")
### * partvar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: partvar
### Title: Variance partitioning accross nested scales
### Aliases: partvar barPartvar piePartvar

### ** Examples


	data(finch.ind)
	
	cond<-seq(1,length(sp.finch)*2, by = 2)
	genus <- as.vector(unlist(strsplit(as.vector(sp.finch),"_"))[cond])

	res.partvar.finch <- partvar(traits = traits.finch, 
	factors = cbind(sites = as.factor(as.vector(ind.plot.finch)), 
	species = as.factor(as.vector(sp.finch)), genus = as.factor(genus)))
	
	res.partvar.finch
	
	oldpar<-par()
	par(mfrow = c(2,2), mai = c(0.2,0.2,0.2,0.2))
	piePartvar(res.partvar.finch)
	par(oldpar)
	
	barPartvar(res.partvar.finch)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("partvar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot.listofindex")
### * plot.listofindex

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.listofindex
### Title: Plot community assembly index
### Aliases: plot.listofindex

### ** Examples
	
	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)

	## Not run: 
##D 		#### Use a different regional pool than the binding of studied communities
##D 		#create a random regional pool for the example
##D 	
##D 		reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
##D 	
##D 		res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch,res.finch2)))
##D     
##D     
##D 	    #### Use a different regional pool for each communities
##D 		#create a random regional pool for each communities for the example
##D 		list.reg.p <- list(
##D 		traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
##D 		traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
##D 		traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )
##D 
##D 		# Warning: the regional pool need to be larger than the observed communities
##D 		table(ind.plot.finch)
##D 		# For exemple, the third community need a regional pool of more than 981 individuals
##D 		
##D 		res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 	    sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)	
##D 	    
##D 	    plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))	
##D 	
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.listofindex", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotCorTstats")
### * plotCorTstats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotCorTstats
### Title: Plot the bivariate relationships between T-statistics
### Aliases: plotCorTstats

### ** Examples

	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)
	
	plotCorTstats(res.finch, bysite = FALSE)
	plotCorTstats(res.finch, bysite = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotCorTstats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotDistri")
### * plotDistri

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotDistri
### Title: Plot function to represent density of trait values
### Aliases: plotDistri

### ** Examples
	
data(finch.ind)

## Not run: 
##D 	#Plot the distribution of trait values for populations, 
##D 	#species, sites and regional scales. 
##D 	
##D 	### First, let try the distribution for all populations 
##D 	#of Darwin finches.
##D 	
##D 	par(mfrow = c(4,4), cex = 0.5)
##D 	plotDistri(traits.finch, sp.finch, ind.plot.finch, ylim.cex = 3, 
##D 	plot.ask = FALSE, multipanel = FALSE, leg = FALSE)
##D 	
##D 	### Then we can inverse the second and the third arguments 
##D 	#to plot the distribution for all finches species. 
##D 	
##D 	par(mfrow = c(4,4), cex = 0.5)
##D 	plotDistri(traits.finch, ind.plot.finch, sp.finch, ylim.cex = 8, 
##D 	plot.ask = FALSE, multipanel = FALSE, leg = FALSE)
##D 	
##D 	### Only one trait to plot using leg = TRUE to plot the legend
##D 
##D 	par(mfrow=c(2,3))
##D 	plotDistri(as.matrix(traits.finch[,1]), ind.plot.finch, sp.finch, 
##D     ylim.cex=8, plot.ask = FALSE, multipanel = FALSE, leg = TRUE, cex.leg=0.5)
##D 
##D 	### You can also plot trait distribution for all species in the region
##D 	
##D 	par(mfrow = c(1,1), cex = 1)
##D 	plotDistri(traits.finch, rep("region", times = dim(traits.finch)[1]), 
##D 	sp.finch, ylim.cex = 6, plot.ask = FALSE, leg = FALSE)
##D 	
## End(Not run)
	
### You can also plot trait distribution for all sites
#without taking into account species identity

plotDistri(traits.finch, rep("toutes_sp", times = dim(traits.finch)[1]), 
ind.plot.finch, ylim.cex = 3, plot.ask = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotDistri", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotRandtest")
### * plotRandtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotRandtest
### Title: Plot result of observed indices values against null distribution
### Aliases: plotRandtest

### ** Examples

	data(finch.ind)
	## Not run: 
##D 		res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
##D 		sp = sp.finch, nperm = 99, print = FALSE)
##D 	
##D 	par(mfrow = c(4,4))
##D 	
##D 	plotRandtest(res.finch)
##D 	plotRandtest(res.finch, alter = "less")
##D 	
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotRandtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotSESvar")
### * plotSESvar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSESvar
### Title: Plot SES values against a variable
### Aliases: plotSESvar

### ** Examples

	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, sp = sp.finch, 
	nperm = 9)

	par(mfrow = c(2,2))
	species.richness <- table(ind.plot.finch)
	plotSESvar(as.listofindex(list(res.finch)), species.richness, 
	multipanel = FALSE)

	#Same plot with resume = TRUE.
	
	par(mfrow = c(2,2))
	plotSESvar(as.listofindex(list(res.finch)), species.richness, 
	resume = TRUE, multipanel = FALSE)
	par(mfrow = c(1,1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSESvar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotSpPop")
### * plotSpPop

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSpPop
### Title: Plot populations values against species values
### Aliases: plotSpPop

### ** Examples

	data(finch.ind)
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE)

	#If we change the value of the threshold 
	#(alpha = 10% instead of 5% 
	#and the minimum individual to represent significativity 
	#fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpPop(traits.finch, ind.plot.finch, sp.finch, p.val = 0.1,  
	min.ind.signif = 3, silent = TRUE)


	#For a more simple figure, add the option resume = TRUE. 
	#Again if we change the value of the threshold 
	#(alpha = 10% instead of 5% 
	#and the minimum individual to represent significativity
	# fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey")
	
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey", col.sp = "black")
	
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey", col.sp = "black", 
	p.val = 0.1,  min.ind.signif = 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSpPop", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotSpVar")
### * plotSpVar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSpVar
### Title: Plot populations values against species values
### Aliases: plotSpVar

### ** Examples

	data(finch.ind)
	
	#Random variable for this example
	variable <- c(1,5,15,6,3,25)
	
	plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
	silent = TRUE)

	#If we change the value of the threshold 
	#(alpha = 10% instead of 5% 
	#and the minimum individual to represent significativity 
	#fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
	p.val = 0.1,  min.ind.signif = 3, silent = TRUE)


	#For a more simple figure, add the option resume = TRUE. 
	#Again if we change the value of the threshold 
	#(alpha = 10% instead of 5% 
	#and the minimum individual to represent significativity
	# fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
	silent = TRUE, resume = TRUE, col.pop = "grey")
	
	plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
	silent = TRUE, resume = TRUE, col.pop = "grey", col.sp = "black")
	
	plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
	silent = TRUE, resume = TRUE, col.pop = "grey", col.sp = "black", 
	p.val = 0.1,  min.ind.signif = 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSpVar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ses")
### * ses

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ses
### Title: Standardized effect size and confidence interval for a matrix of
###   statistics
### Aliases: ses

### ** Examples

	data(finch.ind)
	
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)

	ses(res.finch$Tstats$T_IP.IC, res.finch$Tstats$T_IP.IC_nm)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ses", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
