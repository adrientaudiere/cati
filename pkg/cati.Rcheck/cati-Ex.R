pkgname <- "cati"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('cati')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CVNND")
### * CVNND

flush(stderr()); flush(stdout())

### Name: CVNND
### Title: Coefficient of variation of the nearest neigbourhood distance
###   (CVNND)
### Aliases: CVNND

### ** Examples

	data(finch.ind)
	CVNND(traits.finch[,1])



cleanEx()
nameEx("ComIndex")
### * ComIndex

flush(stderr()); flush(stdout())

### Name: ComIndex
### Title: Computing the moments of the trait distribution and other
###   metrics to test and quantify the non-random assembly of communities
### Aliases: ComIndex

### ** Examples
	
	data(finch.ind)
	oldpar <- par()
	
	####
	#The function ComIndex allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	require(e1071)
	
	funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
	"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )
	
	## Not run: 
##D 		res.finch.sp_mn2 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.sp_mn3 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 	
##D 		####
##D 		#We can represent Standardized Effect Size (ses)
##D 		#using the function plot(as.listofindex(list1, list2, list3))
##D 		
##D 		list.ind2 <- list(res.finch.sp_mn2, res.finch.sp_mn3)
##D 		index.list2 <- as.listofindex(list.ind2)
##D 		
##D 		plot(index.list2, type = "bytraits")
##D 		
##D 		plot(index.list2)
##D 	
##D 		####
##D 		#This allows to calcul index per site 
##D 		#for example using "tapply(x, sites, mean)".
##D 		
##D 		funct <- c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))", 
##D 		"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))", 
##D 		"tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE) - 
##D 		min(x, na.rm = TRUE) )", "tapply(x, ind.plot.finch, function(x) 
##D 		CVNND(x))"  )
##D 		
##D 		
##D 		##Null model 1 is trivial for this function 
##D 		#because randomisation is within community only
##D 		
##D 		res.finch.ind_mn1 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c(1,1,1,1), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.ind_mn2 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 	
##D 		
##D 		####
##D 		#We can calcul metrics with or without intraspecific variance.
##D 		#Calculation of trait averages per population 
##D 		#(name_sp_site is a name of a population) 
##D 		#like in the function ComIndex
##D 		#and determine the site for each population (sites_bypop)
##D 	
##D 		name_sp_sites = paste(sp.finch, ind.plot.finch, sep = "_")
##D 		traits.by.pop <- apply(traits.finch, 2 , function (x) 
##D 		tapply(x, name_sp_sites, mean , na.rm = TRUE))
##D 		
##D 		sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), 
##D 		split = "_"), function(x) x[3])
##D 		
##D 		funct.withoutIV <- c("tapply(x, unlist(sites_bypop), function(x) 
##D 		mean(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), function(x) 
##D 		kurtosis(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), function(x) 
##D 		max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
##D 		
##D 		
##D 		funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
##D 		mean(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		kurtosis(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
##D 		
##D 		
##D 		res.finch.withIV <- ComIndex(traits = traits.finch, index = funct.withIV, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.withoutIV <- ComIndex(traits = traits.finch, index = funct.withoutIV, 
##D 		sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		
##D 		####
##D 		#We can also represent T-statistics and custom index thanks to
##D 		#the plot.listofindex function.
##D 		
##D 		res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, sp = sp.finch, 
##D 		nperm = 9, print = FALSE)
##D 	
##D 		list.ind <- list(res.finch.withIV, res.finch.withoutIV ,res.finch)
##D 		
##D 		index.list1 <- as.listofindex(list.ind, namesindex = c("mean", "kurtosis", 
##D 		"range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
##D 		"T_IP.IC", "T_IC.IR", "T_PC.PR"))
##D 		
##D 		class(index.list1)
##D 		
##D 		par(mfrow = c(2,3))
##D 		plot(index.list1,type = "bytraits", bysite = TRUE)
##D 		
##D 		par(mfrow = c(2,2))
##D 		plot(index.list1,type = "bytraits")
##D 		par(mfrow = c(1,1))
##D 		
##D 		plot(index.list1,type = "simple")
##D 		plot(index.list1,type = "simple_range")
##D 		plot(index.list1,type = "normal")
##D 		plot(index.list1,type = "barplot")
##D 	
## End(Not run)
	
	############################
	####Using ind.value = FALSE and community data matrix if there is no data 
	#available at the individual level.
	## Not run: 
##D 	
##D 		#create traits data at the species level
##D 		traits_by_sp <- apply(traits.finch,2,function(x) tapply(x,sp.finch,
##D 		function(x) mean(x, na.rm = T)))  
##D 		
##D 		#create traits data at the populational level
##D 		names_sp_ind_plot <- as.factor(paste(sp.finch, ind.plot.finch, sep = "@")) 
##D 		traits_by_pop <- apply(traits.finch,2,function(x) tapply(x,names_sp_ind_plot, 
##D 		function(x) mean(x, na.rm = T) ))  
##D 		
##D 		#create community data matrix at the species or populational level
##D 		w1 <- table(sp.finch,ind.plot.finch)
##D 		dim(w1)
##D 		dim(traits_by_sp)
##D 		
##D 		w2 <- table(names_sp_ind_plot,ind.plot.finch)
##D 		dim(w2)
##D 		dim(traits_by_pop)
##D 		
##D 		#Choose indices
##D 		require(e1071)
##D 		
##D 		funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
##D 		"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )
##D 		
##D 		
##D 		#################
##D 		#with species value
##D 		
##D 		res <- AbToInd(traits_by_sp, w1)
##D 		
##D 		ComIndex(traits_by_sp, nullmodels = 2, ind.value = FALSE, index = funct, 
##D 		sp = rownames(traits_by_sp), com = w1, nperm = 9)
##D 		
##D 		
##D 		#################
##D 		#with population value
##D 		res <- AbToInd(traits_by_pop, w2)
##D 		sp.sp <- unlist(strsplit(rownames(traits_by_pop),"@"))[seq(1,39*2,2)]
##D 		
##D 		ComIndex(traits_by_pop, nullmodels = 2, ind.value = FALSE, index = funct, 
##D 		sp = sp.sp, com = w2)
##D 		
##D 	
## End(Not run)
	
	############################
	####Simple example using null model 2sp.prab (species level without taking 
	# into acount for species abundance, prab for presence/absence)
	
	## Not run: 
##D 	traits_by_sp <- apply(traits.finch, 2, function(x) tapply(x, name_sp_sites, mean, na.rm=T))
##D 	sites_bysp<-unlist(strsplit(rownames(traits_by_sp), split="_"))[seq(3,3*dim(traits_by_sp)[1], by=3) ]
##D 	
##D 	funct.withoutIV.prab <- c("tapply(x, unlist(sites_bysp), function(x) mean(x, na.rm = TRUE))", 
##D 		"tapply(x, unlist(sites_bysp), function(x) kurtosis(x, na.rm = TRUE))", 
##D 		"tapply(x, unlist(sites_bysp), function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, unlist(sites_bysp), function(x) CVNND(x))")
##D 		
##D 	res.finch.withoutIV.prab <- ComIndex(traits = traits.finch, index = funct.withoutIV.prab, 
##D 		sp = sp.finch, nullmodels = rep("2sp.prab", times=4), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 	
##D 	list.ind2 <- list(res.finch.withoutIV, res.finch.withoutIV.prab)
##D 	index.list2 <- as.listofindex(list.ind2, namesindex = 
##D 		c("mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
##D 		"mean.prab", "kurtosis.prab", "range.prab", "CVNND.prab"))
##D 	
##D 	plot(index.list2)
##D 	
##D 	
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ComIndexMulti")
### * ComIndexMulti

flush(stderr()); flush(stdout())

### Name: ComIndexMulti
### Title: Computing multitraits metrics to test and quantify the
###   non-random assembly of communities
### Aliases: ComIndexMulti

### ** Examples

	data(finch.ind)
	
	####
	#For most multivariate functions we need to replace (or exclude)
	#NA values.
	
	#For this example, we use the package mice to complete the data.
	
	## Not run: 
##D 		names.sp_ind_plot <- as.factor(paste(sp.finch, ind.plot.finch, sep = "_")) 
##D 		
##D 		comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))
##D 		
##D 		library(mice)
##D 		traits = traits.finch
##D 		mice <- mice(traits.finch)
##D 		traits.finch.mice <- complete(mice)
##D 		
##D 		####
##D 		#A simple example to illustrate the concept of the function 
##D 		#ComIndexMulti
##D 		
##D 		res.sum.1 <- ComIndexMulti(traits.finch, 
##D 		index = c("sum(scale(x), na.rm = TRUE)", "sum(x, na.rm = TRUE)"), 
##D 		by.factor = names.sp_ind_plot, nullmodels = c(2,2), 
##D 		ind.plot = ind.plot.finch, nperm = 50, sp = sp.finch)
##D 		
##D 		attributes(ses.listofindex(as.listofindex(res.sum.1)))
##D 		
##D 		####
##D 		#A more interesting example using the function hypervolume 
##D 		#from the package hypervolume. 
##D 		#We show here several results which differe in there factor 
##D 		#that delimit the group to calculate different hypervolume 
##D 		#(argument by_factor). 
##D 		
##D 		require(hypervolume)
##D 		
##D 		res.hv.1 <- ComIndexMulti(traits.finch.mice, index = 
##D 		c("as.numeric (try(hypervolume(na.omit(x), warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
##D 		by.factor = rep(1,length(names.sp_ind_plot)), nullmodels = c(2,2), 
##D 		ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
##D 		
##D 		res.hv.2 <- ComIndexMulti(traits.finch.mice, index = 
##D 		c("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
##D 		by.factor = names.sp_ind_plot, nullmodels = c(2,2), 
##D 		ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
##D 		
##D 		res.hv.3 <- ComIndexMulti(traits.finch.mice, index = 
##D 		c("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
##D 		by.factor = ind.plot.finch,	nullmodels = c(2,2), 
##D 		ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
##D 		
##D 		res.hv.4 <- ComIndexMulti(traits.finch.mice, index = 
##D 		c("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
##D 		by.factor = sp.finch, nullmodels = c(2,2), ind.plot = ind.plot.finch, 
##D 		nperm = 9, sp = sp.finch)
##D 		
##D 		
##D 		list.ind.multi <- as.listofindex(list(res.hv.2, res.hv.3, res.hv.4))
##D 		
##D 		ses.listofindex(list.ind.multi)
##D 		
##D 		plot(list.ind.multi)
##D 		plot(list.ind.multi, xlim = c(-200,20))
##D 		
##D 		
##D 
##D 		
##D 	
## End(Not run)



cleanEx()
nameEx("RaoRel")
### * RaoRel

flush(stderr()); flush(stdout())

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
##D 		res.rao.w <- RaoRel(sample=as.matrix(comm.sp), dfunc=mat.dist, dphyl=NULL, 
##D 		weight=TRUE, Jost=FALSE, structure=NULL)
##D 
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
##D 		traits=traits.finch
##D 
##D 		mice <- mice(traits.finch)
##D 		traits.finch.mice <- complete(mice)
##D 		
##D 		
##D 		traits.finch.mice.sp <- apply(apply(traits.finch.mice, 2, scale ), 2, 
##D 
##D 		function(x) tapply(x, sp.finch, mean, na.rm = TRUE))
##D 
##D 		function(x) tapply(x, sp.finch, mean, na.rm=TRUE))
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
##D 		}
##D 		
##D 		for(t in 1 : 4){
##D 		  trait.rao.w[[t]] <- RaoRel(sample=as.matrix(comm.sp), 
##D 		  dfunc=dist(traits.finch.mice.sp[,t]), dphyl=NULL, weight=TRUE, 
##D 		  Jost=FALSE, structure=NULL)
##D 
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



cleanEx()
nameEx("Tstats")
### * Tstats

flush(stderr()); flush(stdout())

### Name: Tstats
### Title: Computing observed T-statistics (T for Traits) and null
###   expectations.
### Aliases: Tstats barplot.Tstats plot.Tstats summary_Tstats

### ** Examples

	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)
	
	attributes(res.finch)

	#Tstats class is associated to S3 methods plot, barplot and summary
	
	plot(res.finch)
	
	plot(res.finch, type = "color_cond")
	plot(res.finch, type = "simple")
	plot(res.finch, type = "simple_sd")
	plot(res.finch, type = "barplot")
	
	attributes(summary_Tstats(res.finch))
	head(summary_Tstats(res.finch)$p.value, 10)
	
	summary_Tstats(res.finch, type = "binary")
	summary_Tstats(res.finch, type = "percent")
	summary_Tstats(res.finch, type = "site")
	summary_Tstats(res.finch, type = "p.value")
	summary_Tstats(res.finch, type = "all")
	
	barplot(res.finch)
	
	attributes(summary_Tstats(res.finch))
	head(summary_Tstats(res.finch)$p.value, 10)
	
	#### An other way to see "ses values" of T-statistics
	
	# Custom theme (from rasterVis package)
	require(rasterVis)
	
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col = my.theme$regions$col)
	
	levelplot(t(ses(res.finch$T_IP.IC,res.finch$T_IP.IC_nm)$ses), 
	colorkey = my.ckey, par.settings = my.theme,border = "black")
	
	
	#### Use a different regional pool than the binding of studied communities
	
	#create a random regional pool for the example

	reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])

	res.finch2 <- Tstats(traits.finch, ind_plot = ind.plot.finch, 
    sp = sp.finch, nperm = 9, print = FALSE)	

	reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300),])

	res.finch2 <- Tstats(traits.finch, ind_plot=ind.plot.finch, 
    sp=sp.finch, nperm=9, print=FALSE)	

	



cleanEx()
nameEx("as.listofindex")
### * as.listofindex

flush(stderr()); flush(stdout())

### Name: as.listofindex
### Title: Transform index results in a list of index
### Aliases: as.listofindex

### ** Examples

	data(finch.ind)
	oldpar <- par()
	
	####
	#The function ComIndex allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	require(e1071)
	

	funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)",
	"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)" )
	
	res.finch.sp_mn2 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
	nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, nperm = 9, print = FALSE)
	
	res.finch.sp_mn3 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
	nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

	funct <- c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)" )
	
	res.finch.sp_mn2 <- ComIndex(traits=traits.finch, index=funct, sp=sp.finch,
	nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	res.finch.sp_mn3 <- ComIndex(traits=traits.finch, index=funct, sp=sp.finch,
	nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot=ind.plot.finch, nperm=9, print=FALSE)


	####
	#We can represent Standardized Effect Size (ses) using the 
	#function plot(as.listofindex(list1, list2, list3))

	#The function ComIndex allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	

	funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
	"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )
	
	funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
	"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )
	
	res.finch.sp_mn2 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, 
	nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, nperm = 9, print = FALSE)
	
	res.finch.sp_mn3 <- ComIndex(traits = traits.finch, index = funct, sp = sp.finch, 
	nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

	funct <- c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
	
	funct <- c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
	
	res.finch.sp_mn2 <- ComIndex(traits=traits.finch, index=funct, sp=sp.finch, 
	nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	res.finch.sp_mn3 <- ComIndex(traits=traits.finch, index=funct, sp=sp.finch, 
	nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot=ind.plot.finch, nperm=9, print=FALSE)


	####
	#We can represent Standardized Effect Size (ses) 
	#using the function plot(as.listofindex(list1, list2, list3))
	
	list.ind2 <- list(res.finch.sp_mn2, res.finch.sp_mn3)
	index.list2 <- as.listofindex(list.ind2)
	
	plot(index.list2, type = "bytraits")
	
	plot(index.list2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("decompWithin")
### * decompWithin

flush(stderr()); flush(stdout())

### Name: decompWithin
### Title: Variance partitioning for multiple traits
### Aliases: decompWithin barplot.decompWithin

### ** Examples

	data(finch.ind)
	
	res.decomp <- decompWithin(traits = traits.finch, sp = sp.finch, 
	ind.plot = ind.plot.finch, print = FALSE)
	
	barplot.decompWithin(res.decomp)
	
	par(mfrow = c(2,2))
	barplot.decompWithin(res.decomp, resume = FALSE)
	par(mfrow = c(1,1))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("finch.ind")
### * finch.ind

flush(stderr()); flush(stdout())

### Name: finch.ind
### Title: Finch morphological data
### Aliases: finch.ind ind.plot.finch sp.finch traits.finch .Random.seed

### ** Examples

	data(finch.ind)



cleanEx()
nameEx("partvar")
### * partvar

flush(stderr()); flush(stdout())

### Name: partvar
### Title: Variance partitioning accross nested scales
### Aliases: partvar barPartvar piePartvar

### ** Examples

	data(finch.ind)

	genus <- as.vector(unlist(strsplit(as.vector(sp.finch),"_"))[seq(1,length(sp.finch)*2, by = 2)])

	res.partvar.finch <- partvar(traits = traits.finch, 
	factors = cbind(sites = as.factor(as.vector(ind.plot.finch)), 
	species = as.factor(as.vector(sp.finch)), genus = as.factor(genus)))
	
	res.partvar.finch <- partvar(traits=traits.finch, 
	factors=cbind(sites=as.factor(as.vector(ind.plot.finch)), 
	species=as.factor(as.vector(sp.finch)), genus=as.factor(genus)))


	res.partvar.finch
	
	oldpar <- par()

	par(mfrow = c(2,2), mai = c(0.2,0.2,0.2,0.2))
	piePartvar(res.partvar.finch, col = c("red", "green", "blue", "purple"))
	par(oldpar)

	barPartvar(res.partvar.finch, col = c("red", "green", "blue", "purple"))

	par(mfrow=c(2,2), mai=c(0.2,0.2,0.2,0.2))
	piePartvar(res.partvar.finch, col=c("red", "green", "blue", "purple"))
	par(oldpar)

	barPartvar(res.partvar.finch, col=c("red", "green", "blue", "purple"))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot.listofindex")
### * plot.listofindex

flush(stderr()); flush(stdout())

### Name: plot.listofindex
### Title: Plot community assembly index
### Aliases: plot.listofindex

### ** Examples
	
	data(finch.ind)
	oldpar <- par()
	
	####
	#The function ComIndex allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	require(e1071)
	

	funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
	"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )
	
	## Not run: 
##D 		res.finch.sp_mn2 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.sp_mn3 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 
##D 	funct <- c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
##D 	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
##D 	
## End(Not run)
	
	## Not run: 
##D 		res.finch.sp_mn2 <- ComIndex(traits=traits.finch, index=funct, 
##D 		sp=sp.finch, nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 		
##D 		res.finch.sp_mn3 <- ComIndex(traits=traits.finch, index=funct, 
##D 		sp=sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 
##D 
##D 		####
##D 		#We can represent Standardized Effect Size (ses)
##D 		#using the function plot(as.listofindex(list1, list2, list3))
##D 		
##D 		list.ind2 <- list(res.finch.sp_mn2, res.finch.sp_mn3)
##D 		index.list2 <- as.listofindex(list.ind2)
##D 		
##D 		plot(index.list2, type = "bytraits")
##D 		
##D 		plot(index.list2)
##D 	
## End(Not run)
	
	####
	#This allows to calcul index per site 
	#for example using "tapply(x, sites, mean)".
	

	funct <- c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))", 
	"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))", 
	"tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE) - 
	min(x, na.rm = TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))")

	funct <- c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", 
	"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", 
	"tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - 
	min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))")

	
	
	##Null model 1 is trivial for this function 
	#because randomisation is within community only
	
	## Not run: 
##D 
##D 		res.finch.ind_mn1 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c(1,1,1,1), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.ind_mn2 <- ComIndex(traits = traits.finch, index = funct, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 
##D 		res.finch.ind_mn1 <- ComIndex(traits=traits.finch, index=funct, 
##D 		sp=sp.finch, nullmodels = c(1,1,1,1), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 		
##D 		res.finch.ind_mn2 <- ComIndex(traits=traits.finch, index=funct, 
##D 		sp=sp.finch, nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 
##D 		
##D 		####
##D 		#We can calcul metrics with or without intraspecific variance.
##D 		#Calculation of trait averages per population 
##D 		#(name_sp_site is a name of a population) 
##D 		#like in the function ComIndex
##D 		#and determine the site for each population (sites_bypop)
##D 		
##D 		name_sp_sites = paste(sp.finch, ind.plot.finch, sep = "_")
##D 		
##D 		traits.by.pop <- apply(traits.finch, 2 , function (x) 
##D 
##D 		tapply(x, name_sp_sites, mean , na.rm = TRUE))
##D 		
##D 		sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), 
##D 		split = "_"), function(x) x[3])
##D 		
##D 		sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), 
##D 		split = "_"), function(x) x[3])
##D 		
##D 		funct.withoutIV <- c("tapply(x, unlist(sites_bypop), 
##D 		function(x) mean(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) kurtosis(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
##D 		
##D 		funct.withoutIV <- c("tapply(x, unlist(sites_bypop), 
##D 		function(x) mean(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) kurtosis(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
##D 		
##D 		funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
##D 		mean(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		kurtosis(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
##D 		
##D 		funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
##D 		mean(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		kurtosis(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
##D 		"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
##D 		
##D 		res.finch.withIV <- ComIndex(traits = traits.finch, index = funct.withIV, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.withIV <- ComIndex(traits = traits.finch, index = funct.withIV, 
##D 		sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		res.finch.withoutIV <- ComIndex(traits = traits.finch, index = funct.withoutIV, 
##D 		sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot = ind.plot.finch, 
##D 		nperm = 9, print = FALSE)
##D 
##D 		tapply(x, name_sp_sites, mean , na.rm=TRUE))
##D 		
##D 		sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep="_"), 
##D 		split="_"), function(x) x[3])
##D 		
##D 		sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep="_"), 
##D 		split="_"), function(x) x[3])
##D 		
##D 		funct.withoutIV <- c("tapply(x, unlist(sites_bypop), 
##D 		function(x) mean(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
##D 		"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
##D 		
##D 		funct.withoutIV <- c("tapply(x, unlist(sites_bypop), 
##D 		function(x) mean(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
##D 		function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
##D 		"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
##D 		
##D 		funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
##D 		mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
##D 		"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
##D 		
##D 		funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
##D 		mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
##D 		max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
##D 		"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
##D 		
##D 		res.finch.withIV <- ComIndex(traits=traits.finch, index=funct.withIV, 
##D 		sp=sp.finch, nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 		
##D 		res.finch.withIV <- ComIndex(traits=traits.finch, index=funct.withIV, 
##D 		sp=sp.finch, nullmodels = c("2","2","2","2"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 		
##D 		res.finch.withoutIV <- ComIndex(traits=traits.finch, index=funct.withoutIV, 
##D 		sp=sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), ind.plot=ind.plot.finch, 
##D 		nperm=9, print=FALSE)
##D 
##D 	
##D 	
## End(Not run)
	
	####
	#We can also represent T-statistics and custom index thanks to
	#the plot.listofindex function.
	
	## Not run: 
##D 
##D 		res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, sp = sp.finch, 
##D 		nperm = 9, print = FALSE)
##D 		
##D 		list.ind <- list(res.finch.withIV, res.finch.withoutIV ,res.finch)
##D 		
##D 		index.list1 <- as.listofindex(list.ind, namesindex = c("mean", "kurtosis", 
##D 
##D 		res.finch <- Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, 
##D 		nperm=9, print=FALSE)
##D 		
##D 		list.ind <- list(res.finch.withIV, res.finch.withoutIV ,res.finch)
##D 		
##D 		index.list1 <- as.listofindex(list.ind, namesindex=c("mean", "kurtosis", 
##D 
##D 		"range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
##D 		"T_IP.IC", "T_IC.IR", "T_PC.PR"))
##D 		
##D 		class(index.list1)
##D 		
##D 		par(mfrow = c(2,3))
##D 		plot(index.list1,type = "bytraits", bysite = TRUE)
##D 		
##D 		par(mfrow = c(2,2))
##D 		plot(index.list1,type = "bytraits")
##D 		par(mfrow = c(1,1))
##D 		
##D 		plot(index.list1,type = "simple")
##D 		plot(index.list1,type = "simple_range")
##D 		plot(index.list1,type = "normal")
##D 		plot(index.list1,type = "barplot")
##D 	
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotCorTstats")
### * plotCorTstats

flush(stderr()); flush(stdout())

### Name: plotCorTstats
### Title: Plot the bivariate relationships between T-statistics
### Aliases: plotCorTstats

### ** Examples

	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)
	
	plotCorTstats(res.finch, bysite = FALSE)
	plotCorTstats(res.finch, bysite = TRUE)



cleanEx()
nameEx("plotDistri")
### * plotDistri

flush(stderr()); flush(stdout())

### Name: plotDistri
### Title: Plot function to represent density of trait values
### Aliases: plotDistri

### ** Examples
	
	data(finch.ind)
	
	#Plot the distribution of trait values for populations, 
	#species, sites and regional scales. 
	
	#First, let try the distribution for all populations 
	#of Darwin finches.
	
	par(mfrow = c(4,4), cex = 0.5)
	
	plotDistri(traits.finch, sp.finch, ind.plot.finch, ylim.cex = 3, 
	plot.ask = FALSE, multipanel = FALSE, leg = FALSE)
	
	par(mfrow = c(1,1), cex = 1)
	
	
	#Then we can inverse the second and the third arguments 
	#to plot the distribution for all finches species. 
	
	par(mfrow = c(4,4), cex = 0.5)
	
	plotDistri(traits.finch, ind.plot.finch, sp.finch, ylim.cex = 8, 
	plot.ask = FALSE, multipanel = FALSE, leg = FALSE)
	
	par(mfrow = c(1,1), cex = 1)
	
	
	#You can also plot trait distribution for all species in the region
	
	plotDistri(traits.finch, rep("region", times = dim(traits.finch)[1]), 
	sp.finch, ylim.cex = 6, plot.ask = FALSE, leg = FALSE)
	
	
	#You can also plot trait distribution for all sites
	#without taking into account species identity
	
	plotDistri(traits.finch, rep("toutes_sp", times = dim(traits.finch)[1]), 
	ind.plot.finch, ylim.cex = 3, plot.ask = FALSE)
	
	par(mfrow = c(4,4), cex = 0.5)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotRandtest")
### * plotRandtest

flush(stderr()); flush(stdout())

### Name: plotRandtest
### Title: Plot result of observed indices values against null distribution
### Aliases: plotRandtest

### ** Examples

	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)
	
	par(mfrow = c(3,4))
	
	plotRandtest(res.finch)
	plotRandtest(res.finch, alter = "two-sided")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotSESvar")
### * plotSESvar

flush(stderr()); flush(stdout())

### Name: plotSESvar
### Title: Plot SES values against a variable
### Aliases: plotSESvar

### ** Examples

	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, sp = sp.finch, 
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



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotSpPop")
### * plotSpPop

flush(stderr()); flush(stdout())

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



cleanEx()
nameEx("plotSpVar")
### * plotSpVar

flush(stderr()); flush(stdout())

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



cleanEx()
nameEx("ses")
### * ses

flush(stderr()); flush(stdout())

### Name: ses
### Title: Standardized effect size and confidence interval for a matrix of
###   statistics
### Aliases: ses

### ** Examples

	data(finch.ind)
	

	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)

	res.finch <- Tstats(traits.finch, ind_plot=ind.plot.finch, 
	sp=sp.finch, nperm=9)

	
	ses(res.finch$T_IP.IC, res.finch$T_IP.IC_nm)
	
	ses(t(res.finch$T_IP.IC), res.finch$T_IP.IC_nm)



cleanEx()
nameEx("ses.listofindex")
### * ses.listofindex

flush(stderr()); flush(stdout())

### Name: ses.listofindex
### Title: Standardized effect size for a list of index.
### Aliases: ses.listofindex

### ** Examples

	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind_plot = ind.plot.finch, sp = sp.finch, 
	nperm = 9, print = FALSE)
	
	#Calcul of means by population (name_sp_site is a name of a population) 
	#like in the function ComIndex and determine the site 
	#for each population (sites_bypop)

	name_sp_sites = paste(sp.finch, ind.plot.finch,sep = "_")
	traits.by.pop <- apply(traits.finch, 2 , function (x) 
	tapply(x, name_sp_sites, mean , na.rm = TRUE))
	
	require(e1071)
	
	sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_")
	, split = "_"), function(x) x[3])
	
	
	funct.withoutIV <- c("tapply(x, unlist(sites_bypop), function(x) mean(x, na.rm=TRUE))",
	"tapply(x, unlist(sites_bypop), function(x) kurtosis(x, na.rm=TRUE))",
	"tapply(x, unlist(sites_bypop), function(x)	max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
	"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
	
	
	funct.withIV <- c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))",
	"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))",
	"tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
	"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	

	res.finch.withIV <- ComIndex(traits = traits.finch, index = funct.withIV,
	sp = sp.finch, nullmodels = rep("2", times=4), ind.plot = ind.plot.finch, nperm = 9
	, print = FALSE)
	
	res.finch.withoutIV <- ComIndex(traits = traits.finch, 
	index = funct.withoutIV, sp = sp.finch, nullmodels = rep("2sp", times=4), 
	ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

	res.finch.withIV <- ComIndex(traits=traits.finch, index=funct.withIV, 
	sp=sp.finch, nullmodels = rep("2", times=4), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)
	
	res.finch.withoutIV <- ComIndex(traits=traits.finch, 
	index=funct.withoutIV, sp=sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), 
	ind.plot=ind.plot.finch, nperm=9, print=FALSE)


	##Plot T-statistics and custom metrics thanks to 
	#the plot.listofindex function.
	
	list.ind <- list(res.finch.withIV, res.finch.withoutIV, res.finch)
	index.list <- as.listofindex(list.ind, 
				  namesindex=c("mean", "kurtosis", "range", "CVNND",
				  "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop",
                  "T_IP.IC", "T_IC.IR", "T_PC.PR"))
	
	class(index.list)
	
	plot(index.list, plot.ask = FALSE)

	plot(index.list, plot.ask = FALSE, bysite = FALSE)

	ses.list <- ses.listofindex(index.list)
	ses.list
	attributes(ses.list)
	
	#### An other way to see "ses values" 
	
	# Custom theme (from rasterVis package)
	require(rasterVis)
	
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col = my.theme$regions$col)
	
	levelplot(t(rbind(ses.list[[1]]$ses, ses.list[[2]]$ses, 
	ses.list[[3]]$ses,  ses.list[[4]]$ses)), colorkey = my.ckey, 
	par.settings = my.theme,border = "black")
	
	levelplot(t(rbind(ses.list[[1]]$ses>ses.list[[1]]$ses.sup, 
	ses.list[[2]]$ses>ses.list[[2]]$ses.sup, 
	ses.list[[3]]$ses>ses.list[[3]]$ses.sup,
	ses.list[[4]]$ses>ses.list[[4]]$ses.sup)), 
	colorkey = my.ckey, par.settings = my.theme,border = "black")
	

	#For all metrics of the list of index
	ses.list.levelplot <- c()

	for(i in 1: length(ses.list)){

		ses.list.levelplot <- rbind(ses.list.levelplot, ses.list[[i]]$ses)
	}
	
	levelplot(t(ses.list.levelplot), colorkey = my.ckey, 
	par.settings = my.theme,border = "black")



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
