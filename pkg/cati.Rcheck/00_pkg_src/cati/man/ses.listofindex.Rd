\name{ses.listofindex}
\alias{ses.listofindex}

\title{
	Standardized effect size for a list of index. 
}

\description{
	Standardized effect size and confidence interval for a list of index.
}

\usage{
	ses.listofindex(index.list = NULL, val.quant = c(0.025, 0.975))
}

\arguments{
	\item{index.list}{
		A list of index obtain using the function as.listofindex.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculation confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
}

\value{
	A list which each component correspond to the result of the ses function for an index. Further, each component is a list of three components:
	\item{$ses }{Observed value of standardized effect size.}
	\item{$ses.inf }{Lower limit of the confidence interval.}
	\item{$ses.sup }{Upper limit of the confidence interval.}
}

\author{
	Adrien	Taudiere
}

\seealso{
	\code{\link{as.listofindex}};
	\code{\link{ses}}
}


\examples{
	data(finch.ind)


	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, sp = sp.finch, 
	nperm = 9, print = FALSE)
	
	#calculation of means by population (name_sp_site is a population) 
	#like in the function ComIndex and determine the site 
	#for each population (sites_bypop)

	name_sp_sites = paste(sp.finch, ind.plot.finch,sep = "_")
	traits.by.pop <- apply(traits.finch, 2 , function (x) 
	tapply(x, name_sp_sites, mean , na.rm = TRUE))
	
	require(e1071)
	
	sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_")
	, split = "_"), function(x) x[3])

	funct.withoutIV <- c("tapply(x, unlist(sites_bypop), 
	function(x) mean(x, na.rm=TRUE))",
	"tapply(x, unlist(sites_bypop), function(x) kurtosis(x, na.rm=TRUE))",
	"tapply(x, unlist(sites_bypop), function(x)	max(x, na.rm = TRUE) 
	- min(x, na.rm = TRUE) )", 
	"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
		
	funct.withIV <- c("tapply(x, ind.plot.finch, 
	function(x) mean(x, na.rm = TRUE))",
	"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))",
	"tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE) 
	- min(x, na.rm = TRUE) )", 
	"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )


	res.finch.withIV <- ComIndex(traits = traits.finch, 
	index = funct.withIV, sp = sp.finch, nullmodels = rep("2", times=4), 
	ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

	res.finch.withoutIV <- ComIndex(traits = traits.finch, 
	index = funct.withoutIV, sp = sp.finch, nullmodels = rep("2sp", times=4), 
	ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

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
}
