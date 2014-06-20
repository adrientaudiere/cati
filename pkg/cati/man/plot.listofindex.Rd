\name{plot.listofindex}
\alias{plot.listofindex}

\title{
	Plot community assembly index
}

\description{
	Plot community assembly index and confidence intervals using a list of index. S3 method for class listofindex.
}

\usage{
	\method{plot}{listofindex}(x, type = "normal", 
	col.index = c("red", "purple", "green"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), 
	grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, 
	cex.text = 0.8, plot.ask = FALSE, srt.text = 90, 
	bysite = FALSE, ...)
}

\arguments{

	\item{x}{
		A list of index and related null models obtained from to the as.listofindex function.
	}
	
	\item{type}{
		Type of plot. Possible type = "simple", "simple_range", "normal", "barplot" and "bytraits".
	}
	
	\item{col.index}{
		Vector of colors for index.
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{color.cond}{
		Logical value; If color.cond=TRUE, color points indicate T-statistics values significatively different from the null model and grey points are not different from null model.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculate confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{grid.v}{
		Logical value; print vertical grid or not
	}
	
	\item{grid.h}{
		Logical value; print horizontal grid or not
	}
	
	\item{xlim}{
		Numeric vectors of length 2, giving the x coordinates range 
	}
	
	\item{ylim}{
		Numeric vectors of length 2, giving the y coordinates range 
	}
		
	\item{cex.text}{
		Numeric value; the magnification to be used for text relative to the current setting of cex
	}
	
	\item{plot.ask}{
		Logical value; ask for plotting the next plot or not.
	}
	
	\item{srt.text}{
		Degree of rotation for text.
	}
	
	\item{bysite}{
		Logical value; plot per site or by traits. 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}

\value{
	None; used for the side-effect of producing a plot.
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{as.listofindex}};
	\code{\link{plot.Tstats}};
	\code{\link{ses.listofindex}}
}

\examples{	
	data(finch.ind)
	oldpar<-par()
	
	####
	#The function com.index allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	require(e1071)
	
	funct<-c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
	
	res.finch.sp_mn2<-com.index(traits=traits.finch, index=funct, 
	sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)
	
	res.finch.sp_mn3<-com.index(traits=traits.finch, index=funct, 
	sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)

	####
	#We can represent Standardized Effect Size (ses)
	#using the function plot(as.listofindex(list1, list2, list3))
	
	list.ind2<-list(res.finch.sp_mn2, res.finch.sp_mn3)
	index.list2<-as.listofindex(list.ind2)
	
	plot(index.list2, type="bytraits")
	
	plot(index.list2)

	####
	#This allows to calcul index per site 
	#for example using "tapply(x, sites, mean)".
	
	funct<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", 
	"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", 
	"tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - 
	min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))")
	
	
	##Null model 1 is trivial for this function 
	#because randomisation is within community only
	
	res.finch.ind_mn1<-com.index(traits=traits.finch, index=funct, 
	sp=sp.finch, nullmodels=c(1,1,1,1), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)
	
	res.finch.ind_mn2<-com.index(traits=traits.finch, index=funct, 
	sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)

	
	####
	#We can calcul metrics with or without intraspecific variance.
	#Calculation of trait averages per population 
	#(name_sp_site is a name of a population) 
	#like in the function com.index
	#and determine the site for each population (sites_bypop)

	name_sp_sites=paste(sp.finch, ind.plot.finch, sep="_")
	
	traits.by.pop<-apply(traits.finch, 2 , function (x) 
	tapply(x, name_sp_sites, mean , na.rm=TRUE))
	
	sites_bypop<-lapply(strsplit(paste(rownames(traits.by.pop), sep="_"), 
	split="_"), function(x) x[3])
	
	funct.withoutIV<-c("tapply(x, unlist(sites_bypop), 
	function(x) mean(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
	function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), 
	function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
	"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
	
	
	funct.withIV<-c("tapply(x, ind.plot.finch, function(x) 
	mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
	kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) 
	max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", 
	"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
	
	res.finch.withIV<-com.index(traits=traits.finch, index=funct.withIV, 
	sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)
	
	res.finch.withoutIV<-com.index(traits=traits.finch, index=funct.withoutIV, 
	sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, 
	nperm=9, print=FALSE)
	
	
	####
	#We can also represent T-statistics and custom index thanks to
	#the plot.listofindex function.
	
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, 
	nperm=9, print=FALSE)

	list.ind<-list(res.finch.withIV, res.finch.withoutIV ,res.finch)
	
	index.list1<-as.listofindex(list.ind, namesindex=c("mean", "kurtosis", 
	"range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
	"T_IP.IC", "T_IC.IR", "T_PC.PR"))
	
	class(index.list1)
	
	par(mfrow=c(2,3))
	plot(index.list1,type="bytraits", bysite=TRUE)
	
	par(mfrow=c(2,2))
	plot(index.list1,type="bytraits")
	par(mfrow=c(1,1))
	
	plot(index.list1,type="simple")
	plot(index.list1,type="simple_range")
	plot(index.list1,type="normal")
	plot(index.list1,type="barplot")
}




