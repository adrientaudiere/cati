\name{com.index}
\alias{com.index}

\title{
	Calcul of custom statistics to test community assembly
}

\description{
	Calcul of custom statistics to test community assembly using null models.
}

\usage{
	com.index(traits = NULL, index = NULL, namesindex = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, nperm = 99, printprogress = TRUE)
}

\arguments{
	
	\item{traits}{
		Individual Traits matrix with traits in column.
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)".
	}
	
	\item{namesindex}{
		A vector of names for index.
	}
	
	\item{nullmodels}{
		A vector of values for null models coresponding to index. A value of 1 correspond to randomization of individual values whithin community. A value of 2 correspond to randomization of individual values whithin region. A values of 3 correspond to randomization of population values whithin region.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individu is. 
	}
	
	\item{sp}{
		Factor defining the species which the individu belong to.
	}
	
	\item{nperm}{
		Number of permutation. If NULL, only observed values are returned; 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
}

\details{
	Calcul of statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models.
	For each statistic this function return observed value and correspondant Null distribution.
	This function implement three null models which keep unchanged the number of individual per community.
	Models 1 correspond to randomization of individual values whithin community.
	Models 2 correspond to randomization of individual values whithin region.
	Models 3 correspond to randomization of population values whithin region.
	
	In most case, model 1 and 2 correspond to index at the individual level and the model 3 to index at the species (or any other aggregate variable like genus or family) level.
}

\value{
	An object of class "com.index" corresponding to a list of list:
	
	\item{$obs}{List of observed values for each traits in each community. Each component of the list correspond to a matrix containing the result for each custom function.}
	
	\item{$Null}{List of null values for each traits in each community. Each component of the list correspond to an array containing the result of the permutation for each custom function.}
	
	\item{$list.index}{List of index values and associate null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and associate null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species by sites.}
	\item{$namestraits}{Names of traits.}
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}

\seealso{
	\code{\link{com.index.multi}};
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{
	data(finch.ind)
	oldpar<-par()
	
	####
	#The function \code(com.index) allow to choose your own function (like mean, range, variance...) to calculate customize index.
	
	funct<-c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", "max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
	res.finch.sp_mn2<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	res.finch.sp_mn3<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, nperm=9, print=FALSE)

	####
	#We can represent Standardized Effect Size (ses) using the function \code(plot(as.listofindex(list1, list2, list3)))
	
	list.ind2<-list(res.finch.sp_mn2, res.finch.sp_mn3)
	index.list2<-as.listofindex(list.ind2)
	
	plot(index.list2, type="bytraits")
	
	plot(index.list2)

	####
	#This allows to calcul index by sites for example using \code("tapply(x, sites, mean)").
	
	funct<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
	
	##Null model 1 is trivial for this function because randomisation is whithin community only
	res.finch.ind_mn1<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(1,1,1,1), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	res.finch.ind_mn2<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)

	
	####
	#We can calcul index with or whithout intraspecific variance.
	#Calcul of means by population (name_sp_site is a name of a population) like in the function com.index and determine the site for each population (sites_bypop)

	name_sp_sites=paste(sp.finch, ind.plot.finch,sep="_")
	traits.by.pop<-apply(traits.finch, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=TRUE))
	
	sites_bypop<-lapply(strsplit(paste(rownames(traits.by.pop), sep="_"), split="_"), function(x) x[3])
	
	funct.whithoutIV<-c("tapply(x, unlist(sites_bypop), function(x) mean(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
	
	
	funct.whithIV<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
	
	res.finch.whithIV<-com.index(traits=traits.finch, index=funct.whithIV, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	res.finch.whithoutIV<-com.index(traits=traits.finch, index=funct.whithoutIV, sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	####
	#We can also represent T-statistics and custom index thanks to the \code(plot.listofindex) function.
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=9, print=FALSE)

	list.ind<-list(res.finch.whithIV, res.finch.whithoutIV ,res.finch)
	index.list1<-as.listofindex(list.ind, namesindex=c("mean", "kurtosis", "range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", "T_IP.IC", "T_IC.IR", "T_PC.PR"))
	class(index.list1)
	
	par(mfrow=c(2,3))
	plot(index.list1,type="bytraits", bysites=TRUE)
	
	par(mfrow=c(2,2))
	plot(index.list1,type="bytraits")
	par(mfrow=c(1,1))
	
	plot(index.list1,type="simple")
	plot(index.list1,type="simple_range")
	plot(index.list1,type="normal")
	plot(index.list1,type="barplot")
	
}



\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
