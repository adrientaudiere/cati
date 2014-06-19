\name{com.index}
\alias{com.index}

\title{
	Computing the moments of the trait distribution and other metrics to test and quantify the non-random assembly of communities
}

\description{
	Computing the moments of the trait distribution and other metrics to test and quantify the non-random assembly of communities. 
}

\usage{
	com.index(traits = NULL, index = NULL, namesindex = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, nperm = 99, printprogress = TRUE)
}

\arguments{
	
	\item{traits}{
		Individual Matrix of traits with traits in column.
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)", see examples for more complexe functions.
	}
	
	\item{namesindex}{
		A vector of names for index.
	}
	
	\item{nullmodels}{
		A vector of values corresponding to null models tu use for each index. A value of 1 corresponds to a randomization of individual values within a given community. A value of 2 corresponds to randomization of individual values within region, ie within all the dataset. A value of 3 corresponds to randomization of population values (each individual value are replaced by the mean value of it population) within region. For example, if nullmodels=c(1,2), the first index will be calculated on the null model 1 and the second index on the null model 2.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot (site or community) in which the individual is. 
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
	}
	
	\item{nperm}{
		Number of permutations. If NULL, only observed values are returned; 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
}

\details{
	Compute statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models.
	For each statistic this function returns observed values and the related null distribution.
	This function implement three null models which keep unchanged the number of individual per community.
	Model 1 corresponds to randomization of individual values within community.
	Model 2 corresponds to randomization of individual values within region.
	Model 3 corresponds to randomization of population values within region.
	
	In most cases, models 1 and 2 correspond to index at the individual level and the model 3 to index at the species level (or any other aggregate variable like genus, family or functionnal group).
}

\value{
	An object of class "com.index" corresponding to a list of lists:
	
	\item{$obs}{List of observed values for each trait in each community. Each component of the list corresponds to a matrix containing the result for each custom function.}
	
	\item{$Null}{List of null values for each trait in each community. Each component of the list corresponds to an array containing the result of the permutations for each custom function.}
	
	\item{$list.index}{List of index values and related null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and related null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species per site.}
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
	#The function com.index allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	funct<-c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", "max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)"  )
	res.finch.sp_mn2<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	res.finch.sp_mn3<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, nperm=9, print=FALSE)

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
	
	funct<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
	
	##Null model 1 is trivial for this function 
	#because randomisation is within community only
	res.finch.ind_mn1<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(1,1,1,1), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	res.finch.ind_mn2<-com.index(traits=traits.finch, index=funct, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)

	
	####
	#We can calcul metrics with or without intraspecific variance.
	#Calculation of trait averages per population 
	#(name_sp_site is a name of a population) 
	#like in the function com.index
	#and determine the site for each population (sites_bypop)

	name_sp_sites=paste(sp.finch, ind.plot.finch, sep="_")
	traits.by.pop<-apply(traits.finch, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=TRUE))
	
	sites_bypop<-lapply(strsplit(paste(rownames(traits.by.pop), sep="_"), split="_"), function(x) x[3])
	
	funct.withoutIV<-c("tapply(x, unlist(sites_bypop), function(x) mean(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, unlist(sites_bypop), function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )
	
	
	funct.withIV<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm=TRUE))", "tapply(x, ind.plot.finch, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE) )", "tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
	
	res.finch.withIV<-com.index(traits=traits.finch, index=funct.withIV, sp=sp.finch, nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	res.finch.withoutIV<-com.index(traits=traits.finch, index=funct.withoutIV, sp=sp.finch, nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	
	####
	#We can also represent T-statistics and custom index thanks to
	#the plot.listofindex function.
	
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=9, print=FALSE)

	list.ind<-list(res.finch.withIV, res.finch.withoutIV ,res.finch)
	index.list1<-as.listofindex(list.ind, namesindex=c("mean", "kurtosis", "range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", "T_IP.IC", "T_IC.IR", "T_PC.PR"))
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
