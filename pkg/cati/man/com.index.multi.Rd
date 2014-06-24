\name{com.index.multi}
\alias{com.index.multi}

\title{
	Computing multitraits metrics to test and quantify the non-random assembly of communities
}

\description{
	Computing multitraits metrics to test and quantify the non-random assembly of communities
}

\usage{
	com.index.multi(traits = NULL, index = NULL, by.factor = NULL, 
	namesindex = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, 
	com=NULL, reg.pool=NULL, nperm = 99, printprogress = TRUE, 
	ind.value=TRUE, type="count")
}

\arguments{
 
	\item{traits}{
		Individual Matrix of traits with traits in column.
	}
	
	\item{index}{
		A vector of functions to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)".
	}
	
	\item{by.factor}{
		A factor to split the Matrix of traits and compute index for each subset eg for each site.  
	}
	
	\item{namesindex}{
		A vector of names for metrics.
	}
	
	\item{nullmodels}{
		A vector of values for null models coresponding to index. A value of 1 correspond to randomization of individual values within community. A value of 2 correspond to randomization of individual values within region. A values of 3 correspond to randomization of population values within region.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot (site or community) in which the individual is. 
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
	}
	
	\item{com}{
		Community data matrix with species (or populations) in rows and sites in column. Use only if ind.value=FALSE. 
		"traits" matrix and "com" matrix must have the same number of rows.
	}
	
	\item{reg.pool}{
		Regional pool data for traits. If not informed, traits is considere as the regional pool. This matrix need to be larger than traits. Use only for null model 2.
	}	
	
	\item{nperm}{
		Number of permutations. If NULL, only observed values are returned; 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
	
	\item{ind.value}{
		Do the data are from individual value. If not, an internal function transform the abundance data matrix of com into a individual like matrix to compute the function.
	}
	
	\item{type}{
		Only if ind.value=FALSE. Either "count" or "abundance". Use abundance when all values in the com matrix are not superior to one.
	}
	
	
}

\value{
	A list of lists:
	
	\item{$obs}{List of observed values for each trait in each community. Each component of the list correspond to a matrix containing the result for each custom function.}
	
	\item{$Null}{List of null values for each trait in each community. Each component of the list correspond to an array containing the result of the permutations for each custom function.}
	
	\item{$list.index}{List of index values and associate null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and associate null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species per site.}
	\item{$namestraits}{Names of traits.}
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{com.index}};
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{
	data(finch.ind)
	
	####
	#For most multivariate functions we need to replace (or exclude)
	#NA values.
	
	#For this example, we use the package mice to complete the data.
	
	\dontrun{
		names.sp_ind_plot<-as.factor(paste(sp.finch, ind.plot.finch, sep="_")) 
		
		comm<-t(table(ind.plot.finch,1:length(ind.plot.finch)))
		
		library(mice)
		traits=traits.finch
		mice<-mice(traits.finch)
		traits.finch.mice<-complete(mice)
		
		####
		#A simple example to illustrate the concept of the function 
		#com.index.multi
		
		res.sum.1<-com.index.multi(traits.finch, 
		index=c("sum(scale(x), na.rm=TRUE)", "sum(x, na.rm=TRUE)"), 
		by.factor=names.sp_ind_plot, nullmodels=c(2,2), 
		ind.plot=ind.plot.finch, nperm=50, sp=sp.finch)
		
		attributes(ses.listofindex(as.listofindex(res.sum.1)))
		
		####
		#A more interesting example using the function hypervolume 
		#from the package hypervolume. 
		#We show here several results which differe in there factor 
		#that delimit the group to calculate different hypervolume 
		#(argument by_factor). 
		
		require(hypervolume)
		
		res.hv.1<-com.index.multi(traits.finch.mice, index=c("as.numeric
		(try(hypervolume(na.omit(x), reps=100, bandwidth=0.2, verbose=FALSE, 
		warnings=FALSE)@Volume))"), by.factor=rep(1,length(names.sp_ind_plot)), 
		nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=9, sp=sp.finch)
		
		res.hv.2<-com.index.multi(traits.finch.mice, index=c("as.numeric
		(try(hypervolume(na.omit(x), reps=100, bandwidth=0.2, verbose=FALSE, 
		warnings=FALSE)@Volume))"), by.factor=names.sp_ind_plot, 
		nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=9, sp=sp.finch)
		
		res.hv.3<-com.index.multi(traits.finch.mice, index=c("as.numeric
		(try(hypervolume(na.omit(x), reps=100, bandwidth=0.2, verbose=FALSE, 
		warnings=FALSE)@Volume))"), by.factor=ind.plot.finch, 
		nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=9, sp=sp.finch)
		
		res.hv.4<-com.index.multi(traits.finch.mice, index=c("as.numeric
		(try(hypervolume(na.omit(x), reps=100, bandwidth=0.2, verbose=FALSE, 
		warnings=FALSE)@Volume))"), by.factor=sp.finch, 
		nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=9, sp=sp.finch)
		
		list.ind.multi<-as.listofindex(list(res.hv.2, res.hv.3, res.hv.4))
		
		ses.listofindex(list.ind.multi)
		
		plot(list.ind.multi)
		plot(list.ind.multi, xlim=c(-200,20))
	}
}





