\name{as.listofindex}
\alias{as.listofindex}
\title{
	Transform index results in a list of index
}

\description{
	Transform various results from functions Tstast, com.index or com.index.multi in a list of index. Useful to use the functions plot.listofindex (S3 method) and ses.listofindex.
}

\usage{
	as.listofindex(x, namesindex = NULL)
}

\arguments{
	\item{x}{
		A list of objects of class Tstast, com.index or com.index.multi
	}
	
	\item{namesindex}{
		Optionnal, the names of index in the same order as in x.
	}
}

\details{
}

\value{
	A list of observed values and corresponding "null" values (i.e. produced by null models) in the form "list(index1, null model index1, index2, null model index2 ...)"
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}

\seealso{
	\code{\link{ses.listofindex}};
	\code{\link{plot.listofindex}}
}

\examples{
	data(finch.ind)
	oldpar<-par()
	
	####
	#The function com.index allow to choose your own function 
	#(like mean, range, variance...) to calculate customize index.
	
	funct<-c("mean(x, na.rm=TRUE)", "kurtosis(x, na.rm=TRUE)", 
	"max(x, na.rm=TRUE) - min(x, na.rm=TRUE)", "CVNND(x)" )
	
	res.finch.sp_mn2<-com.index(traits=traits.finch, index=funct, sp=sp.finch,
	nullmodels=c(2,2,2,2), ind.plot=ind.plot.finch, nperm=9, print=FALSE)
	
	res.finch.sp_mn3<-com.index(traits=traits.finch, index=funct, sp=sp.finch,
	nullmodels=c(3,3,3,3), ind.plot=ind.plot.finch, nperm=9, print=FALSE)

	####
	#We can represent Standardized Effect Size (ses) using the 
	#function plot(as.listofindex(list1, list2, list3))
	
	list.ind2<-list(res.finch.sp_mn2, res.finch.sp_mn3)
	index.list2<-as.listofindex(list.ind2)
	
	plot(index.list2, type="bytraits")
	
	plot(index.list2)
}



