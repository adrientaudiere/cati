\name{partvar}
\alias{partvar}

\title{
	Variance partitioning accross nested scales
}

\description{
	Variance partitioning accross nested scales using a decomposition (varcomp function) of variance on restricted maximum likelihood (REML) method (lme function). See Messier et al. 2010 for more information. 
}

\usage{
	partvar(traits, factors, printprogress = TRUE)
}

\arguments{
	\item{traits}{
		Matrix of traits with traits in column
	}
  
	\item{factors}{
		A matrix of factors with the first column corresponds to the higher level factor, the second row the second higher level factor and so on.  
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
}

\details{
}

\value{
	An object of class "partvar" corresponding to a matrix of variance values with traits in rows and nested factors in column.
}

\references{
	Messier, Julie, Brian J. McGill, et Martin J. Lechowicz. 2010. How do traits vary across ecological scales? A case for trait-based ecology: How do traits vary across ecological scales? Ecology Letters 13(7): 838-848. doi:10.1111/j.1461-0248.2010.01476.x.
}

\author{
	Adrien Taudiere
	Julie Messier
}

\note{
}

\seealso{
	\code{\link{pie_partvar}};
	\code{\link{bar_partvar}}
}

\examples{
	data(finch.ind)

	genus<-as.vector(unlist(strsplit(as.vector(sp.finch),"_")) [seq(1,length(sp.finch)*2, by=2)])
	res.partvar.finch<-partvar(traits=traits.finch, factors=cbind(sites=as.factor(as.vector(ind.plot.finch)), species=as.factor(as.vector(sp.finch)), genus=as.factor(genus)))

	res.partvar.finch
	
	oldpar<-par()
	par(mfrow=c(2,2), mai=c(0.2,0.2,0.2,0.2))
	pie_partvar(res.partvar.finch, col=c("red", "green", "blue", "purple"))
	par(oldpar)

	bar_partvar(res.partvar.finch, col=c("red", "green", "blue", "purple"), leg=TRUE)
}





