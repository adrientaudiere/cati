\name{bar_partvar}
\alias{bar_partvar}
\title{
	Barplot of variance partitionning
}

\description{
	Plot bars of the partition of variance calculate with the function partvar.
}

\usage{
	bar_partvar(partvar, col.bar = NA, leg = FALSE)
}

\arguments{
	\item{partvar}{
	The result of the partvar function.
	}
	
	\item{col.bar}{
	Vector of color for bar
	}
	
	\item{leg}{
	Logical value; print legend or not.
	}	
}

\details{
}

\value{
	None; used for the side-effect of producing a plot.
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}



\seealso{
	\code{\link{partvar}};
	\code{\link{pie_partvar}}
}

\examples{
	data(finch.ind)
	genus<-as.vector(unlist(strsplit(as.vector(sp.finch),"_")) [seq(1,length(sp.finch)*2, by=2)])
	res.partvar.finch<-partvar(traits = traits.finch, factors = cbind(sites = as.factor(as.vector(ind.plot.finch)), species = as.factor(as.vector(sp.finch)), genus=as.factor(genus)))
	res.partvar.finch

	oldpar<-par()
	par(mfrow=c(2,2), mai=c(0.2,0.2,0.2,0.2))
	pie_partvar(res.partvar.finch, col=c("red", "green", "blue", "purple"))
	par(oldpar)

	bar_partvar(res.partvar.finch, col=c("red", "green", "blue", "purple"), leg=TRUE)
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
