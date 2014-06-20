\name{pie_partvar}
\alias{pie_partvar}

\title{
	Pie plot of variance partitioning.
}

\description{
	Plot pies of the partitioning of variance calculated using the function partvar.
}

\usage{
	pie_partvar(partvar, col.pie = NA)
}

\arguments{
	\item{partvar}{
		The result of the partvar function.
	}
	
	\item{col.pie}{
		Vector of color for pie.
	}
	
	\item{\dots}{
		Any additional arguments are passed to the pie function.
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
	\code{\link{bar_partvar}}
}

\examples{
	data(finch.ind)
	
	res.partvar.finch<-partvar( traits=traits.finch, 
	factors=cbind(as.factor(as.vector(ind.plot.finch)), 
	as.factor(as.vector(sp.finch)), 
	as.factor(as.vector(unlist(strsplit(as.vector(sp.finch),"_")) 
	[seq(1,length(sp.finch)*2, by=2)])) ) )
	
	pie_partvar(res.partvar.finch, col=c("red", "green", "blue"))

	bar_partvar(res.partvar.finch, col=c("red", "green", "blue"))
}





