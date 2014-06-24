\name{plot_randtest}
\alias{plot_randtest}


\title{
	Plot result of observed indices values against null distribution
}
\description{
	Function to plot result of observed indices values against null distribution.
}

\usage{
	plot_randtest(x, alter=c("greater", "less", "two-sided"), ...)
}

\arguments{
		
	\item{x}{
		An object of class listofindex, com.index, com.index.multi or Tstats.
	}
	
	\item{alter}{
		An object of class listofindex, com.index, com.index.multi or Tstats.
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
	\code{\link{com.index}};
	\code{\link{com.index.multi}};
	\code{\link{Tstats}};
	\code{\link{as.listofindex}};
	\code{\link{plot.listofindex}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, 
	sp=sp.finch, nperm=9, print=FALSE)
	
	par(mfrow=c(3,4))
	
	plot_randtest(res.finch)
	plot_randtest(res.finch, alter="two-sided")
}





