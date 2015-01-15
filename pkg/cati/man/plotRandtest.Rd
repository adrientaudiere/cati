\name{plotRandtest}
\alias{plotRandtest}


\title{
	Plot result of observed indices values against null distribution
}
\description{
	Function to plot result of observed indices values against null distribution.
}

\usage{
	plotRandtest(x, alternative = "two-sided", ...)
}

\arguments{
		
	\item{x}{
		An object of class listofindex, ComIndex, ComIndexMulti or Tstats.
	}
	
	\item{alternative}{
		Indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. "greater" corresponds to positive association, "less" to negative association.
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
	\code{\link{ComIndex}};
	\code{\link{ComIndexMulti}};
	\code{\link{Tstats}};
	\code{\link{as.listofindex}};
	\code{\link{plot.listofindex}}
}

\examples{
	data(finch.ind)
	\dontrun{
		res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
		sp = sp.finch, nperm = 99, print = FALSE)
	
	par(mfrow = c(4,4))
	
	plotRandtest(res.finch)
	plotRandtest(res.finch, alter = "less")
	}
}





