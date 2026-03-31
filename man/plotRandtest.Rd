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
\donttest{
  # Simulate a small NA-free dataset to avoid empty null distributions
  set.seed(42)
  n <- 60
  traits_sim <- matrix(rnorm(n * 2), nrow = n, ncol = 2,
                       dimnames = list(NULL, c("T1", "T2")))
  sp_sim <- factor(rep(paste0("sp", 1:10), each = 6))
  plot_sim <- factor(rep(paste0("P", 1:3), each = 20))
  res_sim <- Tstats(traits_sim, ind.plot = plot_sim,
                   sp = sp_sim, nperm = 9, printprogress = FALSE)
  plotRandtest(res_sim)
}
}





