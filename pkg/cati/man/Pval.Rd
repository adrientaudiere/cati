\name{Pval}
\alias{Pval}

\title{
	Calcul of p-value for object of class Tstats, ComIndex, ComIndexMulti and listofindex
}

\description{
	Calcul of p-value for object of class Tstats, ComIndex, ComIndexMulti and listofindex. 
	This test equates to finding the quantile in exp in which obs would be found (under a one-tailed test).
}

\usage{
	Pval(x, na.rm = TRUE)
}

\arguments{
	\item{x}{
		An object of class Tstats, ComIndex, ComIndexMulti or listofindex. 
	}
	
	\item{na.rm}{
		A logical value indicating whether NA values should be stripped before the computation proceeds. 
	}
}

\value{
	A list of p-value for each metrics, traits and grouping if needed (e.g. sites) 
}

\author{
	Adrien Taudiere
}

\examples{

data(finch.ind)
res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)
  \dontrun{
    Pval(res.finch)
  }
}





