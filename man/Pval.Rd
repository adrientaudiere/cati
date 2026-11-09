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
		An object of class Tstats, ComIndex, ComIndexMulti or listofindex. Objects made by these functions or by \code{\link{as.listofindex}} carry the layout of their null values. For a hand-built list whose numbers of plots, traits or permutations are equal, the layout is guessed from sizes and may be wrong; use \code{\link{ses}} with \code{null.dims} in that case.
	}
	
	\item{na.rm}{
		A logical value indicating whether NA values should be stripped before the computation proceeds. 
	}
}

\value{
	A list of p-value for each metrics, traits and grouping if needed (e.g. sites). \code{pval.inf} is small when the observed value is lower than the null values, \code{pval.sup} when it is higher. Ties count against significance; NA null values are dropped; the p-value is NA when the observed value or all null values are NA.
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





