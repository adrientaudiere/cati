\name{MinMaxMST}
\alias{MinMaxMST}

\title{
	Ratio of the shortest distance to the longest distance in a minimum spanning tree
}

\description{
	Ratio of the shortest distance to the longest distance in a minimum spanning tree.
}

\usage{
	MinMaxMST(traits, gower.dist = TRUE, scale.tr = TRUE, method.dist = "euclidian")
}

\arguments{
	\item{traits}{
		Traits matrix (traits in column)
	}
	
	\item{gower.dist}{
		Calculate gower distance using the function gowdis from package FD.
	}	
	
	\item{scale.tr}{
		Does traits need to be scale before multi-traits metric calculation? Only use when gower.dist = FALSE. Default is yes.
	}
	
	\item{method.dist}{
		Method to calculate the distance in case of multi-traits metric (function dist). Only use when gower.dist = FALSE. Default is euclidian.
	}
}

\value{
	The value of the ratio of the shortest distance to the longest distance in a minimum spanning tree.
}

\references{
	Stubbs, WJ., and Wilson, JB. 2004. Evidence for limiting similarity in a sand dune community. Journal of Ecology 92: 557-567.
	Aiba, M., Katabuchi, M., Takafumi, H., Matsuzaki, S.S., Sasaki, T. & Hiura, T. 2013. Robustness of trait distribution metrics for community assembly studies under the uncertainties of assembly processes. Ecology, 94, 2873-2885.
}

\author{
	Aiba et al., 2013 modified by Adrien Taudiere
}

\examples{
\dontrun{

	data(finch.ind)
	
	MinMaxMST(traits.finch[1:10,])
	MinMaxMST(traits.finch[1:10,], gower.dist = FALSE)
	MinMaxMST(traits.finch[1:10,], gower.dist = FALSE, scale.tr = FALSE)
}

}











