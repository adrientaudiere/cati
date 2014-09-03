\name{SumBL}
\alias{SumBL}

\title{
	Sum of branch length of a classification dendrogram (Petchey and Gaston, 2002) 
}

\description{
	Sum of branch length of a classification dendrogram (Petchey and Gaston, 2002) 
}

\usage{
	SumBL(traits, gower.dist = TRUE, method.hclust = "average", 
	scale.tr = TRUE, method.dist = "euclidian")
}

\arguments{
	\item{traits}{
		Traits matrix (traits in column)
	}
	
	\item{gower.dist}{
		Calculate gower distance using the function gowdis from package FD.
	}	
	
	\item{method.hclust}{
		Define the method for the hclust function (default is "average" i.e. UPGMA).
	}
		
	\item{scale.tr}{
		Does traits need to be scale before multi-traits metric calculation? Only use when gower.dist = FALSE. Default is yes.
	}
	
	\item{method.dist}{
		Method to calculate the distance in case of multi-traits metric (function dist). Only use when gower.dist = FALSE. Default is euclidian.
	}
}

\value{
	The value of the sum of branch length from a classification dendrogram of traits.
}

\references{
Petchey, OL., and Gaston, KJ. 2002. Functional diversity (FD), species richness and community composition. Ecology Letters 5:402-411
}

\author{
	Adrien Taudiere
}

\examples{

\dontrun{

data(finch.ind)
SumBL(traits.finch)
SumBL(traits.finch, gower.dist = FALSE)

}

}





