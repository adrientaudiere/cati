\name{CVNND}
\alias{CVNND}

\title{
	Coefficient of variation of the nearest neigbourhood distance (CVNND)
}

\description{
	Coefficient of variation of the nearest neigbourhood distance (CVNND). Internal use for other function.
}

\usage{
CVNND(trait, div_range= FALSE)
}

\arguments{
	\item{trait}{
		Trait vector.
	}
	
	\item{div_range}{
		Does CVNND need to be divise by the range.
	}
}

\details{
}

\value{
	One value corresponding to the coefficient of variation of the nearest neigbourhood for the trait.
}

\references{
	Jung, Vincent, Cyrille Violle, Cedric Mondy, Lucien Hoffmann, et Serge Muller. 2010. Intraspecific variability and trait-based community assembly: Intraspecific variability and community assembly. Journal of Ecology 98 (5): 1134-1140. doi:10.1111/j.1365-2745.2010.01687.x.
}

\author{
	Adrien Taudiere
}

\note{
}

\seealso{
}

\examples{
	data(finch.ind)
	CVNND(traits.finch[,1])
}





