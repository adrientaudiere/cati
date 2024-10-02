\name{Fred}
\alias{Fred}

\title{
	Functional richness, evenness and divergence following Villeger et al. 2008
}

\description{
	Compute the 3 functional diversity indices (multi-traits) presented in Villeger et al. 2008 (Ecology 89 2290-2301): Functional richness (FRic), Functional evenness (FEve), Functional divergence (FDiv) 
}

\usage{
	Fred(traits, ind.plot)  
}

\arguments{
	\item{traits}{
		Individual Matrix of traits with traits in columns. NA are not allowed .
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individual is.
	}
}

\value{
	list of 4 vectors with values of indices in each sites
    \item{$nbind}{number of individuals}
	\item{$FRic}{functional richness index}
	\item{$FEve}{functional evenness index}
	\item{$FDiv}{functional divergence index}
}

\details{
	For each trait, values are standardized (mean=0 and standard deviation=1)
	For FRic computation, number of individuals must be higher than number of traits
}

\author{
	Sebastien Villeger sligthy modified by Adrien Taudiere
}

\seealso{
	\code{\link{ComIndexMulti}}
	\code{\link{ComIndex}}
}

\examples{
data(finch.ind)
\dontrun{
#For most multivariate functions we need to replace (or exclude) NA values. 
#For this example, we use the package mice to complete the data.

comm<-t(table(ind.plot.finch,1:length(ind.plot.finch)))

library(mice)
traits = traits.finch
mice<-mice(traits.finch)
traits.finch.mice<-complete(mice)
fred<-Fred(traits.finch.mice, ind.plot.finch)
}
}




