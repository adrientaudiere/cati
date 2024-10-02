\name{samplingSubsetData}
\alias{samplingSubsetData}

\title{
	Sampling subset of data.
}

\description{
	Sampling subset of data.
}

\usage{
	samplingSubsetData(d = NULL, sampUnit = NULL, nperm = 9, 
	type = "proportion", prop = seq(10, 100, by = 10), MinSample = 1,
	Size = NULL)
}

\arguments{

	\item{d}{
		Dataframe of data to sample. Each line is an individual.
	}

	\item{sampUnit}{
		A Factor defining the sampling unit to impoverish. For example it can be the species or the plot attributes of each individual. 
	}

	\item{nperm}{
		Number of permutations. 
	}

	\item{type}{
		Type of sampling. Either \strong{proportion}, \strong{count}, \strong{propBySize} or \strong{factorBySize}. See details.
	}

	\item{prop}{
		Integer between 1 and 100. Categorical proportions to sample in percentage. 
	}
	
	\item{MinSample}{
		Minimum number of individual to sample by sample unit. Default is one.
	}
	
	\item{Size}{
		A vector of value for each individual (type \strong{propBySize} and \strong{factorBySize}) or for each level of factor (\strong{factorBySize} only). Determine the rank of individual/factor when using the sampling schemes \strong{propBySize} and \strong{factorBySize}.
	}
}

\details{
	Sampling scheme \strong{count} sample a number of individuals wheras \strong{proportion} sample a proportion of individuals by sample unit. Sampling scheme \strong{propBySize} sample in each sampling unit (sampUnit) a proportion of the individual ranked using the argument \strong{Size}. Consequently, the bigest individuals (higher Size) will be sample before the smaller one. \strong{factorBySize} sample a proportion of sampling unit (sampUnit) ranked using the argument \strong{Size}. For example you can sample only the individuals of the 20\% of the more aboundant species.
}

\value{
	Return a list list of sample dataframe. The first level of the list depicts the permutation and the second level depicts the different proportion/number of individual sampled by factor. 
}

\author{
	Adrien Taudiere
}
