\name{samplingSubsetData}
\alias{samplingSubsetData}

\title{
	Sampling subset of data.
}

\description{
	Sampling subset of data.
}

\usage{
	sampling_subset_data(d = NULL, by.factor = NULL, nperm = 2, 
	type = "proportion_1", prop = seq(10,100, by = 10)) {
}

\arguments{
	
	\item{d}{
		Dataframe of data to sample. Each line is an individual.
	}
	
	\item{by.factor}{
		This factor define the sampling unit to impoverish.
	}
	
	\item{nperm}{
		Number of permutation. 
	}
	
	\item{type}{
		Type of sampling. \strong{count_1} sample a number of individual with one individual by factor at least. \strong{count_0} sample a number of individual without any minimum by factor. \strong{proportion_1} sample a proportion of individual with one individual by factor at least. \strong{proportion_0} sample a proportion of individual without any minimum by factor.
	}
	
	\item{prop}{
		Categorical proportion to sample.
	}
}

\details{
	Sampling subset of data.
}

\value{
}

\author{
	Adrien Taudiere
}

