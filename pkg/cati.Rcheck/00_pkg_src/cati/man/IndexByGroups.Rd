\name{IndexByGroups}
\alias{IndexByGroups}

\title{
	Apply metrics to groups.
}

\description{
	Transforme a list of metrics to apply them to groups, typically to populations.
}

\usage{
	IndexByGroups(metrics, groups)
}

\arguments{
	\item{metrics}{
		A vector of metrics like the argument "index" of function ComIndex 
	}
	
	\item{groups}{
		Name of the factor to apply the metrics to groups in the form "pop", e.g. population
	}
}

\value{
	A vector of transformed metrics
}

\author{
	Adrien Taudiere
}





