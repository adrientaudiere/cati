\name{AbToInd}
\alias{AbToInd}

\title{
	Internal function. Transform abundance data matrix into individual like matrix.
}
\description{
	Transform abundance data matrix into individual like matrix to allows the use of ComIndex and ComIndexMulti on populationnal or specific traits values. 
}

\usage{
	AbToInd(traits, com, type.sp.val = "count")
}

\arguments{
	\item{traits}{
		Individual Matrix of traits with traits in columns. "traits" matrix must have row names (e.g. species or populationnal names).
	}
	
	\item{com}{
		Community data matrix with species in rows and sites in column.
	}
	
	\item{type.sp.val}{
		Either "count" or "abundance". Use abundance when all values in the com matrix are not superior to one. Using abundance is EXPERIMENTAL. This function round abundance to fit count data. 
	}
}

\details{
	Internal function
}

\value{
	A list of objects:
	
	\item{$traits}{Individual traits matrix}
	\item{$sp}{Vector of species attributes}
	\item{$ind.plot}{Vector of sites attributes}
	
}

\author{
	Adrien Taudiere
}
