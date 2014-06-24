\name{ab_to_ind}
\alias{ab_to_ind}

\title{
	Internal function. Transform abundance data matrix into individual like matrix.
}
\description{
	Transform abundance data matrix into individual like matrix to allows the use of com.index and com.index.multi on populationnal or specific traits values. 
}

\usage{
	ab_to_ind(traits, com, type="count")
}

\arguments{
	\item{traits}{
		Individual Matrix of traits with traits in columns. "traits" matrix must have row names (e.g. species or populationnal names).
	}
	
	\item{com}{
		Community data matrix with species in rows and sites in column.
	}
	
	\item{type}{
		Either "count" or "abundance". Use abundance when all values in the com matrix are not superior to one.
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
