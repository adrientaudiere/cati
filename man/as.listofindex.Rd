\name{as.listofindex}
\alias{as.listofindex}
\title{
	Transform index results in a list of index
}

\description{
	Transform index result in a list of index. Useful to use the functions plot.listofindex (S3 method) and ses.listofindex.
}

\usage{
	as.listofindex(x, namesindex = NULL)
}

\arguments{
	\item{x}{
		A list of objects of class Tstast, com.index or com.index.multi
	}
	
	\item{namesindex}{
		Optionnal, the names of index in the same order than in x.
	}
}

\details{
}

\value{
	A list of observed and corresponding nullmodel value in the forme "list(index1, null model index1, index2, null model index2 ...)"
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}

\seealso{
	\code{\link{ses.listofindex}};
	\code{\link{plot.listofindex}}
}

\examples{
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
