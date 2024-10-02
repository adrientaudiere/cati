\name{ses.listofindex}
\alias{ses.listofindex}

\title{
	Standardized effect size for a list of index. 
}

\description{
	Standardized effect size and confidence interval for a list of index.
}

\usage{
	ses.listofindex(index.list = NULL, val.quant = c(0.025, 0.975))
}

\arguments{
	\item{index.list}{
		A list of index obtain using the function as.listofindex.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculation confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
}

\value{
	A list which each component correspond to the result of the ses function for an index. Further, each component is a list of three components:
	\item{$ses }{Observed value of standardized effect size.}
	\item{$ses.inf }{Lower limit of the confidence interval.}
	\item{$ses.sup }{Upper limit of the confidence interval.}
}

\author{
	Adrien	Taudiere
}

\seealso{
	\code{\link{as.listofindex}};
	\code{\link{ses}}
}
