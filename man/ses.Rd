\name{ses}
\alias{ses}
\title{
	Standardized effect size and confidence interval for a matrix of statistics 
}

\description{
	calculation standardized effect size and confidence interval for a matrix of statistics and the related null model expressed as a list or as an array. Internal function use by other functions of the package. You can transpose the observed matrix to represent either the SES by traits or by plots; the result is transposed in the same way, except for a square matrix, whose transposition cannot be detected.
}

\usage{
	ses(obs = NULL, nullmodel = NULL, val.quant = c(0.025, 0.975), null.dims = NULL)
}

\arguments{
	\item{obs}{
		Observed matrix or vector of values.
	}
	
	\item{nullmodel}{
		Either a list or an array of three (two for a vector of observed values) dimensions corresponding to the null model permutations. 
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculation confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}

	\item{null.dims}{
		Integer vector of length 2 giving the dimensions of the \code{nullmodel} array that match the rows and the columns of \code{obs}; the remaining dimension holds the permutations. For example \code{c(3, 2)} for a null array of permutations x traits x plots and an observed matrix of plots x traits. By default (NULL), the layout stored by \code{\link{as.listofindex}} is used, or the dimensions are guessed from their sizes.
	}

}

\details{
	Warning: to detect automatically the correspondence between dimension of observed matrix and null model list or array, observed matrix, null model and permutations need to have different sizes. When two sizes are equal, give the layout with \code{null.dims}. Objects converted with \code{\link{as.listofindex}} (as done by \code{ses.listofindex}, \code{Pval} and the plot functions) carry their layout and need no \code{null.dims}.
}

\value{
	A list of three components:
	\item{$ses}{Observed value of standardized effect size.}
	\item{$ses.inf }{Lower limit of the confidence interval.}
	\item{$ses.sup }{Upper limit of the confidence interval.}
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{plot.listofindex}};
	\code{\link{plotSESvar}};
	\code{\link{ses.listofindex}}
}

\examples{
	data(finch.ind)
	
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)
	\donttest{
		ses(res.finch$Tstats$T_IP.IC, res.finch$Tstats$T_IP.IC_nm)
	}
}





