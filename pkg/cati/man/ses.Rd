\name{ses}
\alias{ses}
\title{
	Standardized effect size and confidence interval for a matrix of statistics 
}

\description{
	Calcul standardized effect size and confidence interval for a matrix of statistics and the related null model expressed as a list or as an array. Internal function use by other functions of the package. You can transpose the observed matrix to represent either the SES by traits or by plots. Warnings, to detect automatically the correspondence between dimension of observed matrix and null model list or array, observed matrix needs to have different numbers of rows and columns.
}

\usage{
	ses(obs = NULL, nullmodel = NULL, val.quant = c(0.025, 0.975))
}

\arguments{
	\item{obs}{
		Observed matrix or vector of values.
	}
	
	\item{nullmodel}{
		Either a list or an array of three (two for a vector of observed values) dimensions corresponding to the null model permutations. 
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}

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

	res.finch <- Tstats(traits.finch, ind.plot=ind.plot.finch, 
	sp=sp.finch, nperm=9)

	
	ses(res.finch$T_IP.IC, res.finch$T_IP.IC_nm)
	
	ses(t(res.finch$T_IP.IC), res.finch$T_IP.IC_nm)
}





