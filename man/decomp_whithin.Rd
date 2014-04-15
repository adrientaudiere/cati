\name{decomp_whithin}
\alias{decomp_whithin}
\alias{barplot.decomp_whithin}

\title{
	Variance decomposition for multiple traits
}

\description{
	This function decomposes the variation in community trait composition into three sources – (i) the intraspecific trait variability, (ii) the
	variability due to species turnover and (iii) their covariation is also separated. This decomposition is computed for the whole variation in the trait
	values and, The formula specified, across the contribution of various explanatory variables considered in the model. Barplot.decomp_whithin allow to plot the result of the decomposition.
}

\usage{
	decomp_whithin(traits = NULL, formula = ~1, ind.plot = NULL, sp = NULL, printprogress = TRUE, ...)
	
	\method{barplot}{decomp_whithin}(height, resume = TRUE, \dots)
}

\arguments{
	\item{traits}{
		Traits matrix with traits in column
	}
	
	\item{height}{
		An object of class decomp_whithin obtain by the function decomp_whithin.
	}
	
	\item{formula}{
		The formula parameter must be a one-sided formula, i.e. starting with a tilde (~) character. The response variable is specified by the next two arguments, specif.avg and const.avg. By default set to ~1.		
	}

	\item{ind.plot}{
		Factor defining the name of the plot in which the individu is. 
	}
	
	\item{sp}{
		Factor defining the species which the individu belong to.
	}

	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
	
	\item{resume}{
		Logical. If resume=FALSE, plot one graphic by traits.
	}
	
	\item{\dots}{
		Optional additional arguments
	}
}

\details{
}

\value{
	An object of class "decomp_whithin".
}

\references{
	Leps, Jan, Francesco de Bello, Petr Smilauer and Jiri Dolezal. 2011. Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects. Ecography 34 (5): 856-863.
}

\author{
	Adrien Taudiere
	Jan Leps
}

\note{
}

\seealso{
	\code{\link{barplot.decomp_whithin}};
	\code{\link{trait.flex.anova}}
}

\examples{
	data(finch.ind)
	
	res.decomp<-decomp_whithin(traits=traits.finch, sp=sp.finch, ind.plot=ind.plot.finch, print=FALSE)
	
	barplot.decomp_whithin(res.decomp)
	
	par(mfrow=c(2,2))
	barplot.decomp_whithin(res.decomp, resume=FALSE)
	par(mfrow=c(1,1))
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
