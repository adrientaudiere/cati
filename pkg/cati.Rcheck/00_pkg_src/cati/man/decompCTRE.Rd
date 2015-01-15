\name{decompCTRE}
\alias{decompCTRE}
\alias{barplot.decompCTRE}

\title{
	Variance partitioning for multiple traits
}

\description{
	This function decomposes the variation in community trait composition into three sources: (i) the intraspecific trait variability, (ii) the
	variability due to species turnover and (iii) their covariation is also separated. This decomposition is computed for the whole variation in the trait
	values and, The formula specified, across the contribution of various explanatory variables considered in the model. Barplot.decompCTRE allow to plot the result of the decomposition.
}

\usage{
	decompCTRE(traits = NULL, formula = ~1, ind.plot = NULL, sp = NULL, 
	printprogress = TRUE, ...)
	
	\method{barplot}{decompCTRE}(height, resume = TRUE, \dots)
}

\arguments{
	\item{traits}{
		Matrix of traits with traits in column
	}
	
	\item{height}{
		An object of class decompCTRE obtain by the function decompCTRE.
	}
	
	\item{formula}{
		The formula parameter must be a one-sided formula, i.e. starting with a tilde (~) character. The response variable is specified by the next two arguments, specif.avg and const.avg. By default set to ~1.		
	}

	\item{ind.plot}{
		Factor defining the name of the plot (site or community) in which the individual is. 
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
	}

	\item{printprogress}{
		Logical value; print progress during the calculation or not.
	}
	
	\item{resume}{
		Logical. If resume = FALSE, plot one graphic by traits.
	}
	
	\item{\dots}{
		Optional additional arguments
	}
}

\value{
	An object of class "decompCTRE".
}

\references{
	Leps, Jan, Francesco de Bello, Petr Smilauer and Jiri Dolezal. 2011. Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects. Ecography 34 (5): 856-863.
}

\author{
	Adrien Taudiere
	Jan Leps
}

\seealso{
	\code{\link{barplot.decompCTRE}};
	\code{\link{traitflex.anova}}
}

\examples{
data(finch.ind)

res.decomp <- decompCTRE(traits = traits.finch, sp = sp.finch, 
ind.plot = ind.plot.finch, print = FALSE)

barplot.decompCTRE(res.decomp)

par(mfrow = c(2,2))
barplot.decompCTRE(res.decomp, resume = FALSE)
par(mfrow = c(1,1))
}



