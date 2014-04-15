\name{plot.listofindex}
\alias{plot.listofindex}

\title{
	Plot community assembly index
}

\description{
	Plot community assembly index and confidence intervals using a list of index. S3 method fot class listofindex.
}

\usage{
	plot.listofindex(x, type = "normal", col.index = c("red", "purple", "green"), add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025, 0.975), grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, cex.text = 0.8, plot.ask = FALSE, srt.text = 90, bysites = FALSE, ...)
}

\arguments{

	\item{x}{
		A list of index and the associate null models obtain thanks to the as.listofindex function.
	}
	
	\item{type}{
		Type of plot. Possible type = "simple", "simple_range", "normal", "barplot" and a very different one: "bytraits".
	}
	
	\item{col.index}{
		Vector of colors for index.
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{color.cond}{
		Logical value; If color.cond=TRUE, color points indicate T-statistics values significatively different from the null model and grey points are not different from null model.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{grid.v}{
		Logical value; print vertical grid or not
	}
	
	\item{grid.h}{
		Logical value; print horizontal grid or not
	}
	
	\item{xlim}{
		Numeric vectors of length 2, giving the x coordinates range 
	}
	
	\item{ylim}{
		Numeric vectors of length 2, giving the y coordinates range 
	}
		
	\item{cex.text}{
		Numeric value; the magnification to be used for text relative to the current setting of cex
	}
	
	\item{plot.ask}{
		Logical value; ask for ploting the next plot or not.
	}
	
	\item{srt.text}{
		Degree of rotation for text.
	}
	
	\item{bysites}{
		Logical value; plot by sites or by traits. 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}

\details{
}

\value{
	None; used for the side-effect of producing a plot.
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}

\seealso{
	\code{\link{as.listofindex}};
	\code{\link{plot.Tstats}};
	\code{\link{ses.listofindex}}
}

\examples{	
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
