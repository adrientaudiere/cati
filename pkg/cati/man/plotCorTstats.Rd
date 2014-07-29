\name{plotCorTstats}
\alias{plotCorTstats}

\title{
	Plot the bivariate relationships between T-statistics
}

\description{
	Plot the bivariate relationships between the three T-statistics namely T_IP.IC, T_IC.IR and T_PC.PR.
}

\usage{
	plotCorTstats(tstats = NULL, val.quant = c(0.025, 0.975), 
	add.text =  FALSE, bysite =  FALSE, col.obj = NULL, plot.ask = TRUE, 
	multipanel = TRUE, ...)
}

\arguments{
	\item{tstats}{
		The list resulting from the function Tstats.
	}
	
	\item{val.quant}{		
		Numeric vector of length 2, giving the quantile to calculate confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{add.text}{
		Logical value; Add text or not.
	}
	
	\item{bysite}{
		Logical value; plot per site or by traits. 
	}
	
	\item{col.obj}{
		Vector of colors for object (either traits or sites).
	}
	
	\item{plot.ask}{
		Logical value; Ask for new plot or not.
	}
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}

\value{
	None; used for the side-effect of producing a plot.
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{Tstats}};
	\code{\link{plot.Tstats}};
	\code{\link{plotSESvar}}
}

\examples{
	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9)
	
	plotCorTstats(res.finch, bysite = FALSE)
	plotCorTstats(res.finch, bysite = TRUE)
}






