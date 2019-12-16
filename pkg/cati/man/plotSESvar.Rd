\name{plotSESvar}
\alias{plotSESvar}

\title{
	Plot SES values against a variable
}

\description{
	Plot standardized effect size values against a variable
}

\usage{
	plotSESvar(index.list, variable = NULL, ylab = "variable", 
	color.traits = NULL, val.quant = c(0.025, 0.975), resume = FALSE, 
	multipanel = TRUE)
}

\arguments{
	\item{index.list}{
		A list of index and the associate null models in the forme: list( index_1 = index_1_observed, index_1_nm = null.model.index_1 ,index_2 = index_2_observed, index_2_nm = null.model.index_2, ...). 
	}
	
	\item{variable}{
		The variable against standardized effect sizes are plotted.
	}
	
	\item{ylab}{
		Label for the variable.
	}
	
	\item{color.traits}{
		A vector of colors corresponding to traits.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculation confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{resume}{
		Logical value; resume = FALSE by default; Simplify the plot by plotting the mean and standard error for index value of multiple traits
	}
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
	}
}

\value{
	None; used for the side-effect of producing a plot.
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{
	data(finch.ind)
	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, sp = sp.finch, 
	nperm = 9)
\dontrun{
	par(mfrow = c(2,2))
	species.richness <- table(ind.plot.finch)
	plotSESvar(as.listofindex(list(res.finch)), species.richness, 
	multipanel = FALSE)

	#Same plot with resume = TRUE.
	
	par(mfrow = c(2,2))
	plotSESvar(as.listofindex(list(res.finch)), species.richness, 
	resume = TRUE, multipanel = FALSE)
	par(mfrow = c(1,1))
	}
}




