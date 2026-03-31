\name{plotSpPop}
\alias{plotSpPop}

\title{
	Plot populations values against species values
}

\description{
	Plot populations values against species values. The objectif is to see the contribution of intra-specific vs inter-specific variation to trait gradient.
}

\usage{
	plotSpPop(traits = NULL, ind.plot = NULL, sp = NULL, 
	col.ind = rgb(0.5, 0.5, 0.5, 0.5), col.pop = NULL, col.sp = NULL, 
	col.site = NULL, resume =  FALSE, p.val = 0.05, min.ind.signif = 10, 
	multipanel = TRUE, col.nonsignif.lm = rgb(0, 0, 0, 0.5), 
	col.signif.lm = rgb(1, 0.1, 0.1, 0.8), silent =  FALSE)
}

\arguments{
	\item{traits}{
		Individual Matrix of traits with traits in columns.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individual is.
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
	}
	
	\item{col.ind}{
		Color for individual values.
	}
	
	\item{col.pop}{
		Color for populational mean values.
	}
	
	\item{col.sp}{
		Color for species mean values.
	}
	
	\item{col.site}{
		Color for sites mean values.
	}
	
	\item{resume}{
		Logical, if TRUE plot a simple form of the plot.
	}
	
	\item{p.val}{
		Choosen p.value to print significant linear relationship using linear model. Argument past to the lm funtion internally.
	}
	
	\item{min.ind.signif}{
		Minimum individual to print significant linear relationship.
	}
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
	}
	
	\item{col.nonsignif.lm}{
		Color for non significant linear relationship.
	}
	
	\item{col.signif.lm}{
		Color for significant linear relationship.
	}
	
	\item{silent}{
		Logical value, if resume = FALSE do not print warning argument.
	}
}

\details{
	Example of utilisation: Cornwell, W.K., Ackerly, D.D., 2009. Community assembly and shifts in plant trait distributions across an environmental gradient in coastal California. Ecological Monographs 79, 109-126.
}

\value{
	None; used for the side-effect of producing a plot.
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{plotDistri}}
}

\examples{
	data(finch.ind)
	\donttest{
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE)

	#If we change the value of the threshold 
	#(alpha = 10\% instead of 5\% 
	#and the minimum individual to represent significativity 
	#fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpPop(traits.finch, ind.plot.finch, sp.finch, p.val = 0.1,  
	min.ind.signif = 3, silent = TRUE)


	#For a more simple figure, add the option resume = TRUE. 
	#Again if we change the value of the threshold 
	#(alpha = 10\% instead of 5\% 
	#and the minimum individual to represent significativity
	# fixed to 3 instead of 10 by default) 
	#we can see some significant relationships.

	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey")
	
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey", col.sp = "black")
	
	plotSpPop(traits.finch, ind.plot.finch, sp.finch, silent = TRUE, 
	resume = TRUE, col.pop = "grey", col.sp = "black", 
	p.val = 0.1,  min.ind.signif = 3)
	}
}




