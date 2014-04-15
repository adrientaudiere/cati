\name{plot_dens}
\alias{plot_dens}

\title{
	Plot function to represent density of traits values
}

\description{
	Plot function to represent density of traits values
}

\usage{
	plot_dens(traits=NULL, var.1=NULL, var.2=NULL, col.dens=NULL, plot.ask=TRUE, ylim.cex=1, cex.leg=0.8, polyg=TRUE, multipanel=TRUE, leg=TRUE)
}

\arguments{
	\item{traits}{
		Traits matrix with traits in column.
	}
	
	\item{var.1}{
		The first variable define the division of each plots, in most case either a vector of species or sites names.
	}
	
	\item{var.2}{
		The second variable define the division by color, in most case either a vector of species or sites names.
	}
	
	\item{col.dens}{
		A vector of colors for the second variable.
	}
	
	\item{plot.ask}{
		Logical value; ask for ploting the next plot or not.
	}
	
	\item{ylim.cex}{
		Numeric value; the magnification to be used for range of y axe
	}
	
	\item{cex.leg}{
		Numeric value; the magnification to be used for legend relative to the current setting of cex
	}
	
	\item{polyg}{
		Logical value; do the mean distribution is full or empty
	}
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
	}
	
	\item{leg}{
		Logical value; if TRUE print the legend.
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
	\code{\link{plot_sp_pop}}
}

\examples{	
	data(finch.ind)
	
	#Plot the distribution of traits values for populations, species, sites and regional scales. First, let try the distribution for all populations of Darwin finches.
	
	par(mfrow=c(4,4), cex=0.5)
	plot_dens(traits.finch, sp.finch, ind.plot.finch, ylim.cex=3, plot.ask=FALSE, multipanel=FALSE, leg=FALSE)
	par(mfrow=c(1,1), cex=1)
	
	#Then we can inverse the second and the third arguments to plot the distribution for all finches species. 
	
	par(mfrow=c(4,4), cex=0.5)
	plot_dens(traits.finch, ind.plot.finch, sp.finch, ylim.cex=8, plot.ask=FALSE, multipanel=FALSE, leg=FALSE)
	par(mfrow=c(1,1), cex=1)
	
	#You can also plot trait distribution for all species in the region
	plot_dens(traits.finch, rep("region", times=dim(traits.finch)[1]), sp.finch, ylim.cex=6, plot.ask=FALSE, leg=FALSE)
	
	
	#You can also plot trait distribution for all sites without regard on species attributes
	plot_dens(traits.finch, rep("toutes_sp", times=dim(traits.finch)[1]), ind.plot.finch, ylim.cex=3, plot.ask=FALSE)
	par(mfrow=c(4,4), cex=0.5)
}




\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
