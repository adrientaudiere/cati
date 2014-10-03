\name{partvar}
\alias{partvar}
\alias{barPartvar}
\alias{piePartvar}

\title{
	Variance partitioning accross nested scales
}

\description{
	Variance partitioning accross nested scales using a decomposition (varcomp function) of variance on restricted maximum likelihood (REML) method (lme function). See Messier et al. 2010 for more information. barPartvar and piePartvar are associated plotting functions.
}

\usage{
	partvar(traits, factors, printprogress = TRUE)
	barPartvar(partvar, col.bar = NA, ...)
	piePartvar(partvar, col.pie = NA, ...)
}

\arguments{
	\item{traits}{
		Matrix of traits with traits in column
	}
  
	\item{factors}{
		A matrix of factors with the first column corresponds to the higher level factor, the second row the second higher level factor and so on.  
	}
	
	\item{printprogress}{
		Logical value; print progress during the calculation or not.
	}
	
	\item{partvar}{
		The result of the partvar function.
	}
	
	\item{col.bar}{
		Vector of colors of bars 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the pie function.
	}
	
	\item{col.pie}{
		Vector of color for pie.
	}	
}

\value{
	An object of class "partvar" corresponding to a matrix of variance values with traits in rows and nested factors in column.
}

\references{
	Messier, Julie, Brian J. McGill, et Martin J. Lechowicz. 2010. How do traits vary across ecological scales? A case for trait-based ecology: How do traits vary across ecological scales? Ecology Letters 13(7): 838-848. doi:10.1111/j.1461-0248.2010.01476.x.
}

\author{
	Adrien Taudiere
	Julie Messier
}

\seealso{
	\code{\link{piePartvar}};
	\code{\link{barPartvar}}
}

\examples{

	data(finch.ind)
	
	cond<-seq(1,length(sp.finch)*2, by = 2)
	genus <- as.vector(unlist(strsplit(as.vector(sp.finch),"_"))[cond])

	res.partvar.finch <- partvar(traits = traits.finch, 
	factors = cbind(sites = as.factor(as.vector(ind.plot.finch)), 
	species = as.factor(as.vector(sp.finch)), genus = as.factor(genus)))
	
	res.partvar.finch
	
	oldpar<-par()
	par(mfrow = c(2,2), mai = c(0.2,0.2,0.2,0.2))
	piePartvar(res.partvar.finch)
	par(oldpar)
	
	barPartvar(res.partvar.finch)

}





