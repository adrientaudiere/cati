\name{plot.listofindex}
\alias{plot.listofindex}

\title{
	Plot community assembly index
}

\description{
	Plot community assembly index and confidence intervals using a list of index. S3 method for class listofindex.
}

\usage{
	\method{plot}{listofindex}(x, type = "normal", 
	col.index = c("red", "purple", "olivedrab3"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), 
	grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, 
	cex.text = 0.8, plot.ask = FALSE, srt.text = 90, alpha = 0.4, ...)
}

\arguments{

	\item{x}{
		A list of index and related null models obtained from to the as.listofindex function.
	}
	
	\item{type}{
		Type of plot. Possible type = "simple", "simple_range", "normal", "barplot" and "bytraits".
	}
	
	\item{col.index}{
		Vector of colors for index.
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{color.cond}{
		Logical value; If color.cond = TRUE, color points indicate T-statistics values significatively different from the null model and grey points are not different from null model.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculate confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
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
		Logical value; ask for plotting the next plot or not.
	}
	
	\item{srt.text}{
		Degree of rotation for text.
	}
	
	\item{alpha}{
		Degree of transparency for null models aera.
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}


\details{
	S3 method plot for class listofindex:
	-Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	-Simple_range type plot means, standard deviations and range of T-statistics
	-Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	-Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
	-Bysites type plot each metrics for each sites
	-Bytraits type plot each metrics for each traits
}

\value{
	None; used for the side-effect of producing a plot.
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{as.listofindex}};
	\code{\link{plot.Tstats}};
	\code{\link{ses.listofindex}}
}

\examples{	
	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)

	\dontrun{
		#### Use a different regional pool than the binding of studied communities
		#create a random regional pool for the example
	
		reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
	
		res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	    sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)	
	    
	    plot(as.listofindex(list(res.finch,res.finch2)))
    
    
	    #### Use a different regional pool for each communities
		#create a random regional pool for each communities for the example
		list.reg.p <- list(
		traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
		traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
		traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )

		# Warning: the regional pool need to be larger than the observed communities
		table(ind.plot.finch)
		# For exemple, the third community need a regional pool of more than 981 individuals
		
		res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	    sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)	
	    
	    plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))	
	}
}




