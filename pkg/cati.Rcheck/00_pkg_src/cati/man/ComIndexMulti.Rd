\name{ComIndexMulti}
\alias{ComIndexMulti}
\alias{plot.ComIndexMulti}
\alias{summary.ComIndexMulti}
\alias{print.ComIndexMulti}

\title{
	Computing multitraits metrics to test and quantify the non-random assembly of communities
}

\description{
	Computing multitraits metrics to test and quantify the non-random assembly of communities
}

\usage{
	
	ComIndexMulti(traits = NULL, index = NULL, by.factor = NULL, 
	nullmodels = NULL, ind.plot = NULL, sp = NULL, com = NULL,
	SE = 0, namesindex = NULL, reg.pool = NULL, SE.reg.pool = NULL,
	nperm = 99, printprogress = TRUE, independantTraits = TRUE, 
	type.sp.val = "count")
	
	\method{plot}{ComIndexMulti}(x, type = "normal", 
	col.index = c("red", "purple", "olivedrab3"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), ...)
		
	\method{print}{ComIndexMulti}(x, \dots)
		
	\method{summary}{ComIndexMulti}(object, \dots)

}

\arguments{
 
	\item{traits}{
		Individual Matrix of traits with traits in column (or species matrix when using "com" instead of "ind.plot").
	}
	
	\item{index}{
		A vector of functions to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)".
	}
	
	\item{by.factor}{
		A factor to split the Matrix of traits and compute index for each subset eg for each site.  
	}
	
	\item{nullmodels}{
		A vector of names corresponding to null models tu use for each index. \strong{local} (or 1) corresponds to a randomization of individual values within a given community. \strong{regional.ind} (or 2) corresponds to randomization of individual values within region, ie within all the dataset. A value of \strong{regional.pop} (or 2sp) corresponds to randomization of population values (each individual value are replaced by the mean value of it population) within region. Finally a value of \strong{regional.pop.prab} (or 2sp.prab)  mirror null model \strong{regional.pop} but without taking indo account species abundance. For example, if nullmodels = c("local", "regional.ind"), the first index will be calculated on the null model \strong{local}  and the second index on the null model \strong{regional.ind}.
		
		If only one value is given, all the the null model will be determined by this value. 
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot (site or community) in which the individual is. 
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
	}
	
	\item{com}{
		Community data matrix with species (or populations) in rows and sites in column. Use only if ind.plot = NULL. 
		"traits" matrix and "com" matrix must have the same number of rows.
	}
	
	\item{SE}{
		A single value or vector of standard errors associated with each traits. Especially allow to handle measurement errors. Not used with populational null model.  
	}
	
	\item{namesindex}{
		A vector of names for metrics.
	}
	
	\item{reg.pool}{
		Regional pool data for traits. If not informed, traits is considere as the regional pool. This matrix need to be larger (more rows) than the matrix "traits". Use only for null model 2.
	}	
	
	\item{SE.reg.pool}{
		A single value or vector of standard errors associated with each traits in each regional pool. Use only if reg.pool is used. Need to have the same dimension as reg.pool.
	}
	
	\item{nperm}{
		Number of permutations. If NULL, only observed values are returned.
	}
	
	\item{printprogress}{
		Logical value; print progress during the calculation or not.
	}
	
	\item{independantTraits}{ 
		Logical value (default: TRUE). If independantTraits is true (default), each traits is sample independently in null models, if not, each lines of the matrix are randomized, keeping the relation (and trade-off) among traits.
	}
	
	\item{type.sp.val}{
		Only if ind.plot = NULL. Either "count" or "abundance". Use abundance when all values in the com matrix are not superior to one.
	}
	
	\item{x}{
		An object of class ComIndexMulti.
	}

	\item{object}{
		An object of class ComIndexMulti.
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
	
	\item{\dots}{
		Any additional arguments are passed to the plot, print or summary function creating the core of the plot and can be used to adjust the look of resulting graph. See \code{\link{plot.listofindex}} for more arguments.
	}
}

\details{
	This function implement four null models which keep unchanged the number of individual per community.
	Model \strong{local} (1) corresponds to randomization of individual values within community.
	Model \strong{regional.ind} (2) corresponds to randomization of individual values within region.
	Model \strong{regional.pop} (2sp) corresponds to randomization of population values within region.
	Model \strong{regional.pop.prab} (2sp.prab) corresponds to randomization of population values within region but whitout taking into account for abundance.
	
	S3 method plot for class listofindex:
	
	-Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	
	-Simple_range type plot means, standard deviations and range of T-statistics
	
	-Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	
	-Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
	
	-Bysites type plot each metrics for each sites
	
	-Bytraits type plot each metrics for each traits
}

\value{
	A list of lists:
	
	\item{$obs}{List of observed values for each trait in each community. Each component of the list correspond to a matrix containing the result for each custom function.}
	
	\item{$null}{List of null values for each trait in each community. Each component of the list correspond to an array containing the result of the permutations for each custom function.}
	
	\item{$sites_richness}{Number of species per site.}
	\item{$namestraits}{Names of traits.}
	
	\item{$traits}{traits data}
 	\item{$ind.plot}{name of the plot in which the individual is}
 	\item{$sp}{groups (e.g. species) which the individual belong to}
 	
 	\item{$nullmodels}{List of null models used for each indices.}
 	 	
 	\item{$call}{call of the function Tstats}
 	
 	\item{$list.index}{List of index values and associate null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and associate null models. Internal use in other function. Traits in rows.}
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{ComIndex}};
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{
data(finch.ind)

\dontrun{
#For most multivariate functions we need to replace (or exclude) NA values. 
#For this example, we use the package mice to complete the data.

comm<-t(table(ind.plot.finch,1:length(ind.plot.finch)))

library(mice)
traits = traits.finch
mice<-mice(traits.finch)
traits.finch.mice<-complete(mice)


#A simple example to illustrate the concept of the function ComIndexMulti

n_sp_plot<-as.factor(paste(sp.finch, ind.plot.finch, sep = "_")) 
res.sum.1<-ComIndexMulti(traits.finch, 
              index = c("sum(scale(x), na.rm = T)", "sum(x, na.rm = T)"), 
              by.factor = n_sp_plot, nullmodels = "regional.ind", 
              ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)
res.sum.1



#A more interesting example using the function hypervolume
library(hypervolume)

hv<-hypervolume(traits.finch.mice, 
        reps = 100,bandwidth = 0.2, 
        verbose = F, warnings = F)
plot(hv)

hv.1<-ComIndexMulti(traits.finch.mice, 
             index = c("as.numeric(try(hypervolume(na.omit(x), reps = 100, 
                 bandwidth = 0.2, verbose = F, warnings = F)@Volume))"),
             by.factor = rep(1,length(n_sp_plot)), nullmodels = "regional.ind",
             ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch) 

hv.1

}
}
