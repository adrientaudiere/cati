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
	namesindex = NULL, reg.pool = NULL, nperm = 99, 
	printprogress = TRUE, type.sp.val = "count")
	
	\method{plot}{ComIndexMulti}(x, type = "normal", 
	col.index = c("red", "purple", "olivedrab3"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), ...)
		
	\method{print}{ComIndexMulti}(x, \dots)
		
	\method{summary}{ComIndexMulti}(object, \dots)

}

\arguments{
 
	\item{traits}{
		Individual Matrix of traits with traits in column.
	}
	
	\item{index}{
		A vector of functions to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)".
	}
	
	\item{by.factor}{
		A factor to split the Matrix of traits and compute index for each subset eg for each site.  
	}
	
	\item{nullmodels}{
		A vector of values corresponding to null models tu use for each index. A value of 1 corresponds to a randomization of individual values within a given community. A value of 2 corresponds to randomization of individual values within region, ie within all the dataset. A value of 2sp corresponds to randomization of population values (each individual value are replaced by the mean value of it population) within region. Finally a value of 2sp.prab mirror null model 2 but without taking indo account species abundance. For example, if nullmodels = c(1,2), the first index will be calculated on the null model 1 and the second index on the null model 2.
		
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
	
	\item{namesindex}{
		A vector of names for metrics.
	}
	
	\item{reg.pool}{
		Regional pool data for traits. If not informed, traits is considere as the regional pool. This matrix need to be larger (more rows) than the matrix "traits". Use only for null model 2.
	}	
	
	\item{nperm}{
		Number of permutations. If NULL, only observed values are returned.
	}
	
	\item{printprogress}{
		Logical value; print progress during the calculation or not.
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

####
#For most multivariate functions we need to replace (or exclude)
#NA values.

#For this example, we use the package mice to complete the data.

\dontrun{
names.sp_ind.plot <- as.factor(paste(sp.finch, ind.plot.finch, sep = "_")) 

comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))

library(mice)
traits = traits.finch
mice <- mice(traits.finch)
traits.finch.mice <- complete(mice)

####
#A simple example to illustrate the concept of the function 
#ComIndexMulti

res.sum.1 <- ComIndexMulti(traits.finch, 
index = c("sum(scale(x), na.rm = TRUE)", "sum(x, na.rm = TRUE)"), 
by.factor = names.sp_ind.plot, nullmodels = c(2,2), 
ind.plot = ind.plot.finch, nperm = 50, sp = sp.finch)

attributes(ses.listofindex(as.listofindex(res.sum.1)))

####
#A more interesting example using the function hypervolume 
#from the package hypervolume. 
#We show here several results which differe in there factor 
#that delimit the group to calculate different hypervolume 
#(argument by_factor). 

require(hypervolume)

res.hv.1 <- ComIndexMulti(traits.finch.mice, index = 
paste("as.numeric (try(hypervolume(na.omit(x), warnings = FALSE,", 
"bandwidth=0.2, verbose=FALSE)@Volume))"), 
by.factor = rep(1,length(names.sp_ind.plot)), nullmodels = c(2,2), 
ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)

res.hv.2 <- ComIndexMulti(traits.finch.mice, index = 
paste("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE,",
"bandwidth=0.2, verbose=FALSE)@Volume))"), 
by.factor = names.sp_ind.plot, nullmodels = c(2,2), 
ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)

res.hv.3 <- ComIndexMulti(traits.finch.mice, index = 
paste("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE,"
"bandwidth=0.2, verbose=FALSE)@Volume))"), 
c("as.numeric (try(hypervolume(na.omit(x), 
warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
by.factor = rep(1,length(names.sp_ind.plot)), nullmodels = c(2,2), 
ind.plot = ind.plot.finch, nperm = 9, sp = sp.finch)

res.hv.4 <- ComIndexMulti(traits.finch.mice, index = 
paste("as.numeric(try(hypervolume(na.omit(x), warnings = FALSE,",
"bandwidth=0.2, verbose=FALSE)@Volume))"), 
c("as.numeric(try(hypervolume(na.omit(x), 
warnings = FALSE, bandwidth=0.2, verbose=FALSE)@Volume))"), 
by.factor = sp.finch, nullmodels = c(2,2), ind.plot = ind.plot.finch, 
nperm = 9, sp = sp.finch)


list.ind.multi <- as.listofindex(list(res.hv.2, res.hv.3, res.hv.4))

ses.listofindex(list.ind.multi)

plot(list.ind.multi)
plot(list.ind.multi, xlim = c(-200,20))
}
}

