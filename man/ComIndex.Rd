\name{ComIndex}
\alias{ComIndex}
\alias{plot.ComIndex}
\alias{summary.ComIndex}
\alias{print.ComIndex}

\title{
	Computing metrics to test and quantify the non-random assembly of communities
}

\description{
	Computing the moments of the trait distribution and other metrics to test and quantify the non-random assembly of communities. 
}

\usage{

	ComIndex(traits = NULL, index = NULL, nullmodels = NULL, 
	ind.plot = NULL, sp = NULL, com = NULL, SE = 0, namesindex = NULL,
	reg.pool = NULL, SE.reg.pool = NULL, nperm = 99, printprogress = TRUE, 
	independantTraits = TRUE, type.sp.val = "count")
	
	\method{plot}{ComIndex}(x, type = "normal", 
	col.index = c("red", "purple", "olivedrab3"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), ...)
		
	\method{print}{ComIndex}(x, \dots)
		
	\method{summary}{ComIndex}(object, \dots)

}

\arguments{
	
	\item{traits}{
		Individual Matrix of traits with traits in column (or species matrix when using "com" instead of "ind.plot"). 
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)", see examples for more complexe functions.
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
		A vector of names for index.
	}
		
	\item{reg.pool}{
		Regional pool data for traits. If not informed, traits is considere as the regional pool. This matrix need to be larger (more rows) than the matrix "traits". Use only for null model \strong{regional.ind}.
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
		Only if ind.plot = NULL. Either "count" or "abundance". Use abundance when one value or more in the com matrix are inferior to one. Using abundance is EXPERIMENTAL. This function round abundance to fit count data.
	}
	
	\item{x}{
		An object of class ComIndex.
	}

	\item{object}{
		An object of class ComIndex.
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
		Any additional arguments are passed to the plot, print or summary function. See \code{\link{plot.listofindex}} for more arguments.
	}	
}

\details{
	Compute statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models.
	For each statistic this function returns observed values and the related null distribution.
	This function implement four null models which keep unchanged the number of individual per community.
	Model \strong{local} (1) corresponds to randomization of individual values within community.
	Model \strong{regional.ind} (2) corresponds to randomization of individual values within region.
	Model \strong{regional.pop} (2sp) corresponds to randomization of population values within region.
	Model \strong{regional.pop.prab} (2sp.prab) corresponds to randomization of population values within region but whitout taking into account for abundance.
	
	In most cases, models \strong{local} and \strong{regional.ind} correspond to index at the individual level and the model \strong{regional.pop} and \strong{regional.pop.prab} to index at the species level (or any other aggregate variable like genus, family or functionnal group).
	
	S3 method plot for class listofindex:
	
	-Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	
	-Simple_range type plot means, standard deviations and range of T-statistics
	
	-Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	
	-Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
	
	-Bysites type plot each metrics for each sites
	
	-Bytraits type plot each metrics for each traits
}

\value{
	An object of class "ComIndex" corresponding to a list of lists:
	
	\item{$obs}{List of observed values for each trait in each community. Each component of the list corresponds to a matrix containing the result for each custom function.}
	
	\item{$null}{List of null values for each trait in each community. Each component of the list corresponds to an array containing the result of the permutations for each custom function.}
	
	\item{$list.index}{List of index values and related null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and related null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species per site.}
	\item{$namestraits}{Names of traits.}
	
	\item{$traits}{traits data}
 	\item{$ind.plot}{name of the plot in which the individual is}
 	\item{$sp}{groups (e.g. species) which the individual belong to}
 	
 	\item{$nullmodels}{List of null models used for each indices.}
 		
 	\item{$call}{call of the function Tstats}
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{ComIndexMulti}};
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{	
data(finch.ind)

\donttest{
#Define the functions that will be calculating
funct<-c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)",
     "max(x, na.rm = TRUE) - min(x, na.rm = TRUE)" )

#Test against the null model regional.ind
res.finch.sp_regional.ind<-ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
                           nullmodels = "regional.ind", ind.plot = ind.plot.finch,
                            nperm = 9, print = FALSE)
 
#Test against the null model regional.pop
#Individuals values are transformed in populational values
res.finch.sp_regional.pop<-ComIndex(traits = traits.finch, index = funct, sp = sp.finch,
               nullmodels = "regional.pop", ind.plot = ind.plot.finch, 
               nperm = 9, print = FALSE)


#We can calculate index with or without intraspecific variance.

#calculate  of means by population (name_sp_site is a name of a population) 
#determine the site for each population (sites_bypop)
 
name_sp_sites = paste(sp.finch, ind.plot.finch,sep = "_")
traits.by.pop<-apply(traits.finch, 2 , 
           function (x) tapply(x, name_sp_sites, mean, na.rm = TRUE))

sites_bypop<-lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), split = "_"), 
          function(x) x[3])

#New list of function "funct"

funct.1<-c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))",
     "tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))",
     "tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE)-min(x, na.rm = TRUE))", 
     "tapply(x, ind.plot.finch, function(x) CVNND(x, na.rm = TRUE))" )

fact<-unlist(sites_bypop)  
funct.2<-c("tapply(x, fact, function(x) mean(x, na.rm = TRUE))",
          "tapply(x, fact, function(x) kurtosis(x, na.rm = TRUE))",
          "tapply(x, fact, function(x) max(x, na.rm = TRUE)-min(x, na.rm = TRUE))", 
          "tapply(x, fact, function(x) CVNND(x, na.rm = TRUE))")


res.finch.withIV<-ComIndex(traits = traits.finch, index = funct.1,
               sp = sp.finch, nullmodels = "regional.ind",
               ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

res.finch.withoutIV<-ComIndex(traits = traits.finch, index = funct.2, 
               sp = sp.finch, nullmodels = "regional.pop", 
               ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

#ComIndex class are associated to S3 methods plot, print and summary.

res.finch.withIV
summary(res.finch.withIV)
plot(res.finch.withIV)
plot(res.finch.withoutIV)

plot(as.listofindex(list(res.finch.withIV, res.finch.withoutIV)))
}

}
