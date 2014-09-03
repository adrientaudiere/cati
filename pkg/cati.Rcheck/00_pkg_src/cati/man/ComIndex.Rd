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
	ind.plot = NULL, sp = NULL, com = NULL, namesindex = NULL,
	reg.pool = NULL, nperm = 99, printprogress = TRUE, 
	ind.value = TRUE, type.sp.val = "count")
	
	\method{plot}{ComIndex}(x, type = "normal", 
	col.index = c("red", "purple", "green"), add.conf = TRUE, 
	color.cond = TRUE, val.quant = c(0.025, 0.975), 
	grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, 
	cex.text = 0.8, plot.ask = FALSE, srt.text = 90, 
	bysite = FALSE, ...)
		
	\method{print}{ComIndex}(x, \dots)
		
	\method{summary}{ComIndex}(object, \dots)

}

\arguments{
	
	\item{traits}{
		Individual Matrix of traits with traits in column.
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)", see examples for more complexe functions.
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
		Community data matrix with species (or populations) in rows and sites in column. Use only if ind.value = FALSE. 
		"traits" matrix and "com" matrix must have the same number of rows.
	}
	
	\item{namesindex}{
		A vector of names for index.
	}
		
	\item{reg.pool}{
		Regional pool data for traits. If not informed, traits is considere as the regional pool. This matrix need to be larger (more rows) than the matrix "traits". Use only for null model 2.
	}	
	
	\item{nperm}{
		Number of permutations. If NULL, only observed values are returned; 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calculation or not.
	}
	
	\item{ind.value}{
		Do the data are from individual value. If not, an internal function transform the abundance data matrix of com into a individual like matrix to compute the function. 
	}
	
	\item{type.sp.val}{
		Only if ind.value = FALSE. Either "count" or "abundance". Use abundance when one value or more in the com matrix are inferior to one. Using abundance is EXPERIMENTAL. This function round abundance to fit count data.
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
	
	\item{bysite}{
		Logical value; plot per site or by traits. 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot, print or summary function.
	}
	
}

\details{
	Compute statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models.
	For each statistic this function returns observed values and the related null distribution.
	This function implement three null models which keep unchanged the number of individual per community.
	Model 1 corresponds to randomization of individual values within community.
	Model 2 corresponds to randomization of individual values within region.
	Model 2sp corresponds to randomization of population values within region.
	Model 2sp.prab corresponds to randomization of population values within region but whitout taking into account for abundance.
	
	In most cases, models 1 and 2 correspond to index at the individual level and the model 2sp and 2sp.prab to index at the species level (or any other aggregate variable like genus, family or functionnal group).
}

\value{
	An object of class "ComIndex" corresponding to a list of lists:
	
	\item{$obs}{List of observed values for each trait in each community. Each component of the list corresponds to a matrix containing the result for each custom function.}
	
	\item{$Null}{List of null values for each trait in each community. Each component of the list corresponds to an array containing the result of the permutations for each custom function.}
	
	\item{$list.index}{List of index values and related null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and related null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species per site.}
	\item{$namestraits}{Names of traits.}
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
oldpar <- par(no.readonly = TRUE)

####
#The function ComIndex allow to choose your own function 
#(like mean, range, variance...) to calculate customize index.

require(e1071)

funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )

\dontrun{

res.finch.sp_mn2 <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = c("2","2","2","2"), 
ind.plot = ind.plot.finch,	nperm = 9, print = FALSE)

res.finch.sp_mn2sp <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), 
ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

####
#We can represent Standardized Effect Size (ses)
#using the function plot(as.listofindex(list1, list2, list3))

list.ind2 <- list(res.finch.sp_mn2, res.finch.sp_mn2sp)
index.list2 <- as.listofindex(list.ind2)

plot(index.list2)
plot(index.list2, type = "bytraits")	

####
#This allows to calcul index per site 
#for example using "tapply(x, sites, mean)".

funct <- c("tapply(x, ind.plot.finch, function(x) mean(x, na.rm = TRUE))", 
"tapply(x, ind.plot.finch, function(x) kurtosis(x, na.rm = TRUE))", 
"tapply(x, ind.plot.finch, function(x) max(x, na.rm = TRUE) - 
min(x, na.rm = TRUE) )", "tapply(x, ind.plot.finch, function(x) 
CVNND(x))"  )


##Null model 1 is trivial for this function 
#because randomisation is within community only

res.finch.ind_mn1 <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = c(1,1,1,1), ind.plot = ind.plot.finch, 
nperm = 9, print = FALSE)

res.finch.ind_mn2 <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
nperm = 9, print = FALSE)

	
####
#We can calcul metrics with or without intraspecific variance.
#Calculation of trait averages per population 
#(name_sp_site is a name of a population) 
#like in the function ComIndex
#and determine the site for each population (sites_bypop)

name_sp_sites = paste(sp.finch, ind.plot.finch, sep = "_")
traits.by.pop <- apply(traits.finch, 2 , function (x) 
tapply(x, name_sp_sites, mean , na.rm = TRUE))

sites_bypop <- lapply(strsplit(paste(rownames(traits.by.pop), sep = "_"), 
split = "_"), function(x) x[3])

funct.withoutIV <- c("tapply(x, unlist(sites_bypop), function(x) 
mean(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), function(x) 
kurtosis(x, na.rm = TRUE))", "tapply(x, unlist(sites_bypop), function(x) 
max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
"tapply(x, unlist(sites_bypop), function(x) CVNND(x))"  )

funct.withIV <- c("tapply(x, ind.plot.finch, function(x) 
mean(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
kurtosis(x, na.rm = TRUE))", "tapply(x, ind.plot.finch, function(x) 
max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
"tapply(x, ind.plot.finch, function(x) CVNND(x))"  )
	
res.finch.withIV <- ComIndex(traits = traits.finch, index = funct.withIV, 
sp = sp.finch, nullmodels = c("2","2","2","2"), ind.plot = ind.plot.finch, 
nperm = 9, print = FALSE)

res.finch.withoutIV <- ComIndex(traits = traits.finch, index = funct.withoutIV, 
sp = sp.finch, nullmodels = c("2sp","2sp","2sp","2sp"), 
ind.plot = ind.plot.finch, nperm = 9, print = FALSE)

	
####
#We can also represent T-statistics and custom index thanks to
#the plot.listofindex function.

res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, sp = sp.finch, 
nperm = 9, print = FALSE)

list.ind <- list(res.finch.withIV, res.finch.withoutIV ,res.finch)

index.list1 <- as.listofindex(list.ind, namesindex = c("mean", "kurtosis", 
"range", "CVNND", "mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
"T_IP.IC", "T_IC.IR", "T_PC.PR"))

class(index.list1)

par(mfrow = c(2,3))
plot(index.list1,type = "bytraits", bysite = TRUE)

par(mfrow = c(2,2))
plot(index.list1,type = "bytraits")
par(mfrow = c(1,1))

plot(index.list1,type = "simple")
plot(index.list1,type = "simple_range")
plot(index.list1,type = "normal")
plot(index.list1,type = "barplot")
}


############################
####Using ind.value = FALSE and community data matrix if there is no data 
#available at the individual level.

\dontrun{

#create traits data at the species level
traits_by_sp <- apply(traits.finch,2,function(x) tapply(x,sp.finch,
function(x) mean(x, na.rm = T)))  

#create traits data at the populational level
names_sp_ind.plot <- as.factor(paste(sp.finch, ind.plot.finch, sep = "@")) 
traits_by_pop <- apply(traits.finch,2,function(x) tapply(x,names_sp_ind.plot, 
function(x) mean(x, na.rm = T) ))  

#create community data matrix at the species or populational level
w1 <- table(sp.finch,ind.plot.finch)
dim(w1)
dim(traits_by_sp)

w2 <- table(names_sp_ind.plot,ind.plot.finch)
dim(w2)
dim(traits_by_pop)

#Choose indices
require(e1071)

funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)", 
"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)", "CVNND(x)"  )


#################
#with species value

res <- AbToInd(traits_by_sp, w1)

ComIndex(traits_by_sp, nullmodels = 2, ind.value = FALSE, index = funct, 
sp = rownames(traits_by_sp), com = w1, nperm = 9)


#################
#with population value
res <- AbToInd(traits_by_pop, w2)
sp.sp <- unlist(strsplit(rownames(traits_by_pop),"@"))[seq(1,39*2,2)]

ComIndex(traits_by_pop, nullmodels = 2, ind.value = FALSE, index = funct, 
sp = sp.sp, com = w2)
	
}

############################
####Simple example using null model 2sp.prab (species level without taking 
# into acount for species abundance, prab for presence/absence)

\dontrun{

traits_by_sp <- apply(traits.finch, 2, function(x) 
tapply(x, name_sp_sites, mean, na.rm=T))

sites_bysp<-unlist(strsplit(rownames(traits_by_sp), 
split="_"))[seq(3,3*dim(traits_by_sp)[1], by=3) ]

funct.withoutIV.prab <- c("tapply(x, unlist(sites_bysp), 
	function(x) mean(x, na.rm = TRUE))", 
	"tapply(x, unlist(sites_bysp), function(x) kurtosis(x, na.rm = TRUE))", 
	"tapply(x, unlist(sites_bysp), function(x) max(x, na.rm = TRUE) 
	- min(x, na.rm = TRUE) )", 
	"tapply(x, unlist(sites_bysp), function(x) CVNND(x))")
	
res.finch.withoutIV.prab <- ComIndex(traits = traits.finch, 
	index = funct.withoutIV.prab, sp = sp.finch, 
	nullmodels = rep("2sp.prab", times=4), ind.plot = ind.plot.finch, 

traits_by_sp <- apply(traits.finch, 2, function(x) tapply(x, 
name_sp_sites, mean, na.rm=T))
sites_bysp<-unlist(strsplit(rownames(traits_by_sp), split="_"))
[seq(3,3*dim(traits_by_sp)[1], by=3) ]

funct.withoutIV.prab <- c("tapply(x, unlist(sites_bysp), 
	function(x) mean(x, na.rm = TRUE))", 
	"tapply(x, unlist(sites_bysp), 
	function(x) kurtosis(x, na.rm = TRUE))", 
	"tapply(x, unlist(sites_bysp), 
	function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )", 
	"tapply(x, unlist(sites_bysp), 
	function(x) CVNND(x))")
	
res.finch.withoutIV.prab <- ComIndex(traits = traits.finch, 
	index = funct.withoutIV.prab, sp = sp.finch, 
	nullmodels = rep("2sp.prab", times=4), ind.plot = ind.plot.finch,
	nperm = 9, print = FALSE)

list.ind2 <- list(res.finch.withoutIV, res.finch.withoutIV.prab)
index.list2 <- as.listofindex(list.ind2, namesindex = 
	c("mean.pop", "kurtosis.pop", "range.pop", "CVNND.pop", 
	"mean.prab", "kurtosis.prab", "range.prab", "CVNND.prab"))

plot(index.list2)
}


}
