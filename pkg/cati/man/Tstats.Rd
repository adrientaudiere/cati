\name{Tstats}
\alias{Tstats}
\alias{barplot.Tstats}
\alias{plot.Tstats}
\alias{summary.Tstats}
\alias{print.Tstats}
\alias{sum_Tstats}

\title{
	Computing observed T-statistics (T for Traits) and null expectations.
}
\description{
	Computing observed T-statistics (T for Traits) as three ratios of variance, namely T_IP.IC, T_IC.IR and T_PC.PR. This function can also return the distribution of this three statistics under null models.
}

\usage{
	Tstats(traits, ind.plot, sp, reg.pool = NULL, nperm = NULL, printprogress = TRUE)
	
	\method{barplot}{Tstats}(height, val.quant = c(0.025, 0.975), 
	col.Tstats = c("red", "purple", "green", "white"), ylim = NULL, \dots)
	
	\method{plot}{Tstats}(x, val.quant = c(0.025, 0.975), 
	col.Tstats = c("red", "purple", "green"), type = "normal", 
	add.conf = TRUE, ylim = NULL, xlim = NULL, \dots)
		
	\method{print}{Tstats}(x, \dots)
		
	\method{summary}{Tstats}(object, \dots)
	
	sum_Tstats(x, val.quant = c(0.025, 0.975), type = "all")
}

\arguments{
	\item{traits}{
		Individual Matrix of traits with traits in columns. For one trait, use as.matrix().
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individual is.
	}
	
	\item{sp}{
		Factor defining the species which the individual belong to.
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
	
	\item{x}{
		An object of class Tstats.
	}
	
	\item{height}{
		An object of class Tstats.
	}
	
	\item{object}{
		An object of class Tstats.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calculation confidence interval. By default val.quant = c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{col.Tstats}{
		A vector of three color correspond to the three T-statistics.
	}
	
	\item{type}{
		For the plot function, type of plot. Either "color_cond", "simple", "simple_sd", "normal" and "barplot". 
		For the summary function, type of summary statistics. Either "binary", "percent", "p.value", "site" or "all". 
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{ylim}{
		limits for the y axis.
	}
	
	\item{xlim}{
		limits for the x axis.
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}

\details{
	S3 method plot:
	-Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	-Color_cond type plot T-statistics for each site and traits with color for significant values and grey for non signifivative ones. 
	-Simple_sd type plot means, standard deviations and confidence intervals of T-statistics
	-Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	-Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
	
	S3 method print: print the structure if the object of class Tstats
	
	S3 method summary: print the summary statistics of the three T-statistics
	
	Method summary sum_Tstats:
	-Binary type only test if a T-statistics is significatively different from the null expectation for each trait.
	-Percent type determine the percentage of sites were the T-statistics is significatively different from the null expectation for each trait. Asterix shows global significance of the test.
	-P-value type determine the p-value (two unilateral tests) of the T-statistics for each trait and sites.
	-Site type allows to know in which sites T-statistics deviate from the null expectation.
	-All type do all the precedent type of summary.
}

\value{
	A list of statistics:
	
		
	\item{$tstats$Tstats$T_IP.IC}{Observed ratio between variance of individuals in populations and individuals in communities}
	\item{$tstats$Tstats$T_IC.IR}{Observed ratio between variance of individuals in communities and individuals in the region}
	\item{$tstats$Tstats$T_PC.PR}{Observed ratio between variance of populations in communities and populations in the region}
	
	\item{$tstats$T_IP.IC_nm}{If nperm is numeric; Result of simulation for T_IP.IC}
	\item{$tstats$T_IC.IR_nm}{If nperm is numeric; Result of simulation for T_IC.IR}
	\item{$tstats$T_PC.PR_nm}{If nperm is numeric; Result of simulation for T_PC.PR}
	
	\item{$variances$var_IP}{variance of individuals within populations}    
	\item{$variances$var_PC}{variance of populations within communities}
	\item{$variances$var_CR}{variance of communities within the region}
	\item{$variances$var_IC}{variance of individuals within communities}
	\item{$variances$var_PR}{variance of populations within the region}
	\item{$variances$var_IR}{variance of individuals within the region}
	    
	\item{$variances$var_IP_nm1}{variance of individuals within populations in null model 1}    
	\item{$variances$var_PC_nm2sp}{variance of populations within communities in null model 2sp}
	\item{$variances$var_IC_nm1}{variance of communities within the region in null model 1}
	\item{$variances$var_IC_nm2}{variance of individuals within communities in null model 2}
	\item{$variances$var_PR_nm2sp}{variance of populations within the region in null model 2sp}
	\item{$variances$var_IR_nm2}{variance of individuals within the region in null model 2}
}

\references{
	Violle, Cyrille, Brian J. Enquist, Brian J. McGill, Lin Jiang, Cecile H. Albert, Catherine Hulshof, Vincent Jung, et Julie Messier. 2012. The return of the variance: intraspecific variability in community ecology. Trends in Ecology & Evolution 27 (4): 244-252. doi:10.1016/j.tree.2011.11.014.
}

\author{
	Adrien Taudiere
	Cyrille Violle
}

\seealso{
	\code{\link{ComIndex}};
	\code{\link{ComIndexMulti}};
	\code{\link{plotCorTstats}};
	\code{\link{plotSESvar}};
	\code{\link{plot.listofindex}}
}

\examples{
	data(finch.ind)

	res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
	sp = sp.finch, nperm = 9, print = FALSE)
	
	attributes(res.finch)

	#Tstats class is associated to S3 methods plot, barplot and summary
	
	plot(res.finch)
	
	plot(res.finch, type = "color_cond")
	plot(res.finch, type = "simple")
	plot(res.finch, type = "simple_sd")
	plot(res.finch, type = "barplot")
	
	attributes(sum_Tstats(res.finch))
	head(sum_Tstats(res.finch)$p.value, 10)
	
	sum_Tstats(res.finch, type = "binary")
	sum_Tstats(res.finch, type = "percent")
	sum_Tstats(res.finch, type = "site")
	sum_Tstats(res.finch, type = "p.value")
	sum_Tstats(res.finch, type = "all")
	
	barplot(res.finch)
	
	attributes(sum_Tstats(res.finch))
	head(sum_Tstats(res.finch)$p.value, 10)
	
	#### An other way to see "ses values" of T-statistics
	
	# Custom theme (from rasterVis package)
	require(rasterVis)
	
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col = my.theme$regions$col)
	
	levelplot(t(ses(res.finch$Tstats$T_IP.IC,res.finch$Tstats$T_IP.IC_nm)$ses), 
	colorkey = my.ckey, par.settings = my.theme,border = "black")
	
	
	#### Use a different regional pool than the binding of studied communities
	
	#create a random regional pool for the example

	reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])

	res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
    sp = sp.finch, nperm = 9, print = FALSE)	
}





