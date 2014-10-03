\name{as.listofindex}
\alias{as.listofindex}
\title{
	Transform index results in a list of index
}

\description{
	Transform various results from functions Tstast, ComIndex or ComIndexMulti in a list of index. Useful to use the functions plot.listofindex (S3 method) and ses.listofindex.
}

\usage{
	as.listofindex(x, namesindex = NULL)
}

\arguments{
	\item{x}{
		A list of objects of class Tstast, ComIndex or ComIndexMulti
	}
	
	\item{namesindex}{
		Optionnal, the names of index in the same order as in x.
	}
}

\value{
	A list of observed values and corresponding "null" values (i.e. produced by null models) in the form "list(index1, null model index1, index2, null model index2 ...)"
}

\author{
	Adrien Taudiere
}

\seealso{
	\code{\link{ses.listofindex}};
	\code{\link{plot.listofindex}}
}

\examples{
data(finch.ind)
oldpar <- par(no.readonly = TRUE)

####
#The function ComIndex allow to choose your own function 
#(like mean, range, variance...) to calculate customize index.

require(e1071)

funct <- c("mean(x, na.rm = TRUE)", "kurtosis(x, na.rm = TRUE)",
"max(x, na.rm = TRUE) - min(x, na.rm = TRUE)")

res.finch.sp_loc <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = "regional.ind", ind.plot = ind.plot.finch, 
nperm = 9, print = FALSE)

res.finch.sp_reg <- ComIndex(traits = traits.finch, index = funct, 
sp = sp.finch, nullmodels = "regional.pop", 
ind.plot = ind.plot.finch, nperm = 9, print = FALSE)


####
#We can represent Standardized Effect Size (ses) 
#using the function plot(as.listofindex(list1, list2, list3))

list.ind2 <- list(res.finch.sp_loc, res.finch.sp_reg)
index.list2 <- as.listofindex(list.ind2)

plot(index.list2, type = "bytraits")

}



