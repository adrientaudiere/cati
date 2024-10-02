\name{RaoRel}
\alias{RaoRel}

\title{
	Alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity
}

\description{
	The Rao function computes alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity with the Rao index. The script integrates two functions: "Qdecomp", by Villeger et Mouillot (J Ecol, 2008) modified by Wilfried Thuiller, and "disc", by S. Pavoine, in the package ade4. For a regional assemblage of C local communities gamma = mean(alpha) + beta, where: gamma is the diversity of the regional pool, alpha is the diversity of the local community and beta is the turn over between local communities diversity is estimated with the Rao quadratic entropy index (Rao 1982)
} 
	
\usage{
	RaoRel(sample, dfunc, dphyl, weight = FALSE, Jost = FALSE, 
	structure = NULL)
}

\arguments{
	\item{sample}{
		 Community matrix of abundance (c x s) of the s species for the c local communities.
	}
	
	\item{dfunc}{
		matrix (s x s) or dist object with pairwise functional trait distances between the s species
	}
	
	\item{dphyl}{
		As dfunct but for phylogenetic distances
	}
	
	\item{weight}{
		Defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
	}
	
	\item{Jost}{
		Defining if the Jost correction is applied (Jost 2007)
	}
	
	\item{structure}{
		A data frame containing the name of the group to which samples belong see de Bello et al, 2011 for more details.
	}
}

\details{
	NA are automatically replaced by 0 in "sample". This function use the function "Qdecomp" by Sebastien Villeger & David Mouillot (J Ecol, 2008) modified by Wilfried Thuiller and the function disc originally proposed by Sandrine Pavoine.
}

\value{
	The results are organized for Taxonomic diversity ($TD), Functional diversity ($FD) and phylogenetical diversity ($PD). Beta and gamma diversities are calculated for the whole data set and for each pair of samples ("Pairwise_samples"):
	
	$Richness_per_plot(number of species per sample)
	
	$Relative_abundance (species relative abundances per plot)
	
	$Pi (species regional relative abundance)
	
	$Wc (weigthing factor),                               
	
	$Mean_Alpha (mean aplpha diversity; for taxonomic diversity the Simpson index is calculated)                                   
	
	$Alpha (alpha diversity for each sample; for taxonomic diversity the Simpson index is calculated)                                       
	
	$Gamma (gamma diversity; for taxonomic diversity the Simpson index is calculated)                                      
	
	$Beta_add (Gamma-Mean_Alpha)                                          
	
	$Beta_prop (Beta_add*100/Gamma)                                                         
	
	$Pairwise_samples$Alpha (mean alpha for each pair of samples)                        
	
	$Pairwise_samples$Gamma (gamma for each pair of samples)
	
	$Pairwise_samples$Beta_add (beta for each pair of samples as Gamma-Mean_Alpha)  
	
	$Pairwise_samples$Beta_prop (beta for each pair of samples as Beta_add*100/Gamma) 
}

\references{
	De Bello, Francesco, Sandra Lavorel, Cecile H. Albert, Wilfried Thuiller, Karl Grigulis, Jiri Dolezal, stepan Janecek, et Jan Leps. 2011. Quantifying the relevance of intraspecific trait variability for functional diversity: Intraspecific variability in functional diversity. Methods in Ecology and Evolution 2: 163-174. 
}

\author{
	Francesco De Bello et al., 2011 modified by Adrien Taudiere
}

\examples{
	data(finch.ind)
	
	\dontrun{
		comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))
		comm.sp <- table(sp.finch, ind.plot.finch)
		class(comm.sp) <- "matrix"
		
		traits.finch.sp <- apply( apply(traits.finch, 2, scale ), 2, 

		function(x) tapply(x, sp.finch, mean, na.rm = TRUE))
		
		mat.dist <- (as.matrix(dist(traits.finch.sp))^2)/2
		
		res.rao <- RaoRel(sample = as.matrix(comm.sp), dfunc = mat.dist, dphyl = NULL, 
		weight = FALSE, Jost = FALSE, structure = NULL)

		function(x) tapply(x, sp.finch, mean, na.rm=TRUE))
		
		mat.dist <- as.matrix(dist(traits.finch.sp))^2
		
		res.rao <- RaoRel(sample=as.matrix(comm.sp), dfunc=mat.dist, dphyl=NULL, 
		weight=FALSE, Jost=FALSE, structure=NULL)

		
		witRao <- res.rao$FD$Mean_Alpha  #overall within species variance
		betRao <- res.rao$FD$Beta_add    #between species variance
		totRao <- res.rao$FD$Gamma       #the total variance
		
		witRao+betRao
		totRao
		
		
		#Now let"s take the abundance to calculate Rao diversity.
		

		res.rao.w <- RaoRel(sample = as.matrix(comm.sp), dfunc = mat.dist, dphyl = NULL, 
		weight = TRUE, Jost = FALSE, structure = NULL)
		
		witRao.w <- res.rao.w$FD$Mean_Alpha  #overall within species variance
		betRao.w <- res.rao.w$FD$Beta_add    #between species variance
		totRao.w <- res.rao.w$FD$Gamma       #the total variance
		
		witRao.w
		betRao.w
		
		
		#Plot the results
		
		barplot(cbind(c(witRao.w, betRao.w), c(witRao, betRao)), 
		names.arg  = c("abundance" ,"presence"), 
		legend.text = c("within species", "between species"), 
		ylab = "Rao", ylim = c(0,10))
		
		
		#We can do this analysis for each trait separately. 
		#First we need to replace (or exclude) NA values. 
		#For this example, we use the package mice to complete the data.
		
		comm <- t(table(ind.plot.finch,1:length(ind.plot.finch)))
		
		library(mice)

		traits = traits.finch

		mice <- mice(traits.finch)
		traits.finch.mice <- complete(mice)
		
		traits.finch.mice.sp <- apply(apply(traits.finch.mice, 2, scale ), 2, 
		function(x) tapply(x, sp.finch, mean, na.rm = TRUE))


		trait.rao.w <- list()
		witRao.w.bytrait <- c()
		betRao.w.bytrait <- c()

		for (t in 1 : 4){
		  trait.rao.w[[t]] <- RaoRel(sample = as.matrix(comm.sp), 
		  dfunc = (dist(traits.finch.mice.sp[,t])^2)/2, dphyl = NULL, weight = TRUE, 
		  Jost = FALSE, structure = NULL)
				  
		  witRao.w.bytrait <- c(witRao.w.bytrait, trait.rao.w[[t]]$FD$Mean_Alpha)
		  betRao.w.bytrait <- c(betRao.w.bytrait, trait.rao.w[[t]]$FD$Beta_add)
		}
		
		#Plot the results by traits.
		
		barplot(t(cbind( witRao.w.bytrait, betRao.w.bytrait)), 
		names.arg = colnames(traits.finch),
		legend.text = c("within species", "between species"), 
		ylab = "Rao", ylim = c(0,1.5))	
	}
}



