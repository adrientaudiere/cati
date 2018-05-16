\name{RandCom}
\alias{RandCom}


\title{
	Toy model to simulate internal and/or external filtering
}

\description{
	Toy model to simulate internal and/or external filtering
}

\usage{
	RandCom(Ncom = 10, Nsp = 20, Nind.com = 100, sdlog = 1.5, 
	min_value_traits = 80, max_value_traits = 200, 
	cv_intra_sp = 1.5, cv_intra_com = 1.5, 
	Int_Filter_Strength = 50, Ext_Filter_Strength = 50, Filter="None")
}

\arguments{
	
	\item{Ncom}{
		Number of communities (or sites).
	}
	
	\item{Nsp}{
		Number of species at the regional scale.
	}
	
	\item{Nind.com}{
		Number of individuals by communities.
	}
	
	\item{sdlog}{
		Parameter of the log normal distribution for species abundances distribution within communities.
	}
		
	\item{min_value_traits}{
		Minimum mean value for traits distributions.
	}
	
	\item{max_value_traits}{
		Maximum mean value for traits distributions.
	}
		
	\item{cv_intra_sp}{
		Coefficient of variation for intra-specific distributions. The more the value is high the less there is internal filtering. Used only for the trait 1 (normally distributed).
	}
	
	\item{cv_intra_com}{
		Coefficient of variation for intra-community distributions. The more the value is high the less there is external filtering. Used only for the trait 1 (normally distributed)
	}
		
	\item{Int_Filter_Strength}{
		Strength of internal filtering in percentage. Use in addition to cv_intra_sp by distributing mean species trait more or less evenly. In the most extreme case (if Int_Filter_Strength==100), species have equally distributed mean values along the trait gradient.
	}
	
	\item{Ext_Filter_Strength}{
		Strength of external filtering in percentage. Use in addition to cv_intra_com by distributing mean communities trait more or less evenly. In the most extreme case (if Ext_Filter_Strength==100), communities have equally distributed mean values along the trait gradient.
	}
				
	\item{Filter}{
		The type of filter to simulate. Either "None", "Internal", "External" or "Both"
	}	
}

\details{
	In this version of the function, the trait 1 follows a normal distribution wheras the trait 2 follows a uniform distribution.
}

\value{
	\item{$data$com}{Vector of simulated communities for each individual.}
	\item{$data$sp}{Vector of simulated species for each individual.}
	\item{$data$trait1}{Vector of simulated value for the trait 1: normally distributed.}
	\item{$data$trait2}{Vector of simulated value for the trait 2: normally distributed.}
	
 	\item{$call}{call of the function Tstats}
}

\author{
	Cecile Albert and Adrien Taudiere
}

\examples{	
	res <- RandCom()
}
