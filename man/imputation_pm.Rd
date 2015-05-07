\name{imputation_pm}
\alias{imputation_pm}
\title{imputation polymorphic markers for F2 populations 
}
\description{coming...
}
\usage{
imputation_pm(progenygeno, genotypes, map,...)
}
\arguments{
  \item{progenygeno}{coming..
}
  \item{genotypes}{coming..
}
  \item{map}{coming...
}
\item{...}{unused}
}
\details{imputation with conditional expectation does not really depend on parent genotypes. It only depends on the flanking markers. But, if you do not know the parents gentoypes, you will not know the linkage phase of the markers(you can impute the genotype of that marker as homozygotes for parent1 or homozygotes for parent2, but you do not know what parent 1 is and what parent2 is, so, you cannot transform the marker genotype back to the original genotypes coded as A, T, C, G.So, if parent is NA, the imputation of this marker is not done with C++
}
\value{coming..

}
\references{coming...
}
\author{
Lian Lian
}



\examples{
	cat("imputation\n")
}
