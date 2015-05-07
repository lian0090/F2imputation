\name{imputation}
\alias{imputation}
\title{
Impute genotypes of F2 populations based on Flanking markers
}
\description{
This is currently only for F2 populations.
}
\usage{
imputation(combinedgeno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combinedgeno}{nmar+1 columns and nLINE+4 rows. First column is LINE names, first four rows are map chr, map pos, P1 genotype and P2 genotype
}
}
\value{
return out:a list of pimputedgeno and genoprob
out$pimputedgeno is the same format as combinedgeno, but with the imputed geotypes
out$genoprob an array of nLINE * nmar* 3. The probabilities of three gentoypes (-1,0,1)
}
\references{
conditional expectation table for F2
}
\author{
Lian Lian 
}




\examples{
data(combinedgeno)
out=imputation(combinedgeno)
  names(out)
}
