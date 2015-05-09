\name{imputation.cmbgeno}
\alias{imputation.cmbgeno}
\title{
Impute genotypes of F2 populations based on flanking markers. 
}
\description{
Impute genotypes of F2 populations based on flanking markers. The input format should be combinedgeno.
}
\usage{
imputation(combinedgeno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combinedgeno}{
  	nmar+1 columns and nLINE+4 rows. First column is LINE names, first four rows are map chr, map pos, P1 genotype and P2 genotype
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
outcmb=imputation.cmbgeno(combinedgeno)
outcmb.geno=outcmb$pimputedgeno[-c(1:4),]
outcmb.prob=outcmb$genoprob
data(progeno)
data(pargeno)
data(map)
out=imputation(progeno,pargeno,map,nColSkip=1,IDcol=NULL)
out.geno=out$pimputedgeno
out.prob=out$genoprob
all(out.geno==outcmb.geno)
all(out.prob==outcmb.prob)
}
