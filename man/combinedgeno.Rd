\name{combinedgeno}
\alias{combinedgeno}
\docType{data}
\title{
map positions cominbed with parent genotypes and progeny genotypes

}
\description{
A dataframe or matrix of characters (Do not use factors!)	
colname 1 is the name of LINEs 
row 1:4 are map chr, map pos, P1 genotypes, P2 genotypes
All genotypes should be coded as -1,0,1
}
\usage{data(combinedgeno)}
\format{
  A data frame with 10 observations on the following 60 variables.
  \describe{
    LINE:a character vector of : chr, pos,P1, P2,the name of progeny LINEs. chr,pos,P1 and P2 can be empty.
    M1-M60: marker genetypes:first two are chr and map position, 3 to 4 are parents genotypes,all others are progeny genotypes
      }
}
\examples{
data(combinedgeno)
combinedgeno[1:6,1:10]
}
