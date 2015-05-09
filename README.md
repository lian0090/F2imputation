#Impute F2 genotypes based on parents' genotypes

This is currently only for F2 populations.

     (1) each chromosome is imputed separately
     (2) If parents are the same for that marker, marker is imputed as the parent gentoype
     (3) if there is only one marker in the chromsome, it is imputed as the mode genotype
     (4) If one parent is NA, the genotype is imputed as the mode genotype
     (5) for (3) and (4), if the mode genotype is NA, marker is not imputed, and the marker is removed from the final imputed data.
     (6) All other situations are imputed based on conditional distribution of flanking markers

##Install
F2imputation is not available on [CRAN](http://cran.r-project.org/) yet. However, it can be installed directly from GitHub using the [devtools](https://github.com/hadley/devtools) package.

1. Install `devtools` package: `install.packages('devtools')`
2. Load `devtools` package: `library(devtools)`
3. Install `BGData` package from GitHub: `install_github('lian0090/F2imputation')`
4. Load `SKAT2` package: `library(F2imputation)`

##Usage
`imputation(progeno,pargeno,map,nColSkip)`
- Arguments
    - `progeno`: progenty genotype, each row is an individual, the 1:nColSkip columns can be pedigree information, marker names must be in column names. genotypes must be coded in -1,0,1
    - `pargeno`: parent genotype, each row is an individual, the 1:nColSkip columns can be  pedigree information, marker names must be in column names. genotypes must be coded in -1,0,1
    - `map`: each row of map describes a single marker and must contain 3 columns: chromosome, SNPid, Genetic distance (cM)
    - `nColSkip`: number of columns to skip before the first marker genotype  
    - `IDcol`: column number for LINE ID, if specified, the LINE IDs will be added to the dimnames of genotype probability array. 

- Return Values
    - `pimputedgeno`: a dataframe of imputed genotypes  
    - `genoprob`: an array of genotype probabilities
##Example
```R
data(map)
data(progeno)
data(pargeno)
out=imputation(progeno,pargeno,map,nColSkip=1)
```
