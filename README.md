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
    - `progeno` and `pargeno`: typed markers of F2 progeny genotypes (`progeno`) and parents genotypes (`pargeno`). Each row is an individual and there should be only two rows for `pargeno` for the two parents.  The `1:nColSkip` columns can be pedigree information. Marker names must be in column names. Genotypes must be coded in `-1,0,1`. For example, if progeny are only genotyped for `M1`, `M3`, and parents are genotyped for `M1`, `M2`, `M3`, `M4`, `M5`. If ID is included in the genotype, then, `progeno` will have three columns with column names `ID`, `M1`, `M3`. `pargeno` will have six columns with column names `ID`, `M1`, `M2`, `M3`, `M4`, `M5`. 
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

##Compatability with old version
`imputation.cmbgeno(combinedgeno)`
- arguments
      - `combinedgeno`: a matrix or data.frame combined the map information, parents genotype and progeny genotype. Entries of combinedgeno are `character` types so that all the information can be put together into a single matrix. Marker genotypes must be coded as -1, 0, 1. 
           - column 1 is LINE names, first two rows of column 1 are NA. 
           - row 1 is chrosomes number
           - row 2 is genetic distance in centimorgan (cM)
           - row 3 and row 4 are marker genotype for the two parents
           - row 5 onwards are the progeny genotypes


- Return Values
    - `pimputedgeno`: a dataframe of imputed genotypes  
    - `genoprob`: an array of genotype probabilities
- R code comparing the old version and new version

```R
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
```


