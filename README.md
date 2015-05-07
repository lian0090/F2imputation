#Impute F2 genotypes based on parents' genotypes

This is currently only for F2 populations.
(1) each chromosome is imputed separately
(2) If parents are the same for that marker, marker is imputed as the parent gentoype
(3) if there is only one marker in the chromsome, it is imputed as the mode genotype
(4) If one parent is NA, the genotype is imputed as the mode genotype
(5) for (3) and (4), if the mode genotype is NA, marker is not imputed, and the marker is removed from the final imputed data.
(6) All other situations are imputed based on conditional distribution of flanking markers


##Example
'''R
data(combinedgeno)
outdat=imputation(combinedgeno)
outdat$pimputedgeno
outdat$genoprob
'''
