calc_genoprob <-
function(l_geno,r_geno,recombinations=c("l","r","total"),genotypes=c(1,2,3),...){
 genotypes=as.integer(genotypes)
 l_geno=as.integer(l_geno)
 r_geno=as.integer(r_geno)
 r_lk=as.double(recombinations[1])
 r_rk=as.double(recombinations[2])
 r=as.double(recombinations[3])
 genoprob=as.double(rep(0,length(genotypes)))
 .C("calc_genoprob_R",genotypes,l_geno,r_geno,r_lk,r_rk,r,genoprob=genoprob,PACKAGE="F2imputation")$genoprob	
 }
