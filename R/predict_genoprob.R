predict_genoprob <-
function(l_geno,r_geno,recombinations=c("l","r","total"),...){
	genotypes=c(1,2,3)
 genotypes=as.integer(genotypes)
 r_prob=l_prob=rep(0,3)
 l_prob[l_geno]=1
 r_prob[r_geno]=1
 l_prob=as.double(l_prob)
 r_prob=as.double(r_prob)
 rlk_rrk_r=as.double(recombinations)
 genoprob=as.double(rep(0,length(genotypes)))
 .C("predict_prob",genotypes,l_prob,r_prob,rlk_rrk_r,genoprob=genoprob,PACKAGE="F2imputation")$genoprob	
 }
