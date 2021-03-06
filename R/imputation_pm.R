imputation_pm<-
function(progenygeno,genotypes,SNPcM,...){
 n_pos=ncol(progenygeno)
 n_ind=nrow(progenygeno)
 progenygeno=as.integer(as.matrix(progenygeno))
 genotypes=as.integer(genotypes) 
 n_genotype=as.integer(length(genotypes))
 SNPcM=as.double(SNPcM)
 genoprob=as.double(array(0,dim=c(n_ind,n_pos,n_genotype)))
 naind=as.integer(rep(NA,n_ind))
output=.C("imputation",progenygeno=progenygeno,genotypes,n_genotype,n_pos,n_ind,SNPcM,genoprob=genoprob,naind=naind,NAOK=T,PACKAGE="F2imputation")
progenygeno=matrix(output$progenygeno,nrow=n_ind,byrow=F)
genoprob=array(output$genoprob,dim=c(n_ind,n_pos,n_genotype))
naind=(output$naind)
return(list(progenygeno=progenygeno,genoprob=genoprob,naind=naind[!is.na(naind)]))
}
