##expand progeny genotypes to the same number of columns as parent genotypes
getExpProgeno=function(pargeno,progeno,nColSkip=1){
  #nColSkip: number of columns to skip before the first marker genotype	
  if(nColSkip>0){
  par.markername=colnames(pargeno)[-c(1:nColSkip)]
  pro.markername=colnames(progeno)[-c(1:nColSkip)]
  skip.name=colnames(progeno)[1:nColSkip]
  }else{
  par.markername=colnames(pargeno)
  pro.markername=colnames(progeno)
  skip.name=NULL
  	}
  nanames=setdiff(par.markername,pro.markername)
  nalength=length(nanames)
  progeno=cbind(progeno,matrix(NA,nrow=nrow(progeno),ncol=nalength,dimnames=list(NULL,nanames)))
  progeno=progeno[,c(skip.name,par.markername)]
  return(progeno)
}