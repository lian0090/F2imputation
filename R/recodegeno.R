
##turn  genotypes into newgeno types based on their parent genotypes genotypes.If geno==parent1, newgeno=P1,if geno==oriheter,newgeno=heter,if geno==parent2,newgeno=P2. Non-homozygotes loci are coded into NA.

recodegeno_IBD=function(genomatrix,parent1,parent2,oriheter=0,newgeno=c(-1,0,1),...){
parent1=as.integer(parent1)
parent2=as.integer(parent2)
geno=as.integer(as.matrix(genomatrix))
n_pos=as.integer(ncol(genomatrix))
n_ind=as.integer(nrow(genomatrix))
oriheter=as.integer(oriheter)
newgeno=as.integer(newgeno)
geno=.C("RECODEGENO_IBD",parent1,parent2,geno,n_pos,n_ind,oriheter,newgeno,NAOK=T,PACKAGE="F2imputation")[[3]]  
matrix(geno,nrow=n_ind,byrow=F)

 }
recodegeno_IBDtoOri=function(genomatrix,parent1,parent2,oriheter=0,IBDgeno=c(-1,0,1),...){
	#description:transform the relative genotypes(relative to the parents) into the absolute genotypes, which is univeral across all populations.
	#example: if in the progeny, all markers equal to parent 1 is coded as 1, all heterozygotes coded as 2, all markers equal to parent2 is coded as 3, then IBDgeno=c(1,2,3), for two markers,matrix(rep(c(1,2,3),2),ncol=2) the absolute values of parent1 is c(1,-1), the absolute values of parent2 is c(-1,1) and the absolute values of heterozygotes is 0, so parent1=c(1,-1),parent2=c(-1,1),oriheter=0, then the marker genotypes is transformed to matrix(c(c(1,0,-1),c(-1,0,1)),ncol=2)
#recodegeno_IBDtoOri(matrix(rep(c(1,2,3),2),ncol=2),parent1=c(1,-1),parent2=c(-1,1),oriheter=0,IBDgeno=c(1,2,3))	
	parent1=as.integer(parent1)
	parent2=as.integer(parent2)
	geno=as.integer(as.matrix(genomatrix))
	n_pos=as.integer(ncol(genomatrix))
	n_ind=as.integer(nrow(genomatrix))
	oriheter=as.integer(oriheter)
	IBDgeno=as.integer(IBDgeno)
	geno=.C("recodegeno_IBDtoOri",parent1,parent2,geno,n_pos,n_ind,oriheter,IBDgeno,NAOK=T,PACKAGE="F2imputation")[[3]]
	matrix(geno,nrow=n_ind,byrow=F)

	 }


 recodegeno<-function(method=c("ChrToInt","IntToChr"),genomatrix,...){
 ##it recodes A, H, B to 1,2,3 this is stupid...	
  if(is.vector(genomatrix)) n.row=1
  else n.row=nrow(genomatrix)
  if (method=="ChrToInt") 
 {cfunc="recodegeno_ChrToInt";
  chrgeno=as.character(genomatrix)	
  intgeno=as.integer(rep(0,length(chrgeno)))
  size=as.integer(length(chrgeno))
 geno=.C(cfunc,intgeno,chrgeno,size,NAOK=T)[[1]]
  }
 if (method=="IntToChr") {cfunc="recodegeno_IntToChr";
 intgeno=as.integer(genomatrix)	
 chrgeno=as.character(rep("N",length(intgeno)))
  size=as.integer(length(intgeno))
  geno=.C(cfunc,intgeno,chrgeno,size,NAOK=T,PACKAGE="F2imputation")[[2]]
  }
  geno[geno==""]=NA
  if (n.row==1) return(geno)
  else return(matrix(geno,nrow=n.row,byrow=F))
   	
 }
 
genoprob_IBDtoOri=function(Genoprob,parent1,parent2,origeno=c(-1,0,1),...)
{#the IBD genoprobabilities is arranged as geno(P1), geno(heter),geno(P2), the raw genoprobabilities should be arranged as origeno[1],origeno[2], origeno[3], in this case, as (-1,0,1)
	npos=dim(Genoprob)[2]
        for (i in 1:npos){
if(parent1[i]==origeno[3]){
	tmp1=Genoprob[,i,3]
	tmp3=Genoprob[,i,1]
	Genoprob[,i,1]=tmp1
	Genoprob[,i,3]=tmp3}
	}
	return(Genoprob)
}
	
