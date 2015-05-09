imputation=function(progeno,pargeno,map,nColSkip,IDcol=NULL){
  #progeno:progenty genotype matrix, the 1:nColSkip columns can be  pedigree information, marker names must be in column names
  #pargeno:parent genotype matrix, the 1:nColSkip columns can be  pedigree information, marker names must be in column names
  ##map: each row of map describes a single marker and must contain 3 columns: chromosome, SNPid, Genetic distance (cM)
  #nColSkip: number of columns to skip before the first marker genotype  
  #IDcol: column number for LINE ID. 
  pimputedgeno=getExpProgeno(pargeno,progeno,nColSkip=nColSkip)
  rm("progeno")
  if(nColSkip>0){
    geno.ped=pimputedgeno[,c(1:nColSkip),drop=F]
    pimputedgeno=pimputedgeno[,-c(1:nColSkip)]
    pargeno=pargeno[,-c(1:nColSkip)]
  }else geno.ped=NULL 
  if(!is.null(IDcol)){
    LINEnames=geno.ped[,IDcol]
  }else{
    LINEnames=NULL
  }
  genoprob=array(NA,dim=c(nrow(pimputedgeno),ncol(pimputedgeno),3))
  sameparent=which(pargeno[1,]==pargeno[2,])
  NAparent=which(is.na(pargeno[1,])|is.na(pargeno[2,]))
  non_polymorph=union(sameparent,NAparent)
  #NApa_pr:
  #(1) only one marker on this chromosme,most progeny are NA for this marker 
  #or (2) parents genotypes NA and most progeny genotypes NA for this marker.
  NApa_pr=vector()
  
  #only the polymorphic geno enters into the imutation based on conditional probability
  pmgeno=pimputedgeno[,-non_polymorph]
  pmpargeno=pargeno[,-non_polymorph]
  pmgenoprob=genoprob[,-non_polymorph,]
  pmmap=map[-non_polymorph,]
  
  ######Begin imputation on polymorphic markers 
  #seperate chromosomes
  nchr=length(unique(as.numeric(pmmap[,1])))
  for (j in c(1:nchr)){
    whichchr=which(as.integer(pmmap[,1])==j)
    if(length(whichchr)>0){
      pmgeno.j=as.matrix(pmgeno[,whichchr])
      p1.j=pmpargeno[1,whichchr]
      p2.j=pmpargeno[2,whichchr]
      map.j=pmmap[whichchr,]
      mknames.j=map.j[,2]
      ###if only one marker, impute with the mode genotype
      if(length(whichchr)==1){
        warning(paste("there is only 1 marker on chr",j,sep=""))
        tmp=pmgeno.j
        tabtmp=table(factor(tmp,levels=c(-1,0,1,NA)),useNA="always")
        modegeno=names(tabtmp[which.max(tabtmp)])
        if(!is.na(modegeno)){
          NAi=which(is.na(tmp))
          pmgenoprob[NAi,whichchr,]=rep(tabtmp[1:3]/length(na.omit(tmp)),each=length(NAi))
          pmgenoprob[-NAi,whichchr,]=model.matrix(~factor(tmp[-NAi],levels=c(-1,0,1))-1)
          tmp[NAi]=modegeno        
          pmgeno[,whichchr]=tmp
        } else{ 
          NApa_pr=c(NApa_pr,which(map[,2]==mknames.j))
        }        
        ### end imputation with mode genotype
      }else {
        #imputation with conditonal probability 
        pmgeno.j=recodegeno_IBD(genomatrix=pmgeno.j,p1.j,p2.j,oriheter=0,newgeno=c(1,2,3)) 
        output=imputation_pm(pmgeno.j,genotypes=c(1,2,3),map.j[,3])
        pmgeno[,whichchr]=recodegeno_IBDtoOri(genomatrix=output$progenygeno,p1.j,p2.j,oriheter=0,IBDgeno=c(1,2,3))
        pmgenoprob[,whichchr,]=genoprob_IBDtoOri(Genoprob=output$genoprob,p1.j,p2.j,origeno=c(-1,0,1))
      }
    }
  }
  #fill the imputed geno with the imputed genotypes for polymorhpic markers
  pimputedgeno[,-non_polymorph]=pmgeno
  genoprob[,-non_polymorph,]=pmgenoprob
  ######end of polymorphic markers
  
  #if parents are the same for the marker, then the marker is imputed as parents marker.
  if(length(sameparent)>0){
    for(i in 1:length(sameparent)){
      sameparenti=sameparent[i]
      pimputedgeno[,sameparenti]=pargeno[1,sameparenti]
      genoprob[,sameparenti,]=model.matrix(~factor(pimputedgeno[,sameparenti],levels=c(-1,0,1))-1)}
  }
  #end of same parents
  
  #Imputation for NAparent:
  #If parents are NA for the marker, then the missing marker of this genotype is imputed as the mode genotype of the other progenies.
  if(length(NAparent)>0){
    for (i in 1:length(NAparent)){
      tmp=pimputedgeno[,NAparent[i]]
      tabtmp=table(factor(tmp,levels=c(-1,0,1,NA)),useNA="always")
      modegeno=names(tabtmp[which.max(tabtmp)])
      if(!is.na(modegeno)){
        NAi=which(is.na(tmp))
        if(length(NAi)>0){pimputedgeno[NAi,NAparent[i]]=modegeno
                          genoprob[NAi,NAparent[i],]=rep(tabtmp[1:3]/length(na.omit(tmp)),each=length(NAi))
                          genoprob[-NAi,NAparent[i],]=model.matrix(~factor(pimputedgeno[-NAi,NAparent[i]],levels=c(-1,0,1))-1)}else{
                          genoprob[,NAparent[i],]=model.matrix(~factor(pimputedgeno[,NAparent[i]],levels=c(-1,0,1))-1)}} else
                          #the mode genotype is NA for this marker, so this marker is noted and will be removed from analysis.	
                          NApa_pr=c(NApa_pr,NAparent[i])
    }
  }
  #end of NAparent
  #remove non-imputed markers
  if(length(NApa_pr)>0){
    pimputedgeno=pimputedgeno[,-NApa_pr]
    genoprob=genoprob[,-NApa_pr,]
  }
  #
  dimnames(genoprob)=list(LINEnames,map[,2],c(-1,0,1))
  colnames(pimputedgeno)=map[,2]
  pimputedgeno=cbind(geno.ped,pimputedgeno)
  rownames(pimputedgeno)=NULL
  out=list(pimputedgeno=pimputedgeno,genoprob=genoprob)
  return(out)
}