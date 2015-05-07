imputation=function(combinedgeno){
LINEnames=combinedgeno[,1]
rowheaders=combinedgeno[c(1:4),-1]
combinedgeno=combinedgeno[,-1]	
pimputedgeno=combinedgeno[-c(1:4),]  ##the header were removed and the column of LINE was removed from the last step. 
sameparent=which(combinedgeno[3,]==combinedgeno[4,])
NAparent=which(is.na(combinedgeno[3,])|is.na(combinedgeno[4,]))
non_polymorph=union(sameparent,NAparent)
#only the polymorphic geno enters into the imutation based on conditional probability
pmgeno=combinedgeno[,-non_polymorph]
genoprob=array(NA,dim=c(nrow(pimputedgeno),ncol(pimputedgeno),3))
pmgenoprob=genoprob[,-non_polymorph,]
#seperate chromosomes
nchr=length(unique(as.numeric(combinedgeno[1,])))
#NApa_pr:only one marker on this chromosme,most progeny are NA for this marker or parents genotypes NA and most progeny genotypes NA for this marker.
NApa_pr=vector()
for (j in c(1:nchr)){
whichchr=which(as.integer(pmgeno[1,])==j)
if(length(whichchr)>0){
chr=as.matrix(pmgeno[,whichchr])
parent1=chr[3,]
parent2=chr[4,]
progenygeno=chr[-c(1:4),]
map=chr[2,]
chrmknames=colnames(chr)
###if only one marker, impute with the mode genotype
if(length(whichchr)==1){
	warning(paste("there is only 1 marker on chr",j,sep=""))
	tmp=progenygeno
tabtmp=table(factor(tmp,levels=c(-1,0,1,NA)),useNA="always")
modegeno=names(tabtmp[which.max(tabtmp)])
if(!is.na(modegeno)){
NAi=which(is.na(tmp))
pmgenoprob[NAi,whichchr,]=rep(tabtmp[1:3]/length(na.omit(tmp)),each=length(NAi))
pmgenoprob[-NAi,whichchr,]=model.matrix(~factor(tmp[-NAi],levels=c(-1,0,1))-1)
tmp[NAi]=modegeno        
pmgeno[-c(1:4),whichchr]=tmp
} else{ 
	NApa_pr=c(NApa_pr,which(colnames(pimputedgeno)==chrmknames))}        
}else {
#imputation with C++
progenygeno=recodegeno_IBD(genomatrix=progenygeno,parent1,parent2,oriheter=0,newgeno=c(1,2,3)) 
output=imputation_pm(progenygeno,genotypes=c(1,2,3),map)
pmgeno[-c(1:4),whichchr]=recodegeno_IBDtoOri(genomatrix=output$progenygeno,parent1,parent2,oriheter=0,IBDgeno=c(1,2,3))
pmgenoprob[,whichchr,]=genoprob_IBDtoOri(Genoprob=output$genoprob,parent1,parent2,origeno=c(-1,0,1))
}
}
}
#fill the imputed geno with the imputed genotypes for polymorhpic markers
pimputedgeno[,-non_polymorph]=pmgeno[-c(1:4),]
genoprob[,-non_polymorph,]=pmgenoprob
#end of polymorphic markers

#if parents are the same for the marker, then the marker is imputed as parents marker.
pimputedgeno[,sameparent]=combinedgeno[3,sameparent]
for(sameparenti in sameparent){
genoprob[,sameparenti,]=model.matrix(~factor(pimputedgeno[,sameparenti],levels=c(-1,0,1))-1)}
#end of same parents

#if parents are NA for the marker, then the missing marker of this genotype is imputed as the mode genotype of the other progenies.
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
if(length(NApa_pr)>0){
pimputedgeno=pimputedgeno[,-NApa_pr]
genoprob=genoprob[,-NApa_pr,]}
#end of NA parent
dimnames(genoprob)=list(LINEnames[-c(1:4)],colnames(pimputedgeno),c(-1,0,1))
pimputedgeno=cbind(LINEnames,rbind(combinedgeno[c(1:4),colnames(pimputedgeno)],pimputedgeno))
out=list(pimputedgeno=pimputedgeno,genoprob=genoprob)
return(out)
}