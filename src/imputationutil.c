#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "imputationutil.h"



//calculate recombination rates from map distance
void map2r(double *map,double *r,int pos1, int pos2){
double d=fabs(map[pos1]-map[pos2])/100;//use Morgan rather than cM
*r=(1-exp(-2*d))/2;}

void calc_genoprob_R(int *genotypes,int *l_geno,int *r_geno,double *r_lk,double *r_rk,double *r,double *genoprob){
calc_genoprob_c(genotypes,*l_geno,*r_geno,*r_lk,*r_rk,*r,genoprob);}


//calculate genoprobs based on predicted flanking marker 
void predict_prob(int *genotypes,double *l_prob,double *r_prob, double *rlk_rrk_r,double *genoprob){
//probilitity is stored as prob[left][right][target}
double r_lk=rlk_rrk_r[0];
double r_rk=rlk_rrk_r[1];
double r=rlk_rrk_r[2];
double prob[3][3][3];
 prob[0][0][0]=pow((1-r_lk)*(1-r_rk)/(1-r),2);
 prob[0][0][1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/pow((1-r),2);
 prob[0][0][2]=pow(r_lk,2)*pow(r_rk,2)/pow((1-r),2);

 prob[0][1][0]=pow((1-r_lk),2)*(1-r_rk)*r_rk/((1-r)*r);
 prob[0][1][1]=r_lk*(1-r_lk)*pow((1-r_rk),2)/((1-r)*r)+r_lk*(1-r_lk)*pow(r_rk,2)/((1-r)*r);
 prob[0][1][2]=pow(r_lk,2)*r_rk*(1-r_rk)/((1-r)*r);

 prob[0][2][0]=pow((1-r_lk)*r_rk/r,2);
 prob[0][2][1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/(r*r);
 prob[0][2][2]=pow(r_lk*(1-r_rk)/r,2);

prob[1][0][0]=(1-r_lk)*r_lk*pow((1-r_rk),2)/((1-r)*r);
prob[1][0][1]=pow((1-r_lk),2)*r_rk*(1-r_rk)/((1-r)*r)+r_lk*r_lk*r_rk*(1-r_rk)/((1-r)*r);
prob[1][0][2]=(1-r_lk)*r_lk*r_rk*r_rk/((1-r)*r);

prob[1][1][0]=2*(1-r_lk)*(r_lk)*(1-r_rk)*r_rk/((1-r)*(1-r)+r*r);
prob[1][1][1]=(r_lk*r_lk*r_rk*r_rk+pow((1-r_lk),2)*pow(r_rk,2)+pow(r_lk,2)*pow((1-r_rk),2)+pow((1-r_lk),2)*pow((1-r_rk),2))/((1-r)*(1-r)+r*r);
prob[1][1][2]=2*(1-r_lk)*r_lk*r_rk*(1-r_rk)/((1-r)*(1-r)+r*r);

prob[1][2][0]=r_lk*(1-r_lk)*pow(r_rk,2)/(r*(1-r));
prob[1][2][1]=pow((1-r_lk),2)*r_rk*(1-r_rk)/(r*(1-r))+pow(r_lk,2)*r_rk*(1-r_rk)/(r*(1-r));
prob[1][2][2]=(1-r_lk)*r_lk*pow((1-r_rk),2)/(r*(1-r));

prob[2][0][0]=pow(r_lk,2)*pow((1-r_rk),2)/(r*r);
prob[2][0][1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/(r*r);
prob[2][0][2]=pow((1-r_lk)*r_rk/r,2);

prob[2][1][0]=r_lk*r_lk*(1-r_rk)*r_rk/((1-r)*r);
prob[2][1][1]=r_lk*(1-r_lk)*pow((1-r_rk),2)/((1-r)*r)+r_lk*(1-r_lk)*pow(r_rk,2)/((1-r)*r);
prob[2][1][2]=pow((1-r_lk),2)*r_rk*(1-r_rk)/((1-r)*r);

prob[2][2][0]=pow(r_lk*r_rk/(1-r),2);
prob[2][2][1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/((1-r)*(1-r));
prob[2][2][2]=pow((1-r_lk)*(1-r_rk)/(1-r),2);

int i; int j;int k;
for (k=0;k<3;k++){
for (i=0;i<3;i++){
for (j=0;j<3;j++){
genoprob[k]+=prob[i][j][k]*l_prob[i]*r_prob[j];
}}}

}

//calculate genoprobs based on really typed flanking marker
void calc_genoprob_c(int *genotypes,int l_geno,int r_geno,double r_lk,double r_rk,double r,double *genoprob){

if (l_geno==genotypes[0]){
if (r_geno==genotypes[0]){
genoprob[0]=pow((1-r_lk)*(1-r_rk)/(1-r),2);
genoprob[1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/pow((1-r),2);
genoprob[2]=pow(r_lk,2)*pow(r_rk,2)/pow((1-r),2);
}
else if (r_geno==genotypes[1]){
genoprob[0]=pow((1-r_lk),2)*(1-r_rk)*r_rk/((1-r)*r);
genoprob[1]=r_lk*(1-r_lk)*pow((1-r_rk),2)/((1-r)*r)+r_lk*(1-r_lk)*pow(r_rk,2)/((1-r)*r);
genoprob[2]=pow(r_lk,2)*r_rk*(1-r_rk)/((1-r)*r);
}
else if (r_geno==genotypes[2]){
genoprob[0]=pow((1-r_lk)*r_rk/r,2);
genoprob[1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/(r*r);
genoprob[2]=pow(r_lk*(1-r_rk)/r,2);
}
}

else if (l_geno==genotypes[1])
{
if (r_geno==genotypes[0]){
genoprob[0]=(1-r_lk)*r_lk*pow((1-r_rk),2)/((1-r)*r);
genoprob[1]=pow((1-r_lk),2)*r_rk*(1-r_rk)/((1-r)*r)+r_lk*r_lk*r_rk*(1-r_rk)/((1-r)*r);
genoprob[2]=(1-r_lk)*r_lk*r_rk*r_rk/((1-r)*r);
}
else if (r_geno==genotypes[1]){
genoprob[0]=2*(1-r_lk)*(r_lk)*(1-r_rk)*r_rk/((1-r)*(1-r)+r*r);
genoprob[1]=(r_lk*r_lk*r_rk*r_rk+pow((1-r_lk),2)*pow(r_rk,2)+pow(r_lk,2)*pow((1-r_rk),2)+pow((1-r_lk),2)*pow((1-r_rk),2))/((1-r)*(1-r)+r*r);
genoprob[2]=2*(1-r_lk)*r_lk*r_rk*(1-r_rk)/((1-r)*(1-r)+r*r);
}
else if (r_geno==genotypes[2]){
genoprob[0]=r_lk*(1-r_lk)*pow(r_rk,2)/(r*(1-r));
genoprob[1]=pow((1-r_lk),2)*r_rk*(1-r_rk)/(r*(1-r))+pow(r_lk,2)*r_rk*(1-r_rk)/(r*(1-r));
genoprob[2]=(1-r_lk)*r_lk*pow((1-r_rk),2)/(r*(1-r));
}
}

else if (l_geno==genotypes[2]){
if (r_geno==genotypes[0]){
genoprob[0]=pow(r_lk,2)*pow((1-r_rk),2)/(r*r);
genoprob[1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/(r*r);
genoprob[2]=pow((1-r_lk)*r_rk/r,2);
}
else if (r_geno==genotypes[1]){
genoprob[0]=r_lk*r_lk*(1-r_rk)*r_rk/((1-r)*r);
genoprob[1]=r_lk*(1-r_lk)*pow((1-r_rk),2)/((1-r)*r)+r_lk*(1-r_lk)*pow(r_rk,2)/((1-r)*r);
genoprob[2]=pow((1-r_lk),2)*r_rk*(1-r_rk)/((1-r)*r);
}
else if (r_geno==genotypes[2]){
genoprob[0]=pow(r_lk*r_rk/(1-r),2);
genoprob[1]=2*r_lk*(1-r_lk)*r_rk*(1-r_rk)/((1-r)*(1-r));
genoprob[2]=pow((1-r_lk)*(1-r_rk)/(1-r),2);
}
}

}



void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
    int i;
   int size=sizeof(int *); 
    *Geno = (int **) R_alloc(n_pos, size);
    (*Geno)[0] = geno;
    for(i=1; i< n_pos; i++) 
        (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

void reorg_genoprob(int n_genotype, int n_pos, int n_ind,
 double *genoprob,double ****Genoprob)
{int i,j;
double **a;
*Genoprob=(double ***)R_alloc(n_genotype,sizeof(double **));
//by R_alloc, you actually defined a pointer *Genoprob that points to a memory block with n_genotype element. However, you did notdefine any element in this memory block. so,the pointer a is actually an intermediate step,you cannot directly assigna specific address to (*Genoprob)[0], so, it has to be relied on the help of another variable.
a=(double **)R_alloc(n_pos*n_genotype, sizeof(double*));
(*Genoprob)[0]=a;
for (i=1;i<n_genotype;i++)
(*Genoprob)[i]=(*Genoprob)[i-1]+n_pos;
for (i=0;i<n_genotype;i++){
for (j=0;j<n_pos;j++)(*Genoprob)[i][j]=genoprob+i*n_ind*n_pos + j*n_ind;
//after doing this, you actually also changed the value
//that a is pointing to. (*Genoprob)[i][j] is a vector with
//n_ind element. 
}
}


//recode genotypes of A,H,B into 1,2,3.
void recodegeno_ChrToInt(int *intgeno, char **chrgeno,int *size)
{int i;
for (i=0;i< *size;i++){
if (chrgeno[i][0]=='N'||chrgeno[i][0]=='\0') intgeno[i]=NA_REAL;
if (*chrgeno[i]=='A') intgeno[i]=1;
else if (*chrgeno[i]=='H') intgeno[i]=2;
else if (*chrgeno[i]=='B') intgeno[i]=3;
}
}
//recode genotypes of 1,2,3, into A,B,H
void recodegeno_IntToChr(int *intgeno, char **chrgeno,int *size)
{int i;
for (i=0;i< *size;i++){
if (intgeno[i]==NA_INTEGER) *chrgeno[i]=NA_INTEGER;
if (intgeno[i]==1) *chrgeno[i]='A';
else if (intgeno[i]==2) *chrgeno[i]='H';
else if (intgeno[i]==3) *chrgeno[i]='B';
}
}

//recode one set of integer genotypes into another set of integer genotypes
void recodegeno_int2int(int *geno, int *codes,int *newcodes,int *size)
{int i;
for (i=0;i< *size;i++){
if (geno[i]==codes[0]) geno[i]=newcodes[0];
else if (geno[i]==codes[1]) geno[i]=newcodes[1];
else if (geno[i]==codes[2]) geno[i]=newcodes[2];
}
}

//recode integer genotypes to new gentoypes such as 1,2,3 according to parent genotypes.(identity by descent:IBD).If parents are not polymorphic, or if the one parent is NA,  code all the genotypes as NA.  
//Geno[n_pos][n_ind]
void RECODEGENO_IBD(int *parent1, int *parent2,int *geno, int *n_pos,int *n_ind,int *oriheter,int *newgeno)
{int **Geno,i,j;
reorg_geno(*n_ind,*n_pos,geno,&Geno);
for (i=0;i<*n_pos; i++){
if(NA_INTEGER!=parent1[i] && parent2[i]!=NA_INTEGER && parent1[i]!=parent2[i] ){
	for(j=0;j<*n_ind;j++){
if (Geno[i][j]==parent1[i])Geno[i][j]=newgeno[0];
else if (Geno[i][j]==parent2[i])Geno[i][j]=newgeno[2];
else if (Geno[i][j]==*oriheter)Geno[i][j]=newgeno[1];
}
}else {
	error("Cannot assign progeny genotype to non-polymorphic loci");

}
}
}


//recode IBD genotypes to their original genotypes according to parent information 
//Geno[n_pos][n_ind]
void recodegeno_IBDtoOri(int *parent1, int *parent2,int *geno, int *n_pos,int *n_ind,int *oriheter,int *IBDgeno)
{int **Geno,i,j;
reorg_geno(*n_ind,*n_pos,geno,&Geno);
for (i=0;i<*n_pos; i++){
if(NA_INTEGER!=parent1[i] && parent2[i]!=NA_INTEGER && parent1[i]!=parent2[i] ){
for(j=0;j<*n_ind;j++){
if (Geno[i][j]==IBDgeno[0])Geno[i][j]=parent1[i];
else if (Geno[i][j]==IBDgeno[2])Geno[i][j]=parent2[i];
else if (Geno[i][j]==IBDgeno[1])Geno[i][j]=*oriheter;
}
}
else {
error("Cannot assign progeny genotype to non-polymorphic loci");
}

}
} 


