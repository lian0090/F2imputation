//note: this is only for imputation for F2 population based on the conditional expectation

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "imputationutil.h"
void ends_impute(int typed, int pos,int ind,int *genotypes,double *map,int **Geno,double ***Genoprob); 
void know_prob(int pos, int ind,int *genotypes,int **Geno, double ***Genoprob );
void fillgeno(int pos,int ind,int *genotypes,double ***Genoprob,int **Geno);



void imputation(int *progenygeno,int *genotypes, int *n_genotype,int *n_pos, int *n_ind, double *map,double *genoprob,int *naind){
int i, j, **Geno;
double ***Genoprob; 
reorg_geno(*n_ind,*n_pos,progenygeno,&Geno);
reorg_genoprob(*n_genotype,*n_pos,*n_ind,genoprob,&Genoprob);
	/* you are actually editing *geno (the value both **Geno  and *geno are poiting to for  in the homozygotes loci,stop. Cannot assign progeny genotype to one particular genotyp.*/

//check whether there is enough marker,if not, stop.
//if(*n_pos<2){Rprintf("error:there is only one marker");abort();}
//check whether there is enough typed marker. (give a warning if a individual has no typed marker or only 1 typed marker.
int genotyped[*n_ind];
for (j=0;j<*n_ind;j++){
genotyped[j]=0;
for (i=0;i<*n_pos;i++){if(Geno[i][j]!=NA_INTEGER) genotyped[j]++;};

if (genotyped[j]<1) {
warning("there is no genotyped marker for  individual %d no imputation was done\n",j+1);
naind[j]=j+1;
}
else if (genotyped[j]==1){
warning("there is only one genotyped marker for individual %d\n",j+1);
//impute with only one marker
//find the typed marker position. 
int typed=0;
while (Geno[typed][j]==NA_INTEGER && typed<*n_pos) typed++;
know_prob(typed,j,genotypes,Geno,Genoprob);
//If the only typed marker is at the first marker
if(typed==0){
for (i=1;i<*n_pos;i++){
ends_impute(0,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);
}
}
//If the only typed marker is at the last marker
else if (typed==*n_pos-1){
for (i=*n_pos-2;i>-1;i--){
ends_impute(*n_pos-2,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);

}
}
//if the only typed marker is at the middle
else {
for (i=typed-1;i>-1;i--){
ends_impute(typed,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);}
for (i=typed+1;i<*n_pos;i++){
ends_impute(typed,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);
}
}
} 
}
// the end markers should be calculated differently//
int typed;
i=0;//the first marker
for (j=0;j<*n_ind;j++){
if (genotyped[j]>1){
 if(Geno[i][j]==NA_INTEGER){
 typed=1;//find the nearest typed markers 
while(Geno[typed][j]==NA_INTEGER) typed++;
ends_impute(typed,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);
}
else know_prob(i,j,genotypes,Geno,Genoprob);
}
}
i=*n_pos-1;
for (j=0;j<*n_ind;j++){
if (genotyped[j]>1){
if(Geno[i][j]==NA_INTEGER){
 typed=*n_pos-2;
while(Geno[typed][j]==NA_INTEGER) typed--;
ends_impute(typed,i,j,genotypes,map,Geno,Genoprob);
fillgeno(i,j,genotypes,Genoprob,Geno);
}
else know_prob(i,j,genotypes,Geno,Genoprob);
}
}
// the middle markers
int l_typed,r_typed;
double d,d_lk,d_rk;//left and right distance to the target marker.
double r,r_lk,r_rk;//left and right recombination to the traget marker.

//give the typed markers a probability first, so that, when a particular marker needs the probability of its right marker to calculate its genotype probability, it will have that information.probability of NA markers are calcualted based on the nearest typed markers.
for (i=1;i<*n_pos-1;i++)
{for (j=0;j<*n_ind;j++){
if (genotyped[j]>1){
if (Geno[i][j]!=NA_INTEGER){know_prob(i,j,genotypes,Geno,Genoprob);}
}
}
}
for (i=1;i<*n_pos-1;i++){
for (j=0;j<*n_ind;j++){
if (genotyped[j]>1){
if (Geno[i][j]==NA_INTEGER){
l_typed=i-1;
r_typed=i+1;
while(Geno[l_typed][j]==NA_INTEGER)l_typed--;
while(Geno[r_typed][j]==NA_INTEGER)r_typed++;
d=fabs(map[l_typed]-map[r_typed])/100;
d_lk=fabs(map[i]-map[l_typed])/100;
d_rk=fabs(map[i]-map[r_typed])/100;
r=(1-exp(-2*d))/2;
r_lk=(1-exp(-2*d_lk))/2;
r_rk=(1-exp(-2*d_rk))/2;
double l_prob[3],r_prob[3],newgenoprob[3];int k;
for(k=0;k<3;k++){
l_prob[k]=Genoprob[k][l_typed][j];
r_prob[k]=Genoprob[k][r_typed][j];
newgenoprob[k]=Genoprob[k][i][j];};
double rlk_rrk_r[3]={r_lk,r_rk,r};
predict_prob(genotypes,l_prob, r_prob,rlk_rrk_r, newgenoprob);

for(k=0;k<3;k++){Genoprob[k][i][j]=newgenoprob[k];}
fillgeno(i,j,genotypes,Genoprob,Geno);
}
}
}
}
}







//a function to calculate the genotype probability of the end markers
void ends_impute(int typed, int pos,int ind,int *genotypes,double *map,int **Geno,double ***Genoprob){
double d; //genetic distance
double r; //recombination frequency
int i=pos;int j=ind;
d=fabs(map[typed]-map[i])/100;//use Morgan rather than cM
r=(1-exp(-2*d))/2;
if (Geno[typed][j]==genotypes[0])//homozygotes for parent1
{Genoprob[0][i][j]=(1-r)*(1-r);
Genoprob[1][i][j]=2*r*(1-r);
Genoprob[2][i][j]=r*r;}
else if (Geno[typed][j]==genotypes[1])//heterozytoes
{
Genoprob[0][i][j]=r*(1-r);
Genoprob[1][i][j]=r*r + (1-r)*(1-r);
Genoprob[2][i][j]=r*(1-r);
}
else if (Geno[typed][j]==genotypes[2]) //homozygotes for parent2
{Genoprob[0][i][j]=r*r;
Genoprob[1][i][j]=2*r*(1-r);
Genoprob[2][i][j]=(1-r)*(1-r);
}
}

//give the genotype probability of know genotypes.
void know_prob(int pos, int ind,int *genotypes,int **Geno, double ***Genoprob )
{int i=pos;int j=ind;
if (Geno[i][j]==genotypes[0]){
Genoprob[0][i][j]=1.0;
Genoprob[1][i][j]=0.0;
Genoprob[2][i][j]=0.0;}
else if (Geno[i][j]==genotypes[1]){
Genoprob[0][i][j]=0.0;
Genoprob[1][i][j]=1.0;
Genoprob[2][i][j]=0.0;}
else if (Geno[i][j]==genotypes[2]){
Genoprob[0][i][j]=0.0;
Genoprob[1][i][j]=0.0;
Genoprob[2][i][j]=1.0;
}
}


void fillgeno(int pos,int ind,int *genotypes,double ***Genoprob,int **Geno)
{
int i=pos;int j=ind; 
int maxgeno=0;
if (Genoprob[maxgeno][i][j]<Genoprob[1][i][j])
maxgeno=1;
if (Genoprob[maxgeno][i][j]<Genoprob[2][i][j])
maxgeno=2;
Geno[i][j]=genotypes[maxgeno];
}






