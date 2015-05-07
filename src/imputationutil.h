#ifndef IMPUTATIONUTIL_H
#define IMPUTATIONUTIL_H
void map2r(double *map,double *r,int pos1, int pos2);
 void calc_genoprob_R(int *genotypes,int *l_geno,int *r_geno,double *r_lk,double *r_rk,double *r,double *genoprob);
void predict_prob(int *genotypes,double *l_prob, double *r_prob,double *rlk_rrk_r,double *genoprob);
void calc_genoprob_c(int *genotypes,int l_geno,int r_geno,double r_lk,double r_rk,double r,double *genoprob);
void reorg_geno(int n_pos, int n_ind, int *geno, int ***Geno);
void reorg_genoprob(int n_genotype, int n_pos, int n_ind,
  double *genoprob,double ****Genoprob);
void recodegeno_ChrToInt(int *intgeno, char **chrgeno,int *size);
void recodegeno_IntToChr(int *intgeno, char **chrgeno,int *size);
 void RECODEGENO_IBD(int *parent1, int *parent2,int *geno, int *n_pos,int *n_ind,int *oriheter,int *newgeno); 
#endif
