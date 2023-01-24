//
//  main.h
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//

//#ifndef main_h
#define main_h

#define TANG_SOFTW  "\nSoftware for calculating iES plus Rsb and iESk plus Rsbk (dividing per individual) statistics." \
"\nbased on initial statistics from Tang, Thornton & Stoneking, PloS Biology 2007." \
"\nversion 20230124"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ran1.h"

int read_row(FILE *plink_file,char *chr_name, double *lox, int **geno, int geno_cols, long int row);

void impute_genotypes(int **geno, long int geno_rows, int geno_cols);

void cut_low_freqs(int **geno, long int geno_rows, int geno_cols, double cut_freq, int npops, int *popsize);

void usage(void);

void calc_iES_iESk_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES, double **iESk, int pop_target);

//void calc_EHHS_Ikij_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *EHH, int **Ikij, int pop_target);

void calc_iRESda_iRESdak_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iESa, double *iESd, double **iRESk,int pop_target);

//void calc_EHHak_EHHdk_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max. double *iEHHa, double *iEHHd, double **iEHHak, double **iEHHdk, int pop_target);

//void calc_iES_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES);

//void calc_EHHS_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, double *EHH);

int compare_(const void *,const void *);

//#endif /* main_h */
