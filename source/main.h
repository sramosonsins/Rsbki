//
//  main.h
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//

//#ifndef main_h
#define main_h

#define TANG_SOFTW  "\nSoftware for calculating iES plus Rsb and iESk plus Rsbk (per individual) statistics." \
"\nfollowing Tang, Thornton & Stoneking, PloS Biology 2007." \
"\nverion 20190718"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ran1.h"

int read_row(FILE *plink_file,char *chr_name, double *lox, int **geno, int geno_cols, long int row);

void impute_genotypes(int **geno, long int geno_rows, int geno_cols);

void usage(void);

void calc_iES_iESk_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES, double **iESk, int pop_target);

void calc_EHHS_Ikij_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *EHH, int **Ikij, int pop_target);

//void calc_iES_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES);

//void calc_EHHS_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, double *EHH);

int compare_(const void *,const void *);

//#endif /* main_h */
