//
//  imputation.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 31/01/2019.
//
#include "ran1.h"

void impute_genotypes(int **geno, long int geno_rows, int geno_cols)
{
    int pj[3][3];
    //int pjt;
    double prob_pj[3];
    double prob_pjt;
    long int i,i1;
    int j,k,l,m;
    int *is_na;
    double ran_val;
    
    is_na = (int *)calloc(geno_cols,sizeof(int));
    
    for(i=0;i<geno_rows;i++) {
        //check if there are missing values in this position
        m = 0;
        for(j=0;j<geno_cols;j++) {
            if(geno[j][i]==9) {
                m=1;
                break;
            }
        }
        if(m==0)  continue;
        
        //init arrays
        for(k=0;k<3;k++) {
            for(l=0;l<3;l++) {
                pj[k][l] = 0;
            }
        }
        //estimate probs
        k=0;
        if(i==0) i1=0;
        else i1 = i-1;
        for(j=0;j<geno_cols;j++) {
            //vector and matrix to estimate genotype
            pj[0][0] += (geno[j][i1]==0 && geno[j][i]==0);
            pj[0][1] += (geno[j][i1]==0 && geno[j][i]==1);
            pj[0][2] += (geno[j][i1]==0 && geno[j][i]==2);
            pj[1][0] += (geno[j][i1]==1 && geno[j][i]==0);
            pj[1][1] += (geno[j][i1]==1 && geno[j][i]==1);
            pj[1][2] += (geno[j][i1]==1 && geno[j][i]==2);
            pj[2][0] += (geno[j][i1]==2 && geno[j][i]==0);
            pj[2][1] += (geno[j][i1]==2 && geno[j][i]==1);
            pj[2][2] += (geno[j][i1]==2 && geno[j][i]==2);
            //the position of missing (j) and the number of missing (k)
            if(geno[j][i]==9) is_na[k++] = j;
            //the real frequency at position i
        }
        //pjt = geno_cols - k;
        
        //assign genotypes
        if(pj[0][0] + pj[0][1] + pj[0][2] + pj[1][0] + pj[1][1] + pj[1][2] + pj[2][0] + pj[2][1] + pj[2][2] == 0) {
            //assign a genotype 0 to all individuals
            geno[is_na[j]][i] = 0;
        }
        else {
            ran_val = ran1();
            for(j=0;j<k;j++) {
                if(geno[k][i1] == 9) {
                    prob_pj[0] = pj[0][0] + pj[1][0] + pj[2][0];
                    prob_pj[1] = pj[0][1] + pj[1][1] + pj[2][1];
                    prob_pj[2] = pj[0][2] + pj[1][2] + pj[2][2];
                }
                else {
                    prob_pj[0] = pj[geno[k][i1]][0];
                    prob_pj[1] = pj[geno[k][i1]][1];
                    prob_pj[2] = pj[geno[k][i1]][2];
                }
                prob_pjt = prob_pj[0] + prob_pj[1] + prob_pj[2];
                if((double)prob_pj[0]/(double)prob_pjt > ran_val) {
                    geno[is_na[j]][i] = 0;
                }
                else {
                    if(((double)prob_pj[0]+(double)prob_pj[1])/(double)prob_pjt > ran_val) {
                        geno[is_na[j]][i] = 1;
                    }
                    else {
                        geno[is_na[j]][i] = 2;
                     }
                }
            }
        }
        for(j=0;j<geno_cols;j++) is_na[j] = 0;
    }
    free(is_na);
}

