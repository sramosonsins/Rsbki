//
//  iRES_functions.c
//  Rsbki_functions
//
//  Created by Sebastian E. Ramos Onsins on 25/07/2019.
//  Copyright © 2019 CRAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

void calc_iRESda_iRESdak_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iESa, double *iESd, double **iRESk,int pop_target)
{
    void calc_EHHak_EHHdk_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *iEHHa, double *iEHHd, double **iEHHak, double **iEHHdk, int pop_target);
    
    double *x;
    double *EHHa;
    double *EHHd;
    //double iESka;
    double iESkd;
    double **iEHHak;
    double **iEHHdk;
    
    long int L;
    long int i,j,k;
    long int min,max;
    
    L = *geno_rows;
    x = (double *)calloc(L-1,sizeof(double));

    iEHHak = 0;
    iEHHdk = 0;
    if(pop_target) {
        iEHHak = (double **)calloc(*geno_cols,sizeof(double *));
        iEHHdk = (double **)calloc(*geno_cols,sizeof(double *));
        for(k=0;k<*geno_cols;k++) {
            iEHHak[k] = (double *)calloc(L,sizeof(double));
            iEHHdk[k] = (double *)calloc(L,sizeof(double));
        }
    }
    EHHa = (double *)calloc(L,sizeof(double));
    EHHd = (double *)calloc(L,sizeof(double));

    for(i=0;i<L-1;i++) {x[i] = lox[i+1] - lox[i];}
    for(i=0;i<L;i++) {
        calc_EHHak_EHHdk_pos(&i,geno,geno_rows, geno_cols,thresh,&min,&max, EHHa,EHHd,iEHHak,iEHHdk,pop_target );
        for(j=min;j<max;j++) { //max value at pos=i
            iESa[i] += (EHHa[j+1] + EHHa[j]) * x[j] / 2.0;
            iESd[i] += (EHHd[j+1] + EHHd[j]) * x[j] / 2.0;
            EHHa[j] = EHHd[j] = 0.; //init for next
        }
        EHHa[max] = EHHd[max] = 0.; //init for next
        if(pop_target) {
            //double sumiESka=0.;
            /**/double sumiESkd=0.;
            for(k=0;k<*geno_cols;k++) {
                //iESka = 0.0;
                iESkd = 0.0;
                for(j=min;j<max;j++) {
                    //iESka += (iEHHak[k][j+1] + iEHHak[k][j]) * x[j] / 2.0;
                    iESkd += (iEHHdk[k][j+1] + iEHHdk[k][j]) * x[j] / 2.0;
                    //iEHHak[k][j] = 0;
                    iEHHdk[k][j] = 0.; //init for next
                }
                //iEHHak[k][max] = 0;
                iEHHdk[k][max] = 0.; //init for next
                //sumiESka += iESka;
                /**/sumiESkd += iESkd;
                //if(iESka)
                //if(iESka + iESkd)
                if(iESa[i])
                    iRESk[k][i] = iESkd / iESa[i];  //(iESka + iESkd); // iESkd / iESka;
                else iRESk[k][i] = -1;
            }
            //if(!(sumiESka >= (iESa[i] - 0.1) && sumiESka <= (iESa[i] + 0.1)))
            //    printf("iESa and sumiESka are different\n");
            if(!(sumiESkd >= (iESd[i] - 0.1) && sumiESkd <= (iESd[i] + 0.1)))
                printf("iESd and sumiESkd are different\n");
        }
    }
    free(x);
    free(EHHa);
    free(EHHd);
    if(pop_target) {
        for(k=0;k<*geno_cols;k++) {free(iEHHak[k]); free(iEHHdk[k]);}
        free(iEHHak);
        free(iEHHdk);
    }
    return;
}

void calc_EHHak_EHHdk_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *iEHHa, double *iEHHd, double **iEHHak, double **iEHHdk, int pop_target)
{
    long int M;
    long int pos;
    long int j,k;
    long int Iia,Iid,Ija,Ijd;
    long int *Ica,*Icd;
    double cur_EHHa,cur_EHHd;
    
    pos = *i;
    M = *geno_rows;
    Iia = Iid = 0;
    for(k=0;k<*geno_cols;k++) {
        Iia += (geno[k][pos]==0);
        Iid += (geno[k][pos]==2);
    }
    if(Iia == 0 || Iid == 0) {
        iEHHa[pos] = -1;
        iEHHd[pos] = -1;
        *min = *max = pos;
        return;
    }
    
    iEHHa[pos] = 1.0;
    iEHHd[pos] = 1.0;

    if(pop_target) {
        for(k=0;k<*geno_cols;k++) {
            iEHHak[k][pos] = (double)(geno[k][pos]==0)/(double)Iia;
            iEHHdk[k][pos] = (double)(geno[k][pos]==2)/(double)Iid;
        }
    }
    
    Ica = (long int *)calloc(*geno_cols,sizeof(long int));
    Icd = (long int *)calloc(*geno_cols,sizeof(long int));

    //left-flank
    for(k=0;k<*geno_cols;k++) {
        Ica[k] = !(geno[k][pos]==0);
        Icd[k] = !(geno[k][pos]==2);
    }
    j = pos - 1;
    while( j >= 0) {
        Ija = Ijd = 0;
        for(k=0;k<*geno_cols;k++) {
            Ica[k] += !(geno[k][j]==0);
            Icd[k] += !(geno[k][j]==2);
            Ija += (Ica[k]==0);
            Ijd += (Icd[k]==0);
        }
        cur_EHHa = (double)Ija / (double)Iia;
        cur_EHHd = (double)Ijd / (double)Iid;
        if(cur_EHHa < *thresh && cur_EHHd < *thresh) {
            break;
        }
        else {
            if(pop_target) {
                for(k=0;k<*geno_cols;k++) {
                    iEHHak[k][j] = (double)(Ica[k]==0)/(double)Iia;
                    iEHHdk[k][j] = (double)(Icd[k]==0)/(double)Iid;
                }
            }
            iEHHa[j] = cur_EHHa;
            iEHHd[j] = cur_EHHd;
        }
        j -= 1;
    }
    *min = j+1;
    
    //right-flank
    for(k=0;k<*geno_cols;k++) {
        Ica[k] = !(geno[k][pos]==0);
        Icd[k] = !(geno[k][pos]==2);
    }
    j = pos + 1;
    while( j < M) {
        Ija = Ijd = 0;
        for(k=0;k<*geno_cols;k++) {
            Ica[k] += !(geno[k][j]==0);
            Icd[k] += !(geno[k][j]==2);
            Ija += (Ica[k]==0);
            Ijd += (Icd[k]==0);
        }
        cur_EHHa = (double)Ija / (double)Iia;
        cur_EHHd = (double)Ijd / (double)Iid;
        if(cur_EHHa < *thresh && cur_EHHd < *thresh) {
            break;
        }
        else {
            if(pop_target) {
                for(k=0;k<*geno_cols;k++) {
                    iEHHak[k][j] = (double)(Ica[k]==0)/(double)Iia;
                    iEHHdk[k][j] = (double)(Icd[k]==0)/(double)Iid;
                }
            }
            iEHHa[j] = cur_EHHa;
            iEHHd[j] = cur_EHHd;
        }
        j += 1;
    }
    *max = j-1;
    
    free(Ica);
    free(Icd);

    return;
}
