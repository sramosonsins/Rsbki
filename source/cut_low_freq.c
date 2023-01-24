//
//  cut_low_freq.c
//  Rsbki_functions
//
//  Created by Sebastian Ramos-Onsins on 18/1/23.
//  Copyright © 2023 CRAG. All rights reserved.
//

void cut_low_freqs(int **geno, long int geno_rows, int geno_cols, double cut_freq, int npops, int *popsize)
{
    int j,k,l;
    long int i;
    double freq,sum_tot;
    
    for(i=0;i<geno_rows;i++) {
        for(k=0;k<npops;k++) {
            l = 0; freq = sum_tot = 0.;
            if(k>0) l += popsize[k-1];
            for(j=0;j<popsize[k];j++) {
                if(geno[j+l][i] != 9) {
                    freq += (double)geno[j+l][i];
                    sum_tot += 2.;
                }
            }
            if(sum_tot > 0. && freq/sum_tot < cut_freq) {
                for(j=0;j<popsize[k];j++) {
                    geno[j+l][i] = 0; /*eliminate variant under cut_freq*/
                }
            }
        }
    }
 }

