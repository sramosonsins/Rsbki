//
//  cut_low_freq.c
//  Rsbki_functions
//
//  Created by Sebastian Ramos-Onsins on 18/1/23.
//  Copyright © 2023 CRAG. All rights reserved.
//

void cut_low_freqs(int **geno, long int geno_rows, int geno_cols, double cut_freq, int npops, int *popsize, long int *erased_rows, long int *n_erased_rows)
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
            if(sum_tot > 0. && (freq/sum_tot < cut_freq || freq/sum_tot > 1.0-cut_freq)) {
                erased_rows[*n_erased_rows] = i;
                *n_erased_rows += 1;
                k=npops;
            }
        }
    }
 }

