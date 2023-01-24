//
//  cut_low_freq.c
//  Rsbki_functions
//
//  Created by Sebastian Ramos-Onsins on 18/1/23.
//  Copyright © 2023 CRAG. All rights reserved.
//

void cut_low_freqs(int **geno, long int geno_rows, int geno_cols, double cut_freq) /*TO DO !!!!!*/
{
    int j;
    long int i;
    double freq,sum_tot;
    
    for(i=0;i<geno_rows;i++) {
        freq = sum_tot = 0.;
        for(j=0;j<geno_cols;j++) {
            if(geno[j][i] != 9) {
                freq += (double)geno[j][i];
                sum_tot += 2.;
            }
        }
        if(sum_tot > 0. && freq/sum_tot < cut_freq) {
            for(j=0;j<geno_cols;j++) {
                geno[j][i] = 0; /*eliminate variant under cut_freq*/
            }
        }
    }
 }

