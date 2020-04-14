//
//  main.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//  Copyright © 2019 Sebastian Ramos-Onsins. All rights reserved.
//

#include "main.h"

int main(int arg, const char *argv[])
{
    double *lox;
    int **geno;
    long int geno_rows;
    int geno_cols;
    int **pop_geno;
    int pop_geno_cols;
    double thresh;
    double **all_iES;
    double **pop_iESk;//new code
    double **all_iESa;//new code
    double **all_iESd;//new code
    double **pop_iRESk;//new code
    long int L;
    int N;
    long int i,h;
    int j,k,l,ll,jj;
    int npops;
    int *popsize;
    int popcum;
    int popn_target; //new code
    
    double *mean_res; //new code
    double *median_res; //new code
    double *med_res; //new code
    double *sd_res; //new code
    long int *nit_res;//new code
    double *mean;
    double *median;
    double *med;
    double *sd;
    double **all_iRES;//new code
    double **all_Rsb;
    double ***pop_Rsbk;//new code
    long int *nit;
    char **pop_name;

    FILE *plink_file = 0;
    FILE *genot_file = 0;
    FILE *results_file = 0;
    int c;
    long int row;
    char chr_name[11]; //the file must contain the same ID for all positions!
    char file_out[1024];
    char file_in[1024];
    int argc = 1;

    if(arg < 9) {
        printf("\nError introducing arguments\n");
        usage();
        exit(1);
    }
    
    printf(TANG_SOFTW);
    printf("\n\nRsbki ");
    while(argc < arg) {
        printf("%s ",argv[argc]);
        argc++;
    }
    printf("\n");
    
    //define variables and arrays:
    strcpy( file_in, argv[1]);
    geno_rows = L = atol(argv[2]);
    geno_cols = N = atoi(argv[3]);
    thresh = atof(argv[4]);
    init_seed1(atol(argv[5]));
    npops = atoi(argv[6]);
    popn_target = atoi(argv[7]) - 1; //new code: note added 1 to next indexes in argv
    if(popn_target<0) popn_target=0;
    popsize = (int *)calloc(npops,sizeof(int));
    //for(i=0;i<npops;i++) {
    //    if(i!=npops-1) popsize[i] = atoi(argv[7+1+i+1]) - atoi(argv[7+1+i]);
    //    else popsize[i] = (N+2+1) - atoi(argv[7+1+i]);
    //}
    for(i=0;i<npops;i++) {
        popsize[i] = atoi(argv[7+1+i]);
    }
    pop_name = (char **)calloc(npops,sizeof(char *));
    for(i=0;i<npops;i++) {
        pop_name[i] = (char *)calloc(100,sizeof(char));
        strcpy(pop_name[i], argv[7+1+npops+i]);
    }
 
    if (!(plink_file = fopen(file_in,"r"))) {
        printf("Error reading the input file %s\n",file_in);
        exit(1);
    }
    
    //geno is transposed to facilitate the analysis
    geno = (int **)calloc(N,sizeof(int *));
    for(j=0;j<N;j++) {geno[j] = (int *)calloc(L,sizeof(int));}
    lox = (double *)calloc(L,sizeof(double));
    all_iES = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iES[i] = (double *)calloc(L,sizeof(double));}
    //new code
    pop_iESk = (double **)calloc(popsize[popn_target],sizeof(double *));
    for(i=0;i<popsize[popn_target];i++) {pop_iESk[i] = (double *)calloc(L,sizeof(double));}
    //end
    //new code
    all_iESa = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iESa[i] = (double *)calloc(L,sizeof(double));}
    all_iESd = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iESd[i] = (double *)calloc(L,sizeof(double));}
    all_iRES = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iRES[i] = (double *)calloc(L,sizeof(double));}
    pop_iRESk = (double **)calloc(popsize[popn_target],sizeof(double *));
    for(i=0;i<popsize[popn_target];i++) {pop_iRESk[i] = (double *)calloc(L,sizeof(double));}

    mean_res = (double *)calloc(npops,sizeof(double));
    median_res = (double *)calloc(npops,sizeof(double));
    sd_res = (double *)calloc(npops,sizeof(double));
    nit_res = (long int *)calloc(npops,sizeof(long int));
    //end

    mean = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    median = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    sd = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    nit = (long int *)calloc(npops*(npops-1)/2,sizeof(long int));

    all_Rsb = (double **)calloc(npops*(npops-1)/2,sizeof(double *));
    for(i=0;i<npops*(npops-1)/2;i++) {all_Rsb[i] = (double *)calloc(L,sizeof(double));}
    //new code
    pop_Rsbk = 0;
    if(popn_target>=0) {
        pop_Rsbk = (double ***)calloc(npops,sizeof(double **));
        for(j=0;j<npops-1;j++) {
            pop_Rsbk[j] = (double **)calloc(popsize[popn_target],sizeof(double *));
            for(i=0;i<popsize[popn_target];i++) {
                pop_Rsbk[j][i] = (double *)calloc(L,sizeof(double));
            }
        }
    }
    //end

    //read input data: skip header
    printf("\nReading input file..."); fflush(stdout);
    while((c=getc(plink_file))!='\n' && c!='\r');
    row = 0;
    for(row=0;row<L;row++) {
        if((c = read_row(plink_file,chr_name,lox,geno,geno_cols,row))==0) {
            printf("\nError: input file has less rows or more cols than defined. nrows: %ld\n",row+1);
            exit(1);
        }
    }
    fclose(plink_file);
    
    //imputation
    printf("\nimputing missing values (previously assigned as genotype=9)...\n"); fflush(stdout);
    impute_genotypes(geno, L, geno_cols);

    //writing imputed file
    memset(file_out, '\0', 1024);
    strcat(file_out,file_in);
    strcat(file_out,"_imputed.txt\0");
    printf("\nWriting imputed genotype file %s...\n",file_out);
        //header
    if (!(genot_file = fopen(file_out,"w"))) {
        printf("Error writing the imputed file %s\n",file_in);
        exit(1);
    }
    fprintf(genot_file,"CHR\tPOS");
    for(j=0;j<npops;j++) {
        for(k=0;k<popsize[j];k++) {
            fprintf(genot_file,"\t%s%d_IND%d",pop_name[j],j+1,k+1);
        }
    }
    fprintf(genot_file,"\n");
        //genotypes
    for(i=0;i<L;i++) {
        fprintf(genot_file,"%s\t%f",chr_name,lox[i]);
        for(j=0;j<N;j++) {
            fprintf(genot_file,"\t%d",geno[j][i]);
        }
        fprintf(genot_file,"\n");
    }
    fclose(genot_file);
    
    //calculation of iRESda and iRESdak
    popcum = 0;
    for(j=0;j<npops;j++) {
        if(j==popn_target) printf("\nComputing iRESdak and iRESda for pop %d of %d ...",j+1,npops);
        else printf("\nComputing iRESda for pop %d of %d ...",j+1,npops);
        fflush(stdout);
        pop_geno_cols = popsize[j];
        pop_geno = &geno[popcum];
        calc_iRESda_iRESdak_slow(pop_geno, lox, &geno_rows, &pop_geno_cols, &thresh, all_iESa[j], all_iESd[j], pop_iRESk,(j==popn_target));
        popcum += popsize[j];
    }

    //calculation iRES: later calculate median, mean and sd:
    printf("\nComputing iRESda and iRESdak..."); fflush(stdout);
    for(i=0;i<L;i++) {
        for(j=0;j<npops;j++) {
            if(all_iESa[j][i] > 0.0 && all_iESd[j][i] > 0.0) {
                all_iRES[j][i] = log(all_iESd[j][i])-log(all_iESa[j][i]);
                mean_res[j] = mean_res[j] + all_iRES[j][i];
                sd_res[j] = sd_res[j] + all_iRES[j][i] * all_iRES[j][i];
                nit_res[j] += 1;
            }
            else {
                all_iRES[j][i] = 1234567890;
            }
        }
    }
    // median ...
    med_res = (double *)calloc(L,sizeof(double));
    for(j=0;j<npops;j++) {
        for(i=0;i<L;i++) {
            if(all_iRES[j][i] != 1234567890) {
                med_res[i] = all_iRES[j][i];
            }
        }
        qsort(med_res,nit_res[j],sizeof(double),compare_);
        if((double)nit_res[j]/2.0 == nit_res[j]/2.0) {
            median_res[j] = (med_res[(long int)((double)nit_res[j]/2.0 - 1.0)] +
                         med_res[(long int)((double)nit_res[j]/2.0)]) / 2.0;
        }
        else median_res[j] = med_res[nit_res[j]/2];
    }
    free(med_res);
    // mean, sd ...
    for(j=0;j<npops;j++) {
            mean_res[j] = mean_res[j]/(double)nit_res[j];
            sd_res[j] = sqrt(sd_res[j]/((double)nit_res[j]-1.0) - mean_res[j]*mean_res[j] * (double)nit_res[j]/((double)nit_res[j]-1.0)); //smapled sd
    }

    //writing output file
    file_out[0] = '\0';
    strcat(file_out,file_in);
    strcat(file_out,"_Results_iESd_iESa_iRES_pop.txt\0");
    printf("\nWriting iESa and iESd and iRES results in the output file %s...",file_out);
    fflush(stdout);
    
    if (!(results_file = fopen(file_out,"w"))) {
        printf("Error writing the input file %s\n",file_out);
        exit(1);
    }
    //header
    fprintf(results_file,"Position\t");
    for(j=0;j<npops;j++) fprintf(results_file,"iESd_%s\tiESa_%s\t",pop_name[j],pop_name[j]);
    for(j=0;j<npops;j++) fprintf(results_file,"log(iESd_%s)\tlog(iESa_%s)\t",pop_name[j],pop_name[j]);
    for(j=0;j<npops;j++) fprintf(results_file,"log_iRESda(%s)\t",pop_name[j]);
    for(j=0;j<npops;j++) fprintf(results_file,"log_iRESdaN(%s)\t",pop_name[j]);

    fprintf(results_file,"\n");
    //data
    for(i=0;i<L;i++) {
        fprintf(results_file,"%f\t",lox[i]);
        for(j=0;j<npops;j++) {
            if(all_iESa[j][i] > 0.0 && all_iESd[j][i] > 0.0)
                fprintf(results_file,"%f\t%f\t",all_iESd[j][i],all_iESa[j][i]);
            else fprintf(results_file,"NA\tNA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iESa[j][i] > 0.0 && all_iESd[j][i] > 0.0)
                fprintf(results_file,"%f\t%f\t",log(all_iESd[j][i]),log(all_iESa[j][i]));
            else fprintf(results_file,"NA\tNA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iRES[j][i] != 1234567890) fprintf(results_file,"%f\t",all_iRES[j][i]);
            else fprintf(results_file,"NA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iRES[j][i] != 1234567890) fprintf(results_file,"%f\t",(all_iRES[j][i]-median_res[j])/sd_res[j]);
            else fprintf(results_file,"NA\t");
        }
        fprintf(results_file,"\n");
    }
    fclose(results_file);
    
    if(popn_target >= 0) {
        //output iRESk
        file_out[0] = '\0';
        strcat(file_out,file_in);
        strcat(file_out,"_Results_iRESdak_pop\0");
        strcat(file_out,pop_name[popn_target]);
        strcat(file_out,".txt\0");
        printf("\nWriting iRESk results (not log, no normalized) in the output file %s...",file_out);
        fflush(stdout);
        if (!(results_file = fopen(file_out,"w"))) {
            printf("Error writing the input file %s\n",file_out);
            exit(1);
        }
        //header
        fprintf(results_file,"Position\t");
        for(k=0;k<popsize[popn_target];k++) fprintf(results_file,"iRESk_ind%d\t",k+1);
        fprintf(results_file,"\n");
        //data
        for(i=0;i<L;i++) {
            fprintf(results_file,"%f\t",lox[i]);
            for(k=0;k<popsize[popn_target];k++) {
                fprintf(results_file,"%f\t",pop_iRESk[k][i]);
            }
            fprintf(results_file,"\n");
        }
        fclose(results_file);
    }

    //converting value 2 (homozygotes alternative) to 0 (homozygotes unpolarized):
    printf("\nConverting all homozygotes to value 0..."); fflush(stdout);
    for(j=0;j<N;j++) {
        for(i=0;i<L;i++) {
            if(geno[j][i]==2) geno[j][i]=0;
        }
    }
    
    //calculation iES and iESk
    popcum = 0;
    for(j=0;j<npops;j++) {
        //new code
        if(j==popn_target) printf("\nComputing iES and iESk for pop %d of %d ...",j+1,npops);
        else printf("\nComputing iES for pop %d of %d ...",j+1,npops);
        fflush(stdout);
        //end
        pop_geno_cols = popsize[j];
        pop_geno = &geno[popcum];
        calc_iES_iESk_slow(pop_geno, lox, &geno_rows, &pop_geno_cols, &thresh, all_iES[j],pop_iESk,(j==popn_target)); //new code
        popcum += popsize[j];
    }
    
    //calculation Rsb: later calculate median, mean and sd:
    printf("\nComputing Rsb and Rsbk..."); fflush(stdout);
    for(i=0;i<L;i++) {
        l=0;ll=0;
        for(j=0;j<npops-1;j++) {
            for(h=j+1;h<npops;h++) {
                if(all_iES[j][i] > 0.0 && all_iES[h][i] > 0.0) {
                    all_Rsb[l][i] = log(all_iES[j][i])-log(all_iES[h][i]);
                    mean[l] = mean[l] + all_Rsb[l][i];
                    sd[l] = sd[l] + all_Rsb[l][i] * all_Rsb[l][i];
                    nit[l] += 1;
                    //new code
                    if((j==popn_target && h != popn_target) ||
                       (h==popn_target && j != popn_target)) {
                        if(j!=popn_target) jj=j;
                        else jj=(int)h;
                        for(k=0;k<popsize[popn_target];k++) {
                            pop_Rsbk[ll][k][i] = (pop_iESk[k][i])/(all_iES[jj][i]);
                        }
                        ll++;
                    }
                    //end
                }
                else {
                    all_Rsb[l][i] = 1234567890;
                }
                l++;
            }
        }
    }
    // median ...
    l=0;
    med = (double *)calloc(L,sizeof(double));
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            h = 0;
            for(i=0;i<L;i++) {
                if(all_Rsb[l][i] != 1234567890) {
                    med[h++] = all_Rsb[l][i];
                }
            }
            qsort(med,nit[l],sizeof(double),compare_);
            if((double)nit[l]/2.0 == nit[l]/2.0) {
                median[l] = (med[(long int)((double)nit[l]/2.0 - 1.0)] +
                             med[(long int)((double)nit[l]/2.0)]) / 2.0;
            }
            else median[l] = med[nit[l]/2];
            l++;
        }
    }
    free(med);
    // mean, sd ...
    l=0;
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            mean[l] = mean[l]/(double)nit[l];
            sd[l] = sqrt(sd[l]/((double)nit[l]-1.0) - mean[l]*mean[l] * (double)nit[l]/((double)nit[l]-1.0)); //smapled sd
            l++;
        }
    }
    
    //writing output file
    file_out[0] = '\0';
    strcat(file_out,file_in);
    strcat(file_out,"_Results_Tang.txt\0");
    printf("\nWriting iES and Rsb results in the output file %s...",file_out);
    fflush(stdout);

    if (!(results_file = fopen(file_out,"w"))) {
        printf("Error writing the input file %s\n",file_out);
        exit(1);
    }
    //header
    fprintf(results_file,"Position\t");
    for(j=0;j<npops;j++) fprintf(results_file,"iES_%s\t",pop_name[j]);
    for(j=0;j<npops;j++) fprintf(results_file,"log(iES_%s)\t",pop_name[j]);
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            fprintf(results_file,"log_Rsb(%s/%s)\t",pop_name[j],pop_name[k]);
        }
    }
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            fprintf(results_file,"log_RsbN(%s/%s)\t",pop_name[j],pop_name[k]);
        }
    }
    fprintf(results_file,"\n");
    //data
    for(i=0;i<L;i++) {
        fprintf(results_file,"%f\t",lox[i]);
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(results_file,"%f\t",all_iES[j][i]);
            else fprintf(results_file,"NA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(results_file,"%f\t",log(all_iES[j][i]));
            else fprintf(results_file,"NA\t");
        }

        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(results_file,"%f\t",all_Rsb[l][i]);
                else fprintf(results_file,"NA\t");
                l++;
            }
        }
        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(results_file,"%f\t",(all_Rsb[l][i]-median[l])/sd[l]);
                else fprintf(results_file,"NA\t");
                l++;
            }
        }
        fprintf(results_file,"\n");
    }
    fclose(results_file);
    
    //new code
    if(popn_target >= 0) {
        //output iESk
        file_out[0] = '\0';
        strcat(file_out,file_in);
        strcat(file_out,"_Results_iESk_pop\0");
        strcat(file_out,pop_name[popn_target]);
        strcat(file_out,".txt\0");
        printf("\nWriting iESk results in the output file %s...",file_out);
        fflush(stdout);
        if (!(results_file = fopen(file_out,"w"))) {
            printf("Error writing the input file %s\n",file_out);
            exit(1);
        }
        //header
        fprintf(results_file,"Position\t");
        for(k=0;k<popsize[popn_target];k++) fprintf(results_file,"iESk_ind%d\t",k+1);
        fprintf(results_file,"\n");
        //data
        for(i=0;i<L;i++) {
            fprintf(results_file,"%f\t",lox[i]);
            for(k=0;k<popsize[popn_target];k++) {
                /*if(pop_iESk[k][i] > 0.0)*/ fprintf(results_file,"%f\t",pop_iESk[k][i]);
                //else fprintf(results_file,"0.0\t");
            }
            fprintf(results_file,"\n");
        }
        fclose(results_file);
        
        //output Rsbk
        l=0;ll=0;
        for(j=0;j<npops-1;j++) {
            for(h=j+1;h<npops;h++) {
                if((j==popn_target && k != popn_target) ||
                   (h==popn_target && j != popn_target)) {
                    file_out[0] = '\0';
                    strcat(file_out,file_in);
                    strcat(file_out,"_Results_Rsbk_TARGETpop\0");
                    strcat(file_out,pop_name[popn_target]);
                    strcat(file_out,"_versus_REFpop\0");
                    if(j!=popn_target) jj=j;
                    else jj=(int)h;
                    strcat(file_out,pop_name[jj]);
                    strcat(file_out,".txt\0");
                    printf("\nWriting Rsbk results in the output file %s...",file_out);
                    fflush(stdout);
                    if (!(results_file = fopen(file_out,"w"))) {
                        printf("Error writing the input file %s\n",file_out);
                        exit(1);
                    }
                    //header
                    fprintf(results_file,"Position\t");
                    for(k=0;k<popsize[popn_target];k++) {
                        fprintf(results_file,"Rsbk[%d]\t",k);
                    }
                    fprintf(results_file,"\n");
                    //data
                    for(i=0;i<L;i++) {
                        fprintf(results_file,"%f\t",lox[i]);
                        for(k=0;k<popsize[popn_target];k++) {
                            /*if(pop_Rsbk[ll][k][i] > 0.0)*/ fprintf(results_file,"%f\t",pop_Rsbk[ll][k][i]);
                            //else fprintf(results_file,"NA\t");
                        }
                        fprintf(results_file,"\n");
                    }
                    ll++;
                    fclose(results_file);
                }
                l++;
            }
        }
    }
    //freed
    if(popn_target>=0) {
        for(j=0;j<npops-1;j++) {
            for(i=0;i<popsize[popn_target];i++) {
                free(pop_Rsbk[j][i]);
            }
            free(pop_Rsbk[j]);
        }
        free(pop_Rsbk);
    }
    for(i=0;i<npops*(npops-1)/2;i++) {
        free(all_Rsb[i]);
    }
    free(all_Rsb);
    free(mean);
    free(sd);
    free(nit);
    for(i=0;i<popsize[popn_target];i++) {
        free(pop_iESk[i]);
        free(pop_iRESk[i]);
    }
    free(pop_iESk);
    free(pop_iRESk);
    for(i=0;i<npops;i++) {
        free(all_iESa[i]);
        free(all_iESd[i]);
        free(all_iES[i]);
    }
    free(mean_res);
    free(sd_res);
    free(nit_res);
    free(all_iESa);
    free(all_iESd);
    free(all_iES);
    for(j=0;j<N;j++) {
        free(geno[j]);
    }
    free(geno);
    printf("\ndone\n");
    exit(0);
}

int read_row(FILE *plink_file, char *chr_name, double *lox, int **geno, int geno_cols, long int row)
{
    int ch,i;
    int ncol;
    int nfield;
    char field[100];
    
    for(i=0;i<11;i++) {chr_name[i]='\0';}
    for(i=0;i<100;i++) {field[i]='\0';}
    ncol = 0; nfield = 0;
    while((ch = fgetc(plink_file)) != '\n' && ch != EOF && ch != '\r') {
        if(ch != 9 && ch != 32) {
            field[nfield++] = ch;
        }
        else {
            if(nfield) {
                switch(ncol) {
                    case 0:
                        for(i=0;i<nfield;i++)
                            chr_name[i] = field[i];
                        break;
                    case 1:
                        lox[row] = atof(field);
                        break;
                    default:
                        geno[ncol - 2][row] = atoi(field);
                        break;
                }
                ncol++;
                for(i=0;i<nfield;i++) {field[i]='\0';}
            }
            nfield = 0;
            if(ncol-2 > geno_cols) {
                printf("\nError: input file having more columns than defined. row: %ld ncols: %d\n",row,ncol);
                exit(1);
            }
        }
    }
    if(nfield) {
        geno[ncol - 2][row] = atoi(field);
        ncol++;
    }
    if(ncol-2 != geno_cols) {
        printf("\nError: input file having different individuals than defined. row: %ld n_indiv: %d \n",row,ncol);
        exit(1);
    }
    if(ch == EOF)
        return(0);
    return(1);
}

/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j)
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

void usage()
{
    printf(TANG_SOFTW);
    printf("\n\nUsage:");
    printf("\nRsbki [Plink-like filename] [number of SNPs] [number of indiv] [threshold value (eg=0.1)] [seed (eg=123456)] [number pops] [number target pop -starting from 1- (eg=1)] [size pop1] [size pop2] ... [size popN]");
    printf("\n\nOutput file is automatically generated using the input filename plus '_Results_Tang.txt'");
    printf("\n\n");
}

