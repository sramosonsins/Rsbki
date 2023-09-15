# Rsbki

Software for calculating iES plus Rsb and iESk plus Rsbk (dividing per individual) statistics.
based on initial statistics from Tang, Thornton & Stoneking, PloS Biology 2007.
version 20230915

##Usage:
Rsbki [genotype filename (one chrom)] [number of SNPs] [number of indiv] [EHH threshold value (eg=0.1)] maf_cut (eg=0.05)] [seed (eg=123456)] [number pops] [number target pop -first is 1- (eg=1)] [size pop1] [size pop2] ... [size popN] [name pop1] ... [name popN]

Output files are automatically generated using the input filename plus an extension
The number of output files are:

 1. file with extension '_imputed.txt' (imputed data)
 2. file with extension 'Results_Tang.txt' (iES statistics)
 3. (if npop>1) files with extension '_Results_Tang.txt_Significant_Results_RsbN' per pop comparison (Rsb (Normalized) statistic)
 4. file with extension '_Results_iESk_pop[POP_name].txt' (iESk statistics per position and individual)
 5. file with extension '_Results_log_iESk_pop[POP_name].txt' (iESk statistics per position and individual)
 6. (if npop>1) files with extension '_Results_Rsbk_TARGETpop[POPname]_versus_REFpop[POPname].txt' per pop comparison versus target pop (Rsbk statistics)
