#run example:
set -x
 
cd ./source
gcc *.c -lm -o ../bin/Rsbki -O3 -Wall

cd ../example
time ../bin/Rsbki ./Window_10K.txt 10000 61 0.1 0.05 123456 6 1 5 10 15 12 15 1 BRG BRM CHCU CHCA CHFR TXL
Rscript --vanilla ../bin/run_plots_Rsb.R ./Window_10K.txt_Results_Tang.txt 20 6 BRG BRM CHCU CHCA CHFR TXL

cd ..

