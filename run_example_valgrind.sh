#run example:
set -x

cd ./source
gcc *.c -lm -o ../bin/Rsbki -O3 -Wall -g

cd ../example
valgrind --leak-check=full --track-origins=yes -v  ../bin/Rsbki ./Window_1.txt 10000 61 0.1 0.05 123456 6 1 5 10 15 12 15 1 BRG BRM CHCU CHCA CHFR TXL 2>./Window_1.txt_valgrind.txt

cd ..
 
