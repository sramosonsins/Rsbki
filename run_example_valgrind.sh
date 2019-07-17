#run example:
set -x

cd ./source
gcc *.c -lm -o ../bin/Tang_stats_plus -O3 -Wall -g

cd ../example
valgrind --leak-check=full --track-origins=yes -v  ../bin/Tang_stats_plus ./Window_1.txt 10000 61 0.1 123456 6 0 3 8 18 33 45 60 BRG BRM CHCU CHCA CHFR TXL 2>./Window_1.txt_valgrind.txt

cd ..
