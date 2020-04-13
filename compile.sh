#To compile:
#enter into the source folder and type: gcc *.c -lm -o Rsbki -O -Wall
set -x

cd ./source
gcc *.c -lm -o ../bin/Rsbki -O -Wall
cd ..
