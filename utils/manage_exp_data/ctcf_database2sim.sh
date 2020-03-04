#! /bin/bash

# this script generates the ctcf file starting from the file in the database
# you need to give the database file, the firs chromosome number
# and the resolution
# The output is the positions of ctcfs in the simulation:
# you have to correct manually :(

input=$1
output=$2
first_chr=$3
resolution=$4

init_mon=$(($first_chr/$resolution))

awk -v x=$init_mon -v y=$resolution '{if ($1 == "chrX") printf("%0.f\t%.d\n", ($2/y)-x, $4)}' $input > $output
