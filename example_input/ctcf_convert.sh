#! /bin/bash

input=$1
option=$2

if [ "$option" = "print_ctcf" ]
then
awk 'NR==2{init_mon=$2/3000} {if ($1 == "chrX") printf("%0.f\t%.d\n", ($2/3000 - init_mon), $4)}' $input 
else
p=$(awk 'NR==2{init_mon=$2/3000} {if ($1 == "chrX") printf("%0.f\n", ($2/3000 - init_mon))}' $input )
v=$(awk 'NR==2{init_mon=$2/3000} {if ($1 == "chrX") printf("%d \n", $4)}' $input )

p=($p)
v=($v)

ctcf_num=${#p[@]}
ctcf_max=${p[$(($ctcf_num-1))]}

for i in $(seq 0 $(($ctcf_max))) 
do
    is_ctcf="0"
    for j in $(seq 0 $(($ctcf_num-1)))
    do
        if [ "$i" = "${p[$j]}" ] 
        then
            echo ${v[$j]}
            is_ctcf="1"
            break
        fi
    done
    if [ "$is_ctcf" = 0 ]
    then 
        echo 0
    fi
done
fi
