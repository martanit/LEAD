#! /bin/bash

input=$1

init_mon=$((100378307/3000))
max_mon=$((100702306/3000))

awk -v x=$init_mon '{if ($1 == "chrX") printf("%0.f\t%.d\n", $2/3000-x, $4)}' $input > ctcf_file.dat

p=$(awk '{print $1}' ctcf_file.dat)
v=$(awk '{print $2}' ctcf_file.dat)

p=($p)
v=($v)

ctcf_num=${#p[@]}
for i in $(seq 0 $(($max_mon-$init_mon))) 
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
