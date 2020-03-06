#! /bin/bash

input=$1
first_chr=$2
last_chr=$3
resolution=$4

init_mon=$(($first_chr/$resolution))
last_mon=$(($last_chr/$resolution))

p=( $(awk '{print $1}' $input) )
v=( $(awk '{print $2}' $input) )

ctcf_num=${#p[@]}
for i in $(seq 0 $(($last_mon-$init_mon))) 
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
