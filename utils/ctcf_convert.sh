#! /bin/bash

input=$1
print_file=$2
output=$3
init_mon=$((100377600/3200))
max_mon=$((100704000/3200))


if [ "$print_file" = "-p" ]
then
    awk -v x=$init_mon '{if ($1 == "chrX") printf("%0.f\t%.d\n", $2/3200-x, $4)}' $input > $output
else
    
p=( $(awk '{print $1}' $output) )
v=( $(awk '{print $2}' $output) )

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
fi
