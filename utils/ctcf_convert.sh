#! /bin/bash

input=$1

awk 'NR==2{init_mon=$2/3000} {if ($1 == "chrX") printf("%0.f\t%.d\n", ($2/3000 - init_mon), $4)}' $input > ctcf_file.dat
awk '(NR>1) && ($1!=p || $2!=s) {print p, s; s=0} {p=$1; s+=$2} END{print p, s}' ctcf_file.dat > ctcf_mod.dat

p=$(awk '{print $1}' ctcf_mod.dat )
v=$(awk '{print $2}' ctcf_mod.dat )

p=($p)
v=($v)

ctcf_num=${#p[@]}
ctcf_max=${p[$(($ctcf_num-1))]}

for i in $(seq 0 $(($ctcf_num-1)))
do
	if [ "$i" != "0" ]
	then
		if [ "${p[$i]}" = "${p[$(($i+1))]}" ]
		then
			if [ "$((${p[$i]}-1))" != "${p[$(($i-1))]}" ]
			then	
				echo "$((${p[$i]}-1))   ${v[$i]}" >> ctcf_fin.dat
			else 
				if [ "$((${v[$((i-1))]}*${v[$i]}))" > "0" ]
				then
					echo "$((${p[$i]}-1))   $((${v[$(($i-1))]}+${v[$i]}))" >> ctcf_fin.dat
				else
					echo "$((${p[$i]}-1))   $((${v[$(($i-1))]}+${v[$(($i+1))]}))" >> ctcf_fin.dat
				fi
			fi	
		else
			echo "${p[$i]}   ${v[$i]}" >> ctcf_fin.dat
		fi
	else
		echo "${p[$i]}   ${v[$i]}" >> ctcf_fin.dat
	fi
done

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
