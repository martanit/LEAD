#! /bin/sh

n_sphere=$(grep -oP "psphere=\K.*" input/parameters.in ) 
m_one=$(grep -n -e "-1" input/ctcf.in | cut -d : -f1 | tr "\n" " ")
p_one=$(grep -n -e "1" input/ctcf.in | grep -v "-" | cut -d : -f1 | tr "\n" " ")

echo "topo clearbonds"
echo "for {set i 0} {\$i<$(($n_sphere-1))} {incr i} {topo addbond \$i [expr \$i+1]}"
echo "mol modstyle 0 0 Bonds 0.300000 17.000000"
echo "mol color Name"
echo "mol selection all"
echo "mol material Opaque"
echo "mol addrep 0"
echo "mol modstyle 1 0 DynamicBonds 3.70000 0.300000 12.000000"
echo "mol modcolor 1 0 ColorID 4"
echo "mol selection all"
echo "mol material Opaque"
echo "mol addrep 0"
echo "mol modselect 2 0 serial $m_one" 
echo "mol modstyle 2 0 VDW 0.500000 12.000000"
echo "mol modcolor 2 0 ColorID 1"
echo "mol material Opaque"
echo "mol addrep 0"
echo "mol modselect 3 0 serial $p_one"
echo "mol modstyle 3 0 VDW 0.500000 12.000000"
echo "mol modcolor 3 0 ColorID 0"
echo "mol material Opaque"
