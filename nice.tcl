topo clearbonds
for {set i 0} {$i<199} {incr i} {topo addbond $i [expr $i+1]}
mol modstyle 0 0 Bonds 0.300000 0.500000
mol color Name
mol selection all
mol material Opaque
mol addrep 0
mol modstyle 1 0 DynamicBonds 0.10000 0.300000 0.500000
mol modcolor 1 0 ColorID 4
mol selection all
mol material Opaque
mol addrep 0
mol modselect 2 0 serial 43 
mol modstyle 2 0 VDW 0.500000 0.800000
mol modcolor 2 0 ColorID 1
mol material Opaque
mol addrep 0
mol modselect 3 0 serial 125 
mol modstyle 3 0 VDW 0.500000 0.800000
mol modcolor 3 0 ColorID 0
mol material Opaque
