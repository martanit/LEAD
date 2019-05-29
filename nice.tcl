topo clearbonds
for {set i 0} {$i<99} {incr i} {topo addbond $i [expr $i+1]}
mol modstyle 0 0 Bonds 0.300000 17.000000
mol color Name
mol selection all
mol material Opaque
mol addrep 0
mol modstyle 1 0 DynamicBonds 3.380000 0.300000 12.000000
mol modcolor 1 0 ColorID 4
mol selection all
mol material Opaque
mol addrep 0
mol modselect 2 0 serial 3 25 43 65 86 
mol modstyle 2 0 VDW 1.000000 12.000000
mol modcolor 2 0 ColorID 1
mol material Opaque
mol addrep 0
mol modselect 3 0 serial 5 21 39 45 61 79 82 100
mol modstyle 3 0 VDW 1.000000 12.000000
mol modcolor 3 0 ColorID 0
mol material Opaque
