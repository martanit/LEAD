topo clearbonds
for {set i 0} {$i<99} {incr i} {topo addbond $i [expr $i+1]}
mol modstyle 0 0 Licorice 0.200000 7.000000 20.000000
mol color Name
mol selection all
mol material Opaque
mol addrep 0
mol modstyle 1 0 DynamicBonds 0.50000 0.100000 12.00000
mol modcolor 1 0 ColorID 4
mol selection all
mol material Opaque
mol addrep 0
mol modselect 2 0 serial 34 43 
mol modstyle 2 0 VDW 0.300000 12.000000
mol modcolor 2 0 ColorID 1
mol material Opaque
mol addrep 0
mol modselect 3 0 serial 9 66 
mol modstyle 3 0 VDW 0.300000 12.000000
mol modcolor 3 0 ColorID 0
mol material Opaque
