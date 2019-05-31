topo clearbonds
for {set i 0} {$i<499} {incr i} {topo addbond $i [expr $i+1]}
mol modstyle 0 0 Bonds 0.300000 17.000000
mol color Name
mol selection all
mol material Opaque
mol addrep 0
mol modstyle 1 0 DynamicBonds 3.40000 0.300000 12.000000
mol modcolor 1 0 ColorID 4
mol selection all
mol material Opaque
mol addrep 0
mol modselect 2 0 serial 3 25 43 65 86 105 127 145 167 188 204 226 244 266 287 303 325 343 365 386 403 425 443 465 486 
mol modstyle 2 0 VDW 0.500000 12.000000
mol modcolor 2 0 ColorID 1
mol material Opaque
mol addrep 0
mol modselect 3 0 serial 5 21 39 45 61 79 82 100 107 123 141 147 163 181 184 202 206 222 240 246 262 280 283 301 305 321 339 345 361 379 382 400 405 421 439 445 461 479 482 500
mol modstyle 3 0 VDW 0.500000 12.000000
mol modcolor 3 0 ColorID 0
mol material Opaque
