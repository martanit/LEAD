topo clearbonds
for {set i 0} {$i<19999} {incr i} {topo addbond $i [expr $i+1]}
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
mol modselect 2 0 serial 3 25 43 65 86 105 127 145 167 188 204 226 244 266 287 303 325 343 365 386 403 425 443 465 486 503 525 543 565 586 605 627 645 667 688 704 726 744 766 787 803 825 843 865 886 903 925 943 965 986 1003 1025 1043 1065 1086 1105 1127 1145 1167 1188 1204 1226 1244 1266 1287 1303 1325 1343 1365 1386 1403 1425 1443 1465 1486 1503 1525 1543 1565 1586 1605 1627 1645 1667 1688 1704 1726 1744 1766 1787 1803 1825 1843 1865 1886 1903 1925 1943 1965 1986 2003 2025 2043 2065 2086 2105 2127 2145 2167 2188 2204 2226 2244 2266 2287 2303 2325 2343 2365 2386 2403 2425 2443 2465 2486 2503 2525 2543 2565 2586 2605 2627 2645 2667 2688 2704 2726 2744 2766 2787 2803 2825 2843 2865 2886 2903 2925 2943 2965 2986 3003 3025 3043 3065 3086 3105 3127 3145 3167 3188 3204 3226 3244 3266 3287 3303 3325 3343 3365 3386 3403 3425 3443 3465 3486 3503 3525 3543 3565 3586 3605 3627 3645 3667 3688 3704 3726 3744 3766 3787 3803 3825 3843 3865 3886 3903 3925 3943 3965 3986 4003 4025 4043 4065 4086 4105 4127 4145 4167 4188 4204 4226 4244 4266 4287 4303 4325 4343 4365 4386 4403 4425 4443 4465 4486 4503 4525 4543 4565 4586 4605 4627 4645 4667 4688 4704 4726 4744 4766 4787 4803 4825 4843 4865 4886 4903 4925 4943 4965 4986 5003 5025 5043 5065 5086 5105 5127 5145 5167 5188 5204 5226 5244 5266 5287 5303 5325 5343 5365 5386 5403 5425 5443 5465 5486 5503 5525 5543 5565 5586 5605 5627 5645 5667 5688 5704 5726 5744 5766 5787 5803 5825 5843 5865 5886 5903 5925 5943 5965 5986 6003 6025 6043 6065 6086 6105 6127 6145 6167 6188 6204 6226 6244 6266 6287 6303 6325 6343 6365 6386 6403 6425 6443 6465 6486 6503 6525 6543 6565 6586 6605 6627 6645 6667 6688 6704 6726 6744 6766 6787 6803 6825 6843 6865 6886 6903 6925 6943 6965 6986 7003 7025 7043 7065 7086 7105 7127 7145 7167 7188 7204 7226 7244 7266 7287 7303 7325 7343 7365 7386 7403 7425 7443 7465 7486 7503 7525 7543 7565 7586 7605 7627 7645 7667 7688 7704 7726 7744 7766 7787 7803 7825 7843 7865 7886 7903 7925 7943 7965 7986 8003 8025 8043 8065 8086 8105 8127 8145 8167 8188 8204 8226 8244 8266 8287 8303 8325 8343 8365 8386 8403 8425 8443 8465 8486 8503 8525 8543 8565 8586 8605 8627 8645 8667 8688 8704 8726 8744 8766 8787 8803 8825 8843 8865 8886 8903 8925 8943 8965 8986 9003 9025 9043 9065 9086 9105 9127 9145 9167 9188 9204 9226 9244 9266 9287 9303 9325 9343 9365 9386 9403 9425 9443 9465 9486 9503 9525 9543 9565 9586 9605 9627 9645 9667 9688 9704 9726 9744 9766 9787 9803 9825 9843 9865 9886 9903 9925 9943 9965 9986 10003 10025 10043 10065 10086 10105 10127 10145 10167 10188 10204 10226 10244 10266 10287 10303 10325 10343 10365 10386 10403 10425 10443 10465 10486 10503 10525 10543 10565 10586 10605 10627 10645 10667 10688 10704 10726 10744 10766 10787 10803 10825 10843 10865 10886 10903 10925 10943 10965 10986 11003 11025 11043 11065 11086 11105 11127 11145 11167 11188 11204 11226 11244 11266 11287 11303 11325 11343 11365 11386 11403 11425 11443 11465 11486 11503 11525 11543 11565 11586 11605 11627 11645 11667 11688 11704 11726 11744 11766 11787 11803 11825 11843 11865 11886 11903 11925 11943 11965 11986 12003 12025 12043 12065 12086 12105 12127 12145 12167 12188 12204 12226 12244 12266 12287 12303 12325 12343 12365 12386 12403 12425 12443 12465 12486 12503 12525 12543 12565 12586 12605 12627 12645 12667 12688 12704 12726 12744 12766 12787 12803 12825 12843 12865 12886 12903 12925 12943 12965 12986 13003 13025 13043 13065 13086 13105 13127 13145 13167 13188 13204 13226 13244 13266 13287 13303 13325 13343 13365 13386 13403 13425 13443 13465 13486 13503 13525 13543 13565 13586 13605 13627 13645 13667 13688 13704 13726 13744 13766 13787 13803 13825 13843 13865 13886 13903 13925 13943 13965 13986 14003 14025 14043 14065 14086 14105 14127 14145 14167 14188 14204 14226 14244 14266 14287 14303 14325 14343 14365 14386 14403 14425 14443 14465 14486 14503 14525 14543 14565 14586 14605 14627 14645 14667 14688 14704 14726 14744 14766 14787 14803 14825 14843 14865 14886 14903 14925 14943 14965 14986 15003 15025 15043 15065 15086 15105 15127 15145 15167 15188 15204 15226 15244 15266 15287 15303 15325 15343 15365 15386 15403 15425 15443 15465 15486 15503 15525 15543 15565 15586 15605 15627 15645 15667 15688 15704 15726 15744 15766 15787 15803 15825 15843 15865 15886 15903 15925 15943 15965 15986 16003 16025 16043 16065 16086 16105 16127 16145 16167 16188 16204 16226 16244 16266 16287 16303 16325 16343 16365 16386 16403 16425 16443 16465 16486 16503 16525 16543 16565 16586 16605 16627 16645 16667 16688 16704 16726 16744 16766 16787 16803 16825 16843 16865 16886 16903 16925 16943 16965 16986 17003 17025 17043 17065 17086 17105 17127 17145 17167 17188 17204 17226 17244 17266 17287 17303 17325 17343 17365 17386 17403 17425 17443 17465 17486 17503 17525 17543 17565 17586 17605 17627 17645 17667 17688 17704 17726 17744 17766 17787 17803 17825 17843 17865 17886 17903 17925 17943 17965 17986 18003 18025 18043 18065 18086 18105 18127 18145 18167 18188 18204 18226 18244 18266 18287 18303 18325 18343 18365 18386 18403 18425 18443 18465 18486 18503 18525 18543 18565 18586 18605 18627 18645 18667 18688 18704 18726 18744 18766 18787 18803 18825 18843 18865 18886 18903 18925 18943 18965 18986 19003 19025 19043 19065 19086 19105 19127 19145 19167 19188 19204 19226 19244 19266 19287 19303 19325 19343 19365 19386 19403 19425 19443 19465 19486 19503 19525 19543 19565 19586 19605 19627 19645 19667 19688 19704 19726 19744 19766 19787 19803 19825 19843 19865 19886 19903 19925 19943 19965 19986 
mol modstyle 2 0 VDW 0.500000 12.000000
mol modcolor 2 0 ColorID 1
mol material Opaque
mol addrep 0
mol modselect 3 0 serial 5 21 39 45 61 79 82 100 107 123 141 147 163 181 184 202 206 222 240 246 262 280 283 301 305 321 339 345 361 379 382 400 405 421 439 445 461 479 482 500 505 521 539 545 561 579 582 600 607 623 641 647 663 681 684 702 706 722 740 746 762 780 783 801 805 821 839 845 861 879 882 900 905 921 939 945 961 979 982 1000 1005 1021 1039 1045 1061 1079 1082 1100 1107 1123 1141 1147 1163 1181 1184 1202 1206 1222 1240 1246 1262 1280 1283 1301 1305 1321 1339 1345 1361 1379 1382 1400 1405 1421 1439 1445 1461 1479 1482 1500 1505 1521 1539 1545 1561 1579 1582 1600 1607 1623 1641 1647 1663 1681 1684 1702 1706 1722 1740 1746 1762 1780 1783 1801 1805 1821 1839 1845 1861 1879 1882 1900 1905 1921 1939 1945 1961 1979 1982 2000 2005 2021 2039 2045 2061 2079 2082 2100 2107 2123 2141 2147 2163 2181 2184 2202 2206 2222 2240 2246 2262 2280 2283 2301 2305 2321 2339 2345 2361 2379 2382 2400 2405 2421 2439 2445 2461 2479 2482 2500 2505 2521 2539 2545 2561 2579 2582 2600 2607 2623 2641 2647 2663 2681 2684 2702 2706 2722 2740 2746 2762 2780 2783 2801 2805 2821 2839 2845 2861 2879 2882 2900 2905 2921 2939 2945 2961 2979 2982 3000 3005 3021 3039 3045 3061 3079 3082 3100 3107 3123 3141 3147 3163 3181 3184 3202 3206 3222 3240 3246 3262 3280 3283 3301 3305 3321 3339 3345 3361 3379 3382 3400 3405 3421 3439 3445 3461 3479 3482 3500 3505 3521 3539 3545 3561 3579 3582 3600 3607 3623 3641 3647 3663 3681 3684 3702 3706 3722 3740 3746 3762 3780 3783 3801 3805 3821 3839 3845 3861 3879 3882 3900 3905 3921 3939 3945 3961 3979 3982 4000 4005 4021 4039 4045 4061 4079 4082 4100 4107 4123 4141 4147 4163 4181 4184 4202 4206 4222 4240 4246 4262 4280 4283 4301 4305 4321 4339 4345 4361 4379 4382 4400 4405 4421 4439 4445 4461 4479 4482 4500 4505 4521 4539 4545 4561 4579 4582 4600 4607 4623 4641 4647 4663 4681 4684 4702 4706 4722 4740 4746 4762 4780 4783 4801 4805 4821 4839 4845 4861 4879 4882 4900 4905 4921 4939 4945 4961 4979 4982 5000 5005 5021 5039 5045 5061 5079 5082 5100 5107 5123 5141 5147 5163 5181 5184 5202 5206 5222 5240 5246 5262 5280 5283 5301 5305 5321 5339 5345 5361 5379 5382 5400 5405 5421 5439 5445 5461 5479 5482 5500 5505 5521 5539 5545 5561 5579 5582 5600 5607 5623 5641 5647 5663 5681 5684 5702 5706 5722 5740 5746 5762 5780 5783 5801 5805 5821 5839 5845 5861 5879 5882 5900 5905 5921 5939 5945 5961 5979 5982 6000 6005 6021 6039 6045 6061 6079 6082 6100 6107 6123 6141 6147 6163 6181 6184 6202 6206 6222 6240 6246 6262 6280 6283 6301 6305 6321 6339 6345 6361 6379 6382 6400 6405 6421 6439 6445 6461 6479 6482 6500 6505 6521 6539 6545 6561 6579 6582 6600 6607 6623 6641 6647 6663 6681 6684 6702 6706 6722 6740 6746 6762 6780 6783 6801 6805 6821 6839 6845 6861 6879 6882 6900 6905 6921 6939 6945 6961 6979 6982 7000 7005 7021 7039 7045 7061 7079 7082 7100 7107 7123 7141 7147 7163 7181 7184 7202 7206 7222 7240 7246 7262 7280 7283 7301 7305 7321 7339 7345 7361 7379 7382 7400 7405 7421 7439 7445 7461 7479 7482 7500 7505 7521 7539 7545 7561 7579 7582 7600 7607 7623 7641 7647 7663 7681 7684 7702 7706 7722 7740 7746 7762 7780 7783 7801 7805 7821 7839 7845 7861 7879 7882 7900 7905 7921 7939 7945 7961 7979 7982 8000 8005 8021 8039 8045 8061 8079 8082 8100 8107 8123 8141 8147 8163 8181 8184 8202 8206 8222 8240 8246 8262 8280 8283 8301 8305 8321 8339 8345 8361 8379 8382 8400 8405 8421 8439 8445 8461 8479 8482 8500 8505 8521 8539 8545 8561 8579 8582 8600 8607 8623 8641 8647 8663 8681 8684 8702 8706 8722 8740 8746 8762 8780 8783 8801 8805 8821 8839 8845 8861 8879 8882 8900 8905 8921 8939 8945 8961 8979 8982 9000 9005 9021 9039 9045 9061 9079 9082 9100 9107 9123 9141 9147 9163 9181 9184 9202 9206 9222 9240 9246 9262 9280 9283 9301 9305 9321 9339 9345 9361 9379 9382 9400 9405 9421 9439 9445 9461 9479 9482 9500 9505 9521 9539 9545 9561 9579 9582 9600 9607 9623 9641 9647 9663 9681 9684 9702 9706 9722 9740 9746 9762 9780 9783 9801 9805 9821 9839 9845 9861 9879 9882 9900 9905 9921 9939 9945 9961 9979 9982 10000 10005 10021 10039 10045 10061 10079 10082 10100 10107 10123 10141 10147 10163 10181 10184 10202 10206 10222 10240 10246 10262 10280 10283 10301 10305 10321 10339 10345 10361 10379 10382 10400 10405 10421 10439 10445 10461 10479 10482 10500 10505 10521 10539 10545 10561 10579 10582 10600 10607 10623 10641 10647 10663 10681 10684 10702 10706 10722 10740 10746 10762 10780 10783 10801 10805 10821 10839 10845 10861 10879 10882 10900 10905 10921 10939 10945 10961 10979 10982 11000 11005 11021 11039 11045 11061 11079 11082 11100 11107 11123 11141 11147 11163 11181 11184 11202 11206 11222 11240 11246 11262 11280 11283 11301 11305 11321 11339 11345 11361 11379 11382 11400 11405 11421 11439 11445 11461 11479 11482 11500 11505 11521 11539 11545 11561 11579 11582 11600 11607 11623 11641 11647 11663 11681 11684 11702 11706 11722 11740 11746 11762 11780 11783 11801 11805 11821 11839 11845 11861 11879 11882 11900 11905 11921 11939 11945 11961 11979 11982 12000 12005 12021 12039 12045 12061 12079 12082 12100 12107 12123 12141 12147 12163 12181 12184 12202 12206 12222 12240 12246 12262 12280 12283 12301 12305 12321 12339 12345 12361 12379 12382 12400 12405 12421 12439 12445 12461 12479 12482 12500 12505 12521 12539 12545 12561 12579 12582 12600 12607 12623 12641 12647 12663 12681 12684 12702 12706 12722 12740 12746 12762 12780 12783 12801 12805 12821 12839 12845 12861 12879 12882 12900 12905 12921 12939 12945 12961 12979 12982 13000 13005 13021 13039 13045 13061 13079 13082 13100 13107 13123 13141 13147 13163 13181 13184 13202 13206 13222 13240 13246 13262 13280 13283 13301 13305 13321 13339 13345 13361 13379 13382 13400 13405 13421 13439 13445 13461 13479 13482 13500 13505 13521 13539 13545 13561 13579 13582 13600 13607 13623 13641 13647 13663 13681 13684 13702 13706 13722 13740 13746 13762 13780 13783 13801 13805 13821 13839 13845 13861 13879 13882 13900 13905 13921 13939 13945 13961 13979 13982 14000 14005 14021 14039 14045 14061 14079 14082 14100 14107 14123 14141 14147 14163 14181 14184 14202 14206 14222 14240 14246 14262 14280 14283 14301 14305 14321 14339 14345 14361 14379 14382 14400 14405 14421 14439 14445 14461 14479 14482 14500 14505 14521 14539 14545 14561 14579 14582 14600 14607 14623 14641 14647 14663 14681 14684 14702 14706 14722 14740 14746 14762 14780 14783 14801 14805 14821 14839 14845 14861 14879 14882 14900 14905 14921 14939 14945 14961 14979 14982 15000 15005 15021 15039 15045 15061 15079 15082 15100 15107 15123 15141 15147 15163 15181 15184 15202 15206 15222 15240 15246 15262 15280 15283 15301 15305 15321 15339 15345 15361 15379 15382 15400 15405 15421 15439 15445 15461 15479 15482 15500 15505 15521 15539 15545 15561 15579 15582 15600 15607 15623 15641 15647 15663 15681 15684 15702 15706 15722 15740 15746 15762 15780 15783 15801 15805 15821 15839 15845 15861 15879 15882 15900 15905 15921 15939 15945 15961 15979 15982 16000 16005 16021 16039 16045 16061 16079 16082 16100 16107 16123 16141 16147 16163 16181 16184 16202 16206 16222 16240 16246 16262 16280 16283 16301 16305 16321 16339 16345 16361 16379 16382 16400 16405 16421 16439 16445 16461 16479 16482 16500 16505 16521 16539 16545 16561 16579 16582 16600 16607 16623 16641 16647 16663 16681 16684 16702 16706 16722 16740 16746 16762 16780 16783 16801 16805 16821 16839 16845 16861 16879 16882 16900 16905 16921 16939 16945 16961 16979 16982 17000 17005 17021 17039 17045 17061 17079 17082 17100 17107 17123 17141 17147 17163 17181 17184 17202 17206 17222 17240 17246 17262 17280 17283 17301 17305 17321 17339 17345 17361 17379 17382 17400 17405 17421 17439 17445 17461 17479 17482 17500 17505 17521 17539 17545 17561 17579 17582 17600 17607 17623 17641 17647 17663 17681 17684 17702 17706 17722 17740 17746 17762 17780 17783 17801 17805 17821 17839 17845 17861 17879 17882 17900 17905 17921 17939 17945 17961 17979 17982 18000 18005 18021 18039 18045 18061 18079 18082 18100 18107 18123 18141 18147 18163 18181 18184 18202 18206 18222 18240 18246 18262 18280 18283 18301 18305 18321 18339 18345 18361 18379 18382 18400 18405 18421 18439 18445 18461 18479 18482 18500 18505 18521 18539 18545 18561 18579 18582 18600 18607 18623 18641 18647 18663 18681 18684 18702 18706 18722 18740 18746 18762 18780 18783 18801 18805 18821 18839 18845 18861 18879 18882 18900 18905 18921 18939 18945 18961 18979 18982 19000 19005 19021 19039 19045 19061 19079 19082 19100 19107 19123 19141 19147 19163 19181 19184 19202 19206 19222 19240 19246 19262 19280 19283 19301 19305 19321 19339 19345 19361 19379 19382 19400 19405 19421 19439 19445 19461 19479 19482 19500 19505 19521 19539 19545 19561 19579 19582 19600 19607 19623 19641 19647 19663 19681 19684 19702 19706 19722 19740 19746 19762 19780 19783 19801 19805 19821 19839 19845 19861 19879 19882 19900 19905 19921 19939 19945 19961 19979 19982 20000 
mol modstyle 3 0 VDW 0.500000 12.000000
mol modcolor 3 0 ColorID 0
mol material Opaque
