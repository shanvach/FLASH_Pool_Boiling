  |   l   k820309              18.0        {i�`                                                                                                          
       gr_hgInterface.F90 GR_HGINTERFACE #         @                                        	               #LEVEL    #IVAR    #NLAYERS    #LEAFONLY    #IOPT    #CALL    #EXTRAP              
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                            #         @                                   	     	               #MYPE 
   #NLAYERS    #EXTRAP              
                                 
                     
                                                      
                                            #         @                                        	                #         @                                        	               #ISRC    #ISOLN              
                                                      
                                            #         @                                        	               #LEVEL    #IVAR1    #IVAR2    #LEAFFLAG              
                                                      
                                                      
                                                      
                                            #         @                                        	               #LEVEL    #IVAR    #SCALAR    #LEAFFLAG              
                                                      
                                                      
                                      
                
                                            #         @                                        	               #LEVEL    #IVAR    #SCALAR    #LEAFFLAG              
                                                      
                                                      
                                      
                
                                            #         @                                         	               #LEVEL !   #IVAR "   #LEAFFLAG #             
                                 !                     
                                 "                     
                                 #           #         @                                   $     	               #BCTYPETOAPPLY %   #BCTYPEFROMGRID &   #VARINDEX '   #GRIDDATASTRUCT (   #AXIS )   #FACE *   #IDEST +                                             %                      
                                 &                     
                                 '                     
                                 (                     
                                 )                     
                                 *                     
                                +           #         @                                   ,     	               #LEVEL -   #NORMTYPE .   #IVAR /   #NORM 0   #LEAFONLY 1             
                                 -                     
                                 .                     
                                 /                                                     0     
                 
                                 1           #         @                                   2     	               #LEVEL 3   #IFROM 4   #ITO 5   #ICHILD 6             
                                 3                     
                                 4                     
                                 5                     
                                 6           #         @                                   7     	               #LEVEL 8   #GR_ISOURCE 9   #GR_ISOLN :   #IRES ;   #DT <   #CHI =   #THETA >             
                                 8                     
                                 9                     
                                 :                     
                                 ;                     
                                <     
                
                                =     
                
                                >     
      #         @                                   ?     	                #         @                                   @     	               #LEVEL A   #ITO B   #IFROM C             
                                 A                     
                                 B                     
                                 C           #         @                                   D     	               #LEVEL E   #GR_ISOLN F             
                                 E                     
                                 F           #         @                                   G     	               #IDIAG H   #IDIR I   #EXTRAP J             
                                 H                     
                                 I                     
                                 J           #         @                                   K     	               #LEVEL L             
                                 L           #         @                                   M     	            	   #GR_ISOURCE N   #GR_ISOLN O   #GR_ISLS P   #GR_ICORR Q   #SOLVEBLOCK R   #BNDTYPES S   #SRC_FACT T   #DT U   #CHI V             
                                 N                     
                                 O                     
                                 P                     
                                 Q                     �                                R               
                                 S                       p          p            p                                    
                                T     
                
                                U     
                
                                V     
      #         @                                   W     	               #LEVEL X   #GR_ISOURCE Y   #GR_ISOLN Z   #SOLVEBLOCK [   #LEAFFLAG \   #DT ]   #CHI ^   #THETA _             
                                 X                     
                                 Y                     
                                 Z                     �                                [               
                                 \                     
                                ]     
                
                                ^     
                
                                _     
      #         @                                   `     	                #         @                                   a     	                #         @                                   b     	               #REFINEMENTLEVEL c   #GRIDCHANGED d   #POISFACT e             
                                 c                     
                                 d                     
                                 e     
      #         @                                   f     	               #ISRC g   #ISOLN h   #LEVEL i             
                                 g                     
                                 h                     
                                 i           #         @                                   j     	                #         @                                   k     	                   �   *      fn#fn    �   �       GR_HGBNDRY !   b  @   a   GR_HGBNDRY%LEVEL     �  @   a   GR_HGBNDRY%IVAR #   �  @   a   GR_HGBNDRY%NLAYERS $   "  @   a   GR_HGBNDRY%LEAFONLY     b  @   a   GR_HGBNDRY%IOPT     �  @   a   GR_HGBNDRY%CALL "   �  @   a   GR_HGBNDRY%EXTRAP    "  k       GR_HGGUARDCELL $   �  @   a   GR_HGGUARDCELL%MYPE '   �  @   a   GR_HGGUARDCELL%NLAYERS &     @   a   GR_HGGUARDCELL%EXTRAP    M  H       GR_HGINIT     �  ]       GR_HGINITSOURCE %   �  @   a   GR_HGINITSOURCE%ISRC &   2  @   a   GR_HGINITSOURCE%ISOLN    r  w       GR_HGLEVELADD $   �  @   a   GR_HGLEVELADD%LEVEL $   )  @   a   GR_HGLEVELADD%IVAR1 $   i  @   a   GR_HGLEVELADD%IVAR2 '   �  @   a   GR_HGLEVELADD%LEAFFLAG $   �  w       GR_HGLEVELADDSCALAR *   `  @   a   GR_HGLEVELADDSCALAR%LEVEL )   �  @   a   GR_HGLEVELADDSCALAR%IVAR +   �  @   a   GR_HGLEVELADDSCALAR%SCALAR -      @   a   GR_HGLEVELADDSCALAR%LEAFFLAG )   `  w       GR_HGLEVELMULTIPLYSCALAR /   �  @   a   GR_HGLEVELMULTIPLYSCALAR%LEVEL .   	  @   a   GR_HGLEVELMULTIPLYSCALAR%IVAR 0   W	  @   a   GR_HGLEVELMULTIPLYSCALAR%SCALAR 2   �	  @   a   GR_HGLEVELMULTIPLYSCALAR%LEAFFLAG    �	  k       GR_HGLEVELZERO %   B
  @   a   GR_HGLEVELZERO%LEVEL $   �
  @   a   GR_HGLEVELZERO%IVAR (   �
  @   a   GR_HGLEVELZERO%LEAFFLAG      �       GR_HGMAPBCTYPE -   �  @   a   GR_HGMAPBCTYPE%BCTYPETOAPPLY .   �  @   a   GR_HGMAPBCTYPE%BCTYPEFROMGRID (   2  @   a   GR_HGMAPBCTYPE%VARINDEX .   r  @   a   GR_HGMAPBCTYPE%GRIDDATASTRUCT $   �  @   a   GR_HGMAPBCTYPE%AXIS $   �  @   a   GR_HGMAPBCTYPE%FACE %   2  @   a   GR_HGMAPBCTYPE%IDEST    r  �       GR_HGNORM     �  @   a   GR_HGNORM%LEVEL #   5  @   a   GR_HGNORM%NORMTYPE    u  @   a   GR_HGNORM%IVAR    �  @   a   GR_HGNORM%NORM #   �  @   a   GR_HGNORM%LEAFONLY $   5  s       GR_HGPROLONGBNDRIES *   �  @   a   GR_HGPROLONGBNDRIES%LEVEL *   �  @   a   GR_HGPROLONGBNDRIES%IFROM (   (  @   a   GR_HGPROLONGBNDRIES%ITO +   h  @   a   GR_HGPROLONGBNDRIES%ICHILD    �  �       GR_HGRESIDUAL $   ?  @   a   GR_HGRESIDUAL%LEVEL )     @   a   GR_HGRESIDUAL%GR_ISOURCE '   �  @   a   GR_HGRESIDUAL%GR_ISOLN #   �  @   a   GR_HGRESIDUAL%IRES !   ?  @   a   GR_HGRESIDUAL%DT "     @   a   GR_HGRESIDUAL%CHI $   �  @   a   GR_HGRESIDUAL%THETA &   �  H       GR_HGRESTORENODETYPES    G  g       GR_HGRESTRICT $   �  @   a   GR_HGRESTRICT%LEVEL "   �  @   a   GR_HGRESTRICT%ITO $   .  @   a   GR_HGRESTRICT%IFROM %   n  a       GR_HGSETZEROBOUNDARY +   �  @   a   GR_HGSETZEROBOUNDARY%LEVEL .     @   a   GR_HGSETZEROBOUNDARY%GR_ISOLN $   O  i       GR_HGSETEXTBOUNDARY *   �  @   a   GR_HGSETEXTBOUNDARY%IDIAG )   �  @   a   GR_HGSETEXTBOUNDARY%IDIR +   8  @   a   GR_HGSETEXTBOUNDARY%EXTRAP !   x  S       GR_HGSETMAXLEVEL '   �  @   a   GR_HGSETMAXLEVEL%LEVEL      �       GR_HGSOLVE &   �  @   a   GR_HGSOLVE%GR_ISOURCE $   	  @   a   GR_HGSOLVE%GR_ISOLN #   I  @   a   GR_HGSOLVE%GR_ISLS $   �  @   a   GR_HGSOLVE%GR_ICORR &   �  8      GR_HGSOLVE%SOLVEBLOCK $     �   a   GR_HGSOLVE%BNDTYPES $   �  @   a   GR_HGSOLVE%SRC_FACT    �  @   a   GR_HGSOLVE%DT      @   a   GR_HGSOLVE%CHI     U  �       GR_HGSOLVELEVEL &      @   a   GR_HGSOLVELEVEL%LEVEL +   @  @   a   GR_HGSOLVELEVEL%GR_ISOURCE )   �  @   a   GR_HGSOLVELEVEL%GR_ISOLN +   �  8      GR_HGSOLVELEVEL%SOLVEBLOCK )   �  @   a   GR_HGSOLVELEVEL%LEAFFLAG #   8  @   a   GR_HGSOLVELEVEL%DT $   x  @   a   GR_HGSOLVELEVEL%CHI &   �  @   a   GR_HGSOLVELEVEL%THETA    �  H       GR_HGFINALIZE    @  H       GR_HGPFFTINIT "   �  |       GR_HGPFFTINITGRID 2     @   a   GR_HGPFFTINITGRID%REFINEMENTLEVEL .   D  @   a   GR_HGPFFTINITGRID%GRIDCHANGED +   �  @   a   GR_HGPFFTINITGRID%POISFACT $   �  h       GR_HGPFFTSOLVELEVEL )   ,  @   a   GR_HGPFFTSOLVELEVEL%ISRC *   l  @   a   GR_HGPFFTSOLVELEVEL%ISOLN *   �  @   a   GR_HGPFFTSOLVELEVEL%LEVEL "   �  H       GR_HGPFFTFINALIZE &   4   H       GR_HGPFFTFINALIZEGRID 