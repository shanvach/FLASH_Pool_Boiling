  I  Y   k820309    w          19.1        ��d                                                                                                          
       gr_mgInterface.F90 GR_MGINTERFACE #         @                                        	               #LEVEL    #IVAR    #NLAYERS    #LEAF_ONLY    #IOPT    #CALL              
                                                      
                                                      
                                                      
                                                      
                                                      
                                            #         @                                        	               #LEVEL 	   #IFROM 
   #ITO    #LEAF_ONLY              
                                 	                     
                                 
                     
                                                      
                                            #         @                                        	               #LEVEL    #ISOLN    #ICORR    #LEAF_ONLY              
                                                      
                                                      
                                                      
                                            #         @                                        	            
   #LEVEL    #IMG_SOLN    #IMG_SRC    #IMG_RES    #IMG_CORR    #IMG_TEMP    #IMG_TEMP2    #MG_SOLVE    #MG_RESIDUAL    #MG_RELAX              
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      �                                               �                                               �                                     #         @                                        	               #MYPE2    #IVAR    #NLAYERS     #SIMTIME !   #IDIAG "   #IDIR #             
                                                      
                                                      
                                                                                       !     
                 
                                 "                     
                                 #           #         @                                   $     	                #         @                                   %     	               #BNDTYPES &             
                                 &                       p          p            p                          #         @                                   '     	                #         @                                   (     	                #         @                                   )     	                #         @                                   *     	               #ISRC_DENS +   #POISFACT ,   #IMG_SRC -   #IMG_SOLN .             
                                 +                     
                                 ,     
                
                                 -                     
                                 .           #         @                                   /     	               #BCTYPETOAPPLY 0   #BCTYPEFROMGRID 1   #VARINDEX 2   #GRIDDATASTRUCT 3   #AXIS 4   #FACE 5   #IDEST 6                                             0                      
                                 1                     
                                 2                     
                                 3                     
                                 4                     
                                 5                     
                                6           #         @                                   7     	               #LEVEL 8   #IVAR 9   #NORM :   #LEAF_ONLY ;             
                                 8                     
                                 9                     
                                :     
                 
                                 ;           #         @                                   <     	               #LEVEL =   #IFROM >   #ITO ?   #ADD @             
                                 =                     
                                 >                     
                                 ?                     
                                 @           #         @                                   A     	               #LEVEL B   #IFROM C   #ITO D             
                                 B                     
                                 C                     
                                 D           #         @                                   E     	               #LEVEL F   #IVAR G   #LEAF_ONLY H             
                                 F                     
                                 G                     
                                 H           #         @                                   I     	               #LEVEL J             
                                 J           #         @                                   K     	               #ISRC_DENS L   #IMG_SOLN M   #POISFACT N   #IMG_SRC O   #IMG_RES P   #IMG_CORR Q   #IMG_TEMP R   #IMG_TEMP2 S   #BC_TYPES T   #MG_SOLVE U   #MG_RESIDUAL V   #MG_RELAX W             
                                 L                     
                                 M                     
                                 N     
                
                                 O                     
                                 P                     
                                 Q                     
                                 R                     
                                 S                     
                                 T                       p          p            p                                    �                                U               �                                V               �                                W     #         @                                   X     	                   �   *      fn#fn    �   �       GR_MGBNDRY !   W  @   a   GR_MGBNDRY%LEVEL     �  @   a   GR_MGBNDRY%IVAR #   �  @   a   GR_MGBNDRY%NLAYERS %     @   a   GR_MGBNDRY%LEAF_ONLY     W  @   a   GR_MGBNDRY%IOPT     �  @   a   GR_MGBNDRY%CALL    �  v       GR_MGCOPY     M  @   a   GR_MGCOPY%LEVEL     �  @   a   GR_MGCOPY%IFROM    �  @   a   GR_MGCOPY%ITO $     @   a   GR_MGCOPY%LEAF_ONLY    M  x       GR_MGCORRECT #   �  @   a   GR_MGCORRECT%LEVEL #     @   a   GR_MGCORRECT%ISOLN #   E  @   a   GR_MGCORRECT%ICORR '   �  @   a   GR_MGCORRECT%LEAF_ONLY    �  �       GR_MGCYCLE !   �  @   a   GR_MGCYCLE%LEVEL $   �  @   a   GR_MGCYCLE%IMG_SOLN #     @   a   GR_MGCYCLE%IMG_SRC #   X  @   a   GR_MGCYCLE%IMG_RES $   �  @   a   GR_MGCYCLE%IMG_CORR $   �  @   a   GR_MGCYCLE%IMG_TEMP %     @   a   GR_MGCYCLE%IMG_TEMP2 $   X  8      GR_MGCYCLE%MG_SOLVE '   �  8      GR_MGCYCLE%MG_RESIDUAL $   �  8      GR_MGCYCLE%MG_RELAX     	  �       GR_MGGUARDCELL %   �	  @   a   GR_MGGUARDCELL%MYPE2 $   �	  @   a   GR_MGGUARDCELL%IVAR '   
  @   a   GR_MGGUARDCELL%NLAYERS '   L
  @   a   GR_MGGUARDCELL%SIMTIME %   �
  @   a   GR_MGGUARDCELL%IDIAG $   �
  @   a   GR_MGGUARDCELL%IDIR      H       GR_MGINIT    T  V       GR_MGINITSLV &   �  �   a   GR_MGINITSLV%BNDTYPES    >  H       GR_MGFINALIZE    �  H       GR_MGPFFTINIT "   �  H       GR_MGPFFTFINALIZE      �       GR_MGINITSRC '   �  @   a   GR_MGINITSRC%ISRC_DENS &   �  @   a   GR_MGINITSRC%POISFACT %     @   a   GR_MGINITSRC%IMG_SRC &   V  @   a   GR_MGINITSRC%IMG_SOLN    �  �       GR_MGMAPBCTYPE -   F  @   a   GR_MGMAPBCTYPE%BCTYPETOAPPLY .   �  @   a   GR_MGMAPBCTYPE%BCTYPEFROMGRID (   �  @   a   GR_MGMAPBCTYPE%VARINDEX .     @   a   GR_MGMAPBCTYPE%GRIDDATASTRUCT $   F  @   a   GR_MGMAPBCTYPE%AXIS $   �  @   a   GR_MGMAPBCTYPE%FACE %   �  @   a   GR_MGMAPBCTYPE%IDEST      v       GR_MGNORM     |  @   a   GR_MGNORM%LEVEL    �  @   a   GR_MGNORM%IVAR    �  @   a   GR_MGNORM%NORM $   <  @   a   GR_MGNORM%LEAF_ONLY    |  p       GR_MGPROLONG #   �  @   a   GR_MGPROLONG%LEVEL #   ,  @   a   GR_MGPROLONG%IFROM !   l  @   a   GR_MGPROLONG%ITO !   �  @   a   GR_MGPROLONG%ADD    �  g       GR_MGRESTRICT $   S  @   a   GR_MGRESTRICT%LEVEL $   �  @   a   GR_MGRESTRICT%IFROM "   �  @   a   GR_MGRESTRICT%ITO      l       GR_MGZERO       @   a   GR_MGZERO%LEVEL    �  @   a   GR_MGZERO%IVAR $   �  @   a   GR_MGZERO%LEAF_ONLY %   ?  S       MG_RESTORE_NODETYPES +   �  @   a   MG_RESTORE_NODETYPES%LEVEL    �  �       MULTIGRID $   �  @   a   MULTIGRID%ISRC_DENS #     @   a   MULTIGRID%IMG_SOLN #   E  @   a   MULTIGRID%POISFACT "   �  @   a   MULTIGRID%IMG_SRC "   �  @   a   MULTIGRID%IMG_RES #     @   a   MULTIGRID%IMG_CORR #   E  @   a   MULTIGRID%IMG_TEMP $   �  @   a   MULTIGRID%IMG_TEMP2 #   �  �   a   MULTIGRID%BC_TYPES #   Y  8      MULTIGRID%MG_SOLVE &   �  8      MULTIGRID%MG_RESIDUAL #   �  8      MULTIGRID%MG_RELAX &     H       GR_MGPFFTFINALIZEGRID 