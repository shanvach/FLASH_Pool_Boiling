  cd  �   k820309    w          19.1        ֚�d                                                                                                          
       gr_pfftInterface.F90 GR_PFFTINTERFACE          @                                        	   u #GR_PFFTTRANSPOSE    #GR_PFFTTRANSPOSE3DARR    #         @                                        	                #         @                                        	                #         @                                        	               #INARRAY    #OUTARRAY    #TRIG    #LEN 	   #LDA 
   #NUMVEC    #TRANSFORMTYPE    #SCALE              
                                                    
              &                                                                                                        
               &                                                     
                                                    
              &                                                     
                                 	                     
                                 
                     
                                                      
                                                      
                                      
      #         @                                        	               #INARRAY    #OUTARRAY    #TRIG    #LEN    #LDA    #NUMVEC    #TRANSFORMTYPE    #SCALE              
                                                    
              &                                                                                                        
               &                                                     
                                                    
              &                                                     
                                                      
                                                      
                                                      
                                                      
                                      
      #         @                                        	               #INARRAY    #OUTARRAY    #NX    #NY    #NX1              
                                                    
              &                                                                                                        
               &                                                     
                                                      
                                                      
                                            #         @                                        	                #         @                                        	               #NDIM              
                                            #         @                                         	               #INARRAY !   #OUTARRAY "   #NX #   #NY $   #NX1 %             
                                 !                   
 	             &                                                                                     "                   
 
              &                                                     
                                 #                     
                                 $                     
                                 %           #         @                                   &     	               #INARRAY '   #OUTARRAY (   #NX )   #NX1 *   #BASEDATTYPE +             
                                 '                   
              &                                                                                     (                   
               &                                                     
                                 )                     
                                 *                     
                                 +           #         @                                   ,     	               #LENGTH -   #TRANSFORMTYPE .   #TRIG /   #SCALE 0   #INDATSIZE 1   #OUTDATSIZE 2   #DATINC 3   #FACTOR 4             
                                 -                     
                                 .                                                   /                   
               &                                                                                     0     
                                                 1                                                      2                                                      3                                                      4            #         @                                       	               #DIR 5   #BASEDATTYPE 6   #INARRAY 7   #OUTARRAY 8   #INLEN 9   #OUTLEN :   #PR ;   #COMM <             
                                 5                     
                                 6                     
                                7                   
               &                                                     
                                8                   
               &                                                     
                                 9                       p          p            p                                    
                                 :                       p          p            p                                    
                                 ;                     
                                 <           #         @                                        	               #DIR =   #BASEDATTYPE >   #INARRAY ?   #OUTARRAY @   #INLEN A   #OUTLEN B   #PR C   #COMM D             
                                 =                     
                                 >                     
                                ?                   
               &                   &                   &                                                     
                                @                   
               &                   &                   &                                                     
                                 A                       p          p            p                                    
                                 B                       p          p            p                                    
                                 C                     
                                 D           #         @                                   E     	               #AXIS1 F   #AXIS2 G             
                                 F                     
                                 G           #         @                                   H     	                #         @                                   I     	               #TRANSARRAY J             
                                J                   
               &                                           #         @                                   K     	               #DIMS L   #PFFTMYPE M   #PFFTNUMPROCS N   #PFFTGLOBALLEN O   #PFFTPROCGRID P             
                                 L                     
                                 M                     
                                 N                     
                                 O                       p          p            p                                                                    P                        p          p            p                          #         @                                   Q     	               #AXIS R   #FRAGMENTPTR S   #MAXSINGLEPROCDATA T             
                                 R                     
                                 S                                 &                                                                                     T            #         @                                   U     	               #DIRECTION V   #AXIS W   #BUFFER X   #PFFTARRAY Y             
                                 V                     
                                 W                                                     X                   
               &                   &                                                                                     Y                   
               &                                           #         @                                   Z     	               #STARTPOS [   #PFFTPROCCOORDS \             
      �                           [                       p          & p        p            p                                         �                           \                        p          & p        p            p                          #         @                                   ]     	               #BLOCKID ^   #BLOCKSTARTPOS _   #BLOCKENDPOS `   #BSIZE a   #AXIS b             
                                 ^                     
      �                           _                       p          & p        p            p                                    
      �                           `                       p          & p        p            p                                    
                                 a                     
                                 b           #         @                                   c     	               #DIRECTION d   #GRIDVAR e             
                                 d                     
                                 e           #         @                                   f     	               #DIRECTION g             
                                 g           #         @                                   h     	               #BUFFER i   #MAP j   #LOGUNIT k             
                                 i                   
               &                   &                                                     
                                 j                    !             &                   &                                                     
                                 k           #         @                                   l     	               #MYPE m   #GLOBALPROCS n   #ORIGINALCOMM o   #NEWCOMM p             
                                 m                     
                                 n                     
                                 o                                                     p            #         @                                   q     	               #AXIS1 r   #AXIS2 s   #MEAXIS t   #CURRENTGRIDSHAPE u   #BASEDATTYPE v   #CURRENTLOCALLIMITS w             
                                 r                     
                                 s                     
                                t                     
                                u                    "   p          p            p                                    
                                v                     
     �                           w                    #    p          & p        p          p            p          p                          #         @                                   x     	               #IDIRECTION y   #SOLVEFLAG z   #INSIZE {   #LOCALSIZE |   #GLOBALSIZE }   #TRANSFORMTYPE ~   #INARRAY    #OUTARRAY �             
                                 y                     
                                 z                     
                                 {                     
                                 |                    $   p          p            p                                    
                                 }                    %   p          p            p                                    
                                 ~                    &   p          p            p                                   
                                                     
 '   p          5 O p            5 O p                                                                   �                    
 (    p          5 O p            5 O p                          #         @                                   �     	               #IDIRECTION �   #SOLVEFLAG �   #INSIZE �   #LOCALSIZE �   #GLOBALSIZE �   #TRANSFORMTYPE �   #INARRAY �   #OUTARRAY �             
                                 �                     
                                 �                     
                                 �                     
                                 �                    )   p          p            p                                    
                                 �                    *   p          p            p                                    
                                 �                    +   p          p            p                                   
                                 �                    
 ,   p          5 O p            5 O p                                                                   �                    
 -    p          5 O p            5 O p                          #         @                                   �     	               #IDIRECTION �   #ISRC �   #INSIZE �   #BCTYPES �   #BCVALUES �   #INARRAY �   #OUTARRAY �             
                                 �                     
                                 �                     
                                 �                     
                                 �                    .   p          p            p                                    
                                 �                   
 /   p          p            p                                   
                                �                    
 0    p          5 O p            5 O p                                                                   �                    
 1    p          5 O p            5 O p                          #         @                                   �     	               #LOWER �   #MAIN �   #UPPER �   #RHS �   #X �   #LENGTH �   #IERR �            
                                 �                    
 2   p          5 O p            5 O p                                   
                                 �                    
 3   p          5 O p            5 O p                                   
                                 �                    
 4   p          5 O p            5 O p                                   
                                 �                    
 5   p          5 O p            5 O p                                                                   �                    
 6    p          5 O p            5 O p                                    
                                �                                                     �            #         @                                   �     	            	   #LOWER �   #MAIN �   #UPPER �   #ALPHA �   #BETA �   #RHS �   #X �   #LENGTH �   #IERR �            
                                 �                    
 7   p          5 O p            5 O p                                   
                                 �                    
 8   p          5 O p            5 O p                                   
                                 �                    
 9   p          5 O p            5 O p                                    
                                 �     
                
                                 �     
               
                                 �                    
 :   p          5 O p            5 O p                                                                   �                    
 ;    p          5 O p            5 O p                                    
                                �                                                     �            #         @                                   �     	               #GRIDVAR �   #PFFTINPUTARRAY �             
                                 �                                                     �                   
 <              &                                           #         @                                   �     	               #GRIDVAR �   #PFFTOUTPUTARRAY �             
                                 �                                                     �                   
 =              &                                           #         @                                   �     	                #         @                                   �     	                #         @                                   �     	               #FLASHPROCID �   #FLASHBLOCKID �   #FLASHSTARTPOS �   #FLASHENDPOS �   #PFFTPROCID �   #PFFTSTARTPOS �   #PFFTENDPOS �             
                                 �                     
                                 �                     
      �                           �                    >   p          & p        p            p                                    
      �                           �                    ?   p          & p        p            p                                    
                                 �                     
      �                           �                    @   p          & p        p            p                                    
      �                           �                    A   p          & p        p            p                          #         @                                   �     	                #         @                                   �     	               #PFFT_INLEN �             
      �                           �                    B   p          & p        p            p                          &         @                                 �                          #GR_PFFTMAKEPENCILIN3DSPACE%POSSIBLEGRID_T �   #PENCILGLOBALLEN �   #TOTALPROCS �   #GR_PFFTFNARGCONSTRAINT �   #GR_PFFTMAKEPENCILIN3DSPACE%POSSIBLEGRID_T �                    @                           �     '                    #JPROCS �   #KPROCS �                �                               �                                �                               �                            
                                 �                    D   p          p            p                                    
                                 �           %         @                                 �                           #PENCILGLOBALLEN �   #TOTALPROCS �   #IPROCS �   #JPROCS �   #KPROCS �             
                                 �                    C   p          p            p                                    
                                 �                     
                                 �                     
                                 �                     
                                 �           %         @                                �                           #PENCILGLOBALLEN �   #TOTALPROCS �   #IPROCS �   #JPROCS �   #KPROCS �             
                                 �                    E   p          p            p                                    
                                 �                     
                                 �                     
                                 �                     
                                 �           %         @                                �                           #PENCILGLOBALLEN �   #TOTALPROCS �   #IPROCS �   #JPROCS �   #KPROCS �             
                                 �                    F   p          p            p                                    
                                 �                     
                                 �                     
                                 �                     
                                 �           %         @                                �                           #PENCILGLOBALLEN �   #TOTALPROCS �   #IPROCS �   #JPROCS �   #KPROCS �             
                                 �                    G   p          p            p                                    
                                 �                     
                                 �                     
                                 �                     
                                 �           #         @                                   �     	               #TRANSFORMTYPE �   #BASEDATTYPE �   #BCTYPES �                                             �                    H    p          p            p                                         �                           �                    I    p           & p         p            p                                    
                                �                    J   p          p            p                          #         @                                   �     	               #SOLVELEVEL �   #LEAFMAPMODE �             
                                 �                     
                                �           #         @                                   �     	               #INOUTLEVEL �             
                                �            #         @                                   �     	               #FLASHPROCID �   #FLASHBLOCKID �   #LBLOCKFRAGSTART �   #LBLOCKFRAGEND �   #LACTUALBLOCKFRAGSTART �   #LACTUALBLOCKFRAGEND �   #BLKTYPE �   #BLKREFLEV �   #SOLVELEVEL �   #PFFTPROC �   #LPENCILFRAGSTART �   #LPENCILFRAGEND �             
                                 �                     
                                 �                     
      �                           �                    K   p          & p        p            p                                    
      �                           �                    L   p          & p        p            p                                    
      �                           �                    M   p          & p        p            p                                    
      �                           �                    N   p          & p        p            p                                    
                                 �                     
                                 �                     
                                 �                     
                                 �                     
      �                           �                    O   p          & p        p            p                                    
      �                           �                    P   p          & p        p            p                             �   .      fn#fn %   �   q       gen@GR_PFFTTRANSPOSE    ?  H       GR_PFFTINIT     �  H       GR_PFFTFINALIZE #   �  �       GR_PFFTDCFTFORWARD +   x  �   a   GR_PFFTDCFTFORWARD%INARRAY ,     �   a   GR_PFFTDCFTFORWARD%OUTARRAY (   �  �   a   GR_PFFTDCFTFORWARD%TRIG '     @   a   GR_PFFTDCFTFORWARD%LEN '   \  @   a   GR_PFFTDCFTFORWARD%LDA *   �  @   a   GR_PFFTDCFTFORWARD%NUMVEC 1   �  @   a   GR_PFFTDCFTFORWARD%TRANSFORMTYPE )     @   a   GR_PFFTDCFTFORWARD%SCALE #   \  �       GR_PFFTDCFTINVERSE +     �   a   GR_PFFTDCFTINVERSE%INARRAY ,   �  �   a   GR_PFFTDCFTINVERSE%OUTARRAY (     �   a   GR_PFFTDCFTINVERSE%TRIG '   �  @   a   GR_PFFTDCFTINVERSE%LEN '   �  @   a   GR_PFFTDCFTINVERSE%LDA *   )  @   a   GR_PFFTDCFTINVERSE%NUMVEC 1   i  @   a   GR_PFFTDCFTINVERSE%TRANSFORMTYPE )   �  @   a   GR_PFFTDCFTINVERSE%SCALE     �  |       GR_PFFTDISPERSE (   e	  �   a   GR_PFFTDISPERSE%INARRAY )   �	  �   a   GR_PFFTDISPERSE%OUTARRAY #   }
  @   a   GR_PFFTDISPERSE%NX #   �
  @   a   GR_PFFTDISPERSE%NY $   �
  @   a   GR_PFFTDISPERSE%NX1    =  H       GR_PFFTGENMAP $   �  R       GR_PFFTINITMETADATA )   �  @   a   GR_PFFTINITMETADATA%NDIM #     |       GR_PFFTINTERSPERSE +   �  �   a   GR_PFFTINTERSPERSE%INARRAY ,     �   a   GR_PFFTINTERSPERSE%OUTARRAY &   �  @   a   GR_PFFTINTERSPERSE%NX &   �  @   a   GR_PFFTINTERSPERSE%NY '   +  @   a   GR_PFFTINTERSPERSE%NX1 &   k  �       GR_PFFTLOCALTRANSPOSE .   �  �   a   GR_PFFTLOCALTRANSPOSE%INARRAY /   |  �   a   GR_PFFTLOCALTRANSPOSE%OUTARRAY )     @   a   GR_PFFTLOCALTRANSPOSE%NX *   H  @   a   GR_PFFTLOCALTRANSPOSE%NX1 2   �  @   a   GR_PFFTLOCALTRANSPOSE%BASEDATTYPE     �  �       GR_PFFTSETUPDIM '   {  @   a   GR_PFFTSETUPDIM%LENGTH .   �  @   a   GR_PFFTSETUPDIM%TRANSFORMTYPE %   �  �   a   GR_PFFTSETUPDIM%TRIG &   �  @   a   GR_PFFTSETUPDIM%SCALE *   �  @   a   GR_PFFTSETUPDIM%INDATSIZE +     @   a   GR_PFFTSETUPDIM%OUTDATSIZE '   G  @   a   GR_PFFTSETUPDIM%DATINC '   �  @   a   GR_PFFTSETUPDIM%FACTOR !   �  �       GR_PFFTTRANSPOSE %   m  @   a   GR_PFFTTRANSPOSE%DIR -   �  @   a   GR_PFFTTRANSPOSE%BASEDATTYPE )   �  �   a   GR_PFFTTRANSPOSE%INARRAY *   y  �   a   GR_PFFTTRANSPOSE%OUTARRAY '     �   a   GR_PFFTTRANSPOSE%INLEN (   �  �   a   GR_PFFTTRANSPOSE%OUTLEN $   -  @   a   GR_PFFTTRANSPOSE%PR &   m  @   a   GR_PFFTTRANSPOSE%COMM &   �  �       GR_PFFTTRANSPOSE3DARR *   S  @   a   GR_PFFTTRANSPOSE3DARR%DIR 2   �  @   a   GR_PFFTTRANSPOSE3DARR%BASEDATTYPE .   �  �   a   GR_PFFTTRANSPOSE3DARR%INARRAY /   �  �   a   GR_PFFTTRANSPOSE3DARR%OUTARRAY ,   K  �   a   GR_PFFTTRANSPOSE3DARR%INLEN -   �  �   a   GR_PFFTTRANSPOSE3DARR%OUTLEN )   s  @   a   GR_PFFTTRANSPOSE3DARR%PR +   �  @   a   GR_PFFTTRANSPOSE3DARR%COMM &   �  ^       GR_PFFTGETLOCALLIMITS ,   Q  @   a   GR_PFFTGETLOCALLIMITS%AXIS1 ,   �  @   a   GR_PFFTGETLOCALLIMITS%AXIS2    �  H       GR_PFFTWAVE      X       GR_PFFTDERIVS )   q  �   a   GR_PFFTDERIVS%TRANSARRAY #   �  �       GR_PFFTGETPROCGRID (   �  @   a   GR_PFFTGETPROCGRID%DIMS ,   �  @   a   GR_PFFTGETPROCGRID%PFFTMYPE 0     @   a   GR_PFFTGETPROCGRID%PFFTNUMPROCS 1   T  �   a   GR_PFFTGETPROCGRID%PFFTGLOBALLEN 0   �  �   a   GR_PFFTGETPROCGRID%PFFTPROCGRID $   |   z       GR_PFFTGENMAPHELPER )   �   @   a   GR_PFFTGENMAPHELPER%AXIS 0   6!  �   a   GR_PFFTGENMAPHELPER%FRAGMENTPTR 6   �!  @   a   GR_PFFTGENMAPHELPER%MAXSINGLEPROCDATA &   "  |       GR_PFFTBUFFERTRANSFER 0   ~"  @   a   GR_PFFTBUFFERTRANSFER%DIRECTION +   �"  @   a   GR_PFFTBUFFERTRANSFER%AXIS -   �"  �   a   GR_PFFTBUFFERTRANSFER%BUFFER 0   �#  �   a   GR_PFFTBUFFERTRANSFER%PFFTARRAY )   .$  j       GR_PFFTGETDESTPFFTCOORDS 2   �$  �   a   GR_PFFTGETDESTPFFTCOORDS%STARTPOS 8   <%  �   a   GR_PFFTGETDESTPFFTCOORDS%PFFTPROCCOORDS %   �%  �       GR_PFFTCOPYTOSENDMAP -   n&  @   a   GR_PFFTCOPYTOSENDMAP%BLOCKID 3   �&  �   a   GR_PFFTCOPYTOSENDMAP%BLOCKSTARTPOS 1   R'  �   a   GR_PFFTCOPYTOSENDMAP%BLOCKENDPOS +   �'  @   a   GR_PFFTCOPYTOSENDMAP%BSIZE *   6(  @   a   GR_PFFTCOPYTOSENDMAP%AXIS ,   v(  d       GR_PFFTHANDLEJAXISFRAGMENTS 6   �(  @   a   GR_PFFTHANDLEJAXISFRAGMENTS%DIRECTION 4   )  @   a   GR_PFFTHANDLEJAXISFRAGMENTS%GRIDVAR ,   Z)  W       GR_PFFTHANDLEKAXISFRAGMENTS 6   �)  @   a   GR_PFFTHANDLEKAXISFRAGMENTS%DIRECTION (   �)  j       GR_PFFTPRINTCOMMBUFFERS /   [*  �   a   GR_PFFTPRINTCOMMBUFFERS%BUFFER ,   �*  �   a   GR_PFFTPRINTCOMMBUFFERS%MAP 0   �+  @   a   GR_PFFTPRINTCOMMBUFFERS%LOGUNIT ,   �+  �       GR_PFFTGROUPUSABLEPROCESSES 1   e,  @   a   GR_PFFTGROUPUSABLEPROCESSES%MYPE 8   �,  @   a   GR_PFFTGROUPUSABLEPROCESSES%GLOBALPROCS 9   �,  @   a   GR_PFFTGROUPUSABLEPROCESSES%ORIGINALCOMM 4   %-  @   a   GR_PFFTGROUPUSABLEPROCESSES%NEWCOMM -   e-  �       GR_PFFTGETLOCALLIMITSANYTIME 3   .  @   a   GR_PFFTGETLOCALLIMITSANYTIME%AXIS1 3   N.  @   a   GR_PFFTGETLOCALLIMITSANYTIME%AXIS2 4   �.  @   a   GR_PFFTGETLOCALLIMITSANYTIME%MEAXIS >   �.  �   a   GR_PFFTGETLOCALLIMITSANYTIME%CURRENTGRIDSHAPE 9   b/  @   a   GR_PFFTGETLOCALLIMITSANYTIME%BASEDATTYPE @   �/  �   a   GR_PFFTGETLOCALLIMITSANYTIME%CURRENTLOCALLIMITS %   f0  �       GR_PFFTPOISSONDIRECT 0   &1  @   a   GR_PFFTPOISSONDIRECT%IDIRECTION /   f1  @   a   GR_PFFTPOISSONDIRECT%SOLVEFLAG ,   �1  @   a   GR_PFFTPOISSONDIRECT%INSIZE /   �1  �   a   GR_PFFTPOISSONDIRECT%LOCALSIZE 0   z2  �   a   GR_PFFTPOISSONDIRECT%GLOBALSIZE 3   3  �   a   GR_PFFTPOISSONDIRECT%TRANSFORMTYPE -   �3  �   a   GR_PFFTPOISSONDIRECT%INARRAY .   F4  �   a   GR_PFFTPOISSONDIRECT%OUTARRAY )   �4  �       GR_PFFTPOISSONTRIGDIRECT 4   �5  @   a   GR_PFFTPOISSONTRIGDIRECT%IDIRECTION 3   �5  @   a   GR_PFFTPOISSONTRIGDIRECT%SOLVEFLAG 0   *6  @   a   GR_PFFTPOISSONTRIGDIRECT%INSIZE 3   j6  �   a   GR_PFFTPOISSONTRIGDIRECT%LOCALSIZE 4   �6  �   a   GR_PFFTPOISSONTRIGDIRECT%GLOBALSIZE 7   �7  �   a   GR_PFFTPOISSONTRIGDIRECT%TRANSFORMTYPE 1   &8  �   a   GR_PFFTPOISSONTRIGDIRECT%INARRAY 2   �8  �   a   GR_PFFTPOISSONTRIGDIRECT%OUTARRAY (   n9  �       GR_PFFTPOISSONHOMBCTRIG 3   :  @   a   GR_PFFTPOISSONHOMBCTRIG%IDIRECTION -   R:  @   a   GR_PFFTPOISSONHOMBCTRIG%ISRC /   �:  @   a   GR_PFFTPOISSONHOMBCTRIG%INSIZE 0   �:  �   a   GR_PFFTPOISSONHOMBCTRIG%BCTYPES 1   f;  �   a   GR_PFFTPOISSONHOMBCTRIG%BCVALUES 0   �;  �   a   GR_PFFTPOISSONHOMBCTRIG%INARRAY 1   �<  �   a   GR_PFFTPOISSONHOMBCTRIG%OUTARRAY    B=  �       GR_PFFTTRIDIAG %   �=  �   a   GR_PFFTTRIDIAG%LOWER $   t>  �   a   GR_PFFTTRIDIAG%MAIN %   ?  �   a   GR_PFFTTRIDIAG%UPPER #   �?  �   a   GR_PFFTTRIDIAG%RHS !   `@  �   a   GR_PFFTTRIDIAG%X &   A  @   a   GR_PFFTTRIDIAG%LENGTH $   DA  @   a   GR_PFFTTRIDIAG%IERR %   �A  �       GR_PFFTCYCLICTRIDIAG +   'B  �   a   GR_PFFTCYCLICTRIDIAG%LOWER *   �B  �   a   GR_PFFTCYCLICTRIDIAG%MAIN +   oC  �   a   GR_PFFTCYCLICTRIDIAG%UPPER +   D  @   a   GR_PFFTCYCLICTRIDIAG%ALPHA *   SD  @   a   GR_PFFTCYCLICTRIDIAG%BETA )   �D  �   a   GR_PFFTCYCLICTRIDIAG%RHS '   7E  �   a   GR_PFFTCYCLICTRIDIAG%X ,   �E  @   a   GR_PFFTCYCLICTRIDIAG%LENGTH *   F  @   a   GR_PFFTCYCLICTRIDIAG%IERR "   [F  i       GR_PFFTMAPTOINPUT *   �F  @   a   GR_PFFTMAPTOINPUT%GRIDVAR 1   G  �   a   GR_PFFTMAPTOINPUT%PFFTINPUTARRAY %   �G  j       GR_PFFTMAPFROMOUTPUT -   �G  @   a   GR_PFFTMAPFROMOUTPUT%GRIDVAR 5   :H  �   a   GR_PFFTMAPFROMOUTPUT%PFFTOUTPUTARRAY )   �H  H       GR_PFFTINITIALISESTORAGE '   I  H       GR_PFFTFINALISESTORAGE &   VI  �       GR_PFFTCREATESENDNODE 2   J  @   a   GR_PFFTCREATESENDNODE%FLASHPROCID 3   WJ  @   a   GR_PFFTCREATESENDNODE%FLASHBLOCKID 4   �J  �   a   GR_PFFTCREATESENDNODE%FLASHSTARTPOS 2   ;K  �   a   GR_PFFTCREATESENDNODE%FLASHENDPOS 1   �K  @   a   GR_PFFTCREATESENDNODE%PFFTPROCID 3   L  �   a   GR_PFFTCREATESENDNODE%PFFTSTARTPOS 1   �L  �   a   GR_PFFTCREATESENDNODE%PFFTENDPOS /   gM  H       GR_PFFTCOMMUNICATENODEMETADATA &   �M  X       GR_PFFTGRIDPOINTTABLE 1   N  �   a   GR_PFFTGRIDPOINTTABLE%PFFT_INLEN +   �N  �       GR_PFFTMAKEPENCILIN3DSPACE S   �O  h      GR_PFFTMAKEPENCILIN3DSPACE%POSSIBLEGRID_T+GR_PFFTINTERFACETYPEDECL Z   P  H   a   GR_PFFTMAKEPENCILIN3DSPACE%POSSIBLEGRID_T%JPROCS+GR_PFFTINTERFACETYPEDECL Z   JP  H   a   GR_PFFTMAKEPENCILIN3DSPACE%POSSIBLEGRID_T%KPROCS+GR_PFFTINTERFACETYPEDECL ;   �P  �   a   GR_PFFTMAKEPENCILIN3DSPACE%PENCILGLOBALLEN 6   &Q  @   a   GR_PFFTMAKEPENCILIN3DSPACE%TOTALPROCS B   fQ  �      GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT R   �Q  �   a   GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT%PENCILGLOBALLEN M   �R  @   a   GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT%TOTALPROCS I   �R  @   a   GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT%IPROCS I   S  @   a   GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT%JPROCS I   SS  @   a   GR_PFFTMAKEPENCILIN3DSPACE%GR_PFFTFNARGCONSTRAINT%KPROCS +   �S  �       GR_PFFTFNARGHARDCONSTRAINT ;   ,T  �   a   GR_PFFTFNARGHARDCONSTRAINT%PENCILGLOBALLEN 6   �T  @   a   GR_PFFTFNARGHARDCONSTRAINT%TOTALPROCS 2    U  @   a   GR_PFFTFNARGHARDCONSTRAINT%IPROCS 2   @U  @   a   GR_PFFTFNARGHARDCONSTRAINT%JPROCS 2   �U  @   a   GR_PFFTFNARGHARDCONSTRAINT%KPROCS -   �U  �       GR_PFFTFNARGMEDIUMCONSTRAINT =   YV  �   a   GR_PFFTFNARGMEDIUMCONSTRAINT%PENCILGLOBALLEN 8   �V  @   a   GR_PFFTFNARGMEDIUMCONSTRAINT%TOTALPROCS 4   -W  @   a   GR_PFFTFNARGMEDIUMCONSTRAINT%IPROCS 4   mW  @   a   GR_PFFTFNARGMEDIUMCONSTRAINT%JPROCS 4   �W  @   a   GR_PFFTFNARGMEDIUMCONSTRAINT%KPROCS +   �W  �       GR_PFFTFNARGEASYCONSTRAINT ;   �X  �   a   GR_PFFTFNARGEASYCONSTRAINT%PENCILGLOBALLEN 6   Y  @   a   GR_PFFTFNARGEASYCONSTRAINT%TOTALPROCS 2   ZY  @   a   GR_PFFTFNARGEASYCONSTRAINT%IPROCS 2   �Y  @   a   GR_PFFTFNARGEASYCONSTRAINT%JPROCS 2   �Y  @   a   GR_PFFTFNARGEASYCONSTRAINT%KPROCS (   Z  y       GR_PFFTSPECIFYTRANSFORM 6   �Z  �   a   GR_PFFTSPECIFYTRANSFORM%TRANSFORMTYPE 4   '[  �   a   GR_PFFTSPECIFYTRANSFORM%BASEDATTYPE 0   �[  �   a   GR_PFFTSPECIFYTRANSFORM%BCTYPES $   _\  i       GR_PFFTGENSINGLEMAP /   �\  @   a   GR_PFFTGENSINGLEMAP%SOLVELEVEL 0   ]  @   a   GR_PFFTGENSINGLEMAP%LEAFMAPMODE -   H]  X       GR_PFFTVALIDATESELECTEDLEVEL 8   �]  @   a   GR_PFFTVALIDATESELECTEDLEVEL%INOUTLEVEL *   �]  +      GR_PFFTCREATESENDFRAGMENT 6   _  @   a   GR_PFFTCREATESENDFRAGMENT%FLASHPROCID 7   K_  @   a   GR_PFFTCREATESENDFRAGMENT%FLASHBLOCKID :   �_  �   a   GR_PFFTCREATESENDFRAGMENT%LBLOCKFRAGSTART 8   /`  �   a   GR_PFFTCREATESENDFRAGMENT%LBLOCKFRAGEND @   �`  �   a   GR_PFFTCREATESENDFRAGMENT%LACTUALBLOCKFRAGSTART >   wa  �   a   GR_PFFTCREATESENDFRAGMENT%LACTUALBLOCKFRAGEND 2   b  @   a   GR_PFFTCREATESENDFRAGMENT%BLKTYPE 4   [b  @   a   GR_PFFTCREATESENDFRAGMENT%BLKREFLEV 5   �b  @   a   GR_PFFTCREATESENDFRAGMENT%SOLVELEVEL 3   �b  @   a   GR_PFFTCREATESENDFRAGMENT%PFFTPROC ;   c  �   a   GR_PFFTCREATESENDFRAGMENT%LPENCILFRAGSTART 9   �c  �   a   GR_PFFTCREATESENDFRAGMENT%LPENCILFRAGEND 