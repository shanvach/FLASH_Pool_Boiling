  �-  Q   k820309              18.0        vi�`                                                                                                          
       gr_bcInterface.F90 GR_BCINTERFACE #         @                                        	               #AXIS    #ISWORK              
                                                      
                                            #         @                                        	            	   #AXIS    #BCTYPE    #GRIDDATASTRUCT    #VARCOUNT    #REGIONTYPE 	   #BLKLIMITS 
   #BLKLIMITSGC    #BLOCKHANDLE    #IDEST              
                                                      
                                                      
                                                      
                                                      
                                 	                       p          p            p                                    
      �                           
                       p          & p        p          p            p          p                                    
      �                                                  p          & p        p          p            p          p                                    
                                                      
                                            #         @                                        	               #GRIDDATASTRUCT    #AXIS    #ENDPOINTS    #REGIONSIZE    #MASK    #REGION    #BLOCKID    #IDEST              
                                                      
                                                      
      �                                                  p          & p        p          p            p          p                                    
                                                        p          p            p                                                                                            p          5 2 O p          n                                           4      5 2 O p          n                                      4                                                                                      
           p        5 2 O p          n                                       3  p        5 2 O p          n                                          2  p        5 2 O p          n                                          1  p          5 2 O p          n                                      1    5 2 O p          n                                          2    5 2 O p          n                                          3    5 2 O p          n                                          4      5 2 O p          n                                      1    5 2 O p          n                                          2    5 2 O p          n                                          3    5 2 O p          n                                          4                                                    
                                                      
                                            #         @                                        	               #GRIDDATASTRUCT    #AXIS    #SECONDDIR    #THIRDDIR    #ENDPOINTS    #REGIONSIZE    #REGIONDATAC    #REGIONDATAFN    #REGIONDATAFT1     #REGIONDATAFT2 !   #BLOCKID "   #IDEST #             
                                                      
                                                      
                                                      
                                                      
      �                                                  p          & p        p          p            p          p                                    
                                                     	   p          p            p                                                                                     
               &                   &                   &                   &                                                                                                      
 
              &                   &                   &                   &                                                                                                       
               &                   &                   &                   &                                                                                   !                   
               &                   &                   &                   &                                                     
                                 "                     
                                 #           #         @                                   $     	               #GRIDDATASTRUCT %   #AXIS &   #ENDPOINTS '   #REGIONSIZE (   #MASK )   #REGION *   #BLOCKID +   #IDEST ,             
                                 %                     
                                 &                     
      �                           '                       p          & p        p          p            p          p                                    
                                 (                       p          p            p                                   
                                 )                        p          5 2 O p          n                                           4      5 2 O p          n                                      4                                  
                                 *                    
          p        5 2 O p          n                                       3  p        5 2 O p          n                                          2  p        5 2 O p          n                                          1  p          5 2 O p          n                                      1    5 2 O p          n                                          2    5 2 O p          n                                          3    5 2 O p          n                                          4      5 2 O p          n                                      1    5 2 O p          n                                          2    5 2 O p          n                                          3    5 2 O p          n                                          4                                                    
                                 +                     
                                 ,           #         @                                   -     	               #GRIDDATASTRUCT .   #AXIS /   #SECONDDIR 0   #THIRDDIR 1   #ENDPOINTS 2   #REGIONSIZE 3   #REGIONDATAC 4   #REGIONDATAFN 5   #REGIONDATAFT1 6   #REGIONDATAFT2 7   #BLOCKID 8   #IDEST 9             
                                 .                     
                                 /                     
                                 0                     
                                 1                     
      �                           2                       p          & p        p          p            p          p                                    
                                 3                       p          p            p                                                                  4                   
               &                   &                   &                   &                                                                                   5                   
               &                   &                   &                   &                                                                                   6                   
               &                   &                   &                   &                                                                                   7                   
               &                   &                   &                   &                                                     
                                 8                     
                                 9           #         @                                   :     	               #BCTYPETOAPPLY ;   #BCTYPEFROMGRID <   #VARINDEX =   #GRIDDATASTRUCT >   #AXIS ?   #FACE @   #IDEST A                                             ;                      
                                 <                     
                                 =                     
                                 >                     
                                 ?                     
                                 @                     
                                A           #         @                                   B     	                #         @                                   C     	                #         @                                   D     	               #DENS E   #TEMP F   #YE G   #SUMY H   #N I   #INPUTG J   #DELTA K   #DIRECTION L   #ORDER M   #MODE N   #MASSFRAC O             
                                E                   
               &                                                     
                                F                   
               &                                                     
                                 G                   
              &                                                     
                                 H                   
              &                                                     
                                 I                     
                                 J     
                
                                 K     
                
                                 L                     
                                 M                     
                                 N                     
                                O                   
              &                                           #         @                                   P     	                   �   *      fn#fn $   �   ^       GR_BCAPPLYTOALLBLKS )   (  @   a   GR_BCAPPLYTOALLBLKS%AXIS +   h  @   a   GR_BCAPPLYTOALLBLKS%ISWORK $   �  �       GR_BCAPPLYTOONEFACE )   t  @   a   GR_BCAPPLYTOONEFACE%AXIS +   �  @   a   GR_BCAPPLYTOONEFACE%BCTYPE 3   �  @   a   GR_BCAPPLYTOONEFACE%GRIDDATASTRUCT -   4  @   a   GR_BCAPPLYTOONEFACE%VARCOUNT /   t  �   a   GR_BCAPPLYTOONEFACE%REGIONTYPE .     �   a   GR_BCAPPLYTOONEFACE%BLKLIMITS 0   �  �   a   GR_BCAPPLYTOONEFACE%BLKLIMITSGC 0   �  @   a   GR_BCAPPLYTOONEFACE%BLOCKHANDLE *   �  @   a   GR_BCAPPLYTOONEFACE%IDEST      �       GR_BCGETREGION .   �  @   a   GR_BCGETREGION%GRIDDATASTRUCT $     @   a   GR_BCGETREGION%AXIS )   C  �   a   GR_BCGETREGION%ENDPOINTS *     �   a   GR_BCGETREGION%REGIONSIZE $   �  &  a   GR_BCGETREGION%MASK &   �	  w  a   GR_BCGETREGION%REGION '   8  @   a   GR_BCGETREGION%BLOCKID %   x  @   a   GR_BCGETREGION%IDEST (   �        GR_BCGETREGIONSMIXEDGDS 7   �  @   a   GR_BCGETREGIONSMIXEDGDS%GRIDDATASTRUCT -   �  @   a   GR_BCGETREGIONSMIXEDGDS%AXIS 2   ;  @   a   GR_BCGETREGIONSMIXEDGDS%SECONDDIR 1   {  @   a   GR_BCGETREGIONSMIXEDGDS%THIRDDIR 2   �  �   a   GR_BCGETREGIONSMIXEDGDS%ENDPOINTS 3     �   a   GR_BCGETREGIONSMIXEDGDS%REGIONSIZE 4     �   a   GR_BCGETREGIONSMIXEDGDS%REGIONDATAC 5   �  �   a   GR_BCGETREGIONSMIXEDGDS%REGIONDATAFN 6   �  �   a   GR_BCGETREGIONSMIXEDGDS%REGIONDATAFT1 6   �  �   a   GR_BCGETREGIONSMIXEDGDS%REGIONDATAFT2 0   c  @   a   GR_BCGETREGIONSMIXEDGDS%BLOCKID .   �  @   a   GR_BCGETREGIONSMIXEDGDS%IDEST    �  �       GR_BCPUTREGION .   �  @   a   GR_BCPUTREGION%GRIDDATASTRUCT $   �  @   a   GR_BCPUTREGION%AXIS )     �   a   GR_BCPUTREGION%ENDPOINTS *   �  �   a   GR_BCPUTREGION%REGIONSIZE $   n  &  a   GR_BCPUTREGION%MASK &   �  w  a   GR_BCPUTREGION%REGION '     @   a   GR_BCPUTREGION%BLOCKID %   K  @   a   GR_BCPUTREGION%IDEST (   �        GR_BCPUTREGIONSMIXEDGDS 7   �  @   a   GR_BCPUTREGIONSMIXEDGDS%GRIDDATASTRUCT -   �  @   a   GR_BCPUTREGIONSMIXEDGDS%AXIS 2      @   a   GR_BCPUTREGIONSMIXEDGDS%SECONDDIR 1   N   @   a   GR_BCPUTREGIONSMIXEDGDS%THIRDDIR 2   �   �   a   GR_BCPUTREGIONSMIXEDGDS%ENDPOINTS 3   R!  �   a   GR_BCPUTREGIONSMIXEDGDS%REGIONSIZE 4   �!  �   a   GR_BCPUTREGIONSMIXEDGDS%REGIONDATAC 5   �"  �   a   GR_BCPUTREGIONSMIXEDGDS%REGIONDATAFN 6   �#  �   a   GR_BCPUTREGIONSMIXEDGDS%REGIONDATAFT1 6   b$  �   a   GR_BCPUTREGIONSMIXEDGDS%REGIONDATAFT2 0   6%  @   a   GR_BCPUTREGIONSMIXEDGDS%BLOCKID .   v%  @   a   GR_BCPUTREGIONSMIXEDGDS%IDEST    �%  �       GR_BCMAPBCTYPE -   f&  @   a   GR_BCMAPBCTYPE%BCTYPETOAPPLY .   �&  @   a   GR_BCMAPBCTYPE%BCTYPEFROMGRID (   �&  @   a   GR_BCMAPBCTYPE%VARINDEX .   &'  @   a   GR_BCMAPBCTYPE%GRIDDATASTRUCT $   f'  @   a   GR_BCMAPBCTYPE%AXIS $   �'  @   a   GR_BCMAPBCTYPE%FACE %   �'  @   a   GR_BCMAPBCTYPE%IDEST    &(  H       GR_BCINIT    n(  H       GR_BCHSEINIT    �(  �       GR_HSESTEP     t)  �   a   GR_HSESTEP%DENS      *  �   a   GR_HSESTEP%TEMP    �*  �   a   GR_HSESTEP%YE     +  �   a   GR_HSESTEP%SUMY    �+  @   a   GR_HSESTEP%N "   �+  @   a   GR_HSESTEP%INPUTG !   $,  @   a   GR_HSESTEP%DELTA %   d,  @   a   GR_HSESTEP%DIRECTION !   �,  @   a   GR_HSESTEP%ORDER     �,  @   a   GR_HSESTEP%MODE $   $-  �   a   GR_HSESTEP%MASSFRAC    �-  H       GR_BCFINALIZE 