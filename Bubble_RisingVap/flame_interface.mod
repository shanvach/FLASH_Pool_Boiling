  �  0   k820309    w          19.1        6Gb                                                                                                          
       Flame_interface.F90 FLAME_INTERFACE #         @                                        	               #NUM_BLOCKS    #BLOCKLIST    #DT              
                                                     
                                                         p          5 O p            5 O p                                    
                                      
      #         @                                        	               #X    #F              
                                      
                                                     
       #         @                                        	               #LAMINARWIDTH 	                                             	     
       #         @                                   
     	                #         @                                        	               #DENS    #PRES    #S              
                                      
                
                                      
                                                     
       #         @                                        	               #Q    #FLAG                                                   
                 
                                            #         @                                        	               #DENS_U    #PRES_U    #TEMP_U    #ENER_U    #YE_U    #SUMY_U    #DENS_B    #PRES_B    #TEMP_B    #ENER_B    #YE_B    #SUMY_B    #Q    #S     #FLAG !             
                                     
                 
                                     
                 
                                     
                 
                                     
                 
                                      
                
                                      
                                                     
                                                      
                                                      
                                                      
                 
                                      
                
                                      
                
                                      
                
                                       
                
                                 !           #         @                                   "     	               #EOSDATA_U #   #QBAR_U $   #EOSDATA_B %   #QBAR_B &   #EOS_MODE '             
                                #                   
     p          p            p                                    
                                 $     
                                                %                   
     p          p            p                                                                    &     
                 
                                 '           #         @                                   (     	               #SOLNDATA )   #PHI1DOT *   #BLKLIMITS +   #BLKLIMITSGC ,   #TIME -   #DT .   #BLOCKID /                                           )                   
               &                   &                   &                   &                                                     
      �                           *     �             
    p �        & p        p          & p        p          & p        p            p          p          p                                    
      �                           +                       p          & p        p          p            p          p                                    
      �                           ,                       p          & p        p          p            p          p                                    
                                 -     
                
                                 .     
                
                                 /              �   ,      fn#fn    �   o       FLAME_STEP &   ;  @   a   FLAME_STEP%NUM_BLOCKS %   {  �   a   FLAME_STEP%BLOCKLIST      @   a   FLAME_STEP%DT !   _  V       FLAME_GETPROFILE #   �  @   a   FLAME_GETPROFILE%X #   �  @   a   FLAME_GETPROFILE%F    5  Z       FLAME_GETWIDTH ,   �  @   a   FLAME_GETWIDTH%LAMINARWIDTH    �  H       FLAME_INIT #     c       FLAME_LAMINARSPEED (   z  @   a   FLAME_LAMINARSPEED%DENS (   �  @   a   FLAME_LAMINARSPEED%PRES %   �  @   a   FLAME_LAMINARSPEED%S "   :  Y       FLAME_HEATRELEASE $   �  @   a   FLAME_HEATRELEASE%Q '   �  @   a   FLAME_HEATRELEASE%FLAG      �       FLAME_RHJUMP $   �  @   a   FLAME_RHJUMP%DENS_U $   ?  @   a   FLAME_RHJUMP%PRES_U $     @   a   FLAME_RHJUMP%TEMP_U $   �  @   a   FLAME_RHJUMP%ENER_U "   �  @   a   FLAME_RHJUMP%YE_U $   ?  @   a   FLAME_RHJUMP%SUMY_U $     @   a   FLAME_RHJUMP%DENS_B $   �  @   a   FLAME_RHJUMP%PRES_B $   �  @   a   FLAME_RHJUMP%TEMP_B $   ?	  @   a   FLAME_RHJUMP%ENER_B "   	  @   a   FLAME_RHJUMP%YE_B $   �	  @   a   FLAME_RHJUMP%SUMY_B    �	  @   a   FLAME_RHJUMP%Q    ?
  @   a   FLAME_RHJUMP%S "   
  @   a   FLAME_RHJUMP%FLAG %   �
  �       FLAME_RHJUMPREACTIVE /   K  �   a   FLAME_RHJUMPREACTIVE%EOSDATA_U ,   �  @   a   FLAME_RHJUMPREACTIVE%QBAR_U /     �   a   FLAME_RHJUMPREACTIVE%EOSDATA_B ,   �  @   a   FLAME_RHJUMPREACTIVE%QBAR_B .   �  @   a   FLAME_RHJUMPREACTIVE%EOS_MODE    3  �       FLAME_EFFECTS '   �  �   a   FLAME_EFFECTS%SOLNDATA &   �    a   FLAME_EFFECTS%PHI1DOT (   �  �   a   FLAME_EFFECTS%BLKLIMITS *   q  �   a   FLAME_EFFECTS%BLKLIMITSGC #   5  @   a   FLAME_EFFECTS%TIME !   u  @   a   FLAME_EFFECTS%DT &   �  @   a   FLAME_EFFECTS%BLOCKID 