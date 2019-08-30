subroutine Heat_calGradT_3D_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,nx,ny,nz,mflg,lambda)
        
     use Heat_AD_data
     use Heat_AD_interface, only: Heat_GFMstencil_o1, Heat_GFMstencil_o2

     implicit none
     real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
     real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg,lambda
     real, intent(in) :: dx,dy,dz
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

     real :: th,tol
     integer :: i,j,k
     real :: Tx_plus,Tx_mins,Ty_plus,Ty_mins,Tz_plus,Tz_mins
     real :: Tx,Ty,Tz,Tij
     real :: thxp1, thxm1, thyp1, thym1, thzp1, thzm1
     real :: thxp2, thxm2, thyp2, thym2, thzp2, thzm2

     tol = 0.01

   do k=kz1,kz2
    do j=jy1,jy2
     do i=ix1,ix2

     Tx_plus = T(i+1,j,k)
     Tx_mins = T(i-1,j,k)
     Ty_plus = T(i,j+1,k)
     Ty_mins = T(i,j-1,k)
     Tz_plus = T(i,j,k+1)
     Tz_mins = T(i,j,k-1)

     Tij  = T(i,j,k)

     thxp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
     thxm1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))
     thyp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))
     thym1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))
     thzp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1)))
     thzm1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1)))
        
     thxp2 = abs(s(i-1,j,k))/(abs(s(i-1,j,k))+abs(s(i+1,j,k)))
     thxm2 = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i-1,j,k)))
     thyp2 = abs(s(i,j-1,k))/(abs(s(i,j-1,k))+abs(s(i,j+1,k)))
     thym2 = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j-1,k)))
     thzp2 = abs(s(i,j,k-1))/(abs(s(i,j,k-1))+abs(s(i,j,k+1)))
     thzm2 = abs(s(i,j,k+1))/(abs(s(i,j,k+1))+abs(s(i,j,k-1)))

      !______________________Diffusion Terms_______________________!

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i+1,j,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Tx_plus,Tij,ht_Tsat,max(tol,thxp1))
     ! End of Case 1 !

     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i-1,j,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Tx_mins,Tij,ht_Tsat,max(tol,thxm1))
     ! End of Case 2 !

     ! Case 3 !
     if(s(i,j,k)*s(i,j+1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j+1,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Ty_plus,Tij,ht_Tsat,max(tol,thyp1))
     ! End of Case 3 !

     ! Case 4 !
     if(s(i,j,k)*s(i,j-1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j-1,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Ty_mins,Tij,ht_Tsat,max(tol,thym1))
     ! End of Case 4 ! 
 
     ! Case 5 !
     if(s(i,j,k)*s(i,j,k+1) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j,k+1) .lt. 0.0) &
     call Heat_GFMstencil_o1(Tz_plus,Tij,ht_Tsat,max(tol,thzp1))
     ! End of Case 5 !

     ! Case 6 !
     if(s(i,j,k)*s(i,j,k-1) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j,k-1) .lt. 0.0) &
     call Heat_GFMstencil_o1(Tz_mins,Tij,ht_Tsat,max(tol,thzm1))
     ! End of Case 6 ! 


     Tx = (Tx_plus - Tx_mins)/(2*dx)
     Ty = (Ty_plus - Ty_mins)/(2*dy)
     Tz = (Tz_plus - Tz_mins)/(2*dz)

     !--------------------------------------------------------------------------!
     !----------------------------Calculate Fluxes------------------------------!
     !--------------------------------------------------------------------------!

     if(lambda(i,j,k) .lt. 0.0) then
     if (pf(i,j,k) .eq. 0.) then
          
         Tnl(i,j,k) = ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty + nz(i,j,k)*Tz)
         
     else
         
         Tnv(i,j,k) = ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty - nz(i,j,k)*Tz)
         
     end if
     end if                          

     end do
    end do
   end do

end subroutine Heat_calGradT_3D_central
