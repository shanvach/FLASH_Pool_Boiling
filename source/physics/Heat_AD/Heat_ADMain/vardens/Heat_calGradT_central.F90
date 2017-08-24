subroutine Heat_calGradT_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny,mflg)

         use Heat_AD_data

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2

    real :: th,tol
    integer :: i,j,k
    real :: Tij,Tx_plus,Tx_mins,Ty_plus,Ty_mins,Tx,Ty,Tz

    tol = 0.01

    k = 1

    do i=ix1,ix2
      do j=jy1,jy2

         Tx_plus = T(i+1,j,k)
         Tx_mins = T(i-1,j,k)
         Ty_plus = T(i,j+1,k)
         Ty_mins = T(i,j-1,k)
         Tij     = T(i,j,k)

         ! Case 1 !
         if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then

             th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))))

             if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then
             Tx_plus = (ht_Tsat-Tij)/th + Tij

             else
             Tx_plus = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T(i-1,j,k))/(th + th*th)

             end if

         end if
         ! End of Case 1 !


         ! Case 2 !
         if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then

             th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))

             if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then
             Tx_mins = (ht_Tsat-Tij)/th + Tij

             else
             Tx_mins = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T(i+1,j,k))/(th + th*th)

             end if

         end if
         ! End of Case 2 !

         ! Case 3 !
         if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then

             th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))

             if(s(i,j,k)*s(i,j-1,k) .le. 0.0d0) then
             Ty_plus = (ht_Tsat-Tij)/th + Tij

             else
             Ty_plus = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T(i,j-1,k))/(th + th*th)

             end if

         end if
         ! End of Case 3 !

         ! Case 4 !
         if(s(i,j,k)*s(i,j-1,k) .le. 0.d0) then

             th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))))

             if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then
             Ty_mins = (ht_Tsat-Tij)/th + Tij

             else
             Ty_mins = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T(i,j+1,k))/(th + th*th)

             end if

         end if
         ! End of Case 4 ! 

         Tx = (Tx_plus - Tx_mins)/(2*dx)
         Ty = (Ty_plus - Ty_mins)/(2*dy)

         if (pf(i,j,k) .eq. 0.) then
          
            Tnl(i,j,k) = ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty) 
         
         else
         
            Tnv(i,j,k) = ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty)
         
         end if
   
      end do
    end do
#endif

end subroutine Heat_calGradT_central
