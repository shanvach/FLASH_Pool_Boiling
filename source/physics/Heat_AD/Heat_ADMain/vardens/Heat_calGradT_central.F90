subroutine Heat_calGradT_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny,mflg)
        
    use Heat_AD_data

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2

    real :: th,tol,th2
    integer :: i,j,k
    real :: Tij,Tipj,Timj,Tijp,Tijm,dxp,dxm,dyp,dym,Tx,Ty,Tax,Tbx,Tay,Tby,dyc,dxc,Tpx,Tmx,Tpy,Tmy
    logical :: int_xp,int_xm,int_yp,int_ym

    tol = 0.01

    k = 1

    do i=ix1,ix2
      do j=jy1,jy2

         Tipj = T(i+1,j,k)
         Timj = T(i-1,j,k)
         Tijp = T(i,j+1,k)
         Tijm = T(i,j-1,k)
         Tij  = T(i,j,k)
         Tmx  = T(i,j,k)
         Tpx  = T(i,j,k)
         Tmy  = T(i,j,k)
         Tpy  = T(i,j,k)

         dxp = dx
         dxm = dx
         dym = dy
         dyp = dy

         int_xp = .FALSE.
         int_xm = .FALSE.
         int_yp = .FALSE.
         int_ym = .FALSE.

         if(s(i,j,k)*s(i+1,j,k) .le. 0.) then

            Tpx = T(i,j,k)
            th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))))
      
            dxp = th*dx
            Tipj = ht_Tsat

            int_xp = .TRUE.

         end if

         if(s(i,j,k)*s(i-1,j,k) .le. 0.) then

            Tmx = T(i,j,k)
            th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))

            dxm = th*dx
            Timj = ht_Tsat

            int_xm = .TRUE.

         end if

         if(s(i,j,k)*s(i,j+1,k) .le. 0.) then

            Tpy = T(i,j,k)
            th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))

            dyp = th*dy
            Tijp = ht_Tsat

            int_yp = .TRUE.

         end if

         if(s(i,j,k)*s(i,j-1,k) .le. 0.) then

            Tmy = T(i,j,k)
            th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))))

            dym = th*dy
            Tijm = ht_Tsat

            int_ym = .TRUE.

         end if

         Tx = ((-dxp**2)*(Timj) + (dxp**2)*(Tmx) - (dxm**2)*(Tpx) + (dxm**2)*(Tipj))/(dxp*dxm*(dxp+dxm))
         Ty = ((-dyp**2)*(Tijm) + (dyp**2)*(Tmy) - (dym**2)*(Tpy) + (dym**2)*(Tijp))/(dyp*dym*(dyp+dym))

         !if(int_xm .or. int_xp .or. int_ym .or. int_yp) then

         if (pf(i,j,k) .eq. 0.) then
          
            Tnl(i,j,k) = ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty) 
         
         else
         
            Tnv(i,j,k) = ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty)
         
         end if

         !end if
                           
      end do
    end do

end subroutine Heat_calGradT_central
