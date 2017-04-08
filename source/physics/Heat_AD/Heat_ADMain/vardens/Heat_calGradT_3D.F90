subroutine Heat_calGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,nx,ny,nz,mflg)
        
    use Heat_AD_data

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

    real :: th,tol
    integer :: i,j,k
    real :: Tij,Tipj,Timj,Tijp,Tijm,dxp,dxm,dyp,dym,Tx,Ty,Tax,Tbx,Tay,Tby,dyc,dxc,Tpx,Tmx,Tpy,Tmy
    real :: Tmz,Tpz,dzp,dzm,Tkm,Tkp,Tz
    logical :: int_xp,int_xm,int_yp,int_ym,int_zp,int_zm


   tol = 0.01

   do k=kz1,kz2
    do j=jy1,jy2
      do i=ix1,ix2

!__________METHOD 1_____________!

         Tipj = T(i+1,j,k)
         Timj = T(i-1,j,k)
         Tijp = T(i,j+1,k)
         Tijm = T(i,j-1,k)
         Tkm  = T(i,j,k-1)
         Tkp  = T(i,j,k+1)

         Tij  = T(i,j,k)
         Tmx  = T(i,j,k)
         Tpx  = T(i,j,k)
         Tmy  = T(i,j,k)
         Tpy  = T(i,j,k)
         Tmz  = T(i,j,k)
         Tpz  = T(i,j,k)

         dxp = dx
         dxm = dx
         dym = dy
         dyp = dy
         dzm = dz
         dzp = dz

         int_xp = .FALSE.
         int_xm = .FALSE.
         int_yp = .FALSE.
         int_ym = .FALSE.
         int_zp = .FALSE.
         int_zm = .FALSE.

         if(s(i,j,k)*s(i+1,j,k) .le. 0.) then

            Tpx = T(i,j,k)

            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
      
            if(th .lt. tol) then

               th = tol

            end if

            dxp = th*dx
            Tipj = ht_Tsat

            int_xp = .TRUE.

         end if

         if(s(i,j,k)*s(i-1,j,k) .le. 0.) then

            Tmx = T(i,j,k)
      
            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))

            if(th .lt. tol) then

               th = tol

            end if

            dxm = th*dx
            Timj = ht_Tsat

            int_xm = .TRUE.
         
         end if

         if(s(i,j,k)*s(i,j+1,k) .le. 0.) then

            Tpy = T(i,j,k)

            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))

            if(th .lt. tol) then
               
                th = tol

            end if

            dyp = th*dy
            Tijp = ht_Tsat

            int_yp = .TRUE.

         end if

         if(s(i,j,k)*s(i,j-1,k) .le. 0.) then

            Tmy = T(i,j,k)

            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

            if(th .lt. tol) then

                th = tol

            end if

            dym = th*dy
            Tijm = ht_Tsat

            int_ym = .TRUE.

         end if

        if(s(i,j,k)*s(i,j,k+1) .le. 0.) then

            Tpz = T(i,j,k)

            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1)))

            if(th .lt. tol) then
               
                th = tol

            end if

            dzp = th*dz
            Tkp = ht_Tsat

            int_zp = .TRUE.

         end if

         if(s(i,j,k)*s(i,j,k-1) .le. 0.) then

            Tmz = T(i,j,k)

            th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1)))

            if(th .lt. tol) then

                th = tol

            end if

            dzm = th*dz
            Tkm = ht_Tsat

            int_zm = .TRUE.

         end if

         Tx = ((-dxp**2)*(Timj) + (dxp**2)*(Tmx) - (dxm**2)*(Tpx) + (dxm**2)*(Tipj))/(dxp*dxm*(dxp+dxm))
         Ty = ((-dyp**2)*(Tijm) + (dyp**2)*(Tmy) - (dym**2)*(Tpy) + (dym**2)*(Tijp))/(dyp*dym*(dyp+dym))
         Tz = ((-dzp**2)*(Tkm)  + (dzp**2)*(Tmz) - (dzm**2)*(Tpz) + (dzm**2)*(Tkp))/(dzp*dzm*(dzp+dzm))


         !if(int_xm .or. int_xp .or. int_ym .or. int_yp .or. int_zm .or. int_zp) then

         if (pf(i,j,k) .eq. 0.) then
          
            Tnl(i,j,k) = mflg(i,j,k)* ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty + nz(i,j,k)*Tz)
         
         else
         
            Tnv(i,j,k) = mflg(i,j,k)* ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty - nz(i,j,k)*Tz)
         
         end if

         !end if
                           
      end do
    end do
  end do

end subroutine Heat_calGradT_3D
