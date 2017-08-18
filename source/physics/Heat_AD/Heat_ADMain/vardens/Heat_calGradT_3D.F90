subroutine Heat_calGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,nx,ny,nz,mflg)
        
    use Heat_AD_data

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

    real :: th,tol,th2
    integer :: i,j,k
    real :: Tij,Tipj,Timj,Tijp,Tijm,dxp,dxm,dyp,dym,Tx,Ty,Tax,Tbx,Tay,Tby,dyc,dxc,Tpx,Tmx,Tpy,Tmy
    real :: Tmz,Tpz,dzp,dzm,Tkm,Tkp,Tz
    logical :: int_xp,int_xm,int_yp,int_ym,int_zp,int_zm


   tol = 0.01

   do k=kz1,kz2
    do j=jy1,jy2
      do i=ix1,ix2

        !--------------------------------------------------------------------------!
        !-----------------------------X - Direction--------------------------------!
        !--------------------------------------------------------------------------!

        if((s(i+1,j,k)*s(i,j,k) .le. 0.0) .and. (s(i-1,j,k)*s(i,j,k) .le. 0.0)) then

            th  =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))
            th2 =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))))
            Tx  = 0.5*((T(i,j,k) - ht_Tsat)/(th*dx) + (ht_Tsat-T(i,j,k))/(th2*dx))

        else if((abs(s(i-1,j,k)) .le. abs(s(i+1,j,k))) .or. (s(i-1,j,k)*s(i,j,k) .le. 0.0)) then


            if(s(i-1,j,k)*s(i,j,k) .le. 0.0) then

                if((abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))) .gt. tol) then

                    th = (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))
                    Tx = (T(i,j,k)-ht_Tsat)/(th*dx)

                else

                    th = (abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i-1,j,k))))
                    Tx = (T(i+1,j,k)-ht_Tsat)/(th*dx)

                end if
            else

                Tx = (T(i,j,k)-T(i-1,j,k))/dx

            end if

        else if((abs(s(i-1,j,k)) .gt. abs(s(i+1,j,k))) .or. (s(i+1,j,k)*s(i,j,k) .le. 0.0)) then

            if(s(i+1,j,k)*s(i,j,k) .le. 0.0) then

                if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))) .gt. tol) then

                  th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
                  Tx = (ht_Tsat-T(i,j,k))/(th*dx)

                else

                  th = abs(s(i-1,j,k))/(abs(s(i-1,j,k))+abs(s(i+1,j,k)))
                  Tx = (ht_Tsat-T(i-1,j,k))/(th*dx)

                end if
            else

                Tx = (T(i+1,j,k)-T(i,j,k))/dx


            end if

        endif

        !--------------------------------------------------------------------------!
        !-----------------------------Y - Direction--------------------------------!
        !--------------------------------------------------------------------------!

        if((s(i,j+1,k)*s(i,j,k) .le. 0.0) .and. (s(i,j-1,k)*s(i,j,k) .le. 0.0)) then

            th  =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))))
            th2 =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))
            Ty  = 0.5*((T(i,j,k) - ht_Tsat)/(th*dy) + (ht_Tsat-T(i,j,k))/(th2*dy))

        else if((abs(s(i,j-1,k)) .le. abs(s(i,j+1,k))) .or. (s(i,j-1,k)*s(i,j,k) .le. 0.0)) then

            if(s(i,j-1,k)*s(i,j,k) .le. 0.0) then

                if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))) .gt. tol) then

                  th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))
                  Ty = (T(i,j,k)-ht_Tsat)/(th*dy)

                else

                  th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j-1,k)))
                  Ty = (T(i,j+1,k)-ht_Tsat)/(th*dy)

                end if
            else

                Ty = (T(i,j,k)-T(i,j-1,k))/dy


            end if

        else if((abs(s(i,j-1,k)) .gt. abs(s(i,j+1,k))) .or. (s(i,j+1,k)*s(i,j,k) .le. 0.0)) then

            if(s(i,j+1,k)*s(i,j,k) .le. 0.0) then


                if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) .gt. tol) then

                  th = (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))
                  Ty = (ht_Tsat-T(i,j,k))/(th*dy)

                else

                  th = (abs(s(i,j-1,k))/(abs(s(i,j-1,k))+abs(s(i,j+1,k))))
                  Ty = (ht_Tsat-T(i,j-1,k))/(th*dy)

                end if
            else

                Ty = (T(i,j+1,k)-T(i,j,k))/dy

            end if

        endif

        !--------------------------------------------------------------------------!
        !-----------------------------Z - Direction--------------------------------!
        !--------------------------------------------------------------------------!

        if((s(i,j,k+1)*s(i,j,k) .le. 0.0) .and. (s(i,j,k-1)*s(i,j,k) .le. 0.0)) then

            th  =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1))))
            th2 =  max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1))))
            Tz  = 0.5*((T(i,j,k) - ht_Tsat)/(th*dz) + (ht_Tsat-T(i,j,k))/(th2*dz))

        else if((abs(s(i,j,k-1)) .le. abs(s(i,j,k+1))) .or. (s(i,j,k-1)*s(i,j,k) .le. 0.0)) then

            if(s(i,j,k-1)*s(i,j,k) .le. 0.0) then

                if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1))) .gt. tol) then

                  th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1)))
                  Tz = (T(i,j,k)-ht_Tsat)/(th*dz)

                else

                  th = abs(s(i,j,k+1))/(abs(s(i,j,k+1))+abs(s(i,j,k-1)))
                  Tz = (T(i,j,k+1)-ht_Tsat)/(th*dz)

                end if
            else

                Tz = (T(i,j,k)-T(i,j,k-1))/dz


            end if

        else if((abs(s(i,j,k-1)) .gt. abs(s(i,j,k+1))) .or. (s(i,j,k+1)*s(i,j,k) .le. 0.0)) then

            if(s(i,j,k+1)*s(i,j,k) .le. 0.0) then

                if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1))) .gt. tol) then

                  th = (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1))))
                  Tz = (ht_Tsat-T(i,j,k))/(th*dz)

                else

                  th = (abs(s(i,j,k-1))/(abs(s(i,j,k-1))+abs(s(i,j,k+1))))
                  Tz = (ht_Tsat-T(i,j,k-1))/(th*dz)

                end if
            else

                Tz = (T(i,j,k+1)-T(i,j,k))/dz

            end if

        endif

        !--------------------------------------------------------------------------!
        !----------------------------Calculate Fluxes------------------------------!
        !--------------------------------------------------------------------------!

        if (pf(i,j,k) .eq. 0.) then
          
            Tnl(i,j,k) = ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty + nz(i,j,k)*Tz)
         
        else
         
            Tnv(i,j,k) = ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty - nz(i,j,k)*Tz)
         
        end if
                          
      end do
    end do
  end do

end subroutine Heat_calGradT_3D
