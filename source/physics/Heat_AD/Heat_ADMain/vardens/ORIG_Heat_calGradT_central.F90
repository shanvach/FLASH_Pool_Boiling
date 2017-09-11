subroutine Heat_calGradT_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny,mflg)

    use Heat_AD_data
    use Heat_AD_interface, only: Heat_GFMstencil_o1, Heat_GFMstencil_o2

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2

    real :: th,tol
    integer :: i,j,k
    real :: Tij,Tx_plus,Tx_mins,Ty_plus,Ty_mins,Tx,Ty,Tz

    real :: thxp1,thxm1,thyp1,thym1
    real :: thxp2,thxm2,thyp2,thym2

    tol = 0.01

    k = 1

    do i=ix1,ix2
      do j=jy1,jy2

        Tx_plus = T(i+1,j,k)
        Tx_mins = T(i-1,j,k)
        Ty_plus = T(i,j+1,k)
        Ty_mins = T(i,j-1,k)
        Tij     = T(i,j,k)

        thxp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
        thxm1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))
        thyp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))
        thym1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

        thxp2 = abs(s(i-1,j,k))/(abs(s(i-1,j,k))+abs(s(i+1,j,k)))
        thxm2 = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i-1,j,k)))
        thyp2 = abs(s(i,j-1,k))/(abs(s(i,j-1,k))+abs(s(i,j+1,k)))
        thym2 = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j-1,k)))

        ! Case 1 !
        if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then

           if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then

                call Heat_GFMstencil_o1(Tx_plus,Tij,ht_Tsat,max(tol,thxp1))

           else

                if(thxp1 .gt. tol) then
                call Heat_GFMstencil_o2(Tx_plus,Tij,T(i-1,j,k),ht_Tsat,thxp1)

                else 

                        if(s(i,j,k)*s(i-2,j,k) .ge. 0.d0) then
                        call Heat_GFMstencil_o2(Tx_plus,T(i-1,j,k),T(i-2,j,k),ht_Tsat,thxp2)

                        else
                        call Heat_GFMstencil_o2(Tx_plus,Tij,T(i-1,j,k),ht_Tsat,tol)

                        end if
                end if

           end if

        end if
        ! End of Case 1 !


        ! Case 2 !
        if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then

            if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then

                call Heat_GFMstencil_o1(Tx_mins,Tij,ht_Tsat,max(tol,thxm1))

           else

                if(thxm1 .gt. tol) then
                call Heat_GFMstencil_o2(Tx_mins,Tij,T(i+1,j,k),ht_Tsat,thxm1)

                else 

                        if(s(i,j,k)*s(i+2,j,k) .ge. 0.d0) then
                        call Heat_GFMstencil_o2(Tx_mins,T(i+1,j,k),T(i+2,j,k),ht_Tsat,thxm2)

                        else
                        call Heat_GFMstencil_o2(Tx_mins,Tij,T(i+1,j,k),ht_Tsat,tol)

                        end if
                end if

           end if

        end if
        ! End of Case 2 !

        ! Case 3 !
        if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then

           if(s(i,j,k)*s(i,j-1,k) .le. 0.d0) then

                call Heat_GFMstencil_o1(Ty_plus,Tij,ht_Tsat,max(tol,thyp1))

           else

                if(thyp1 .gt. tol) then
                call Heat_GFMstencil_o2(Ty_plus,Tij,T(i,j-1,k),ht_Tsat,thyp1)

                else 

                        if(s(i,j,k)*s(i,j-2,k) .ge. 0.d0) then
                        call Heat_GFMstencil_o2(Ty_plus,T(i,j-1,k),T(i,j-2,k),ht_Tsat,thyp2)

                        else
                        call Heat_GFMstencil_o2(Ty_plus,Tij,T(i,j-1,k),ht_Tsat,tol)

                        end if

                end if

           end if

        end if
        ! End of Case 3 !

        ! Case 4 !
        if(s(i,j,k)*s(i,j-1,k) .le. 0.d0) then

            if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then

                call Heat_GFMstencil_o1(Ty_mins,Tij,ht_Tsat,max(tol,thym1))

           else

                if(thym1 .gt. tol) then
                call Heat_GFMstencil_o2(Ty_mins,Tij,T(i,j+1,k),ht_Tsat,thym1)

                else 

                        if(s(i,j,k)*s(i,j+2,k) .ge. 0.d0) then
                        call Heat_GFMstencil_o2(Ty_mins,T(i,j+1,k),T(i,j+2,k),ht_Tsat,thym2)

                        else
                        call Heat_GFMstencil_o2(Ty_mins,Tij,T(i,j+1,k),ht_Tsat,tol)

                        end if

                end if

           end if

         end if
         ! End of Case 4 ! 

         Tx = (Tx_plus - Tx_mins)/(2*dx)
         Ty = (Ty_plus - Ty_mins)/(2*dy)

         !--------------------------------------------------------------------------!
         !----------------------------Calculate Fluxes------------------------------!
         !--------------------------------------------------------------------------!

         if (pf(i,j,k) .eq. 0.) then
          
            Tnl(i,j,k) = ( + nx(i,j,k)*Tx + ny(i,j,k)*Ty) 
         
         else
         
            Tnv(i,j,k) = ( - nx(i,j,k)*Tx - ny(i,j,k)*Ty)
         
         end if
   
      end do
    end do

end subroutine Heat_calGradT_central
