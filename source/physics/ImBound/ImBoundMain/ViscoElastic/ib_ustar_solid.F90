!===============================================================================
!!
!! subroutine ib_ustar_solid
!!
!! corrects precursor velocity ustar by adding solid stress effect
!! take ustar, add solid stress term and return to a new ustar
!===============================================================================      

                     
      subroutine ib_ustar_solid(ustr, vstr, xms, Tau1,Tau2,Tau3,Tau4,& 
                                ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz,dt)

        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: ustr, vstr
        real, dimension(:,:,:), intent(in)    :: xms,Tau1,Tau2,Tau3,Tau4

        real, intent(in)    :: dx,dy,dz,dt
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        integer :: i,j,k
        !real, dimension(ix2-ix1+1,jy2-jy1+1,kz2-kz1+1)  :: ustrB, vstrB
        real :: ustrB, vstrB
        real :: xmsccc, xmsrcc, xmslcc, xmscrc
        real :: xmsrrc, xmsrlc, xmslrc, xmsclc
        real :: re_s

        
        re_s = 10.d0 !set re_s
        !------------calculate ustarB and vstarB and add them to ustar vstar-------------

        ustrB = 0.d0
        vstrB = 0.d0

           k = 1
           do j = jy1, jy2
              do i = ix1, ix2
              
                 !--- elastic modulus                  
                 xmsccc = xms(i,  j,  k)
                 xmsrcc = xms(i+1,j,  k)
                 xmslcc = xms(i-1,j,  k)
                 xmscrc = xms(i,  j+1,k)
                 xmsclc = xms(i,  j-1,k)
                 xmsrrc = xms(i+1,j+1,k)
                 xmsrlc = xms(i+1,j-1,k)
                 xmslrc = xms(i-1,j+1,k)


                !---Update ustar by adding elastic solid force term to solid region 

                 !ustrB(i,j,1) =  &                                          !
                 ustrB =  & 
                                                                            !
                          (xmsrcc  * Tau1(i+1,j,k) -             &           !
                           xmsccc  * Tau1(i,j,k)) / 1.d0 / dx/ re_s   &      !  
                                                                            !
                        + (xmscrc*Tau2(i,j+1,k) +xmsrrc*Tau2(i+1,j+1,k) &     ! 
                        - xmsclc*Tau2(i,j-1,k) -xmsrlc*Tau2(i+1,j-1,k)) & 
                            /2.d0/ 2.d0/dy/re_s
                                                                        
                 !ustr(i,j,1) = ustr(i,j,1) + ustrB(i,j,1)
                 ustr(i,j,1) = ustr(i,j,1) + ustrB*dt

                                                                            !
                !---Update vstar by adding elastic solid force term to solid region 

                 !vstrB(i,j,1) =  &                                          !
                 vstrB =  &
                                                                            !
                         (xmsrcc*Tau3(i+1,j,k) + xmsrrc*Tau3(i+1,j+1,k)  & 
                         -xmslcc*Tau3(i-1,j,k)-xmslrc*Tau3(i-1,j+1,k))   &    !  
                          / 2.d0 / 2.d0 / dx/ re_s                     &
 
                        +(xmscrc  * Tau4(i,j+1,k) -                     &    ! 
                              xmsccc  * Tau4(i,j,k)) / 1.d0 / dy/ re_s    
                                                                        

                 !vstr(i,j,1) = vstr(i,j,1) + vstrB(i,j,1)
                 vstr(i,j,1) = vstr(i,j,1) + vstrB*dt



              end do
           end do

      end subroutine ib_ustar_solid
