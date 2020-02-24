!===============================================================================
!!
!! subroutine ib_ustar_solid
!!
!! corrects precursor velocity ustar by adding solid stress effect
!! take ustar, add solid stress term and return to a new ustar
!===============================================================================      

                     
      subroutine ib_ustar_solid(ustr, vstr, xms, Tau,& 
                                ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)

        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: ustr, vstr, Tau
        real, dimension(:,:,:), intent(in)    :: xms

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        integer :: i,j,k
        real :: ustrB, vstrB
        real :: xmsccc, xmsrcc, xmslcc, xmscrc
        real :: xmsrrc, xmsrlc, xmslrc, xmsclc
        real :: re_s

        
        re_s = 10.d0 !set re_s
        !------------calculate ustarB and vstarB and add them to ustar vstar-------------

        ustrB = 0.d0
        vstrB = 0.d0

        do k = 1
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

                 ustrB(i,j,1) =  &                                          !
                                                                            !
                          (xmsrcc  * Tau(i+1,j,1) -             &           !
                           xmsccc  * Tau(i,j,1)) / 1.d0 / dx/ re_s   &      !  
                                                                            !
                        + (xmscrc*Tau(i,j+1,2) +xmsrrc*Tau(i+1,j+1,2) &     ! 
                        - xmsclc*Tau(i,j-1,2) -xmsrlc*Tau(i+1,j-1,2)) & 
                            /2.d0/ 2.d0/dy/re_s
                                                                        
                 ustr(i,j,1) = ustr(i,j,1) + ustrB(i,j,1)

                                                                            !
                !---Update vstar by adding elastic solid force term to solid region 

                 vstrB(i,j,1) =  &                                          !
                                                                            !
                         (xmsrcc*Tau(i+1,j,3) + xmsrrc*Tau(i+1,j+1,3)  & 
                         -xmslcc*Tau(i-1,j,3)-xmslrc*Tau(i-1,j+1,3))   &    !  
                          / 2.d0 / 2.d0 / dx/ re_s                     &
 
                        +(xmscrc  * Tau(i,j+1,4) -                     &    ! 
                              xmsccc  * Tau(i,j,4)) / 1.d0 / dy/ re_s    
                                                                        

                 vstr(i,j,1) = vstr(i,j,1) + vstrB(i,j,1)



              end do
           end do
        end do

      end subroutine ib_ustar_solid
