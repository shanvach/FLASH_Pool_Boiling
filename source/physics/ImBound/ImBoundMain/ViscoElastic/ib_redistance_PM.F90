!===============================================================================
!!
!! subroutine distance function; projection method 
!!
!=============================================================================== 

      subroutine ib_redistance_PM(s,rho1, rho2, xmu1, xmu2, xmus,&
                                     rho, xmu, xms, blockID,     &
                                     ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
#include "constants.h"
#include "Flash.h"

        use Grid_interface, only: Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

          implicit none

          real, dimension(:,:,:), intent(inout) :: s, rho, xmu, xms
          real, intent(in)    :: rho1, rho2, xmu1, xmu2, xmus
          real, intent(in)    :: dx,dy,dz
          integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
          integer, intent(in) :: blockID

          real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: phi, s1, dns

          integer :: nx, ny, nz
          integer :: i,j,k
          integer :: nip,nip_max,inv,n1, n3, nn1, nn2, n, m, n2 
          real    :: ul,ur,vl,vr
          real    :: th, xx, yy, xx1, yy1, xx2, yy2
          real    :: eps, psi
          real    :: r1, r2, r3
          real    :: xx01, yy01, xx02, yy02, xx03, yy03
          real    :: xcell1,ycell1,zcell1,xcell2,ycell2,zcell2

          !nx = ix2-ix1+1
          !ny = jy2-jy1+1
          !nz = kz2-kz1+1
          !nip_max = max(nx,ny) * min(nx,ny)
          !real, dimension(ix2-ix1+1,jy2-jy1+1,kz2-kz1+1)  :: phi, s1, dns
          real, dimension((ix2-ix1+1)*(jy2-jy1+1),6)  :: xip
          !real, dimension(nip_max,6)  :: xip, xip1
          !end header
          real, dimension(2,MDIM) :: boundBox
          real bsize(MDIM),coord(MDIM)
          real del(MDIM)

        !implicit none
        !include 'mpif.h'
        !integer i, j, k, npid, nop, nopx, nopy, nx, ny, nz, n
        !integer itr, itr_max
        !integer nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !integer nip, nip1, nip_max, inv, n1, n3, nn1, nn2, m, n2
        !integer ierr
        !integer cpt
        !integer iont, iont2, iwnt
        !integer nx11, ny11, nx21, ny21
        !integer send_request, recv_request
        !double precision s, phi, x, y, z, rho, xmu, xms
        !double precision xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz
        !double precision xip, xip1
        !double precision th, xx, yy, xx1, yy1, xx2, yy2
        !double precision r1, r2, r3
        !double precision xx01, yy01, xx02, yy02, xx03, yy03
        !double precision s1, dns
        !double precision eps, psi
        !double precision send_buff, recv_buff
        !double precision rho1, rho2, rho3, xmu1, xmu2, xmu3, xmus, re, re_s, g, st, stv
        !double precision umx, vmx, prms, res, div, tp, tpi, vf, vfi
        !double precision iwt
        !double precision t, dt, ti, tf, cfl
        !dimension :: s(nx,ny,nz), phi(nx,ny,nz)
        !dimension :: rho(nx,ny,nz), xmu(nx,ny,nz) ,xms(nx,ny,nz)
        !dimension :: cpt(2)
        !dimension :: x(nx), y(ny), z(nz)
        !dimension :: iwt(1000)
        !allocatable :: s1(:,:,:), dns(:,:,:), xip(:,:), xip1(:,:)
        !allocatable :: send_buff(:), recv_buff(:)
        !common/grid/xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        !common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !common/param/     rho1, rho2, rho3, xmu1, xmu2, xmu3, xmus, re, re_s, g, st, stv
        !common/mpi/npid, nop, nopx, nopy
        !common/max/umx, vmx, prms, res, div, tp, tpi, vf, vfi
        !common/intw/iont, iont2, iwnt
        !common/time/t, dt, ti, tf, cfl

!        call Grid_getBlkBoundBox(blockID,boundBox)
!        bsize(:) = boundBox(2,:) - boundBox(1,:)
!
!        call Grid_getBlkCenterCoords(blockID,coord)
!
!        call Grid_getDeltas(blockID,del)
!
!        eps = 1.E-10
!        !1 is material property inside the interface, 2 is material property outside the interface
!        !rho1 = 1.d0
!        !rho2 = 1.d0
!        !xmu1 = 1.d0
!        !xmu2 = 1.d0 
!        !xmus = 1.d0
!
!        !nip_max = max(nx,ny) * min(nx,ny)
!        !allocate(s1(nx,ny,nz), dns(nx,ny,nz))
!        !allocate(xip(nip_max,6), xip1(nip_max,6))
!
!        !nx11 = nx1
!        !ny11 = ny1
!        !nx21 = nx2
!        !ny21 = ny2
!
!        !if(cpt(1).eq.1) then
!        !   nx11 = 2
!        !end if
!        !if(cpt(1).eq.nopx) then
!        !   nx21 = nx-1
!        !end if                 
!        !if(cpt(2).eq.1) then
!        !   ny11 = 2
!        !end if
!        !-------------------------------------------------------#
!        !find interface location                                #
!        !-------------------------------------------------------# 
!        nip = 0                                                 !
!        xip = 0.d0                                              !
!                                                                !
!        do j = jy1, jy2                                         !
!           do i = ix1, ix2                                      !
!              !calculate x(i),x(i+1),y(j),y(j+1)                !                               
!              xcell1 = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!                    real(i - NGUARD - 1)*del(IAXIS) +   &
!                    0.5*del(IAXIS)
! 
!              ycell1  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!                    real(j - NGUARD - 1)*del(JAXIS)  +  &
!                    0.5*del(JAXIS)
!
!              xcell2 = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!                    real(i+1 - NGUARD - 1)*del(IAXIS) +   &
!                    0.5*del(IAXIS)
! 
!              ycell2  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!                    real(j+1 - NGUARD - 1)*del(JAXIS)  +  &
!                    0.5*del(JAXIS)
!
!              if(s(i,j,1).ge.0..and.s(i+1,j,1).le.0.) then      !
!                                                                !
!                 th = s(i+1,j,1)/(s(i+1,j,1)-s(i,j,1))          !
!                 xx = xcell2 - th*(xcell2-xcell1)                 !
!                 nip = nip + 1                                  !
!                 xip(nip,1) = xx   !interfacial points          !   
!                 xip(nip,2) = ycell1                              !
!                                                                !
!              end if                                            !
!                                                                !
!              if(s(i,j,1).le.0..and.s(i+1,j,1).ge.0.) then      !
!                                                                !
!                 th = s(i,j,1)/(s(i,j,1)-s(i+1,j,1))            ! 
!                 nip = nip + 1                                  !
!                 xx = xcell1 + th*(xcell2-xcell1)                   ! 
!                 xip(nip,1) = xx   !interfacial points          !
!                 xip(nip,2) = ycell1                              !
!                                                                !   
!              end if                                            !
!                                                                !
!              if(s(i,j,1).ge.0..and.s(i,j+1,1).le.0.) then      !
!                                                                ! 
!                 th = s(i,j+1,1)/(s(i,j+1,1)-s(i,j,1))          !
!                 nip = nip + 1                                  !
!                 yy = ycell2 - th*(ycell2-ycell1)                 ! 
!                 xip(nip,1) = xcell1   !interfacial points        !
!                 xip(nip,2) = yy                                !
!                                                                !
!              end if                                            !
!                                                                !  
!              if(s(i,j,1).le.0..and.s(i,j+1,1).ge.0.) then      !
!                                                                !
!                 th = s(i,j,1)/(s(i,j,1)-s(i,j+1,1))            !
!                 nip = nip + 1                                  !
!                 yy = ycell1 + th*(ycell2-ycell1)                   ! 
!                 xip(nip,1) = xcell1   !interfacial points        !
!                 xip(nip,2) = yy                                !
!                                                                !
!              end if                                            !
!                                                                !
!           end do                                               !
!        end do                                                  ! 
!        !-------------------------------------------------------#
!
!
!        !!write out interface -----------------------------------#
!        !if(iwnt.eq.1.and.iwt(iont2).le.t) then                   !
!        !   iont2 = iont2 + 1                                      !
!        !   !sort                                                !
!        !   xip1(1,1) = xip(1,1)
!        !   xip1(1,2) = xip(1,2)
!        !   xip(1,3)  = 1.d0 
!        !   do i = 1, nip-1
!        !      xx1 = xip1(i,1)
!        !      yy1 = xip1(i,2)
!        !      xx = 1E9
!        !      do j = 1, nip
!        !         xx2 = xip(j,1)
!        !         yy2 = xip(j,2)
!        !         yy = sqrt((xx1-xx2)**2 + (yy1-yy2)**2)
!        !         if(xip(j,3).eq.0.d0.and.yy.lt.xx) then
!        !            xx = yy
!        !            k = j
!        !         end if
!        !      end do
!        !      xip1(i+1,1) = xip(k,1)
!        !      xip1(i+1,2) = xip(k,2)
!        !      xip(k,3)  = 1.d0
!        !   end do
!        !   !write                      
!        !   do i = 0,nop-1                                       !
!        !      if(npid.eq.i) then                                !
!        !         open(42,file='interface2.dat',position='append')!
!        !         if(npid.eq.0) then                             !
!        !            write(42,*) "ZONE T = ",'"zone"'            !
!        !         end if                                         !
!        !         do k = 1, nip                                  !
!        !            write(42,'(2f10.5)') xip1(k,1),xip1(k,2)    !
!        !         end do                                         !
!        !         close(42)                                      !
!
!        !          open(45,file='isolines.dat',position='append')!
!        !         !  if(npid.eq.0) then                             !
!        !         !     write(41,*) "ZONE T = ",'"zone"'            !
!        !         !  end if                                         !
!        !           do k = 1, nip                                  !
!        !              write(45,'(2f10.5)') xip1(k,1),xip1(k,2)    !
!        !           end do                                         !
!        !           close(45) 
!
!        !      end if                                            !         
!        !      call MPI_Barrier(MPI_COMM_WORLD,ierr)             !
!        !   end do                                               !
!        !end if                                                  !
!        !!-------------------------------------------------------#
!        
!       
!        !!-------------------------------------------------------!
!        !!share interface location                               !
!        !!-------------------------------------------------------!
!        !allocate(send_buff(nip_max))                            !
!        !allocate(recv_buff(nip_max))                            !
!                                                                !
!        !send_buff = 0.d0                                        !
!        !recv_buff = 0.d0                                        !
!        !                                                        !
!        !nip1 = 0                                                !
!        !xip1 = 0                                                !        
!                                                                ! 
!        !do i = 0, nop-1                                         !   
!                                                                ! 
!        !   if(npid.eq.i) then                                   !
!        !      send_buff(1) = dble(nip)                          !
!        !      j = 1                                             !
!        !      do k = 1, nip                                     !
!        !         send_buff(j+1) = xip(k,1)                      !
!        !         send_buff(j+2) = xip(k,2)                      !
!        !         j = j + 2                                      !
!        !      end do                                            !
!        !   end if                                               !
!                                                                !
!        !   call MPI_Bcast( send_buff, nip_max, MPI_REAL8, i, &  ! 
!        !        MPI_COMM_WORLD, ierr )                          !
!                                                                !
!        !   if(npid.ne.i) then                                   !
!        !      j = 1                                             !
!        !      do k = 1, int(send_buff(1))                       !
!        !         xip1(nip1+k,1) = send_buff(j+1)                !
!        !         xip1(nip1+k,2) = send_buff(j+2)                !
!        !         j = j + 2                                      !
!        !      end do                                            !
!        !      nip1 = nip1 + int(send_buff(1))                   !
!        !   end if                                               ! 
!                                                                !
!        !end do                                                  !
!                                                                !
!        !do i = 1, nip1                                          !
!        !   xip(nip+i,1) = xip1(i,1)                             !
!        !   xip(nip+i,2) = xip1(i,2)                             !  
!        !end do                                                  ! 
!                                                                !
!        !nip = nip + nip1                                        !
!        !!-------------------------------------------------------!
!        
!        
!        !-----------------------------------------!
!        !mark nodes for selective redistance------!
!        dns = 1.d0                                !
!        do j = jy1, jy2                           !
!           do i = ix1, ix2                        !    
!              if(abs(s(i,j,1)).le.1.5d0*dx) then  !
!                 dns(i,j,1) = 0.d0                !
!              end if                              ! 
!           end do                                 !
!        end do                                    !
!        !-----------------------------------------!  
!        
!
!        !------------distance function-----------------------------!
!        !-----------level 1 (P1) projection------------------------!
!        inv = nip                                                  ! 
!        phi = 0.d0                                                 !
!        if(nip.gt.0) then                                          !
!           n1 = 1                                                  !
!           do j = jy1, jy2                                         !
!              do i = ix1, ix2                                      !
!                 r1 = 1E9                                          !
!                 n = 1                                             !
!
!              !calculate x(i),x(i+1),y(j),y(j+1)                   !                               
!              xcell1 = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!                    real(i - NGUARD - 1)*del(IAXIS) +   &
!                    0.5*del(IAXIS)
! 
!              ycell1  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!                    real(j - NGUARD - 1)*del(JAXIS)  +  &
!                    0.5*del(JAXIS)
!
!              xcell2 = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!                    real(i+1 - NGUARD - 1)*del(IAXIS) +   &
!                    0.5*del(IAXIS)
! 
!              ycell2  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!                    real(j+1 - NGUARD - 1)*del(JAXIS)  +  &
!                    0.5*del(JAXIS)
!
!                 do k = 1,inv                                      !
!                    r2 = (xcell1-xip(k,1))**2 + (ycell1-xip(k,2))**2   !          
!                    if(r2.lt.r1) then                              !
!                       n1 = k                                      !
!                       r1 = r2                                     !
!                    end if                                         !
!                 end do                                            !
!                 r1 = 1E9                                          !       
!                 n2 = 1                                            !
!                 do k = 1,inv                                      !
!                    r2 = (xcell1-xip(k,1))**2 + (ycell1-xip(k,2))**2   !          
!                    if(r2.lt.r1.and.k.ne.n1) then                  !
!                       n2 = k                                      !
!                       r1 = r2                                     !
!                    end if                                         !
!                 end do                                            !
!                 r1 = 1E9                                          !
!                 n3 = n1                                           !
!                 do k = 1,inv                                      !
!                    r2 = (xcell1-xip(k,1))**2 + (ycell1-xip(k,2))**2   !          
!                    if(r2.lt.r1.and.k.ne.n1.and.k.ne.n2) then      !
!                       n3 = k                                      !
!                       r1 = r2                                     !
!                    end if                                         !
!                 end do                                            !
!                                                                   !                
!                 !default values                                   !
!                 xx01 = xip(n1,1)                                  !
!                 yy01 = xip(n1,2)                                  !
!                                                                   ! 
!                 !update level set to sdf away from interface------!
!                 s1(i,j,1) = sqrt( (xcell1 - xx01)**2   +   &        !
!                      (ycell1 - yy01)**2 ) *   &                     !
!                      s(i,j,1)/max(eps,abs(s(i,j,1))) * &          !     
!                      dns(i,j,1)           +   &                   !
!                      s(i,j,1) * (1.d0-dns(i,j,1))                 !
!                                                                   !
!                 !-------------------------------------------------!               
!                 !--------level 2 (P2) projection------------------!
!                 !-------------------------------------------------!
!                 !check if 1st pt in middle                        !
!                 xx1 = sqrt((xip(n1,1)-xip(n3,1))**2 + &           !
!                      (xip(n1,2)-xip(n3,2))**2)                    !
!                 xx2 = sqrt((xip(n2,1)-xip(n3,1))**2 + &           !
!                      (xip(n2,2)-xip(n3,2))**2)                    !
!                 if(xx1.gt.xx2) then                               !
!                    nn1 = n1                                       !
!                    nn2 = n2                                       !
!                    n1 = nn2                                       !
!                    n2 = nn1                                       !
!                 end if                                            !                
!                 xx02 = xx01                                       !
!                 yy02 = yy01                                       !
!                 r1 = sqrt((xcell1-xip(n1,1))**2+(ycell1-xip(n1,2))**2)!
!                 m = 11 ! note this parameter is adjustable        ! 
!                 do k = 1,m                                        ! 
!                    xx1 = xip(n1,1)                                ! 
!                    yy1 = xip(n1,2)                                !
!                    xx2 = xip(n2,1)                                !
!                    yy2 = xip(n2,2)                                !
!                    xx = (xx2-xx1)*dble(k-1)/dble(m-1) + xx1       !
!                    yy = (xx-xx2)/(xx1-xx2)*yy1 + &                !
!                         (xx-xx1)/(xx2-xx1)*yy2                    !
!                                                                   !
!                    r2 = sqrt((xcell1-xx)**2+(ycell1-yy)**2)           !
!                    if(r2.lt.r1) then                              !
!                       xx02 = xx                                   ! 
!                       yy02 = yy                                   !
!                       r1 = r2                                     !
!                    end if                                         !
!                 end do                                            !
!                                                                   !
!                 xx03 = xx01                                       !
!                 yy03 = yy01                                       !
!                 r1 = sqrt((xcell1-xip(n1,1))**2+(ycell1-xip(n1,2))**2)!
!                 m = 11 ! note this parameter is adjustable        !
!                 do k = 1,m                                        ! 
!                    xx1 = xip(n1,1)                                !
!                    yy1 = xip(n1,2)                                !
!                    xx2 = xip(n3,1)                                !
!                    yy2 = xip(n3,2)                                !
!                    xx = (xx2-xx1)*dble(k-1)/dble(m-1) + xx1       !
!                    yy = (xx-xx2)/(xx1-xx2)*yy1 + &                !
!                         (xx-xx1)/(xx2-xx1)*yy2                    ! 
!                                                                   !
!                    r2 = sqrt((xcell1-xx)**2+(ycell1-yy)**2)           !
!                    if(r2.lt.r1) then                              !
!                       xx03 = xx                                   !
!                       yy03 = yy                                   !
!                       r1 = r2                                     !
!                    end if                                         !
!                 end do                                            !
!                                                                   !
!                 r1 = sqrt((xcell1-xx01)**2+(ycell1-yy01)**2)      !     
!                 r2 = sqrt((xcell1-xx02)**2+(ycell1-yy02)**2)      !
!                 r3 = sqrt((xcell1-xx03)**2+(ycell1-yy03)**2)      !
!                 if(r2.lt.r1.and.r2.lt.r3) then                    ! 
!                    xx01 = xx02                                    ! 
!                    yy01 = yy02                                    ! 
!                 elseif(r3.lt.r1.and.r3.lt.r2) then                !
!                    xx01 = xx03                                    !
!                    yy01 = yy03                                    ! 
!                 end if                                            !
!                                                                   !
!                 phi(i,j,1) = sqrt((xcell1-xx01)**2   + &          !
!                      (ycell1-yy01)**2)  * &                       !
!                      s(i,j,1)/max(eps,abs(s(i,j,1)))              !
!                                                                   !
!                 !----update density and viscosity-----------------!
!                 psi = (1.d0 + derf(-phi(i,j,1)/(2.d0*dx)))/2.d0   !
!                 xmu(i,j,1) = psi*(xmu1-xmu2) + xmu2               !
!                 !set xmu inside solid to be xmu1, outside to be xmu2
!                 !set xms inside solid to be xmus, outside to be 0
!                 xms(i,j,1) = psi*(xmus-0.d0) + 0.d0
!                 rho(i,j,1) = psi*(rho1-rho2) + rho2               !
!                 !-------------------------------------------------!
!                                                                   !
!              end do                                               !
!           end do                                                  !
!        end if                                                     !
!        !----------------------------------------------------------!
!
!        !convert level set to distance function with P1
!        s = s1
!        
!        !!----volume------------------------------------------------!
!        !vf = 0.d0                                                  !
!        !do j = ny1, ny2                                            !
!        !   do i = nx1, nx2                                         !
!        !      if(phi(i,j,1).le.0.d0) then                          !
!        !         vf = vf + 1                                       ! 
!        !      end if                                               !
!        !   end do                                                  !
!        !end do                                                     !
!        !!----------------------------------------------------------!
         !     phi = s
          do j = jy1, jy2
             do i = ix1, ix2
                  !----update density and viscosity-----------------!
                  psi = (1.d0 + derf(-s(i,j,1)/(2.d0*dx)))/2.d0   !
                  xmu(i,j,1) = psi*(xmu1-xmu2) + xmu2               !
                  !set xmu inside solid to be 0, outside to be xmu2 !
                  !xmu(i,j,1) = psi*(0.d0-xmu2) + xmu2              !
                  !set xms inside solid to be xmus, outside to be 0 !
                  xms(i,j,1) = psi*(xmus-0.d0) + 0.d0               !
                  rho(i,j,1) = psi*(rho1-rho2) + rho2               !
             end do
          end do        

      end subroutine ib_redistance_PM
