! Subroutine interpvars_face:
! Subroutine that interpolates values to a given given block
!
! Written by Marcos Vanella, December 2006
! --------------------------------------------------------------

!!#define RESCALE_LOCATION 1

subroutine interpvars_face(ng,nxb,nyb,nzb,coord,bsize,&
                           nxvar,nyvar,nzvar,xvar,yvar,zvar,&
                           Var0,nx1,ny1,nz1,blkvar,gridflag)


  use Grid_data, only : gr_imax,gr_imin,gr_jmax,gr_jmin,gr_kmax,gr_kmin

  implicit none
#include "constants.h"
      
  integer ng,nxb,nyb,nzb,nxvar,nyvar,nzvar, &
         nx1,ny1,nz1,gridflag

  real coord(3),bsize(3),Var0(nxvar,nyvar,nzvar),&
      xvar(nxvar),yvar(nyvar),zvar(nzvar),&
      blkvar(nx1,ny1,nz1)


  ! Local Variables:
  integer, parameter :: stencil    = 8
  integer, parameter :: interptype = 2  ! =1 MLS, =2 8 node fem interp. functions
  integer, parameter :: interp     = 1
  integer, parameter :: npol       = 4
  integer, parameter :: nderiv     = 4


  integer buildflag
  real dsx,dsy,dsz
  real A(npol,npol,nderiv), B(npol,stencil,nderiv)
  real p(npol,nderiv), phi(stencil,nderiv), &
      gamma(npol,nderiv), indx(npol)

  real xs(stencil),ys(stencil),zs(stencil)
  integer ielem(stencil),jelem(stencil),kelem(stencil)
  real dx,dy,dz,dxaux,dyaux,dzaux,xi,yj,zk
  integer i,j,k,ii,jj,kk,icell,jcell,kcell
  integer is

  real*8 zp,d
  real*8 sitan,etan,xin

  logical :: wrxflg 
  logical :: wryflg 
  logical :: wrzflg 

  logical :: inject = .true.

  real*8 distx,disty,distz

  real, parameter :: eps = 1.e-10

  !------------------------------------------------------------------------------ 
  wrxflg = .true.
  wryflg = .true.
  wrzflg = .true.
      
  ! stencil to use in the interpolation
  dx = bsize(1)/real(nxb)
  dy = bsize(2)/real(nyb)
  dz = bsize(3)/real(nzb)

!!$  write(*,*) 'dx,dy,dz=',dx,dy,dz

  select case (gridflag)
  case(1)
     dxaux = 0.0
     dyaux = 0.5*dy
     dzaux = 0.5*dz
  case(2)
     dxaux = 0.5*dx
     dyaux = 0.0
     dzaux = 0.5*dz
  case(3)
     dxaux = 0.5*dx
     dyaux = 0.5*dy
     dzaux = 0.0
  end select
      

  ! Initialize the variable to zero:
  blkvar = 0.


  do k = ng+1,nz1-ng
     do j = ng+1,ny1-ng
        do i = ng+1,nx1-ng

#ifndef RESCALE_LOCATION
           ! Position of block point i,j,k
           xi = coord(1) - 0.5*bsize(1) + &
                real(i - ng - 1)*dx + dxaux

           yj = coord(2) - 0.5*bsize(2) + &
                real(j - ng - 1)*dy + dyaux

           zk = coord(3) - 0.5*bsize(3) + &
                real(k - ng - 1)*dz + dzaux                  
#else
           ! Rescale:
           ! Position of block point i,j,k
           xi = (coord(1) - 0.5*bsize(1) + &
                real(i - ng - 1)*dx + dxaux)*2.*PI/(gr_imax-gr_imin)

           yj = (coord(2) - 0.5*bsize(2) + &
                real(j - ng - 1)*dy + dyaux)*2.*PI/(gr_jmax-gr_jmin)

           zk = (coord(3) - 0.5*bsize(3) + &
                real(k - ng - 1)*dz + dzaux - PI/2.)*2.*PI/(gr_kmax-gr_kmin)
#endif



!!$           if (gridflag .eq. 1 .and. &
!!$               k .eq. ng+2     .and. &
!!$               j .eq. ng+2     .and. &
!!$               i .eq. ng+1+2    ) then
!!$             
!!$           write(*,*) 'Guardcells',ng
!!$           
!!$           write(*,*) 'Uposition, 2,2,2=', &
!!$                       xi,yj,zk
!!$
!!$           endif

               
           if (inject) then

               ! Now Find the enclosing cell from the Var0 grid:
               icell = 0; jcell = 0; kcell = 0;
               distx = 1.e30; disty = 1.e30;  distz = 1.e30; 

               if ( (xi .gt. (PI+eps)) ) then
                  xi = xi - 2.*PI 
               elseif ( (xi .lt. -(PI+eps)) ) then
                  xi = xi + 2.*PI
               endif

               do ii = 1,nxvar-1
                  if( abs(xvar(ii)-xi) .le. distx) then
                    icell = ii
                    distx = abs(xvar(ii)-xi)
                  endif
               enddo

               if ( (yj .gt. (PI+eps)) ) then
                  yj = yj - 2.*PI
               elseif ( (yj .lt. -(PI+eps)) ) then
                  yj = yj + 2.*PI
               endif               

               do jj = 1,nyvar-1
                  if( abs(yvar(jj)-yj) .le. disty) then
                    jcell = jj
                    disty = abs(yvar(jj)-yj)
                  endif
               enddo

               do kk = 1,nzvar-1
                  if( abs(zvar(kk)-zk) .le. distz) then
                    kcell = kk
                    distz = abs(zvar(kk)-zk)
                  endif
               enddo


           else

               ! Now Find the enclosing cell from the Var0 grid:
               icell = 0; jcell = 0; kcell = 0;
               do ii = 1,nxvar-1
                  if( xvar(ii) .le. xi .and. xvar(ii+1) .ge. xi) then
                    icell = ii
                    goto 10
                  endif
               enddo
 10            continue

               do jj = 1,nyvar-1
                  if( yvar(jj) .le. yj .and. yvar(jj+1) .ge. yj) then
                    jcell = jj
                    goto 11
                  endif
               enddo
 11            continue              

               do kk = 1,nzvar-1
                  if( zvar(kk) .le. zk .and. zvar(kk+1) .ge. zk) then
                    kcell = kk
                    goto 12
                  endif
               enddo
 12            continue  

           endif

           if (icell .eq. 0 .or. jcell  .eq. 0 .or. kcell .eq. 0) then

               write(*,*) 'Error in interpolation of values:'
               write(*,*) 'Gridflag = ',gridflag
               write(*,*) 'Icell,Jcell,Kcell =',icell,jcell,kcell
               write(*,*) 'Block Point coordinates:'
               write(*,*) xi,yj,zk
               stop
                
           endif

           ielem = (/ icell+1,icell+1,icell,icell,&
                     icell+1,icell+1,icell,icell /)               

           jelem = (/ jcell,jcell+1,jcell+1,jcell,&
                     jcell,jcell+1,jcell+1,jcell /)

           kelem = (/ kcell,kcell,kcell,kcell,&
                     kcell+1,kcell+1,kcell+1,kcell+1/)

           do is=1,stencil
              xs(is) = xvar(ielem(is))
              ys(is) = yvar(jelem(is))
              zs(is) = zvar(kelem(is))
           enddo

!!$           xs = xvar(ielem)
!!$           ys = yvar(jelem)
!!$           zs = zvar(kelem)

               
           if (inject) then

              if((zk .lt. -(PI+eps)) .or. (zk .gt. (PI+eps))) cycle 

              zp = Var0(icell,jcell,kcell)

           else

              if (interptype .eq. 1) then

                 ! Having found the cell interpolate:
                 dsx = 1.2*(xvar(icell+1)-xvar(icell))
                 dsy = 1.2*(yvar(jcell+1)-yvar(jcell))
                 dsz = 1.2*(zvar(kcell+1)-zvar(kcell))


                 ! Build A and B matrices:
!                 call buildABLan(stencil,interp,npol,nderiv,dsx,&
!                dsy,dsz,xi,yj,zk,xs,ys,zs,A,B,buildflag, 0 );


                 p(1,1) = 1.; p(2,1) = xi; p(3,1) = yj; p(4,1) = zk; !  p


                 ! Solve for gamma:
                 gamma(:,1) = p(:,1)

!                 call ludcmp(A(:,:,1),npol,npol,indx,d)
!                 call lubksb(A(:,:,1),npol,npol,indx,gamma(:,1))

                 ! Obtain Shape functions:
                 phi(:,1) =0.
                 do ii = 1 , stencil
                  do jj = 1, npol
                     phi(ii,1) = phi(ii,1) + gamma(jj,1)*B(jj,ii,1)
                  enddo
                 enddo                            

                 ! Value of the function in xp and yp:
                 zp = 0.;
                 do ii = 1 , stencil      
                  zp = zp + phi(ii,1)*&
                      Var0(ielem(ii),jelem(ii),kelem(ii));   
                 enddo   


               else

                 ! Having found the cell interpolate:
                 dsx = (xvar(icell+1)-xvar(icell))
                 dsy = (yvar(jcell+1)-yvar(jcell))
                 dsz = (zvar(kcell+1)-zvar(kcell))


                 ! Natural coordinates (from -1 to 1) for rectangle:
                 sitan = -1. + (xi - xs(4))/dsx * 2.
                 etan  = -1. + (yj - ys(4))/dsy * 2.
                 xin   = -1. + (zk - zs(4))/dsz * 2.   

                 ! shape functions:
                 phi(1,1) = (1.+sitan)*(1.-etan)*(1-xin)/8. 
                 phi(2,1) = (1.+sitan)*(1.+etan)*(1-xin)/8.
                 phi(3,1) = (1.-sitan)*(1.+etan)*(1-xin)/8.
                 phi(4,1) = (1.-sitan)*(1.-etan)*(1-xin)/8.
                 phi(5,1) = (1.+sitan)*(1.-etan)*(1+xin)/8.
                 phi(6,1) = (1.+sitan)*(1.+etan)*(1+xin)/8.
                 phi(7,1) = (1.-sitan)*(1.+etan)*(1+xin)/8.
                 phi(8,1) = (1.-sitan)*(1.-etan)*(1+xin)/8.
              

                 ! value of the function in xi, yj, zk:         
                 zp = 0.;
                 do ii = 1 , stencil      
                  zp = zp + phi(ii,1)*&
                      Var0(ielem(ii),jelem(ii),kelem(ii));   
                 enddo    

              endif

           endif

           blkvar(i,j,k) = zp; 

        enddo
     enddo
  enddo

  return

end subroutine interpvars_face
