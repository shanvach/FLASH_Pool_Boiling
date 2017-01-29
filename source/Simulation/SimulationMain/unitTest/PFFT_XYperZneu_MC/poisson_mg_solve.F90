!*******************************************************************************

!  Routine:     mg_solve()

!  Description: Coarse-grid solver for the multigrid module.  This version
!               simply relaxes "to convergence" using mg_relax().

!  Parameters:  level       Level to solve on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.  Receives
!                           the solution.


subroutine poisson_mg_solve (level, irhs, ilhs )

!===============================================================================
#include "Flash.h"
#include "constants.h"

use mg_common
!!$use runtime_parameters
use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use Simulation_data

use Grid_data, only : gr_nBlockX,gr_nBlockY,gr_nBlockZ

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getListOfBlocks,   &
                             Grid_getLocalNumBlks
use tree, only : lrefine,nodetype

implicit none

integer            :: level, irhs, ilhs

integer, save      :: mgrid_solve_max_iter
logical, save      :: first_call = .true.

integer, save :: NX,NY,NZ
real, save, allocatable, dimension(:,:,:,:) :: DP3

real, save :: dx,dy,dz

real, pointer, dimension(:,:,:,:) :: solndata

integer lnblocks2,lb
!===============================================================================

if (first_call) then
  !call get_parm_from_context("mgrid_solve_max_iter", mgrid_solve_max_iter)
  call RuntimeParameters_get('mgrid_solve_max_iter', mgrid_solve_max_iter)
  first_call = .false.

  NX = NXB + 2
  NY = NYB + 2
  NZ = NZB + 2

  allocate(DP3(NX,NY,NZ,2))

  dx = (sim_xMax - sim_xMin)/real(NXB*gr_nBlockX)
  dy = (sim_yMax - sim_yMin)/real(NYB*gr_nBlockY)
  dz = (sim_zMax - sim_zMin)/real(NZB*gr_nBlockZ) 

end if


! load Source to DP3
call Grid_getLocalNumBlks(lnblocks2)
do lb = 1,lnblocks2
   if (lrefine(lb) .eq. level) then

     call Grid_getBlkPtr(lb,solnData,CENTER)

     DP3(2:NX-1,2:NY-1,2:NZ-1,1) = solnData(irhs,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+1:NGUARD+NYB)

     call Grid_releaseBlkPtr(lb,solnData,CENTER)

     goto 111
   endif

enddo

111 continue

call PresDp3P(nx,ny,nz,dx,dy,dz,DP3(:,:,:,1),DP3(:,:,:,2))

!call poisson_mg_relax (level, irhs, ilhs, mgrid_solve_max_iter)

!===============================================================================

do lb = 1,lnblocks2
   if (lrefine(lb) .eq. level) then

     call Grid_getBlkPtr(lb,solnData,CENTER)

     solnData(ilhs,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+1:NGUARD+NZB) = &
     DP3(2:NX-1,2:NY-1,2:NZ-1,2)

     call Grid_releaseBlkPtr(lb,solnData,CENTER)

     goto 112
   endif

enddo
112 continue


return
end

!----------------------------------------------------------------------
! Subroutine : PresDP3P
! Solve the 3D Poisson equation, with periodic boundary conditions
! and 3D ffts. Used for isotropic turbulence runs.
!
! Written by Marcos Vanella, July 2007.
! --------------------------------------------------------------------


      subroutine PresDp3P(nx,ny,nz,dx,dy,dz,fsrc,usol)

      
      implicit none

      integer nx,ny,nz
      real dx,dy,dz
      real fsrc(nx,ny,nz),usol(nx,ny,nz)


      ! Local variables

      real, parameter :: pi = 3.14159265358979

      integer l, n2mh, i , j , k
      real, save :: dx2q,dy2q,dz2q
      real, save, allocatable, dimension(:) :: kx,ky,kz, &
                                   wsavex,wsavey,wsavez     

      complex, save, allocatable, dimension(:,:,:) :: var_c

      logical, save :: firstcall = .true.

      ! Fill used data:
      if (firstcall) then

         allocate(var_c(nx,ny,nz))
         allocate(kx(nx),ky(ny),kz(nz))
         kx = 0.; ky = 0.; kz = 0.;


         !...............................................................
         !     Modified wave numbers
         !...............................................................


         ! X wavenumbers:
         dx2q = 1.0/(dx*dx);
         n2mh = (nx-2)/2;

         do l=1 , n2mh
            kx(l+1) = cos(2.*pi*real(l)/real(nx-2));
            kx(l+1) = 2.*(1.-kx(l+1))*dx2q;
         enddo
         do l=n2mh+1 , (nx-2)-1;
            kx(l+1) = kx(nx-1-l);
         enddo 


         ! Y wavenumbers:
         dy2q = 1.0/(dy*dy);
         n2mh = (ny-2)/2;

         do l=1 , n2mh
            ky(l+1) = cos(2.*pi*real(l)/real(ny-2));
            ky(l+1) = 2.*(1.-ky(l+1))*dy2q;
         enddo
         do l=n2mh+1 , (ny-2)-1
            ky(l+1) = ky(ny-1-l);
         enddo 

         ! Z wavenumbers:
         dz2q = 1.0/(dz*dz);
         n2mh = (nz-2)/2;

         do l=1 , n2mh
            kz(l+1) = cos(2.*pi*real(l)/real(nz-2));
            kz(l+1) = 2.*(1.-kz(l+1))*dz2q;
         enddo
         do l=n2mh+1 , (nz-2)-1
            kz(l+1) = kz(nz-1-l);
         enddo 


         ! Allocate wsave
         allocate(wsavex(4*(nx-2)+15))
         allocate(wsavey(4*(ny-2)+15))
         allocate(wsavez(4*(nz-2)+15))


         ! Initialize FFTS in each dir:
         call CFFTI(NX-2,WSAVEX)
         call CFFTI(NY-2,WSAVEY)
         call CFFTI(NZ-2,WSAVEZ)

         
         firstcall = .false.

      endif


      ! Dump source into complex var_c
      var_c = cmplx(fsrc) 

      ! From real to wave space ...
      ! First fft in the x direction:

      do k = 2 , nz-1
        do j = 2 , ny-1
           
           call CFFTF(NX-2,var_c(2:nx-1,j,k),WSAVEX)

        enddo
      enddo
      var_c = var_c * cmplx(nx-2)**(-1.)

   

      ! Second fft in the y direction:

      do k = 2 , nz-1
         do i = 2 ,nx-1

           call CFFTF(NY-2,var_c(i,2:ny-1,k),WSAVEY)

        enddo
      enddo        
      var_c = var_c * cmplx(ny-2)**(-1.)


      ! Third fft in the z direction:
      do j = 2 , ny-1
         do i = 2 ,nx-1

           call CFFTF(NZ-2,var_c(i,j,2:nz-1),WSAVEZ)

        enddo
      enddo                       
      var_c = var_c * cmplx(nz-2)**(-1.)    

      ! solution
      do k = 2 , nz-1
         do j = 2 , ny-1
            do i = 2 ,nx-1

       var_c(i,j,k)= -var_c(i,j,k)*cmplx((kx(i-1)+ky(j-1)+kz(k-1))**(-1.))            
           
            enddo
         enddo
      enddo

      var_c(2,2,2) = cmplx(0.);


      ! From wave space to real space...
      ! Inverse ffts:
      ! First ifft in the x direction:
      do k = 2 , nz-1
        do j = 2 , ny-1
           
           call CFFTB(NX-2,var_c(2:nx-1,j,k),WSAVEX)

        enddo
      enddo

      ! Second ifft in the y direction:
      do k = 2 , nz-1
         do i = 2 ,nx-1

           call CFFTB(NY-2,var_c(i,2:ny-1,k),WSAVEY)

        enddo
      enddo        

      ! Third ifft in the z direction:
      do j = 2 , ny-1
         do i = 2 ,nx-1

           call CFFTB(NZ-2,var_c(i,j,2:nz-1),WSAVEZ)

        enddo
      enddo                       

      
      ! Dump solution into Usol:
      usol = real(var_c)


      return

      end subroutine




