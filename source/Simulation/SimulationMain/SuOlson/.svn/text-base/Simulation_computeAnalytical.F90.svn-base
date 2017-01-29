!!****if* source/Simulation/SimulationMain/SuOlson/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up a linear conduction problem
!!  to test conduction in a medium with constant isochoric conduction coefficient.
!!  The exact three-dimensional solution for the initial condition
!!  
!!                 T (r, t = 0) = Q delta (0)
!!
!!  is 
!!
!!       T (r, t) = Q / (4 pi \chi t)^(3/2) exp [-r^2 / 4 \chi t]
!!
!!  Here we set up the initial condition with the exact solution
!!  slightly offset from t = 0.
!!
!!  Reference: 
!!
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Simulation_data
  use Conductivity_interface, ONLY : Conductivity
  use Eos_interface, ONLY : Eos
  !use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
  !  Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY: Driver_abortFlash
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  real :: xx, yy, zz, x0, y0, z0, r2, toff, pi
  real :: nexpo, xi, fOfXi
  integer :: ivelx, ively, ivelz

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,vecLen
  integer, dimension(MDIM) :: axis

  real, dimension(:), allocatable :: eosData
  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone, eintZone, tempZone
  real :: cond_zone, diff_coeff, cond_chi
  real, parameter :: bogusUnusedTemperature = -1e20  ! Set to 3.e-5 for Power Law 
  real, save :: bogusUnusedMassFrac(NSPECIES)
  
  logical :: gcell = .true.


  ! Parameters for Analytical Su-Oslon
  real :: clight, ssol, asol, alpha
  real :: erad,trad,trad_ev,tmat,tmat_ev, sim_erad,Ref_Energy

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec



  pi = 4. * atan (1.d0)

  ! WE COULD dump some output to stdout listing the parameters.
  ! But currently we don't.
  if (sim_meshMe == MASTER_PE) then

1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!        03    format(1x,i4,1p8e14.6)
   open(unit=3,file="compare.dat",status='unknown')

   40    format(1x,i4,1p8e14.6)

  endif




  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
  vecLen = 1
  allocate(eosData(vecLen * EOS_NUM))

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
!  the state using an exact solution.

  call Grid_getBlkPtr(1, solnVec)

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCenter(k)

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCenter(j)
                
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx  = xCenter(i)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! Lets try to incorporate Su-oslon solution here.
           ! trad_bc_ev = 1.0E3, opac = 1.0, alpha = 4.0d0*asol

           call PhysicalConstants_get("speed of light",clight)         
           call PhysicalConstants_get("Stefan-Boltzmann",ssol)

           asol    = 4.0 * ssol / clight
           alpha   = 4.0 * asol           ! Relation between C_V and TELE, CV = alpha * TELE^3

           !call Conductivity(bogusUnusedTemperature, sim_rhoInit, bogusUnusedMassFrac, cond_zone, diff_coeff, 3)                     
          
           !call so_wave(tcurr,xx,1.0E3,cond_zone,alpha,erad,trad,trad_ev,tmat,tmat_ev)          

           !write(*,*) i,xx, tmat, trad
           !write (*,*) trad, tmat
           !write(6,40) i,xx, solnVec(TRAD_VAR,i,j,k), trad, solnVec(TELE_VAR,i,j,k), tmat
           !write(3,40) i,xx, solnVec(TRAD_VAR,i,j,k), trad, solnVec(TELE_VAR,i,j,k), tmat

           if (xx .ge. 0.0) then

               Ref_Energy = asol*(1.0E3*1.1604505E4)**4
               write(6,40) i,xx, solnVec(ERAD_VAR,i,j,k)/(Ref_Energy), & 
                                 solnVec(EELE_VAR,i,j,k)/(Ref_Energy), &
                                 solnVec(EION_VAR,i,j,k)/(Ref_Energy)

               write(3,40) i,xx, solnVec(ERAD_VAR,i,j,k)/(Ref_Energy), & 
                                 solnVec(EELE_VAR,i,j,k)/(Ref_Energy), &
                                 solnVec(EION_VAR,i,j,k)/(Ref_Energy)
           endif

            

           !call Grid_putPointData(blockId, CENTER, TMPA_VAR, EXTERIOR, axis, tmat_ev)

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(1, solnVec)

   close(unit=3)


  deallocate(eosData)
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)


  
  return

end subroutine Simulation_computeAnalytical
