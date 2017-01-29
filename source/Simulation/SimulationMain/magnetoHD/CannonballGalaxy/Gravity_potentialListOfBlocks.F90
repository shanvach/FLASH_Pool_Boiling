!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/Gravity_potentialListOfBlocks
!!
!!  NAME 
!!
!!     Gravity_potentialListOfBlocks
!!
!!  SYNOPSIS
!!
!!  Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                integer(IN) :: blockList(blockCount))
!!
!!  DESCRIPTION 
!!      This routine computes the gravitational potential for the gravity
!!      implementations (i.e., various Poisson implementations) which make
!!      use of it in computing the gravitational acceleration.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!!
!!
!!***

subroutine Gravity_potentialListOfBlocks(blockCount,blockList)

  use Gravity_data, ONLY : useGravity, grv_ptdirn
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords
  use IO_interface, ONLY: IO_setScalar

  ! Hacks (shouldn't use this directly)

  use Simulation_data, ONLY : sim_testSingleGalaxy, sim_vInit, &
       sim_rCtr, sim_vrCtr, sim_restart, sim_arCtr, sim_oarCtr, &
       r1, r2, grav1, grav2, gpot1, gpot2, numPoints1, numPoints2, &
       sim_testAtmosphere, sim_xCtr, sim_yCtr, sim_zCtr
  
  use Driver_data, ONLY : dr_dt, dr_dtOld

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, dimension(:), allocatable :: x, y, z

  real, parameter :: onethird = 1./3.
  real, parameter :: onesixth = 1./6.

  integer       :: ierr

  integer       :: size(3)
  integer       :: lb, i, j, k
  real          :: rr1, rr2, simTime, dtOld, dtNew
  real          :: xc, yc, zc
  integer :: error
  character(len=124) :: str_buffer
  real :: wterm, woldterm
  logical, save :: first_call = .true.
  real, external :: interpolate

  !=========================================================================

  if(.not.useGravity) return
  
  call Timers_start("gravity Barrier")
  call MPI_Barrier (MPI_COMM_WORLD, ierr)
  call Timers_stop("gravity Barrier")

  call Timers_start("gravity")
 
  call Driver_getSimTime(simTime)

  dtOld = dr_dtOld
  dtNew = dr_dt

  ! First update the position of the galaxy, unless we are only
  ! testing the main cluster atmosphere

  if (.not. sim_testAtmosphere) then

     ! The galaxy falls as a point mass within the gravitational 
     ! field of the galaxy cluster, unless we are doing a 
     ! translation test
               
     sim_arCtr = -interpolate(grav1,numPoints1,r1,sim_rCtr)
     
     if (first_call .and. (.not. sim_restart)) then
        wterm = 0.5*dtNew
        woldterm = 0.
        first_call = .false.
     else
        wterm = 0.5*dtNew + onethird*dtOld + onesixth*dtNew**2/dtOld
        woldterm = onesixth*(dtOld**2 - dtNew**2)/dtOld
     endif
     
     if (.not. sim_testSingleGalaxy) then
        sim_vrCtr = sim_vrCtr + wterm*sim_arCtr + woldterm*sim_oarCtr
     endif

     sim_rCtr = sim_rCtr + dtOld*sim_vrCtr 
     
     call IO_setScalar("subcluster r", sim_rCtr)
     call IO_setScalar("subcluster vr", sim_vrCtr)
     call IO_setScalar("subcluster ar", sim_arCtr)
     
     if (grv_ptdirn == 1) then
        sim_xCtr = sim_rCtr
     else if (grv_ptdirn == 2) then
        sim_yCtr = sim_rCtr
     else if (grv_ptdirn == 3) then
        sim_zCtr = sim_rCtr
     endif
     
     sim_oarCtr = sim_arCtr
     
  endif

  do lb = 1, blockCount

     call Grid_getBlkPtr(blocklist(lb), solnVec)

     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
     size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

#ifdef GPOL_VAR
     solnVec(GPOL_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:)
#endif
     solnVec(GPOT_VAR,:,:,:) = 0.0

     allocate(x(size(1)), y(size(2)), z(size(3)))

     call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., z, size(3))
     call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., y, size(2))
     call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., x, size(1))

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)

        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)

           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

              if (.not. sim_testSingleGalaxy) then

                 if (grv_ptdirn == 1) then
                    rr1 = x(i)
                 else if (grv_ptdirn == 2) then
                    rr1 = y(j)
                 else if (grv_ptdirn == 3) then
                    rr1 = z(k)
                 endif

                 solnVec(GPOT_VAR,i,j,k) = -interpolate(gpot1, numPoints1, r1, rr1)
                 
              endif

              if (.not. sim_testAtmosphere) then

                 xc = x(i) - sim_xCtr
                 yc = y(j) - sim_yCtr
                 zc = z(k) - sim_zCtr

                 rr2 = sqrt(xc*xc + yc*yc + zc*zc)

                 solnVec(GPOT_VAR,i,j,k) = solnVec(GPOT_VAR,i,j,k) - &
                      interpolate(gpot2, numPoints2, r2, rr2)

              endif
                 
           enddo

        enddo

     enddo

     call Grid_releaseBlkPtr(blocklist(lb), solnVec)

     deallocate(x)
     deallocate(y)
     deallocate(z)
     
  enddo

#ifdef USEBARS
  call MPI_Barrier (MPI_Comm_World, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_potentialListOfBlocks
