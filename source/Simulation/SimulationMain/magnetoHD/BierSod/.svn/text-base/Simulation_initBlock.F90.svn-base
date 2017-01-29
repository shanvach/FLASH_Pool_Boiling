!!****if* source/Simulation/SimulationMain/magnetoHD/BierSod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_putPointData,      &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getDeltas
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: dens, tele, trad, tion
  integer :: species
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData

  integer :: ndx, ndy, ndz, nsctot
  integer :: nsct, nscc, l, m, n
  real :: xsc, ysc, zsc, targFrac, chamFrac, norm
  real :: delta(MDIM)
  real :: r, phi

  ndx = sim_ndiv
  ndy = 1
  ndz = 1
  if(NDIM >= 2) ndy = sim_ndiv
  if(NDIM == 3) ndz = sim_ndiv
  nsctot = ndx * ndy * ndz

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))

  call Grid_getDeltas(blockId, delta)

  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! ****************************
           ! *                          *
           ! *     APPLY FEATHERING     *
           ! *                          *
           ! ****************************
           nsct = 0
           do l = 1, ndx
              do m = 1, ndy
                 do n = 1, ndz
                    xsc = (xcent(i)-delta(IAXIS)/2) + (delta(IAXIS)*(l-1))/ndx
                    ysc = (ycent(j)-delta(JAXIS)/2) + (delta(JAXIS)*(m-1))/ndy
                    zsc = (zcent(k)-delta(KAXIS)/2) + (delta(KAXIS)*(n-1))/ndz

                    r = sqrt(xsc**2 + ysc**2)
                    if (r <= sim_radius) then
                       nsct = nsct + 1
                    end if

                 end do
              end do
           end do

           nscc = nsctot - nsct

           ! Average Cell Density:                   
           dens = (sim_densTarg * nsct + sim_densCham * nscc)/nsctot

           ! Set temperatures:
           tele = (sim_teleTarg * nsct + sim_teleCham * nscc)/nsctot
           tion = (sim_tionTarg * nsct + sim_tionCham * nscc)/nsctot

           ! Set initial condition:
           solnData(DENS_VAR,i,j,k) = dens
           solnData(TEMP_VAR,i,j,k) = tele
#ifdef TION_VAR
           solnData(TION_VAR,i,j,k) = tion
#endif
#ifdef TELE_VAR
           solnData(TELE_VAR,i,j,k) = tele
#endif
        enddo
     enddo
  enddo

  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Simulation_initBlock
