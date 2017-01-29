!!****if* source/Simulation/SimulationMain/magnetoHD/BierShock/Simulation_initBlock
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
                             Grid_getCellCoords, &
                             Grid_getBlkPtr, &
                             Grid_releaseBlkPtr
  
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  real, allocatable :: xcent(:), ycent(:)
  real :: tradActual
  real :: rho, tele, trad, tion, zbar, abar
  integer :: species
  real, pointer :: U(:,:,:,:)

  real :: eele_min
  real :: eele_max
  real :: eele

  ! Get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! Get cell center coordinates:
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, blkLimitsGC(HIGH, IAXIS))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, blkLimitsGC(HIGH, JAXIS))

  call Grid_getBlkPtr(blockID, U, CENTER)

  zbar = sim_singleSpeciesZ

  eele = sim_eint * zbar/(1.0+zbar)
  eele_min = eele * (1.0-sim_skewFactor)
  eele_max = eele + (eele - eele_min)

  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           eele = (eele_max-eele_min)/(sim_ymax-sim_ymin)*ycent(j) + eele_min

           U(ERAD_VAR,i,j,k) = sim_erad
           ! U(EELE_VAR,i,j,k) = sim_eint * zbar / (1.0+zbar)
           U(EELE_VAR,i,j,k) = eele

           U(EION_VAR,i,j,k) = sim_eint - sim_erad - eele

           U(VELX_VAR,i,j,k) = sim_velx
           U(DENS_VAR,i,j,k) = sim_rho
           
        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,U,CENTER)
  deallocate(xcent)
  deallocate(ycent)

  return

end subroutine Simulation_initBlock
