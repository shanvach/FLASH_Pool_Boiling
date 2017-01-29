!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/ConstBier/Simulation_initBlock
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
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr
  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)

  real :: rho  ! Density
  real :: tion ! Ion temperature
  real :: tele ! Electron temperature 
  real :: nele
  real :: pele
  real :: abar
  real :: zbar
  real :: KB
  real :: NA

  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, pointer     :: soln(:,:,:,:),facex(:,:,:,:),facey(:,:,:,:)

  abar = sim_singleSpeciesA
  zbar = sim_singleSpeciesZ
  KB   = sim_boltzmann
  NA   = sim_avogadro

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))

  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))

  call Grid_getBlkPtr(blockId, soln)
  call Grid_getBlkPtr(blockId, facex,FACEX)
  call Grid_getBlkPtr(blockId, facey,FACEY)

  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! Compute mass density in this cell:
           nele = (ycent(j)-sim_ymin) * (sim_nele2-sim_nele1)/(sim_ymax-sim_ymin) + sim_nele1
           rho = abar/(NA*zbar) * nele

           ! Compute electron temperature:
           pele = (xcent(i)-sim_xmin) * (sim_pele2-sim_pele1)/(sim_xmax-sim_xmin) + sim_pele1
           tele = pele/(nele*KB)

           ! Compute ion temperature:
           tion = zbar*(sim_ptot/(nele*KB) - tele)

           soln(DENS_VAR,i,j,k) = rho
           soln(TELE_VAR,i,j,k) = tele
           soln(TION_VAR,i,j,k) = tion
           soln(TEMP_VAR,i,j,k) = tele
        enddo
     enddo
  enddo


  facex = 0.
  facey = 0.
  call Grid_releaseBlkPtr(blockId, soln)
  call Grid_releaseBlkPtr(blockId, facex,FACEX)
  call Grid_releaseBlkPtr(blockId, facey,FACEY)

  deallocate(xcent)
  deallocate(ycent)
  
  return

end subroutine Simulation_initBlock
