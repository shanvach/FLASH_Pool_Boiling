!!****if* source/Simulation/SimulationMain/DiffuseCtC/Simulation_initBlock
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
       Grid_getCellCoords, Grid_putPointData, Grid_getPointData
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)

  real :: tion  ! Ion temperature
  real :: tele  ! Electron temperature 
  real :: trad  ! Radiation temperature
  real :: tradActual
  real :: position

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

  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           if(sim_initGeom == "planar") position = xcent(i)
           if(sim_initGeom == "polar") then
              if (NDIM == 2) then
                 position = sqrt(xcent(i)*xcent(i) + ycent(j)*ycent(j))
              elseif (NDIM == 3) then
                 position = sqrt( xcent(i)*xcent(i) + &
                      ycent(j)*ycent(j) + &
                      zcent(k)*zcent(k) )
              end if
             end if

           call tanh_func(sim_teleOffset, sim_teleSteepness, &
                sim_tele1, sim_tele2, &
                position, tele)

           call tanh_func(sim_tionOffset, sim_tionSteepness, &
                sim_tion1, sim_tion2, &
                position, tion)

           call tanh_func(sim_tradOffset, sim_tradSteepness, &
                sim_trad1, sim_trad2, &
                position, trad)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rho)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, HE_SPEC,  EXTERIOR, axis, 1.0)
           
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
        enddo
     enddo
  enddo

  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)
  return

contains

  subroutine tanh_func(offset, steep, high, low, x, val)
    implicit none

    real, intent(in)  :: offset
    real, intent(in)  :: steep
    real, intent(in)  :: high
    real, intent(in)  :: low
    real, intent(in)  :: x
    real, intent(out) :: val

    val = tanh(-(x-offset)*steep) * 0.5*(high-low) + 0.5*(high+low)
  end subroutine tanh_func

end subroutine Simulation_initBlock
