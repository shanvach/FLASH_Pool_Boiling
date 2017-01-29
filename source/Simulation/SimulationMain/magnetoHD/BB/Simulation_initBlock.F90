!!****if* source/Simulation/SimulationMain/magnetoHD/BB/Simulation_initBlock
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
       Grid_getCellCoords, Grid_putPointData, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData

  integer, intent(in) :: blockId
  
  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: rho, tele, trad, tion, zbar, abar
  integer :: species

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

  
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  #if NFACE_VARS > 0
    if (sim_killdivb) then
       call Grid_getBlkPtr(blockID,facexData,FACEX)
       call Grid_getBlkPtr(blockID,faceyData,FACEY)
    endif
  #endif
  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           species = CHAM_SPEC

           if ( xcent(i) <= sim_targetRadius .and. &
                ycent(j) <= 0.0 .and. &
                ycent(j) >= -sim_targetHeight ) then
              species = TARG_SPEC
           end if

           if(species == TARG_SPEC) then
              rho = sim_rhoTarg
              tele = sim_teleTarg
              tion = sim_tionTarg
              trad = sim_tradTarg
              zbar = sim_zbarTarg
              abar = sim_abarTarg
           else
              rho = sim_rhoCham
              tele = sim_teleCham
              tion = sim_tionCham
              trad = sim_tradCham
              zbar = sim_zbarCham
              abar = sim_abarCham
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, YE_MSCALAR, EXTERIOR, axis, zbar/abar)
           call Grid_putPointData(blockId, CENTER, SUMY_MSCALAR, EXTERIOR, axis, 1/abar)

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)

           if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined.
              ! We put nearly all the mass into either the Xe material if XE_SPEC is defined,
              ! or else into the first species.
              do n = SPECIES_BEGIN,SPECIES_END
                 if (n==species) then
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)
                 else
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
                 end if
              enddo
           end if
           
           solnData(MAGX_VAR,i,j,k) = 0.0
           solnData(MAGY_VAR,i,j,k) = 0.0
           solnData(MAGZ_VAR,i,j,k) = 0.0

           #if NFACE_VARS > 0
              if (sim_killdivb) then
                    facexData(MAG_FACE_VAR,i,j,k)= 0.0
                    faceyData(MAG_FACE_VAR,i,j,k)= 0.0
              endif
           #endif
           
           #ifdef BDRY_VAR
                solnData(BDRY_VAR, i, j, k) = -1.0
           #endif

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  #if NFACE_VARS > 0
    if (sim_killdivb) then
       call Grid_releaseBlkPtr(blockID,facexData,FACEX)
       call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
    endif
  #endif

  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  return

end subroutine Simulation_initBlock
