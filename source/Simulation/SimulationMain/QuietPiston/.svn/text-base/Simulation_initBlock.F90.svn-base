!!****if* source/Simulation/SimulationMain/QuietPiston/Simulation_initBlock
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
!!  a specified block.  This version sets up the QuietPiston problem.
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_pAmbient, sim_rhoAmbient, sim_windVel, sim_gamma, &
     &  sim_smallP, sim_smallX, sim_temp
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords
  use Grid_data, ONLY: gr_quietStartBnd, gr_quietStartDens, gr_pistonDens, &
       gr_pistonVelx, gr_pistonVely, gr_pistonVelz, gr_pistonBnd

  use Logfile_interface, ONLY : Logfile_stamp


  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID
  real,pointer :: solnData(:,:,:,:)

  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: i,j,k
  logical :: gcell = .true.
  integer,dimension(MDIM) :: dSizes
  integer :: maxDim
  real,allocatable, dimension(:,:) :: coords

  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone, temp_zone

  logical, dimension(MDIM) :: inRegion


!===============================================================================


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  dSizes = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  maxDim=maxval(dSizes)
  allocate(coords(maxDim,NDIM))

  do k = 1, NDIM 
     call Grid_getCellCoords(k,blockID,CENTER,gcell,coords(:,k),dsizes(k))
  end do

  call Grid_getBlkPtr(blockID, solnData, CENTER)

!  print *, gr_customRegionBnd(LOW,IAXIS),gr_customRegionBnd(HIGH,IAXIS)
!  print *, gr_customRegionBnd(LOW,JAXIS),gr_customRegionBnd(HIGH,JAXIS)
!  stop

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           rho_zone = sim_rhoAmbient
           pres_zone = sim_pAmbient
           temp_zone = sim_temp

           velx_zone = 0.0
           vely_zone = 0.0
           velz_zone = 0.0
           
              
  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
           ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)
           
           eint_zone = pres_zone / (sim_gamma-1.)
           eint_zone = eint_zone / rho_zone
           ener_zone = eint_zone + ekin_zone
           ener_zone = max(ener_zone, sim_smallP)

! Now check if in quietStart region

           inRegion(IAXIS) = coords(i,IAXIS) > gr_quietStartBnd(LOW,IAXIS) .AND. &
                coords(i,IAXIS) < gr_quietStartBnd(HIGH,IAXIS)

           if (NDIM >= 2) then
              inRegion(JAXIS) = coords(j,JAXIS) > gr_quietStartBnd(LOW,JAXIS) .AND. &
                coords(j,JAXIS) < gr_quietStartBnd(HIGH,JAXIS)
           else
              inRegion(JAXIS) = .TRUE.
           end if

           if (NDIM == 3) then
              inRegion(KAXIS) = coords(k,KAXIS) > gr_quietStartBnd(LOW,KAXIS) .AND. &
                coords(k,KAXIS) < gr_quietStartBnd(HIGH,KAXIS)
           else
              inRegion(KAXIS) = .TRUE.
           end if
              
           if ( ALL(inRegion .eqv. .TRUE.)) then

              rho_zone  = gr_quietStartDens
              pres_zone = sim_pAmbient
              temp_zone = sim_temp

              velx_zone = 0.0
              vely_zone = 0.0
              velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
              ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)
              
              eint_zone = pres_zone / (sim_gamma-1.)
              eint_zone = eint_zone / rho_zone
              ener_zone = eint_zone + ekin_zone
              ener_zone = max(ener_zone, sim_smallP)
              
           end if

!Now check if in piston region

           inRegion(IAXIS) = coords(i,IAXIS) > gr_pistonBnd(LOW,IAXIS) .AND. &
                coords(i,IAXIS) < gr_pistonBnd(HIGH,IAXIS)

           if (NDIM >= 2) then
              inRegion(JAXIS) = coords(j,JAXIS) > gr_pistonBnd(LOW,JAXIS) .AND. &
                coords(j,JAXIS) < gr_pistonBnd(HIGH,JAXIS)
           else
              inRegion(JAXIS) = .TRUE.
           end if

           if (NDIM == 3) then
              inRegion(KAXIS) = coords(k,KAXIS) > gr_pistonBnd(LOW,KAXIS) .AND. &
                coords(k,KAXIS) < gr_pistonBnd(HIGH,KAXIS)
           else
              inRegion(KAXIS) = .TRUE.
           end if
              
           if ( ALL(inRegion .eqv. .TRUE.)) then

              rho_zone = gr_pistonDens
              pres_zone = sim_pAmbient
              
              velx_zone = gr_pistonVelx
              vely_zone = gr_pistonVely
              velz_zone = gr_pistonVelz

  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
              ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)
              
              eint_zone = pres_zone / (sim_gamma-1.)
              eint_zone = eint_zone / rho_zone
              ener_zone = eint_zone + ekin_zone
              ener_zone = max(ener_zone, sim_smallP)
              
           end if


#if NSPECIES > 0
           solnData(SPECIES_BEGIN,i,j,k) =  1.0-(NSPECIES-1)*sim_smallX
           solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k) =     sim_smallX
#endif
  
  ! store the variables in the block's unk data
           solnData(DENS_VAR,i,j,k) = rho_zone
           solnData(PRES_VAR,i,j,k) = pres_zone
           solnData(ENER_VAR,i,j,k) = ener_zone
           solnData(TEMP_VAR,i,j,k) = temp_zone
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k) = eint_zone
#endif
           solnData(GAMC_VAR,i,j,k) = sim_gamma
           solnData(GAME_VAR,i,j,k) = sim_gamma
           
           
           solnData(VELX_VAR,i,j,k) = velx_zone
           solnData(VELY_VAR,i,j,k) = vely_zone
           solnData(VELZ_VAR,i,j,k) = velz_zone
           
        end do
     end do
  end do

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
 

  return
end subroutine Simulation_initBlock



