!!****if* source/Simulation/SimulationMain/Advect/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockId, 
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up a plane-wave advection
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! ARGUMENTS
!!
!!  blockId -          the number of the block to update
!!
!!  PARAMETERS
!!
!!     smallP              smallest pressure allowed
!!     smallX              smallest abundance allowed 
!!     gamma               Gamma value from the EOS
!!     sim_rhoIn           Density inside the pulse
!!     sim_rhoOut          Density outside the pulse
!!     sim_pressure        Pressure
!!     sim_velocity        Fluid velocity
!!     sim_posn            Position of the pulse center at x-axis (y=z=0)
!!     sim_width           Width of the pulse along x-axis
!!     sim_xAngle          Angle made by diaphragm normal w/x-axis (deg)
!!     sim_yAngle          Angle made by diaphragm normal w/y-axis (deg)
!!     sim_pulseFunctn     Which pulse shape function to use
!!***

subroutine Simulation_initBlock(blockId)
  
#include "constants.h"
#include "Flash.h"
  
  use Simulation_data, ONLY: sim_gamma, sim_smallP, &
       sim_smallX, sim_velocity, sim_pressure, sim_rhoIn, &
       sim_rhoOut, sim_posn, sim_xCos, sim_yCos, sim_zCos, sim_pulseFunctn, &
       sim_width
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData

  implicit none

  
  integer, intent(in) :: blockId
  

  integer :: i, j, k, n


  integer :: ii, jj, kk, nInt
  real :: nIntInv1
  real :: sumRho
  real :: delX, xx, delY, yy, delZ, zz
  
  real :: lPosn0, lPosn
  real :: xPos
  
  real :: xMin, xMax, yMin, yMax, zMin, zMax
  
  real :: rhoZone, velXZone, velYZone, velZZone, presZone, & 
       eintZone, enerZone, ekinZone

  real :: wFac
  real :: sim_pulseshape

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(2,MDIM)    :: boundBox

  integer, dimension(MDIM) :: axis,guard
  real, dimension(MDIM) ::  delta





!==============================================================================


  !! Grid_getBlkIndexLimits returns the correct indices for describing 
  !! the size of a block:
  !! blkLimits contains the starting and end indices of the block interior;
  !! blkLimitsGC contains the starting and end indices of the block including
  !! guard cells.

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Grid_getBlkBoundBox(blockId,boundBox)
  
  guard = blkLimits(LOW,:)-blkLimitsGC(LOW,:)

  ! compute the extrema of the current blocks
  xMax = boundBox(HIGH,IAXIS)
  xMin = boundBox(LOW,IAXIS)

  yMax = boundBox(HIGH,JAXIS)
  yMin = boundBox(LOW,JAXIS)

  zMax = boundBox(HIGH,KAXIS)
  zMin = boundBox(LOW,KAXIS)

  call Grid_getDeltas(blockId,delta)

  delX = delta(IAXIS)
  delY = delta(JAXIS)
  delZ = delta(KAXIS)


  !=============================================================================

  ! Loop over cells in the block.  For each, compute the physical position of 
  ! its left and right edge and its center as well as its physical width.  Then 
  ! subcycle to initialize density according to the given shape function, then
  ! initialize the remaining hydro variables appropriately.
  
  nInt = 5
  nIntInv1 = 1.0/(real(nInt) - 1.0)
  
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           
           ! For each zone, find the average 
           ! density by computing the pulse shape at a 
           ! number of points inside the zone.
           sumRho = 0.
           
           do kk = 0, nInt-1
              
              ! compute the subzone z coordinate
              zz = zMin + delZ*(real(k-guard(KAXIS)-1)+kk*nIntInv1)
              lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
              
              do jj = 0, nInt-1
                 
                 ! compute the subzone y coordinate
                 yy = yMin + delY*(real(j-guard(JAXIS)-1)+jj*nIntInv1)
                 lPosn = lPosn0 - yy*sim_yCos/sim_xCos
                 
                 do ii = 0, nInt-1
                    
                    ! compute the subzone x coordinate
                    xx  = xMin + delX*(real(i-guard(IAXIS)-1)+ii*nIntInv1)
                    xPos = (xx - lPosn) * sim_xCos / sim_width
                    
                    ! compute the weighting fraction for the pulse shape, and use this to
                    ! weight the current subzone's contribution to the current zone
                    wFac = sim_pulseshape(xPos, sim_pulseFunctn)
                    
                    sumRho = sumRho + sim_rhoIn*wFac +  & 
                         &                       sim_rhoOut*(1.-wFac) 

                 enddo
                 
              enddo
              
           enddo
           
           ! Initialize the hydro quantities.
           rhoZone  = sumRho / (real(nInt)**3)
           presZone = sim_pressure
           
           velXZone = sim_velocity * sim_xCos
           velYZone = sim_velocity * sim_yCos
           velZZone = sim_velocity * sim_zCos

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           ! initialize the nuclear abundances
           if (NSPECIES > 0) then
              !put in value of default species
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if only 1 species then loop will not execute
              do n=SPECIES_BEGIN+1,SPECIES_END

                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)

              enddo
           end if
           
           
           ! compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           eKinZone = 0.5*(velXZone**2 + & 
                &                          velYZone**2 + & 
                &                          velZZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           
           
           ! update the variables in the current zone via a Grid put method
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)   
#endif
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

        enddo
     enddo
  enddo
  
  return
  
end subroutine Simulation_initBlock






