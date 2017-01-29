!!****if* source/Simulation/SimulationMain/Noh/Simulation_initBlock
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
!!  a specified block.  This version sets up the Noh problem
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_rhoInit, sim_pInit, sim_uInit, sim_gamma, smallp
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, Grid_getGeometry


  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  integer :: geo  

  real :: xx, yy, zz, rr
  
  real,allocatable, dimension(:) ::xCenter,yCenter,zCenter

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone
  
  logical :: gcell = .true.

  integer :: istat 
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCenter(sizeY),stat=istat)
  allocate(zCenter(sizeZ),stat=istat)
  xCenter = 0.0
  yCenter = 0.0
  zCenter = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)

  call Grid_getGeometry(geo)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCenter(k)
          
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCenter(j)
                
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i) 
           
           rhoZone = sim_rhoInit
           presZone = max(sim_pInit,smallp)
           
           if ((NDIM==1) .or. (geo==POLAR) .or. (geo==SPHERICAL)) then
              velxZone = -sim_uInit
              velyZone = 0.0
              velzZone = 0.0
           else if (NDIM==2) then
              rr = sqrt(xx**2 + yy**2)
              velxZone = -sim_uInit * xx / rr
              velyZone = -sim_uInit * yy / rr
              velzZone = 0.0
           else if (NDIM==3) then
              rr = sqrt(xx**2 + yy**2 + zz*2)
              velxZone = -sim_uInit * xx / rr
              velyZone = -sim_uInit * yy / rr
              velzZone = -sim_uInit * zz / rr
           endif

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

#if NSPECIES > 0
           !put in value of default species
           call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                axis, 1.0e0-(NSPECIES-1)*0.0)

           !if there is only 1 species, this loop will not execute
           do n = SPECIES_BEGIN+1,SPECIES_END

              call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                   axis, 0.0)
           enddo
#endif

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
           enerZone = presZone / (sim_gamma-1.)
           enerZone = enerZone / rhoZone
           enerZone = enerZone + ekinZone
           enerZone = max(enerZone, smallp)
           
           ! store the variables in the current zone via the database put methods
           ! data is put stored one cell at a time with this call to Grid_putData           

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif

        enddo
     enddo
  enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)

 
  return
end subroutine Simulation_initBlock










