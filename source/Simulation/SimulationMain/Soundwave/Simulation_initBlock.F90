!!****if* source/Simulation/SimulationMain/Soundwave/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) : blockID)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!!  Parameters:  blockID      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockID          the number of the block to update
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz
  
  real :: lPosn0, lPosn
  real :: wavenumber, pi, phase, ciso
  integer :: ivelx, ively, ivelz

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,istat
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone
  
  logical :: gcell = .true.
  

  pi = 3.1415926535897931
  ciso = 1.
  
  


  ! dump some output to stdout listing the paramters
!!$  if (sim_meshMe == MASTER_PE) then
!!$     
!!$     
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$     
!!$  endif
  
  

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCenter(sizeY),stat=istat)
  allocate(zCenter(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

!  lPosn = sim_posn !- zz*sim_zCos/sim_xCos

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCenter(k)

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCenter(j)
                
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx  = xCenter(i)

           if (sim_orientation .eq. 2) then
              xx = yy
           else if (sim_orientation .eq. 3) then
              xx = zz   
           end if

           wavenumber = 2.0 * pi / sim_wavelength
           phase = wavenumber * xx

           rhoZone = sim_rhoInit * (1 + sim_perturbAmp * sin (phase) )
           velxZone = sim_perturbAmp * sim_cs * sin (phase)
           presZone   = sim_cs**2. * rhoZone / sim_gamma
           velyZone = 0.
           velzZone = 0.

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

#if NSPECIES > 0
           call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                axis, 1.0e0-(NSPECIES-12)*sim_smallX)

           !if there is only 1 species, this loop will not execute
           do n = SPECIES_BEGIN+1,SPECIES_END
              call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                   axis, sim_smallX)
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
           enerZone = max(enerZone, sim_smallP)
           
           ! store the variables in the current zone via the database put methods
           ! data is put stored one cell at a time with this call to Grid_putData           

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, 0.)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, 0.)

           if (sim_orientation .eq. 2) then

              call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, 0.)
              call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velxZone)
              call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, 0.)

           else if (sim_orientation .eq. 3) then

              call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, 0.)
              call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, 0.)
              call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velxZone)

           end if

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
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
 
  return
end subroutine Simulation_initBlock









