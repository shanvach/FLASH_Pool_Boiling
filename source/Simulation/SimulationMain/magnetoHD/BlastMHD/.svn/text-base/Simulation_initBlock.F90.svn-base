!!****if* source/Simulation/SimulationMain/magnetoHD/BlastMHD/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference:  Gardiner & Stone JCP 205(2005),509-539
!!
!!  Parameters:  blockID      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! PARAMETERS
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell, sim_radius, sim_beta_inner, sim_beta_outer,&
                              sim_velx, sim_vely, sim_velz, sim_magx, sim_magy, sim_magz, &
                              sim_gamma, sim_smallP, sim_killdivb

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
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax, istat
  integer :: sizeX, sizeY, sizeZ
  real    :: radius, beta,enerZone,ekinZone,eintZone
  integer, dimension(2,MDIM)    :: blkLimits, blkLimitsGC
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  
  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  if (NDIM >= 2) allocate(yCoord(sizeY),stat=istat)
  if (NDIM == 3) allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  if (NDIM >= 2) yCoord = 0.0
  if (NDIM == 3) zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockId, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockId, CENTER,sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, sim_gCell, xCoord, sizeX)
!------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! Note that we assume that the center of the computational is origin
           if (NDIM == 2) then
              radius = sqrt(xCoord(i)**2 + yCoord(j)**2)
           elseif (NDIM == 3) then
              radius = sqrt(xCoord(i)**2 + yCoord(j)**2 + zCoord(k)**2)
           endif

           if (radius < sim_radius) then
              beta = sim_beta_inner
           else
              beta = sim_beta_outer
           endif

           ! initialize cells.
           solnData(DENS_VAR,i,j,k)  = 1.
           solnData(VELX_VAR,i,j,k) = sim_velx
           solnData(VELY_VAR,i,j,k) = sim_vely
           solnData(VELZ_VAR,i,j,k) = sim_velz

           ! initialize pressure depending on hydro or mhd setup
           solnData(PRES_VAR,i,j,k) = .5*beta*(sim_magx**2 + sim_magy**2 + sim_magz**2)

           solnData(MAGX_VAR,i,j,k) = sim_magx
           solnData(MAGY_VAR,i,j,k) = sim_magy
           solnData(MAGZ_VAR,i,j,k) = sim_magz
           solnData(MAGP_VAR,i,j,k) = .5*(sim_magx**2 + sim_magy**2 + sim_magz**2)
           solnData(DIVB_VAR,i,j,k) = 0.

           solnData(GAMC_VAR,i,j,k) = sim_gamma
           solnData(GAME_VAR,i,j,k) = sim_gamma

           ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                        solnData(VELX_VAR:VELZ_VAR,i,j,k))

           ! specific internal energy
           eintZone = solnData(PRES_VAR,i,j,k)/(solnData(GAME_VAR,i,j,k)-1.)/solnData(DENS_VAR,i,j,k)

           ! total specific gas energy
           enerZone = eintZone + ekinZone

           ! Take a limit value
           enerZone = max(enerZone, sim_smallP)

           solnData(ENER_VAR,i,j,k) = enerZone
           solnData(EINT_VAR,i,j,k) = eintZone
           solnData(TEMP_VAR,i,j,k) = 1.

        enddo
     enddo
  enddo

  deallocate(xCoord)
  if (NDIM >= 2) deallocate(yCoord)
  if (NDIM == 3) deallocate(zCoord)


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)


#if NFACE_VARS > 0
  if (sim_killdivb) then

     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)

     facexData(MAG_FACE_VAR,:,:,:) = sim_magx
     faceyData(MAG_FACE_VAR,:,:,:) = sim_magy
     if (NDIM==3) facezData(MAG_FACE_VAR,:,:,:) = sim_magz

     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


end subroutine Simulation_initBlock










