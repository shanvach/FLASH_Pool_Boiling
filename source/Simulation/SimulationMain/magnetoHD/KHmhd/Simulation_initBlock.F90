!!****if* source/Simulation/SimulationMain/magnetoHD/KHmhd/Simulation_initBlock
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
!!  a specified block.
!!
!!
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell, sim_gamma,     &
                              sim_smallX, sim_smallP,   &
                              sim_killdivb, sim_Vx0,    &
                              sim_Bx0, sim_By0, sim_Bz0,&
                              sim_epsilon, sim_xMin,    &
                              sim_xMax, sim_yMin, sim_yMax

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real :: ymid, yvsc, coff
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData, facezData


  ! Make sure enough fluids have been defined.
  if (NSPECIES /=2) then
     call Driver_abortFlash("Simulation_initBlock needs two fluids!")
  endif

  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)


  ymid = 0.5*(sim_yMin+sim_yMax)
  yvsc = (sim_xMax-sim_xMin)/20.
  coff = (sim_xMax-sim_xMin)/10.


#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! Multiple species
           if (yCoord(j) < ymid) then
              solnData(HEAV_SPEC,i,j,k)=1.-sim_smallX
              solnData(LIGH_SPEC,i,j,k)=sim_smallX
           else
              solnData(HEAV_SPEC,i,j,k)=sim_smallX
              solnData(LIGH_SPEC,i,j,k)=1.-sim_smallX
           endif

           ! Cell-centered values
           solnData(DENS_VAR,i,j,k)=  1.
           solnData(VELX_VAR,i,j,k)=  sim_Vx0*tanh((yCoord(j)-ymid)/yvsc)
           solnData(VELY_VAR,i,j,k)=  sim_epsilon*sim_Vx0*sin(2.0*PI*xCoord(i))*exp(-((yCoord(j)-ymid)/coff)**2)
           solnData(VELZ_VAR,i,j,k)=  0.
           solnData(PRES_VAR,i,j,k)=  1./sim_gamma

           solnData(MAGX_VAR,i,j,k)= sim_Bx0
           solnData(MAGY_VAR,i,j,k)= sim_By0
           solnData(MAGZ_VAR,i,j,k)= sim_Bz0
           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           solnData(DIVB_VAR,i,j,k)= 0.

           ! Compute the gas energy and set the gamma-values needed for the EOS
           ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                        solnData(VELX_VAR:VELZ_VAR,i,j,k))

           ! specific internal energy
           eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.)/solnData(DENS_VAR,i,j,k)

           ! total specific gas energy
           enerZone = eintZone + ekinZone

           ! Take a limit value
           enerZone = max(enerZone, sim_smallP)

           solnData(ENER_VAR,i,j,k)=enerZone
           solnData(EINT_VAR,i,j,k)=eintZone
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(TEMP_VAR,i,j,k)=1.

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)=sim_Bx0
              faceyData(MAG_FACE_VAR,i,j,k)=sim_By0
              if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = sim_Bz0
           endif
#endif
        enddo
     enddo
  enddo



!!$#if NFACE_VARS > 0
!!$           if (sim_killdivb) then
!!$
!!$              do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
!!$                 do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
!!$                    do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
!!$                       if ((j < blkLimitsGC(HIGH,JAXIS)+1) .and. (k < blkLimitsGC(HIGH,KAXIS)+1)) then
!!$                          facexData(MAG_FACE_VAR,i,j,k) = sim_Bx0
!!$                       elseif ((i < blkLimitsGC(HIGH,IAXIS)+1) .and. (k < blkLimitsGC(HIGH,KAXIS)+1)) then
!!$                          faceyData(MAG_FACE_VAR,i,j,k) = sim_By0
!!$                       elseif ((i < blkLimitsGC(HIGH,IAXIS)+1) .and. (j < blkLimitsGC(HIGH,JAXIS)+1)) then
!!$                          facezData(MAG_FACE_VAR,i,j,k) = sim_Bz0
!!$                       endif
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$
!!$           endif
!!$#endif


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
