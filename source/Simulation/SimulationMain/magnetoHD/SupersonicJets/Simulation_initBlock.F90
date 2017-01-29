!!****if* source/Simulation/SimulationMain/magnetoHD/SupersonicJets/Simulation_initBlock
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
!!  Jones, Ryu, Tregillis, ApJ, 473:365-382, 1996
!!  Ryu, Miniati, Jones, Frank, ApJ, 509:244--255, 1998
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  real,allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, taper, radius, r0
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData


  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord (sizeX),stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)
  allocate(yCoord (sizeY),stat=istat)
  allocate(zCoord (sizeZ),stat=istat)

  xCoord = 0.0
  xCoordL= 0.0
  xCoordR= 0.0
  yCoord = 0.0
  zCoord = 0.0

  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE,  sim_gCell, xCoordL, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, sim_gCell, xCoordR, sizeX)
  call Grid_getCellCoords(JAXIS, blockID, CENTER,     sim_gCell, yCoord,  sizeY)
  call Grid_getCellCoords(KAXIS, blockID, CENTER,     sim_gCell, zCoord,  sizeZ)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

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
           !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k)=sim_smallX
           enddo

           
           solnData(DENS_VAR,i,j,k)=sim_gamma*sim_dAmbient
           solnData(VELX_VAR,i,j,k)=sim_uAmbient
           solnData(VELY_VAR,i,j,k)=sim_vAmbient
           solnData(VELZ_VAR,i,j,k)=sim_wAmbient
           solnData(MAGX_VAR,i,j,k)=sim_bxAmbient
           solnData(MAGY_VAR,i,j,k)=sim_byAmbient
           solnData(MAGZ_VAR,i,j,k)=sim_bzAmbient
           solnData(PRES_VAR,i,j,k)=1.

           solnData(HEAV_SPEC,i,j,k)=sim_smallX
           solnData(LIGH_SPEC,i,j,k)=1.-sim_smallX

           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

#ifdef VECZ_VAR
           solnData(VECZ_VAR,i,j,k) = -solnData(MAGY_VAR,i,j,k)*xCoord(i) + &
                                       solnData(MAGX_VAR,i,j,k)*yCoord(j)
#endif

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

        enddo
     enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     facexData(MAG_FACE_VAR,:,:,:) = sim_bxAmbient
     faceyData(MAG_FACE_VAR,:,:,:) = sim_byAmbient
     if (NDIM == 3) facezData(MAG_FACE_VAR,:,:,:) = sim_bzAmbient

     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(xCoordL)
  deallocate(xCoordR)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
