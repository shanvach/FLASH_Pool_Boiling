!!****if* source/Simulation/SimulationMain/magnetoHD/CloudShock/Simulation_initBlock
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
!!  Dai & Woodward, JCP, 142:331--369, 1998
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

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell, sim_gamma,   &
                              sim_smallX, sim_smallP, &
                              sim_killdivb,           &
                              sim_dL,sim_pL,sim_uL,sim_vL,sim_wL,sim_bxL,sim_byL,sim_bzL,&
                              sim_dR,sim_pR,sim_uR,sim_vR,sim_wR,sim_bxR,sim_byR,sim_bzR,&
                              sim_lposn,sim_cloudRadius, sim_cloudXCtr, sim_cloudYCtr,   &
                              sim_cloudZCtr,sim_cloudDensity, sim_Type

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

  r0 = sim_cloudRadius**.98
  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! Multiple species
           !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k)=sim_smallX
           enddo

           ! Cell-centered values
           if (xCoordR(i) < sim_lposn) then
              ! Left states
              solnData(DENS_VAR,i,j,k)=sim_dL
              solnData(VELX_VAR,i,j,k)=sim_uL
              solnData(VELY_VAR,i,j,k)=sim_vL
              solnData(VELZ_VAR,i,j,k)=sim_wL
#if defined(MAGX_VAR)
              solnData(MAGX_VAR,i,j,k)=sim_bxL
#endif
#if defined(MAGY_VAR)
              solnData(MAGY_VAR,i,j,k)=sim_byL
#endif
#if defined(MAGZ_VAR)
              solnData(MAGZ_VAR,i,j,k)=sim_bzL
#endif
              solnData(PRES_VAR,i,j,k)=sim_pL

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_bxL
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_byL
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k)= sim_bzL
              endif
#endif

           elseif ((xCoordL(i) < sim_lposn) .and. (xCoordR(i) > sim_lposn)) then
              solnData(DENS_VAR,i,j,k)=.5*(sim_dL  +  sim_dR)
              solnData(VELX_VAR,i,j,k)=.5*(sim_uL  +  sim_uR)
              solnData(VELY_VAR,i,j,k)=.5*(sim_vL  +  sim_vR)
              solnData(VELZ_VAR,i,j,k)=.5*(sim_wL  +  sim_wR)
#if defined(MAGX_VAR)
              solnData(MAGX_VAR,i,j,k)=.5*(sim_bxL +  sim_bxR)
#endif
#if defined(MAGY_VAR)
              solnData(MAGY_VAR,i,j,k)=.5*(sim_byL +  sim_byR)
#endif
#if defined(MAGZ_VAR)
              solnData(MAGZ_VAR,i,j,k)=.5*(sim_bzL +  sim_bzR)
#endif
              solnData(PRES_VAR,i,j,k)=.5*(sim_pL  +  sim_pR)

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_bxL
                 if (xCoord(i) < sim_lposn) then
                    faceyData(MAG_FACE_VAR,i,j,k)= sim_byL
                 elseif (xCoord(i) == sim_lposn) then
                    faceyData(MAG_FACE_VAR,i,j,k)= .5*(sim_byL+sim_byR)
                 else
                    faceyData(MAG_FACE_VAR,i,j,k)= sim_byR
                 endif
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k)= sim_bzL
              endif
#endif

           else
              ! Right states
              solnData(DENS_VAR,i,j,k)=sim_dR
              solnData(VELX_VAR,i,j,k)=sim_uR
              solnData(VELY_VAR,i,j,k)=sim_vR
              solnData(VELZ_VAR,i,j,k)=sim_wR
#if defined(MAGX_VAR)
              solnData(MAGX_VAR,i,j,k)=sim_bxR
#endif
#if defined(MAGY_VAR)
              solnData(MAGY_VAR,i,j,k)=sim_byR
#endif
#if defined(MAGZ_VAR)
              solnData(MAGZ_VAR,i,j,k)=sim_bzR
#endif
              solnData(PRES_VAR,i,j,k)=sim_pR

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_bxR
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_byR
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k)= sim_bzR
              endif
#endif
           endif

           ! Cloud density
           if (sim_Type == 1) then
              ! sim_Type=1: The original circular cloud
              if (NDIM == 2) then
                 radius = sqrt( (xCoord(i)-sim_cloudXCtr)**2 &
                               +(yCoord(j)-sim_cloudYCtr)**2)
              elseif (NDIM == 3) then
                 radius = sqrt( (xCoord(i)-sim_cloudXCtr)**2 &
                               +(yCoord(j)-sim_cloudYCtr)**2 &
                               +(zCoord(k)-sim_cloudZCtr)**2)
              endif

              taper=(sim_cloudRadius-radius)/(sim_cloudRadius-r0)

              if (radius <= r0) then
                 solnData(DENS_VAR,i,j,k)=sim_cloudDensity
              elseif ((radius > r0).and.(radius < sim_cloudRadius)) then
                 solnData(DENS_VAR,i,j,k)=sim_dR + (sim_cloudDensity-sim_dR)*taper
              endif

           elseif (sim_Type == 2) then
              ! sim_Type=2: A rectangular cloud to test slow-moving shock
              if (NDIM == 2) then
                 if ((abs(xCoord(i)-0.2) .le. 0.1) .and. &
                     (abs(yCoord(j))     .le. sim_cloudRadius)) then

                    solnData(DENS_VAR,i,j,k) = sim_cloudDensity
                 endif
              elseif (NDIM == 3) then
                 if ((abs(xCoord(i)-sim_cloudXCtr) .le. sim_cloudRadius) .and. &
                     (abs(yCoord(j)-sim_cloudYCtr) .le. sim_cloudRadius) .and. &
                     (abs(zCoord(k)-sim_cloudZCtr) .le. sim_cloudRadius)) then
                    solnData(DENS_VAR,i,j,k) = sim_cloudDensity
                 endif
              endif
           endif


#if defined(VECZ_VAR) && defined(MAGX_VAR) && defined(MAGY_VAR)
           solnData(VECZ_VAR,i,j,k) = solnData(MAGX_VAR,i,j,k)*yCoord(j)&
                                     -solnData(MAGY_VAR,i,j,k)*xCoord(i)
#endif
#if defined(MAGP_VAR)
           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
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
