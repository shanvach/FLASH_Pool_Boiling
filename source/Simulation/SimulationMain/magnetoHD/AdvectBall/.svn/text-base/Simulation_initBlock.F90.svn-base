!!****if* source/Simulation/SimulationMain/magnetoHD/AdvectBall/Simulation_initBlock
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
!!  Parameters:  blockID      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!!
!! NOTE
!!
!!  This problem simply advects a 3D high density and pressure ball in a 3D
!!  periodicall box, and especially useful for testing given solver's CFL limit.
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell,sim_xCtr,sim_yCtr,sim_ballRadius,&
                              sim_rx,sim_ry,&
                              sim_Ux_initial,sim_Uy_initial,sim_Uz_initial,&
                              sim_smallP,sim_gamma

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getDeltas, &
                             Grid_getBlkPtr, &
                             Grid_releaseBlkPtr

  use Simulation_data, ONLY : sim_killdivb

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Arguments ------------------------
  integer, intent(in) :: blockID
  !! ----------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, rot, radius, dx, dy, dz
  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: x1,x2,x3,cos_ang,sin_ang,lambda
  real :: xx,yy,zz

  ! dump some output to stdout listing the paramters
!!$!!$   if (sim_meshMe == MASTER_PE) then
!!$!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)

  allocate(yCoord(sizeY), stat=istat)
  allocate(yCoordL(sizeY),stat=istat)
  allocate(yCoordR(sizeY),stat=istat)

  allocate(zCoord(sizeZ), stat=istat)
  allocate(zCoordL(sizeZ),stat=istat)
  allocate(zCoordR(sizeZ),stat=istat)

  xCoord  = 0.0
  xCoordL = 0.0
  xCoordR = 0.0

  yCoord  = 0.0
  yCoordL = 0.0
  yCoordR = 0.0

  zCoord  = 0.0
  zCoordL = 0.0
  zCoordR = 0.0


  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS,blockId,CENTER,    sim_gCell,zCoord, sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,LEFT_EDGE, sim_gCell,zCoordL,sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,RIGHT_EDGE,sim_gCell,zCoordR,sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS,blockId,CENTER,    sim_gCell,yCoord, sizeY)
     call Grid_getCellCoords(JAXIS,blockId,LEFT_EDGE, sim_gCell,yCoordL,sizeY)
     call Grid_getCellCoords(JAXIS,blockId,RIGHT_EDGE,sim_gCell,yCoordR,sizeY)
  endif

  call Grid_getCellCoords(IAXIS,blockId,CENTER,    sim_gCell,xCoord, sizeX)
  call Grid_getCellCoords(IAXIS,blockId,LEFT_EDGE, sim_gCell,xCoordL,sizeX)
  call Grid_getCellCoords(IAXIS,blockId,RIGHT_EDGE,sim_gCell,xCoordR,sizeX)

  rot = atan(sim_rx/sim_ry)
  cos_ang = cos(rot)
  sin_ang = sin(rot)
  lambda  = cos_ang


  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)


  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)


           radius = sqrt(xCoord(i)**2 + yCoord(j)**2 + zCoord(k)**2)

           if (radius <= sim_ballRadius) then
              solnData(DENS_VAR,i,j,k) = 10. !2. !10.
           else
              solnData(DENS_VAR,i,j,k) = 1. !1.
           endif

           solnData(PRES_VAR,i,j,k) = 1. !1.
           solnData(TEMP_VAR,i,j,k) = 1.
           solnData(VELX_VAR,i,j,k) = sim_Ux_initial
           solnData(VELY_VAR,i,j,k) = sim_Uy_initial
           solnData(VELZ_VAR,i,j,k) = sim_Uz_initial

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR) && defined(DIVB_VAR)
           if (NDIM == 2) then
              solnData(MAGX_VAR,i,j,k) = 0.
              solnData(MAGY_VAR,i,j,k) = 0.
           elseif (NDIM == 3) then
              solnData(MAGX_VAR,i,j,k) = 0.
              solnData(MAGY_VAR,i,j,k) = 0.
           endif
           solnData(MAGZ_VAR,i,j,k) = 0.

           ! Update the magnetic pressure
           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           ! Initialize with zero divB
           solnData(DIVB_VAR,i,j,k) = 0.

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


#if NFACE_VARS > 0
           !! In this case we initialized Az using the cell-cornered coordinates.
           if (sim_killdivb) then
              if (NDIM == 2) then
                 facexData(MAG_FACE_VAR,i,j,k)= 0.
                 faceyData(MAG_FACE_VAR,i,j,k)= 0.
              elseif (NDIM == 3) then
                 facexData(MAG_FACE_VAR,i,j,k)= 0.
                 faceyData(MAG_FACE_VAR,i,j,k)= 0.
                 facezData(MAG_FACE_VAR,i,j,k)= 0.
              endif
           endif
#endif

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
  deallocate(yCoordL)
  deallocate(yCoordR)

  deallocate(zCoord)
  deallocate(zCoordL)
  deallocate(zCoordR)


end subroutine Simulation_initBlock



