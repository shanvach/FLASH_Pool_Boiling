!!****if* source/Simulation/SimulationMain/magnetoHD/LinearWave/Simulation_initBlock
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
!!  Reference:
!!
!!  Parameters:  blockID      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***

#define ALFVEN_WAVE  1
#define FAST_WAVE    2
#define SLOW_WAVE    3
#define ENTROPY_WAVE 4

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell,sim_nx,sim_ny,sim_steady,&
                              sim_smallP,sim_gamma,sim_xMax,sim_yMax,&
                              sim_dens0,sim_B0,sim_pres0,sim_choice,&
                              sim_delperturb

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, Grid_getDeltas, &
    Grid_getPointData,  Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_fillGuardCells

  use Simulation_data, ONLY : sim_killdivb

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Arguments ------------------------
  integer, intent(in) :: blockID
  !! ----------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, dx, dy
  real :: a2,ca2,cs2,cf2,in_perturb

  real, allocatable,dimension(:) :: xCoord,xLCoord,xRCoord,yCoord,yLCoord,yRCoord,zCoord
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: x1,x2,rot
  real :: B1,B2,B3,V1,V2,V3,p,d

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az
#else
  real, allocatable, dimension(:,:,:) :: Az
#endif

  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord (sizeX),stat=istat)
  allocate(xLCoord(sizeX),stat=istat)
  allocate(xRCoord(sizeX),stat=istat)
  allocate(yCoord (sizeY),stat=istat)
  allocate(yLCoord(sizeY),stat=istat)
  allocate(yRCoord(sizeY),stat=istat)
  allocate(zCoord (sizeZ),stat=istat)

  xCoord = 0.0
  xLCoord= 0.0
  xRCoord= 0.0
  yCoord = 0.0
  yLCoord= 0.0
  yRCoord= 0.0
  zCoord = 0.0

#ifndef FIXEDBLOCKSIZE
  allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
#endif


  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockId, CENTER,     sim_gCell, yCoord,  sizeY)
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE,  sim_gCell, yLCoord, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, sim_gCell, yRCoord, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockId, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE,  sim_gCell, xLCoord, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, sim_gCell, xRCoord, sizeX)

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  rot = atan(sim_nx/sim_ny)

!!$  !! Initialize wave speeds
!!$  a2  = sim_gamma*sim_pres0/sim_dens0
!!$  ca2 = sim_B0**2/sim_dens0
!!$  cf2 = 0.5*(a2+ca2+sqrt(ca2**2+a2**2))
!!$  cs2 = 0.5*(a2+ca2-sqrt(ca2**2+a2**2))
!!$
!!$  nx = sim_nx/sqrt(sim_nx**2+sim_ny**2)
!!$  ny = sim_ny/sqrt(sim_nx**2+sim_ny**2)
!!$
!!$  if (sim_choice .eq. ALFVEN_WAVE) then
!!$     bx  = 1./sqrt(2.)
!!$     by  = bx
!!$  else
!!$     bx = cos(PI/4.)*nx - sin(PI/4.)*ny
!!$     by = sin(PI/4.)*nx + cos(PI/4.)*ny
!!$  endif

!!$  sim_magx0 = sim_B0
!!$  sim_magy0 = sim_B0
!!$  sim_magz0 = 0.


  p = 1./sim_gamma
  d = 1.

  B1 = 1.
  B2 = sqrt(2.)
  B3 = 0.5
  V1 = 0.
  V2 = 0.
  V3 = 0.

  a2  = sim_gamma*p/d
  ca1_2 = B1**2/d
  ca2   = sqrt(B1**2+B2**2+B3**2)/sqrt(d)
  cf2 = 0.5*(a2+ca2+sqrt((a2+ca2)**2-4.*a2*ca1_2))
  cs2 = 0.5*(a2+ca2-sqrt((a2+ca2)**2-4.*a2*ca1_2))

  !! / x1 \    / cos  sin  0 \
  !! | x2 |  = |-sin  cos  0 |
  !! \ x3 /    \  0    0   1 /

  !! Note Az = A3


  if (sim_choice .eq. ALFVEN_WAVE) then !alfven test
     Rk(1) = 0.
     Rk(2) = 0.
     Rk(3) = 1.
     Rk(4) =-2.*sqrt(2.)
     Rk(5) = 0.
     Rk(6) =-1.
     Rk(7) = 2.*sqrt(2.)
     Rk(8) = 0.
     Rk    = Rk/3.

  elseif (sim_choice .eq. FAST_WAVE) then !fast
     Rk(1) = 6.
     Rk(2) = 12.
     Rk(3) =-4.*sqrt(2.)
     Rk(4) =-2.
     Rk(5) = 0.
     Rk(6) = 8.*sqrt(2.)
     Rk(7) = 4.
     Rk(8) = 27.
     Rk    = Rk/(6.*sqrt(5.))

  elseif (sim_choice .eq. SLOW_WAVE) then !slow
     Rk(1) = 12.
     Rk(2) = 6.
     Rk(3) = 8.*sqrt(2.)
     Rk(4) = 4.
     Rk(5) = 0.
     Rk(6) =-4.*sqrt(2.)
     Rk(7) =-2.
     Rk(8) = 9.
     Rk    = Rk/(6.*sqrt(5.))

  elseif (sim_choice .eq. ENTROPY_WAVE) then !slow
     Rk(1) = 2.
     Rk(2) = 2.
     Rk(3) = 0.
     Rk(4) = 0.
     Rk(5) = 0.
     Rk(6) = 0.
     Rk(7) = 0.
     Rk(8) = 1.
     Rk    = Rk/2.

  endif


  do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
     do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
        if (i <= blkLimitsGC(HIGH,IAXIS) .and. j <= blkLimitsGC(HIGH,JAXIS)) then
           x1 = xLCoord(i)*cos(rot)+yLCoord(j)*sin(rot)
           x2 =-xLCoord(i)*sin(rot)+yLCoord(j)*cos(rot)
        elseif (i == blkLimitsGC(HIGH,IAXIS)+1 .and. j <= blkLimitsGC(HIGH,JAXIS)) then
           x1 = xRCoord(i-1)*cos(rot)+yLCoord(j)*sin(rot)
           x2 =-xRCoord(i-1)*sin(rot)+yLCoord(j)*cos(rot)
        elseif (i <= blkLimitsGC(HIGH,IAXIS) .and. j == blkLimitsGC(HIGH,JAXIS)+1) then
           x1 = xLCoord(i)*cos(rot)+yRCoord(j-1)*sin(rot)
           x2 =-xLCoord(i)*sin(rot)+yRCoord(j-1)*cos(rot)
        elseif (i == blkLimitsGC(HIGH,IAXIS)+1 .and. j == blkLimitsGC(HIGH,JAXIS)+1) then
           x1 = xRCoord(i-1)*cos(rot)+yRCoord(j-1)*sin(rot)
           x2 =-xRCoord(i-1)*sin(rot)+yRCoord(j-1)*cos(rot)
        endif
        Az(i,j) = B1*x2 - B2*x1
     enddo
  enddo


  if (sim_steady) then !standing
     if (sim_choice .eq. ALFVEN_WAVE) then !alfven
        V1 =-sqrt(ca1_2)
        V2 = 0.
        V3 = 0.
     elseif (sim_choice .eq. FAST_WAVE) then !fast
        V1 =-sqrt(cf2)
        V2 = 0.
        V3 = 0.
     elseif (sim_choice .eq. SLOW_WAVE) then !slow
        V1 =-sqrt(cs2)
        V2 = 0.
        V3 = 0.
     elseif (sim_choice .eq. ENTROPY_WAVE) then !entropy
        V1 = 0.
        V2 = 0.
        V3 = 0.
     endif
  else !traveling
     if (sim_choice .eq. ENTROPY_WAVE) then !alfven
        V1 = 1.
        V2 = 0.
        V3 = 0.        
     else
        V1 = 0.
        V2 = 0.
        V3 = 0.
     endif
  endif





  k=1
  do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
     do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

        x1 = xLCoord(i)*cos(rot)+yLCoord(j)*sin(rot)
        cos_perturb = sim_delperturb*cos(2.*PI*x1)

        
!!        totalEner= 

        q(1) = d    + Rk(1)*cos_perturb
        q(2) = d*V1 + Rk(2)*cos_perturb
        q(3) = d*V2 + Rk(3)*cos_perturb
        q(4) = d*V3 + Rk(4)*cos_perturb
        q(5) = B1   + Rk(5)*cos_perturb
        q(6) = B2   + Rk(6)*cos_perturb
        q(7) = B3   + Rk(7)*cos_perturb

        q(8) = V1 + Rk(8)*cos_perturb




        solnData(DENS_VAR,i,j,k)= d + Rk(1)*cos_perturb

        solnData(VELX_VAR,i,j,k)= (V1+ Rk(2)*cos_perturb)/solnData(DENS_VAR,i,j,k)
        solnData(VELY_VAR,i,j,k)= (V2+ Rk(3)*cos_perturb)/solnData(DENS_VAR,i,j,k)
        solnData(VELZ_VAR,i,j,k)= (V3+ Rk(4)*cos_perturb)/solnData(DENS_VAR,i,j,k)

        solnData(PRES_VAR,i,j,k)= p + Rk(1)*cos_perturb
!!        solnData(TEMP_VAR,i,j,k)= 

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


#if NFACE_VARS > 0
  k=1
  do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
        if (sim_killdivb) then
           facexData(MAG_FACE_VAR,i,j,k)= (Az(i,j+1,k)-Az(i,j,k))/dy
           faceyData(MAG_FACE_VAR,i,j,k)=-(Az(i+1,j,k)-Az(i,j,k))/dx

           facexData(MAG_FACE_VAR,i,j,k)= sim_B0*bx
           faceyData(MAG_FACE_VAR,i,j,k)= sim_B0*by
        endif
     enddo
  enddo

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
           solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))

           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
#if NDIM == 1
           solnData(DIVB_VAR,i,j,k) = 0.
#elif NDIM >= 2
           solnData(DIVB_VAR,i,j,k)= &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy
#if NDIM == 3
           solnData(DIVB_VAR,i,j,k)= solnData(DIVB_VAR,i,j,k) &
                   + (facezData(MAG_FACE_VAR,i,  j,  k+1) - facezData(MAG_FACE_VAR,i,j,k))/dz
#endif
#endif
        enddo
     enddo
  enddo
#endif


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  deallocate(xCoord)
  deallocate(xLCoord)
  deallocate(xRCoord)
  deallocate(yCoord)
  deallocate(yLCoord)
  deallocate(yRCoord)
  deallocate(zCoord)

#ifndef FIXEDBLOCKSIZE
  deallocate(Az)
#endif


end subroutine Simulation_initBlock



