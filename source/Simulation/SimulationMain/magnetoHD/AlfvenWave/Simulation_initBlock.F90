!!****if* source/Simulation/SimulationMain/magnetoHD/AlfvenWave/Simulation_initBlock
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
!!  Reference: Gardiner & Stone JCP 205(2005),509-539
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

subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, Grid_getDeltas, &
    Grid_getPointData,  Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_fillGuardCells

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Arguments ------------------------
  integer, intent(in) :: blockID
  !! ----------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, rot, dx, dy,dz
  real :: x1, x2, x3, xx, yy, zz, B1, B2, B3, v1, v2, v3
  real :: a1,a2,a3
  real :: cos_ang2,sin_ang2,cos_ang3,sin_ang3,lambda
  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
  real :: k_perp,k_par,b_perp,b_par,sn,cs

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Ax,Ay,Az
#else
  real, allocatable, dimension(:,:,:) :: Ax,Ay,Az
#endif

  ! dump some output to stdout listing the paramters
!!$!!$   if (sim_meshMe == MASTER_PE) then
!!$!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(yCoordL(sizeY),stat=istat)
  allocate(yCoordR(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)  
  allocate(zCoordL(sizeZ),stat=istat)
  allocate(zCoordR(sizeZ),stat=istat)


  xCoord = 0.0
  xCoordL= 0.0
  xCoordR= 0.0
  yCoord = 0.0
  yCoordL= 0.0
  yCoordR= 0.0
  zCoord = 0.0
  zCoordL= 0.0
  zCoordR= 0.0

#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  elseif (NDIM == 3) then
     allocate(Ax(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
  endif
#endif

  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS, blockId, CENTER,    sim_gCell, zCoord,  sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, sim_gCell, zCoordL, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE,sim_gCell, zCoordR, sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockId, CENTER,     sim_gCell, yCoord,  sizeY)
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE,  sim_gCell, yCoordL, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, sim_gCell, yCoordR, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockId, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE,  sim_gCell, xCoordL, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, sim_gCell, xCoordR, sizeX)

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  if (NDIM == 3) dz=del(3)

  !! Initialize parameters
  rot = atan((sim_xMax-sim_xMin)/(sim_yMax-sim_yMin))
  cos_ang3 = cos(rot)
  sin_ang3 = sin(rot)

  rot =((sim_xMax-sim_xMin)*cos_ang3 + (sim_yMax-sim_yMin)*sin_ang3)/(sim_zMax-sim_zMin)
  rot =atan(0.5*rot)
  cos_ang2 = cos(rot)
  sin_ang2 = sin(rot)

  x1=(sim_xMax-sim_xMin)*cos_ang2*cos_ang3
  x2=(sim_yMax-sim_yMin)*cos_ang2*sin_ang3
  x3=(sim_zMax-sim_zMin)*sin_ang2

  ! choose the smaller of the three values for lambda
  lambda = x1
  lambda = min(lambda,x2)
  lambda = min(lambda,x3)
  k_par = 2.*PI/lambda


  Az = 0.
  if (NDIM == 3) then
     Ax = 0.
     Ay = 0.
  endif

  B1 = 1. !parallel B component
  B2 = sim_B0 !magnitude of the Alfven wave in perpendicular direction
  b_par = B1
  b_perp= B2

  if (sim_steady) then
     v1 = B1
  else
     v1 = 0.
  endif


  if (NDIM == 2) then
     !! Construct Az at each cell corner
     !! Bx=dAz/dy & By=-dAz/dx 
     k=1
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

#if NFACE_VARS > 1
           ! x Coord at cell corner
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              xx = xCoordL(i)
           else
              xx = xCoordR(i-1)
           endif

           ! y Coord at cell corner
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              yy = yCoordL(j)
           else
              yy = yCoordR(j-1)
           endif
#else
           ! x Coord at cell center
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              xx = xCoord(i)
           else
              xx = xCoord(i-1) + dx
           endif

           ! y Coord at cell center
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              yy = yCoord(j)
           else
              yy = yCoord(j-1) + dy
           endif
#endif

           x1 = xx*cos_ang3 + yy*sin_ang3
           x2 =-xx*sin_ang3 + yy*cos_ang3
           Az(i,j,k) = b_par*x2+b_perp/k_par*cos(k_par*x1)
        enddo
     enddo
  elseif (NDIM == 3) then
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

#if NFACE_VARS > 0
              ! x Coord at cell corner
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoordL(i)
              else
                 xx = xCoordR(i-1)
              endif

              ! y Coord at cell corner
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoordL(j)
              else
                 yy = yCoordR(j-1)
              endif

              ! z Coord at cell corner
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoordL(k)
              else
                 zz = zCoordR(k-1)
              endif

#else
              ! x Coord at cell center
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoord(i)
              else
                 xx = xCoord(i-1) + dx
              endif

              ! y Coord at cell center
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoord(j)
              else
                 yy = yCoord(j-1) + dy
              endif

              ! z Coord at cell center
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoord(k)
              else
                 zz = zCoord(k-1) + dz
              endif
#endif
              !! Compute vector potentials Ax,Ay,Az
              !(1) for Ax
              !First compute rotated coordinate system
              x1= xCoord(i)*cos_ang2*cos_ang3+yy*cos_ang2*sin_ang3+zz*sin_ang2
              x2=-xCoord(i)*sin_ang3+yy*cos_ang3
              a2= b_perp/k_par*sin(k_par*x1)
              a3= b_perp/k_par*cos(k_par*x1)+b_par*x2
              Ax(i,j,k)=-a2*sin_ang3-a3*sin_ang2*cos_ang3

              !(2) for Ay
              !First compute rotated coordinate system
              x1= xx*cos_ang2*cos_ang3+yCoord(j)*cos_ang2*sin_ang3+zz*sin_ang2
              x2=-xx*sin_ang3+yCoord(j)*cos_ang3
              a2= b_perp/k_par*sin(k_par*x1)
              a3= b_perp/k_par*cos(k_par*x1)+b_par*x2
              Ay(i,j,k)= a2*cos_ang3-a3*sin_ang2*sin_ang3

              !(3) for Az
              !First compute rotated coordinate system
              x1= xx*cos_ang2*cos_ang3+yy*cos_ang2*sin_ang3+zCoord(k)*sin_ang2
              x2=-xx*sin_ang3+yy*cos_ang3
              a3= b_perp/k_par*cos(k_par*x1)+b_par*x2
              Az(i,j,k)=a3*cos_ang2

           enddo
        enddo
     enddo
  endif


  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  if (NDIM == 2) then
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           x1 = xCoord(i)*cos_ang3+yCoord(j)*sin_ang3
           x2 =-xCoord(i)*sin_ang3+yCoord(j)*cos_ang3

           v2 = sim_U0*sin(2.*PI*x1)
           v3 = sim_U0*cos(2.*PI*x1)
           B3 = sim_B0*cos(2.*PI*x1)

           solnData(DENS_VAR,i,j,k)=  1.
           solnData(VELX_VAR,i,j,k)=  v1*cos_ang3 - v2*sin_ang3
           solnData(VELY_VAR,i,j,k)=  v1*sin_ang3 + v2*cos_ang3
           solnData(VELZ_VAR,i,j,k)=  v3
           solnData(PRES_VAR,i,j,k)=  sim_P0
           solnData(TEMP_VAR,i,j,k)=  1.

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

#ifdef XNOR_VAR
           solnData(XNOR_VAR,i,j,k) = x1
#endif
!!$#ifdef BTRA_VAR
!!$        solnData(BTRA_VAR,i,j,k) =-solnData(MAGX_VAR,i,j,k)*sin_ang3&
!!$                                  +solnData(MAGY_VAR,i,j,k)*cos_ang3
!!$#endif

#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)= (Az(i,j+1,k)-Az(i,j,k))/dy
              faceyData(MAG_FACE_VAR,i,j,k)=-(Az(i+1,j,k)-Az(i,j,k))/dx
           endif
#else
           solnData(MAGX_VAR,i,j,k)= 0.5*(Az(i,j+1,k)-Az(i,j-1,k))/dy
           solnData(MAGY_VAR,i,j,k)=-0.5*(Az(i+1,j,k)-Az(i-1,j,k))/dx
           solnData(DIVB_VAR,i,j,k)= 0.0
#endif
           solnData(MAGZ_VAR,i,j,k)=  B3
     enddo
  enddo
enddo


  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
#if NFACE_VARS > 0
           solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
           solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))

           solnData(DIVB_VAR,i,j,k)= &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy
#endif
           ! Update the magnetic pressure
           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

        enddo
     enddo
  enddo



  elseif (NDIM == 3) then
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           x1 = xCoord(i)*cos_ang2*cos_ang3 + yCoord(j)*cos_ang2*sin_ang3+zCoord(k)*sin_ang2
           sn = sin(k_par*x1)
           cs = cos(k_par*x1)

           ! v1 = either B1 or zero
           v2 = sim_U0*sn
           v3 = sim_U0*cs

           solnData(DENS_VAR,i,j,k)=  1.
           ! this seems wrong - following athena new version
!!$           solnData(VELX_VAR,i,j,k)=  v1*cos_ang2*cos_ang3 + v2*sin_ang3 + v3*sin_ang2*cos_ang3
!!$           solnData(VELY_VAR,i,j,k)=  v1*cos_ang2*sin_ang3 - v2*cos_ang3 + v3*sin_ang2*sin_ang3
!!$           solnData(VELZ_VAR,i,j,k)=  v1*sin_ang2                        - v3*cos_ang2
           ! this seems correct - following athena old version + paper eqn 54
           solnData(VELX_VAR,i,j,k)=  v1*cos_ang2*cos_ang3 - v2*sin_ang3 - v3*sin_ang2*cos_ang3
           solnData(VELY_VAR,i,j,k)=  v1*cos_ang2*sin_ang3 + v2*cos_ang3 - v3*sin_ang2*sin_ang3
           solnData(VELZ_VAR,i,j,k)=  v1*sin_ang2                        + v3*cos_ang2

           solnData(PRES_VAR,i,j,k)=  sim_P0
           solnData(TEMP_VAR,i,j,k)=  1.

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

#ifdef XNOR_VAR
           solnData(XNOR_VAR,i,j,k) = x1
#endif
!!$#ifdef BTRA_VAR
!!$        solnData(BTRA_VAR,i,j,k) =-solnData(MAGX_VAR,i,j,k)*sin_ang3&
!!$                                  +solnData(MAGY_VAR,i,j,k)*cos_ang3
!!$#endif

#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)= &
                   (Az(i,j+1,k)-Az(i,j,k))/dy - (Ay(i,j,k+1)-Ay(i,j,k))/dz

              faceyData(MAG_FACE_VAR,i,j,k)=&
                   (Ax(i,j,k+1)-Ax(i,j,k))/dz - (Az(i+1,j,k)-Az(i,j,k))/dx

              facezData(MAG_FACE_VAR,i,j,k)=&
                   (Ay(i+1,j,k)-Ay(i,j,k))/dx - (Ax(i,j+1,k)-Ax(i,j,k))/dy
           endif
#else
           solnData(MAGX_VAR,i,j,k)= &
                   0.5*((Az(i,j+1,k)-Az(i,j,k-1))/dy - (Ay(i,j,k+1)-Ay(i,j,k-1))/dz)

           solnData(MAGY_VAR,i,j,k)=&
                   0.5*((Ax(i,j,k+1)-Ax(i,j,k-1))/dz - (Az(i+1,j,k)-Az(i-1,j,k))/dx)

           solnData(MAGZ_VAR,i,j,k)=&
                   0.5*((Ay(i+1,j,k)-Ay(i-1,j,k))/dx - (Ax(i,j+1,k)-Ax(i,j-1,k))/dy)

           solnData(DIVB_VAR,i,j,k)= 0.0
#endif
     enddo
  enddo
enddo


  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
#if NFACE_VARS > 0
           solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
           solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))
           solnData(MAGZ_VAR,i,j,k) = 0.5*(facezData(MAG_FACE_VAR,i,j,k)+facezData(MAG_FACE_VAR,i,j,k+1))

           solnData(DIVB_VAR,i,j,k)= &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy &
                   + (facezData(MAG_FACE_VAR,i,  j,  k+1) - facezData(MAG_FACE_VAR,i,j,k))/dz
#endif
           ! Update the magnetic pressure
           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

        enddo
     enddo
  enddo
  endif


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

#ifndef FIXEDBLOCKSIZE
  deallocate(Az)
  if (NDIM == 3) then
     deallocate(Ax)
     deallocate(Ay)
  endif
!!$  deallocate(AzVel)
#endif

end subroutine Simulation_initBlock



