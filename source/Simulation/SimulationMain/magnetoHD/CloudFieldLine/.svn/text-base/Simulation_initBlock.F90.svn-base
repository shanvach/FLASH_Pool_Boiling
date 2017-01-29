!!****if* source/Simulation/SimulationMain/magnetoHD/CloudFieldLine/Simulation_initBlock
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
!!  
!!  Initializes fluid data for a specified block.
!!  This version sets up the Brio-Wu MHD problem.
!!
!! 
!!  Reference:   Brio, M. and Wu, C. C., J. Comput. Phys., 75, 400, 1988
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

  use Simulation_data, ONLY : sim_gCell,   sim_posn,   sim_xcos,   sim_ycos,   sim_zcos,    &
                              sim_dLeft, sim_uLeft,  sim_vLeft,  sim_wLeft,  sim_pLeft,   &
                              sim_dRight,sim_uRight, sim_vRight, sim_wRight, sim_pRight,  &
                              sim_gamma,  sim_smallX, sim_smallP, sim_killdivb,&
                              sim_rx,      sim_ry,     sim_xmax,   sim_xmin,   sim_ymin, sim_ymax,&
                              sim_cloudRadius, sim_cloudXCtr, sim_cloudYCtr,   &
                              sim_cloudZCtr,sim_cloudDensity

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, &
                             Grid_getBlkPtr,     &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"


  !!$ Arguments -------------------
  integer, intent(in) :: blockID
  !!$------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: xx, xxL, xxR, yy, zz, lposn0, lposn
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord, xLCoord, xRCoord,yLCoord,yRCoord
  real :: enerZone, ekinZone, eintZone
  real :: rot,x1,x2,dx,dy
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  real :: r0, taper, radius
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1) :: Az
#else
  real, allocatable, dimension(:,:) :: Az
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

  allocate(xCoord(sizeX), stat=istat)
  allocate(xLCoord(sizeX),stat=istat)
  allocate(xRCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY), stat=istat)
  allocate(yLCoord(sizeY),stat=istat)
  allocate(yRCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ), stat=istat)

  xCoord  = 0.0
  xLCoord = 0.0
  xRCoord = 0.0
  yCoord  = 0.0
  yLCoord  = 0.0
  yRCoord  = 0.0
  zCoord  = 0.0

#ifndef FIXEDBLOCKSIZE
  allocate(Az(sizeX+1,sizeY+1),stat=istat)
#endif

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yCoord,  sizeY)
     call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yLCoord, sizeY)
     call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yRCoord, sizeY)
  endif

  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE,  sim_gCell, xLCoord, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, sim_gCell, xRCoord, sizeX)
!------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)


#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  !! Given a coordinate transformation:
  !! x1 = x*cos_ang + y*sin_ang
  !! x2 =-x*sin_ang + y*cos_ang
  !!
  !! B1 = dAz/dx2
  !! B2 =-dAz/dx1
  !!
  !! Bx = dAz/dy = dAz/dx2*dx2/dy = dAz/dx2*cos_ang
  !! By =-dAz/dx =-dAz/dx1*dx1/dx = dAz/dx1*cos_ang


  if (NDIM >= 2) then
     rot = atan(sim_rx/sim_ry)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS) +1
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS) +1

#if NFACE_VARS > 1
           if (i <= blkLimitsGC(HIGH,IAXIS) .and. j <= blkLimitsGC(HIGH,JAXIS)) then
              x1 = xLCoord(i)*sim_xcos+yLCoord(j)*sim_ycos
              x2 =-xLCoord(i)*sim_ycos+yLCoord(j)*sim_xcos
              yy = yCoord(j)
              xx = xLCoord(i)

           elseif (i == blkLimitsGC(HIGH,IAXIS)+1 .and. j <= blkLimitsGC(HIGH,JAXIS)) then
              x1 = xRCoord(i-1)*sim_xcos+yLCoord(j)*sim_ycos
              x2 =-xRCoord(i-1)*sim_ycos+yLCoord(j)*sim_xcos
              yy = yCoord(j)
              xx = xRCoord(i-1)

           elseif (i <= blkLimitsGC(HIGH,IAXIS) .and. j == blkLimitsGC(HIGH,JAXIS)+1) then
              x1 = xLCoord(i)*sim_xcos+yRCoord(j-1)*sim_ycos
              x2 =-xLCoord(i)*sim_ycos+yRCoord(j-1)*sim_xcos
              yy = yRCoord(j-1)
              xx = xLCoord(i)

           elseif (i == blkLimitsGC(HIGH,IAXIS)+1 .and. j == blkLimitsGC(HIGH,JAXIS)+1) then
              x1 = xRCoord(i-1)*sim_xcos+yRCoord(j-1)*sim_ycos
              x2 =-xRCoord(i-1)*sim_ycos+yRCoord(j-1)*sim_xcos
              yy = yRCoord(j-1)
              xx = xRCoord(i-1)
           endif
#else
           yy = yCoord(j)
           xx = xCoord(i)
#endif

           Az(i,j) = -exp(xx**4.)

        enddo
     enddo
  endif
  Az = Az * sim_xcos


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

           xx  = xCoord (i)
           xxL = xLCoord(i)
           xxR = xRCoord(i)

           ! Cell-centered values
           if (xxR < sim_posn) then
              ! Left states
              solnData(DENS_VAR,i,j,k)=sim_dLeft
              solnData(VELX_VAR,i,j,k)=sim_uLeft
              solnData(VELY_VAR,i,j,k)=sim_vLeft
              solnData(VELZ_VAR,i,j,k)=sim_wLeft
              solnData(PRES_VAR,i,j,k)=sim_pLeft

           elseif ((xxL < sim_posn) .and. (xxR > sim_posn)) then
              solnData(DENS_VAR,i,j,k)=.5*(sim_dLeft  +  sim_dRight)
              solnData(VELX_VAR,i,j,k)=.5*(sim_uLeft  +  sim_uRight)
              solnData(VELY_VAR,i,j,k)=.5*(sim_vLeft  +  sim_vRight)
              solnData(VELZ_VAR,i,j,k)=.5*(sim_wLeft  +  sim_wRight)
              solnData(PRES_VAR,i,j,k)=.5*(sim_pLeft  +  sim_pRight)

           else
              ! Right states
              solnData(DENS_VAR,i,j,k)=sim_dRight
              solnData(VELX_VAR,i,j,k)=sim_uRight
              solnData(VELY_VAR,i,j,k)=sim_vRight
              solnData(VELZ_VAR,i,j,k)=sim_wRight
              solnData(PRES_VAR,i,j,k)=sim_pRight

           endif

           ! Cloud density
           if (NDIM == 2) then
              radius = sqrt( (xCoord(i)-sim_cloudXCtr)**2 &
                            +(yCoord(j)-sim_cloudYCtr)**2)
           elseif (NDIM == 3) then
              radius = sqrt( (xCoord(i)-sim_cloudXCtr)**2 &
                            +(yCoord(j)-sim_cloudYCtr)**2 &
                            +(zCoord(k)-sim_cloudZCtr)**2)
           endif

           taper=(sim_cloudRadius-radius)/(sim_cloudRadius-r0)

           if (radius < r0) then
              solnData(DENS_VAR,i,j,k)=sim_cloudDensity
           elseif ((radius > r0).and.(radius < sim_cloudRadius)) then
              solnData(DENS_VAR,i,j,k)=sim_dLeft + (sim_cloudDensity-sim_dLeft)*taper
           endif

           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))


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



#ifdef VECZ_VAR
           ! vector potential Az
           solnData(VECZ_VAR,i,j,k) = .25*(Az(i,j)+Az(i+1,j)+Az(i,j+1)+Az(i+1,j+1))
#endif

           if (NDIM >= 2) then
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= (Az(i,j+1)-Az(i,j))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=-(Az(i+1,j)-Az(i,j))/dx

                 solnData(DIVB_VAR,i,j,k)= &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy

              else
                 solnData(MAGX_VAR,i,j,k)= 0.5*(Az(i,j+1)-Az(i,j) + Az(i+1,j+1)-Az(i+1,j))/dy
                 solnData(MAGY_VAR,i,j,k)=-0.5*(Az(i+1,j)-Az(i,j) + Az(i+1,j+1)-Az(i,j+1))/dx
                 solnData(MAGZ_VAR,i,j,k)= 0.0
                 solnData(DIVB_VAR,i,j,k)= 0.
              endif
#else
              solnData(MAGX_VAR,i,j,k)= (Az(i,j+1)-Az(i,j))/dy
              solnData(MAGY_VAR,i,j,k)=-(Az(i+1,j)-Az(i,j))/dx
              solnData(MAGZ_VAR,i,j,k)= 0.0
              solnData(DIVB_VAR,i,j,k)= 0.
              solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                        solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
#endif
           endif

        enddo
     enddo
  enddo

#if NFACE_VARS > 0
  if ((NDIM >= 2).and.(sim_killdivb)) then
     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
              solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))
              solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                        solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           enddo
        enddo
     enddo
  endif
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
