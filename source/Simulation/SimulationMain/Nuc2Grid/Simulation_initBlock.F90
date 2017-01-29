!!****if* source/Simulation/SimulationMain/Nuc2Grid/Simulation_initBlock
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
!!  Initializes the variable for convolution kernel in UNK.
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  convoSmearWidI  width of smearing function in I direction
!!  convoSmearWidJ  width of smearing function in J direction
!!  convoSmearWidK  width of smearing function in K direction
!!  convoSmearShapeI  type smearing function in I direction
!!  convoSmearShapeJ  type smearing function in J direction
!!  convoSmearShapeK  type smearing function in K direction
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_convoSmearWidI, sim_convoSmearWidJ, sim_convoSmearWidK, &
        sim_convoSmearShapeI, sim_convoSmearShapeJ, sim_convoSmearShapeK,sim_meshMe
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, Grid_getGlobalIndexLimits, &
    Grid_getBlkCornerID, Grid_getDeltas
  use Grid_data, ONLY : gr_imin, gr_jmin, gr_kmin,&
                        gr_imax, gr_jmax, gr_kmax

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  
  integer :: globalIndexLimits(MDIM), cornerID(MDIM), stride(MDIM)
  real :: del(MDIM)
  real,save :: taccu

  real :: xx, yy,  zz, xxL, xxR
  
  real :: widI,widJ,widK
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,&
       yCoord,yLeft,yRight,&
       zCoord,zLeft,zRight,&
       xSkern,ySkern,zSkern

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: skernZone, accu
  
  real :: wfac
  real :: wfacI, wfacJ, wfacK
  logical :: gcell = .true.

  
  ! dump some output to stdout listing the paramters
  if (sim_meshMe == MASTER_PE) then
     
     
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
     
  endif
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  allocate(xSkern(sizeX))
  allocate(ySkern(sizeY))
  allocate(zSkern(sizeZ))
  xCenter = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

  call Grid_getGlobalIndexLimits(globalIndexLimits)
  call Grid_getBlkCornerID(blockId, cornerID, stride)
  call Grid_getDeltas(blockId, del)

  wfac = 2*sqrt(-alog(0.5))
  if (sim_convoSmearWidI < 0) then
     widI = - sim_convoSmearWidI * del(IAXIS) 
  else
     widI = sim_convoSmearWidI
  end if
  select case (sim_convoSmearShapeI)
  case(2,3)
     wfacI = 1.0 / widI
  case(1)
     wfacI = wfac / widI
  end select
  if (NDIM > 1) then
     if (sim_convoSmearWidJ < 0) then
        widJ = - sim_convoSmearWidJ * del(JAXIS) 
     else
        widJ = sim_convoSmearWidJ
     end if
     select case (sim_convoSmearShapeJ)
     case(2,3)
        wfacJ = 1.0 / widJ
     case(1) 
        wfacJ = wfac / widJ
     end select
  end if
  if (NDIM > 2) then
     if (sim_convoSmearWidK < 0) then
        widK = - sim_convoSmearWidK * del(KAXIS) 
     else
        widK = sim_convoSmearWidK
     end if
     select case (sim_convoSmearShapeK)
     case(2,3)
        wfacK = 1.0 / widK
     case(1) 
        wfacK = wfac / widK
     end select
  end if


  accu = 0.0
  do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
     ! get the cell center position in x
     xx  = xCenter(i) - gr_imin 
     if (gr_imax-xCenter(i) < xx) xx  = xCenter(i) - gr_imax
     xx = xx - 0.5 * del(IAXIS)
     xx = xx * wfacI
     if (sim_convoSmearShapeI==1) then
        xSkern(i) = exp(-xx**2)
     elseif (sim_convoSmearShapeI==2) then
        xSkern(i) = exp(-abs(xx))
     elseif (sim_convoSmearShapeI==3) then
        xSkern(i) = xx**2
     end if
     if (sim_convoSmearShapeI .NE. 3) then
        xx = abs(xx) * 40 / ( (globalIndexLimits(IAXIS))*0.5 )
        xSkern(i) = max(xSkern(i),1.0e-8*exp(-xx))
     end if
     accu = accu + xSkern(i)
  end do
  accu = 1.0/accu
!!$  do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!!$     xSkern(i) = xSkern(i) * accu !Norm so that sum is 1.
!!$  end do
!!$  print*,'wfac,widI,wfacI:',wfac,widI,wfacI
!!$  print*,'xSkern:',xSkern

  if (NDIM > 1) then
     accu = 1.0
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j) - gr_jmin
        if (gr_jmax-yCoord(j) < yy) yy  = yCoord(j) - gr_jmax
        yy = yy - 0.5 * del(JAXIS)
        yy = yy * wfacJ
        if (sim_convoSmearShapeJ==1) then
           ySkern(j) = exp(-yy**2) 
        elseif (sim_convoSmearShapeJ==2) then
           ySkern(j) = exp(-abs(yy)) 
        elseif (sim_convoSmearShapeJ==3) then
           ySkern(j) = yy**2
        end if
        if (sim_convoSmearShapeJ .NE. 3) then
           yy = abs(yy) * 40 / ( (globalIndexLimits(JAXIS))*0.5 )
           ySkern(j) = max(ySkern(j),1.0e-8*exp(-yy))
        end if
        accu = accu + ySkern(j)
     end do
     accu = 1.0/accu
!!$     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!!$        ySkern(j) = ySkern(j) * accu !Norm so that sum is 1.
!!$     end do
  else
     select case (sim_convoSmearShapeJ)
     case(3)
        ySkern(1) = 0.0
     case default
        ySkern(1) = 1.0
     end select
  end if
!!$  print*,'ySkern:',ySkern

  if (NDIM > 2) then
     accu = 1.0
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        ! get the coordinates of the cell center in the z-direction
        zz = zCoord(k) - gr_kmin
        if (gr_kmax-zCoord(k) < zz) zz  = zCoord(k) - gr_kmax
        zz = zz - 0.5 * del(KAXIS)
        zz = zz * wfacK
        if (sim_convoSmearShapeK==1) then
           zSkern(k) = exp(-zz**2)
        elseif (sim_convoSmearShapeK==2) then
           zSkern(k) = exp(-abs(zz))
        elseif (sim_convoSmearShapeK==3) then
           zSkern(k) = zz**2
        end if
        if (sim_convoSmearShapeK .NE. 3) then
           zz = abs(zz) * 40 / ( (globalIndexLimits(KAXIS))*0.5 )
           zSkern(k) = max(zSkern(k),1.0e-8*exp(-zz)) 
        end if
        accu = accu + zSkern(k)
     end do
     accu = 1.0/accu
!!$     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!$        zSkern(k) = zSkern(k) * accu !Norm so that sum is 1.
!!$     end do
  else
     select case (sim_convoSmearShapeK)
     case(3)
        zSkern(1) = 0.0
     case default
        zSkern(1) = 1.0
     end select
  end if
!!$  print*,'zSkern:',zSkern

     !------------------------------------------------------------------------------

     ! Loop over cells in the block.  For each, compute the physical position of 
     ! its left and right edge and its center as well as its physical width.  
     ! Then decide which side of the initial discontinuity it is on and initialize 
     ! the hydro variables appropriately.


  if (sim_convoSmearShapeI==3) then
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        axis(KAXIS) = k
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           axis(JAXIS) = j
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              axis(IAXIS) = i
              
!!$              xx = 160.0 * alog(10.0) / PI
!!$              skernZone = exp( - sin(  &
!!$                           min(0.5*PI, sqrt(xSkern(i)+ySkern(j)+zSkern(k))**1.2/xx) &
!!$                                     ) * xx)
              xx = 16 * alog(10.0)
              skernZone = exp( - fff(  &
                           sqrt(xSkern(i)+ySkern(j)+zSkern(k))**1.5/xx &
                                     ) * xx)
#if(0)
              zz = abs(sqrt(xSkern(i)+ySkern(j)+zSkern(k))) * 40 / ( maxval(globalIndexLimits(IAXIS:NDIM))*0.5 )
              skernZone = max(skernZone,1.0e-8*exp(-zz)) 
#endif
              skernZone = max(skernZone, 1e-16)
              ! store the variables in the current zone via Grid put methods
              ! data is put stored one cell at a time with these calls to Grid_putData           
#ifdef GAUS_VAR
              call Grid_putPointData(blockId, CENTER, GAUS_VAR, EXTERIOR, axis, skernZone)   
#endif


           enddo
        enddo
     enddo
  else                          !handle directions as separate factors
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        axis(KAXIS) = k
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           axis(JAXIS) = j
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              axis(IAXIS) = i

              skernZone = xSkern(i) * ySkern(j) * zSkern(k)
              skernZone = max(skernZone, 1e-24)
              ! store the variables in the current zone via Grid put methods
              ! data is put stored one cell at a time with these calls to Grid_putData           
#ifdef GAUS_VAR
              call Grid_putPointData(blockId, CENTER, GAUS_VAR, EXTERIOR, axis, skernZone)   
#endif


           enddo
        enddo
     enddo
  end if
!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
contains
  real function fff(x)
    real,intent(IN) :: x

    fff = 1.0 - exp(-x)
  end function fff

end subroutine Simulation_initBlock










