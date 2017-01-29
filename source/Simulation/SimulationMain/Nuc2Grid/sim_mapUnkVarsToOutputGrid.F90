#include "constants.h"
#include "Flash.h"

subroutine sim_mapUnkVarsToOutputGrid
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
    Grid_getSingleCellVol, Grid_releaseBlkPtr
  use sim_outputGridData
  implicit none

  integer :: blockCount
  integer :: blockList(MAXBLOCKS), blockID
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer       :: ivar, lb, i,j,k
  real    :: dvol
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer :: point(MDIM), error, sizeX,sizeY,sizeZ
  real :: coords(MDIM)
  ! real :: xOgStep ! in module
  logical,parameter :: gcell=.TRUE.


  real,allocatable :: xCenter(:), yCoord(:), zCoord(:)
  real,allocatable :: xOgFaces(:),yOgFaces(:),zOgFaces(:)

  


  ! prepare output grid
  allocate(xOgFaces(nIOg+1))
  allocate(yOgFaces(nJOg+kOg2d))
  allocate(zOgFaces(nKOg+kOg3d))

  if (xOgMin .LE. 0.0) then     !linear spacing...

     xOgStep = (xOgMax - xOgMin) / nIOg
     do i=1,nIOg
        xOgFaces(i) = xOgMin + (i-1) * xOgStep 
     end do
     xOgFaces(nIOg+1) = xOgMax

  else                          !log spacing...

     xOgFact = 10.0**((alog10(xOgMax) - alog10(xOgMin)) / nIOg)
     do i=1,nIOg
        xOgFaces(i) = xOgMin * xOgFact**(real(i)-1.0)
     end do
     xOgFaces(nIOg+1) = xOgMax
  end if

  if (ndimOg>1) then
     yOgStep = (yOgMax - yOgMin) / nJOg
     do i=1,nJog
        yOgFaces(i) = yOgMin + (i-1) * yOgStep 
     end do
     yOgFaces(nJog+1) = yOgMax
  else
     yOgFaces(1) = 0.0
  end if

  if (ndimOg>2) then
     zOgStep = (zOgMax - zOgMin) / nKOg
     do i=1,nKog
        zOgFaces(i) = zOgMin + (i-1) * zOgStep 
     end do
     zOgFaces(nKog+1) = zOgMax
  else
     zOgFaces(1) = 0.0
  end if

  print*,'LBOUND(outGridData)=',LBOUND(outGridData)
  print*,'UBOUND(outGridData)=',UBOUND(outGridData)
  print*,'LBOUND(outGridNumCount)=',LBOUND(outGridNumCount)
  print*,'UBOUND(outGridNumCount)=',UBOUND(outGridNumCount)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do lb = 1, blockCount
     blockID = blockList(lb)
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
     sizeX = blkLimitsGC(HIGH,IAXIS)
     sizeY = blkLimitsGC(HIGH,JAXIS)
     sizeZ = blkLimitsGC(HIGH,KAXIS)
     allocate(xCenter(sizeX))
     allocate(yCoord(sizeY))
     allocate(zCoord(sizeZ))
     xCenter = 0.0
     yCoord = 0.0
     zCoord = 0.0
     if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
     if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

     call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(IAXIS) = i; coords(IAXIS) = xCenter(i)
              point(JAXIS) = j; coords(JAXIS) = yCoord(j)
              point(KAXIS) = k; coords(KAXIS) = zCoord(k)
              !! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
!!$              dvol = 1.0

              call throwOntoGrid(coords,dvol,solnData(:,i,j,k),outGridData,&
                   outGridNumCount,outGridMappedVol)
              
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)
     deallocate(xCenter)
     deallocate(yCoord)
     deallocate(zCoord)
  end do


contains
  subroutine throwOntoGrid(coords,dvol,solnVect,outGrid,outGridNumCount,outGridMappedVol)
    use ut_interpolationInterface

    real,intent(IN) :: coords(MDIM)
    real,intent(IN) :: dvol
    real,intent(IN) :: solnVect(:)
    real,intent(INOUT) :: outGrid(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0,*)
    integer,intent(INOUT) :: outGridNumCount(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0)
    real,intent(INOUT) :: outGridMappedVol(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0)

    integer,save :: iOut=0,jOut=0,kOut=0

    call ut_hunt(xOgFaces,nIOg+1,coords(IAXIS),iOut)
    if (ndimOg > 1) then
       call ut_hunt(yOgFaces,nJOg+1,coords(JAXIS),jOut)
    else
       jOut = 1
    end if
    if (ndimOg > 2) then
       call ut_hunt(zOgFaces,nKOg+1,coords(KAXIS),kOut)
    else
       kOut = 1
    end if

#ifdef DEBUG_ALL
888 format(3(1xG10.3),A,I4,I4,' dens',1PG15.4)
    print 888,coords,'->',iOut,jOut,solnVect(DENS_VAR)
#endif

    ivar = DENS_VAR
    outGrid(iOut,jOut,kOut,ivar) = outGrid(iOut,jOut,kOut,ivar) + solnVect(ivar)*dvol
    outGridNumCount(iOut,jOut,kOut) = outGridNumCount(iOut,jOut,kOut) + 1
    outGridMappedVol(iOut,jOut,kOut) = outGridMappedVol(iOut,jOut,kOut) + dvol
    

    do ivar = 1,NUNK_VARS
       if (ivar .NE. DENS_VAR .AND. &
            ivar .NE. NUMC_VAR) &
            outGrid(iOut,jOut,kOut,ivar) = outGrid(iOut,jOut,kOut,ivar) + solnVect(ivar)*dvol
    end do

  end subroutine throwOntoGrid
end subroutine sim_mapUnkVarsToOutputGrid
