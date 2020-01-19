subroutine Heat_AD( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Grid_interface, only: Grid_getCellMetrics,    &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_fillGuardCells
   use Grid_data, only : gr_domainBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use Grid_interface, only: Grid_getCellCoords, Grid_solvePoisson, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"   

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld

   integer ::  blockID,lb
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData

   real, dimension(GRID_IHI_GC,3,blockCount) :: iMetrics
   real, dimension(GRID_JHI_GC,3,blockCount) :: jMetrics
   real, dimension(GRID_KHI_GC,3,blockCount) :: kMetrics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real, dimension(GRID_IHI_GC) :: xcell, xedge
   real, dimension(GRID_JHI_GC) :: ycell, yedge
   real, dimension(GRID_KHI_GC) :: zcell, zedge
   real :: pi, fact
   integer, dimension(6) :: bc_types 
   real, dimension(2,6)  :: bc_values = 0.
   integer :: i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,iMetrics(:,LEFT_EDGE,lb), GRID_IHI_GC)
     call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,iMetrics(:,CENTER,lb),    GRID_IHI_GC)
     call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,iMetrics(:,RIGHT_EDGE,lb),GRID_IHI_GC)

     call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,jMetrics(:,LEFT_EDGE,lb), GRID_JHI_GC)
     call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,jMetrics(:,CENTER,lb),    GRID_JHI_GC)
     call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,jMetrics(:,RIGHT_EDGE,lb),GRID_JHI_GC)

     call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,kMetrics(:,LEFT_EDGE,lb), GRID_KHI_GC)
     call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,kMetrics(:,CENTER,lb),    GRID_KHI_GC)
     call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,kMetrics(:,RIGHT_EDGE,lb),GRID_KHI_GC)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

  !gcMask = .FALSE.
  !gcMask(TEMP_VAR)=.TRUE.
  !call Grid_fillGuardCells(CENTER,ALLDIR,&
  !     maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   pi = 3.141592654
   fact = 1.0
   bc_types = reshape(gr_domainBC, (/ 2*MDIM /) ) 

   
   !do i=1,6
   !  write(*,*) i, bc_types(i)
   !end do

   do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     
     call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcell, GRID_IHI_GC)
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycell, GRID_JHI_GC)
     call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcell, GRID_KHI_GC)

     call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, .true., xedge, GRID_IHI_GC)
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, .true., yedge, GRID_JHI_GC)
     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, .true., zedge, GRID_KHI_GC)

     solnData(NSRC_VAR,:,:,:) = 0.0
     solnData(NFLD_VAR,:,:,:) = 0.0
     
         
     !write(*,*) "Neumann Source Term should be ", OUTFLOW 
     if (bc_types(1) == OUTFLOW) then
       !write(*,*) "Neumann Source Term"
       do k=1, blkLimitsGC(HIGH,KAXIS)
         do j=1, blkLimitsGC(HIGH,JAXIS)
           do i=1, blkLimitsGC(HIGH,IAXIS)
#if NDIM == 2
             solnData(NSRC_VAR,i,j,k) = -8.00 * pi**2 * cos(real(2.0 * pi * xcell(i))) * cos(real(2.0 * pi * ycell(j)))
#else
             solnData(NSRC_VAR,i,j,k) = -12.0 * pi**2 * cos(real(2.0 * pi * xcell(i))) * cos(real(2.0 * pi * ycell(j))) * cos(real(2.0 * pi * zcell(k)))
#endif
           enddo
         enddo
       enddo

     else

       do k=1, blkLimitsGC(HIGH,KAXIS)
         do j=1, blkLimitsGC(HIGH,JAXIS)
           do i=1, blkLimitsGC(HIGH,IAXIS)
#if NDIM == 2
             solnData(NSRC_VAR,i,j,k) = -8.00 * pi**2 * sin(real(2.0 * pi * xcell(i))) * sin(real(2.0 * pi * ycell(j)))
#else
             solnData(NSRC_VAR,i,j,k) = -12.0 * pi**2 * sin(real(2.0 * pi * xcell(i))) * sin(real(2.0 * pi * ycell(j))) * sin(real(2.0 * pi * zcell(k)))
#endif
           enddo
         enddo
       enddo

     endif

     call Grid_solvePoisson(NFLD_VAR, NSRC_VAR, bc_types, bc_values, fact)

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

  gcMask = .FALSE.
  gcMask(NFLD_VAR)=.TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

end subroutine Heat_AD
