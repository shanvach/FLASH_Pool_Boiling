subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

#include "Flash.h"
#include "constants.h"

      use mph_interface, only: mph_advect, mph_evolve
      use Multiphase_data, only: mph_thco1, mph_cp1,mph_rho1,mph_meshMe

      ! Used only with analytical equation of interface location

      use Grid_interface, ONLY : Grid_getDeltas,         &
                                 Grid_getBlkIndexLimits, &
                                 Grid_getCellCoords,     &
                                 Grid_getBlkPtr,         &
                                 Grid_releaseBlkPtr,     &
                                 Grid_getBlkBoundBox,    &
                                 Grid_getBlkCenterCoords
      ! End of Case

      use Driver_Data, ONLY: dr_nstep

      implicit none

      include "Flash_mpi.h"

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      integer, intent(in) :: mph_flag

      ! Used only with analytical equation of interface location

      !integer :: blockID,lb,i,j,k
      !real :: del(MDIM),bsize(MDIM),coord(MDIM),xcell,ycell
      !real, dimension(2,MDIM) :: boundBox
      !real, pointer, dimension(:,:,:,:) :: solnData
      !integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

      !real :: solnX

      !if(dr_nstep > 1 .and. mph_flag == 1) then

      !call mph_advect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      !end if

      call mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

end subroutine
