subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

      use mph_interface, only: mph_advect, mph_evolve, mph_iblset, &
                               mph_ibadvect, mph_iblset_3D

      use Driver_Data, only: dr_nstep

#include "Flash.h"

      implicit none
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      integer, intent(in) :: mph_flag

#if NDIM == 2
      if(mph_flag == 1 .and. dr_nstep == 1) &
      call mph_iblset(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
#endif

#if NDIM == 3
      if(mph_flag == 1 .and. dr_nstep == 1) &
      call mph_iblset_3D(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
#endif

      if(dr_nstep > 1 .and. mph_flag == 1) then
        call mph_ibadvect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
        call mph_advect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
      end if

      call mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

end subroutine
