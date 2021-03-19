subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

      use mph_interface, only: mph_advect, mph_evolve

      use ib_lset_interface, only: ib_lset, ib_advect, ib_lset_3D

      use Driver_Data, only: dr_nstep

      use Multiphase_data, only: mph_meshMe

#include "Flash.h"

      implicit none
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      integer, intent(in) :: mph_flag

#if NDIM == 2
      if(mph_flag == 1 .and. dr_nstep == 1) then
        !call ib_lset(blockCount,blockList,dt)

      !else if (mph_flag == 1 .and. mod(dr_nstep,1000) == 0) then
      !  if(mph_meshMe .eq. 0) print *,"Doing IB Level Set Reconstruction" 
      !  call ib_lset(blockCount,blockList,dt)

      !else if(mph_flag == 1 .and. dr_nstep > 1) then
      !  call ib_advect(blockCount,blockList,dt)

      end if
#endif

#if NDIM == 3
      if(mph_flag == 1 .and. dr_nstep == 1) then
        !call ib_lset_3D(blockCount,blockList,dt)

      !else if (mph_flag == 1 .and. mod(dr_nstep,1000) == 0) then
      !  if(mph_meshMe .eq. 0) print *,"Doing IB Level Set Reconstruction" 
      !  call ib_lset_3D(blockCount,blockList,dt)

      !else if (mph_flag == 1 .and. dr_nstep > 1) then
      !  call ib_advect(blockCount,blockList,dt)

      end if
#endif

      if(dr_nstep > 1 .and. mph_flag == 1) &
      call mph_advect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      call mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

end subroutine
