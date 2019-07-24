subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

      use mph_interface, only: mph_advect, mph_evolve

      use Driver_Data, only: dr_nstep

      implicit none
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      integer, intent(in) :: mph_flag

      if(dr_nstep > 1 .and. mph_flag == 1) then

      call mph_advect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      end if

      call mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

end subroutine
