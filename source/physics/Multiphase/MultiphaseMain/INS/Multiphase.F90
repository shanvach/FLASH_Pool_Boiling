subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      use mph_interface, only: mph_advect, mph_evolve

      implicit none
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld

      call mph_advect(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
      call mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

end subroutine
