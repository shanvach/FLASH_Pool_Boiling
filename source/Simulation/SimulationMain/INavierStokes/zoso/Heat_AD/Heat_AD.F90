subroutine Heat_AD( blockCount,  blockList, timeEndAdv, dt,dtOld, sweepOrder)
   implicit none
   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld
end subroutine Heat_AD
