module Plasma_interface

        implicit none

        interface
        subroutine Plasma_init(blockCount,blockList,restart)
                implicit none
                integer, INTENT(INOUT) :: blockCount
                integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
                logical, intent(in)    :: restart
        end subroutine Plasma_init
        end interface

        interface
        subroutine Plasma_finalize()
                implicit none
        end subroutine
        end interface

        interface
        subroutine Plasma(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
                implicit none
                integer, INTENT(INOUT) :: blockCount
                integer, INTENT(IN) :: sweepOrder
                integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
                real,    INTENT(IN) :: timeEndAdv, dt, dtOld
        end subroutine Plasma
        end interface


end module Plasma_interface
