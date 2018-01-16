module Plasma_interface

        implicit none

        interface
        subroutine Plasma_init(blockCount,blockList,restart)
                implicit none
                integer, INTENT(INOUT) :: blockCount
                integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
                logical, INTENT(IN)    :: restart
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


        interface
        subroutine Plasma_Solve(T_p, T_o, dt, dx, dy, ix1,ix2, jy1, jy2)
                implicit none
                real, dimension(:,:,:), intent(inout) :: T_p
                real, dimension(:,:,:), intent(in) :: T_o
                real, intent(in) :: dt, dx, dy
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine
        end interface

        interface
        subroutine Plasma_get_reactions(Te,Th)
                implicit none
                real, intent(in) :: Te, Th
        end subroutine Plasma_get_reactions
        end interface

end module Plasma_interface
