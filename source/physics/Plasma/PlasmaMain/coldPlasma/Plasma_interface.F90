module Plasma_interface

        implicit none

#include "constants.h"
#include "Flash.h"

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
        subroutine Plasma_Solve(T_p, T_o, dfun, dt, dx, dy, ix1,ix2, jy1, jy2,T_res)
                implicit none
                real, dimension(:,:,:), intent(inout) :: T_p
                real, dimension(:,:,:), intent(in) :: T_o, dfun
                real, intent(in) :: dt, dx, dy
                integer, intent(in) :: ix1, ix2, jy1, jy2
                real, intent(out) :: T_res
        end subroutine
        end interface

        interface
        subroutine Plasma_get_reactions(Te,Th)
                implicit none
                real, intent(in) :: Te, Th
        end subroutine Plasma_get_reactions
        end interface

        interface
        subroutine Plasma_computeDt(pls_mindt,pls_minloc)
                implicit none
                real, intent(INOUT) :: pls_mindt
                integer, intent(INOUT) :: pls_minloc(5)
        end subroutine Plasma_computeDt
        end interface

        interface
        subroutine Plasma_computeDtLocal(blockID,   &
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc )
                implicit none
                integer, intent(IN) :: blockID
                integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
                integer, intent(IN) :: isize,jsize,ksize
                real, intent(IN) :: dx, dy, dz
                real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
                real, intent(INOUT) :: dtLocal
                integer, intent(INOUT) :: lminloc(5)

        end subroutine Plasma_computeDtLocal
        end interface
        
end module Plasma_interface
