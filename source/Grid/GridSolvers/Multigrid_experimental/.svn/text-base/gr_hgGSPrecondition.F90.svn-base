subroutine gr_hgGSPrecondition(iSrc, iSoln, iterations, iter_gcFill)

#include "constants.h"
#include "Flash.h"
#include "Multigrid.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use workspace, ONLY : work
  use Grid_data, ONLY: gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getDeltas, Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_fillGuardCells

  implicit none

  include "Flash_mpi.h"

  integer, intent(IN) :: iSrc, iSoln, iterations, iter_gcFill
  integer :: cur_iter, pass
  integer :: nBlocks, lb, b, i, j, k
  integer :: is, js, ks
  integer :: blockList(MAXBLOCKS)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  real, dimension(MDIM) :: deltas
  real :: c, cx, cy, cz
  real, pointer, dimension(:,:,:,:) :: solnData

  call Timers_start("gr_cgGSPrecondition")
  call Grid_getListOfBlocks(LEAF, blockList, nBlocks)
  do cur_iter = 1, iterations
!     write (*,*) "BEGINNING AN ITERATION"
!     if (mod(cur_iter, iter_gcFill) == 0) then
!        call Grid_fillGuardCells(gr_meshMe, CENTER, ALLDIR)
!        call gr_cgSetBoundaries(iSoln)
!     end if

        if (mod(cur_iter, iter_gcFill) == 1 .or. iter_gcFill == 1 .or. cur_iter == 1) then
           call gr_hgBndry(0, iSoln, 1, 0, MG_EXCHANGE_WORK, MG_BEGIN_SERIES, .false.)
        endif

     do pass = 1,2


        do lb = 1, nBlocks
           b = blockList(lb)
           call Grid_getDeltas(b, deltas)
           call Grid_getBlkPtr(b, solnData)
           call Grid_getBlkIndexLimits(b, blkLimits, blkLimitsGC)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                    if (mod(i + j + k, 2) == mod(pass,2)) then
                       if (NDIM == 1) then
                          cx = 0.5
                          c = 0.5 * deltas(IAXIS)**2
                          solnData(iSoln,i,j,k) = &
                               cx*(work(i-1,j,k,b,1)+work(i+1,j,k,b,1))&
                               -c*solnData(iSrc,i,j,k)
                          work(i,j,k,b,1) = solnData(iSoln,i,j,k)
                       end if
                       if (NDIM == 2) then
                          c = 0.5 &
                               / (deltas(IAXIS)**2 + deltas(JAXIS)**2)
                          cx = deltas(JAXIS)**2 * c
                          cy = deltas(IAXIS)**2 * c
                          c = c * deltas(IAXIS)**2 * deltas(JAXIS)**2
                          solnData(iSoln,i,j,k) = & 
                               cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                               cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1))&
                               - c*solnData(iSrc,i,j,k)
                          work(i,j,k,b,1) = solnData(iSoln,i,j,k)
                       end if
                       if (NDIM == 3) then
                          c = 0.5 / (deltas(JAXIS)**2 * deltas(KAXIS)**2 + &
                               deltas(IAXIS)**2 * deltas(KAXIS)**2 + &
                               deltas(IAXIS)**2 * deltas(JAXIS)**2)
                          cx = deltas(JAXIS)**2 * deltas(KAXIS)**2 * c
                          cy = deltas(IAXIS)**2 * deltas(KAXIS)**2 * c
                          cz = deltas(IAXIS)**2 * deltas(JAXIS)**2 * c
                          c = c * deltas(IAXIS)**2 * deltas(JAXIS)**2 * deltas(KAXIS)**2
                          solnData(iSoln,i,j,k) = -c*solnData(iSrc,i,j,k) + &
                               cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                               cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                               cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1))
                          work(i,j,k,b,1) = solnData(iSoln,i,j,k)
                       end if
                    end if
                 end do
              end do
           end do
           call Grid_releaseBlkPtr(b, solnData)
        end do
     end do
  end do
  call Timers_stop("gr_cgGSPrecondition")

end subroutine gr_hgGSPrecondition
