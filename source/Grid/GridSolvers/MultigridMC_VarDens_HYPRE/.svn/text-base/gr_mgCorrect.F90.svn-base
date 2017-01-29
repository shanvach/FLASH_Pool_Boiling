!*******************************************************************************

!  Routine:     mg_correct()

!  Description: Correct the current guess to the solution on a given level
!               using the residual on that level.

!  Parameters:  level       Level to compute the residual on.
!               isoln       Index of variable containing solution.  Receives
!                           corrected solution.
!               icorr       Index of variable containing correction.
!               leaf_only   If /= 0, only correct on leaf nodes on this level.


      subroutine gr_mgCorrect (level, isoln, icorr, leaf_only)

!===============================================================================

      use mg_common, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui

      use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                    Grid_getBlkPtr,       &
                                    Grid_releaseBlkPtr
      use tree, only : lrefine

      implicit none
#include "constants.h"
      integer, intent(in) :: level, isoln, icorr, leaf_only

      integer :: lnblocks
      integer :: lb, i, j, k
      logical :: update
      real, pointer, dimension(:,:,:,:) :: unk

!===============================================================================

      call Grid_getLocalNumBlks(lnblocks)


      !- kpd - unk(isoln,i,j,k) is the pressure solution from the previous
      !           cycle or time step.
      do lb = 1, lnblocks
        update = (lrefine(lb) == level)
        if (leaf_only /= 0) update = update .and. (nodetype_save(lb) == 1)
        if (update) then

          !print*,"Corrector",lb

          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unk,CENTER)
          !print*,"CORR",lb
          do k = kli, kui
            do j = jli, jui
              do i = ili, iui

!if (lb .eq. 58 .or. lb .eq. 57 ) then
!if (lb .eq. 8 .or. lb .eq. 9 ) then
!   print*,"CORR",isoln,lb,i,j,unk(isoln,i,j,k),unk(icorr,i,j,k),unk(isoln,i,j,k) +unk(icorr,i,j,k)
!end if

                unk(isoln,i,j,k) = unk(isoln,i,j,k) + & 
     &                             unk(icorr,i,j,k)

                !if (lb .eq. 58 .or. lb .eq. 57 ) then
                !if (lb .eq. 4 .or. lb .eq. 5 ) then
                !   print*,"CORR",lb,i,j,unk(isoln,i,j,k)
                !end if


              enddo
            enddo
          enddo
          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unk,CENTER)

        endif
      enddo

!===============================================================================

      return
      end
