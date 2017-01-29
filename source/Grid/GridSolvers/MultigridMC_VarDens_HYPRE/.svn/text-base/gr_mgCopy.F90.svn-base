!*******************************************************************************

!  Routine:     mg_copy()

!  Description: Copy one multigrid array to another on a given level.

!  Parameters:  level       Level to copy on.
!               ifrom       Index of variable to copy from.
!               ito         Index of variable to copy to.
!               leaf_only   If /= 0, only copy leaf nodes on this level.


      subroutine gr_mgCopy (level, ifrom, ito, leaf_only)

!===============================================================================

      use mg_common, ONLY: nodetype_save, &
                           ili, iui, jli, jui, kli, kui

      use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                    Grid_getBlkPtr,       &
                                    Grid_releaseBlkPtr
      use tree, only : lrefine

      implicit none
#include "constants.h"
      integer, intent(in) :: level, ifrom, ito, leaf_only
      integer :: lnblocks
      integer :: lb
      logical :: update

      real, pointer, dimension(:,:,:,:) :: unk

!===============================================================================

      call Grid_getLocalNumBlks(lnblocks)

      do lb = 1, lnblocks
        update = (lrefine(lb) == level)
        if (leaf_only /= 0) update = update .and. (nodetype_save(lb) == 1)

        if (update) then

        ! Point to blocks center vars:
        call Grid_getBlkPtr(lb,unk,CENTER)

        !print*,"COPY",lb
        unk(ito,ili:iui,jli:jui,kli:kui) = &
                    unk(ifrom,ili:iui,jli:jui,kli:kui)

        !if (lb .eq. 4 .or. lb .eq. 5 ) then
        !   print*,"COPY",lb,unk(ito,14,6,1),unk(ito,5,6,1)
        !end if

        ! Release pointers:
        call Grid_releaseBlkPtr(lb,unk,CENTER)
         
        endif

      enddo

!===============================================================================

      return
      end
