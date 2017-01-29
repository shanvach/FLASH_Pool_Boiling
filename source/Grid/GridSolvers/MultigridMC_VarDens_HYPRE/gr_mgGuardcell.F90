!*******************************************************************************

! Routine:      mg_guardcell

! Description:  Version of amr_guardcell() for the multigrid solver.
!               The only differences are (1) it only operates on the
!               work array and (2) it calls the special routine
!               mg_set_ext_bndry() to handle exterior guard cells.

!               See documentation for amr_guardcell() for more information.


subroutine gr_mgGuardcell(mype2, ivar, nlayers, simtime, idiag, idir)

!===============================================================================

use mg_common, only :  gr_mgDiffOpDiscretize

use Driver_interface, ONLY : Driver_abortFlash
use paramesh_interfaces, only: amr_guardcell
use paramesh_dimensions, only : nguard_work

use workspace, only : interp_mask_work,interp_mask_work_res

use Grid_interface, ONLY : Grid_fillGuardCells

implicit none

!include 'mpif.h'
#include "Flash.h"
#include "constants.h"

integer, intent(in) :: mype2, nlayers, idiag, idir
integer, intent(in) :: ivar
real    :: simtime
!===============================================================================

if (nlayers > nguard_work) &
  call Driver_abortFlash("mg_guardcell: you must set nlayers <= nguard_work")


select case(gr_mgDiffOpDiscretize) 
case(2)
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
call Grid_fillGuardCells( WORK, ALLDIR, minLayers=1) 
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
case(4)
!interp_mask_work(:) = 4
!interp_mask_work_res(:) = 2
call Grid_fillGuardCells( WORK, ALLDIR, minLayers=3)
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
end select 

!===============================================================================

return
end

