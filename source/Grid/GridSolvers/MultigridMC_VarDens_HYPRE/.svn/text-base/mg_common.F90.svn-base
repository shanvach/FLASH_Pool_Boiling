!*******************************************************************************

!  Module:      mg_common()

!  Description: Global type and variable declarations for the multigrid solver.

module mg_common



implicit none
save

#include "constants.h"  ! for MDIM

!===============================================================================


!logical, save ::  writeflag_mg

! Coarse Solution level
integer :: Solvelevel

! Largest and smallest levels of refinement in the mesh.
integer :: mesh_lrefmax, mesh_lrefmin

! A place to save nodetypes as paramesh has them before we change
! them to perform various one level operations.
integer, allocatable, dimension(:) :: nodetype_save !(maxblocks_tr)

! A place to save child data before the prolongation messes with it
logical, allocatable, dimension(:) :: newchild_save !(maxblocks_tr)

! Ranges of interior indices for blocks.

integer :: ili, iui, jli, jui, kli, kui

! Ranges of exterior indices for blocks.

integer :: ile, iue, jle, jue, kle, kue

! Boundary conditions, flag to substract mean on the source.
integer, dimension(2*MDIM) :: gr_mgBndTypes
integer :: mg_bnd_cond

! Grid geometry.

integer :: mg_geometry

! Just doing a quadrant?

logical :: quadrant

!     Second order interpolation for the work array:
integer, save :: interp_work = 2
integer, save, allocatable, dimension(:) :: interp_mask_work_mg
integer, save, allocatable, dimension(:) :: interp_mask_work_save

integer, save :: gr_mgDiffOpDiscretize


!===============================================================================

end module mg_common
