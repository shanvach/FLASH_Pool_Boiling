!*******************************************************************************

!  Module:      mg_common()

!  Description: Global type and variable declarations for the multigrid solver.

module mg_common

use  dBase, ONLY: maxblocks_tr
implicit none
save

!===============================================================================

! Constants for mg_bndry calls.
 
integer :: BEGIN_SERIES     = 1
integer :: CONTINUE_SERIES  = 2
integer :: END_SERIES       = 3
integer :: STANDALONE       = 4

integer :: COPY_UNK_TO_WORK = 1
integer :: EXCHANGE_WORK    = 2
integer :: UPDATE_UNK       = 3

! Largest and smallest levels of refinement in the mesh.

integer :: mesh_lrefmax, mesh_lrefmin

! Average value of source term (subtracted off when doing periodic boundaries).

real    :: src_avg

! A place to save nodetypes as PARAMESH has them before we change
! them to perform various single-level operations.

integer :: nodetype_save(maxblocks_tr)

! A place to save child data before the prolongation messes with it.

logical :: newchild_save(maxblocks_tr)

! Ranges of interior indices for blocks.

integer :: ili, iui, jli, jui, kli, kui

! Ranges of exterior indices for blocks.

integer :: ile, iue, jle, jue, kle, kue

! Boundary conditions.

integer :: mg_bnd_cond

! Saved variable indices.

integer :: mg_soln_index

! Supported geometry constants.

integer, parameter :: MG_GEOM_1DCARTESIAN = 1, MG_GEOM_1DCYLINDRICAL = 2, &
                      MG_GEOM_1DSPHERICAL = 3, MG_GEOM_2DCARTESIAN   = 4, &
                      MG_GEOM_2DCYLAXISYM = 5, MG_GEOM_2DCYLPOLAR    = 6, &
                      MG_GEOM_3DCARTESIAN = 7, MG_GEOM_3DCYLINDRICAL = 8, &
                      MG_GEOM_3DSPHERICAL = 9, MG_GEOM_INVALID       = 0, &
                      MG_GEOM_2DSPHAXISYM = 10

! Grid geometry.

integer :: mg_geometry

! Supported boundary constants.

integer, parameter :: MG_BND_PERIODIC  = 1, MG_BND_DIRICHLET = 2, &
                      MG_BND_NEUMANN   = 3, MG_BND_GIVENVAL  = 4, &
                      MG_BND_GIVENGRAD = 5, MG_BND_ISOLATED  = 0

! Just doing a quadrant?

logical :: quadrant

! Operations that work on leaf nodes only, parent nodes only, or all nodes.

logical :: MG_NODES_ALL_NODES   = 0
logical :: MG_NODES_LEAF_ONLY   = 1
logical :: MG_NODES_PARENT_ONLY = 2

!===============================================================================

end module mg_common
