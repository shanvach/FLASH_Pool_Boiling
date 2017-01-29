!*******************************************************************************

!  Routine:     poisson()

!  Description: Driver routine for the BiPCGStab Poisson solver.  

!  Parameters:  isrc            Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               isoln           Index for solution array.  The solution is
!                               written directly into this variable.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: The following have been added to the interface only to get this
! version of the paramesh3.x poisson solver to work.  These have no effect. 
!               bc_types(6)     Boundary condition types array:
!                                 GRID_PDE_BND_PERIODIC
!                                 GRID_PDE_BND_DIRICHLET
!                                 GRID_PDE_BND_NEUMANN
!
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z  6 = +z
!               bc_values(2,6)  Values for dirichlet and neumann boundary
!                               conditions.  If Robins, then bc_values(1,*) holds
!                               dirichlet and bc_values(2,*) neumann conditions; if
!                               neumann or dirichlet, then bc_values(1,*) holds
!                               neumann or dirichlet values and bc_values(2,*)goes unused
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine Grid_SolvePoisson (isoln, isrc, bc_types, bc_values, poisfact)

!===============================================================================
#include "Flash.h"

use bicg_common

use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED, &
                              GRID_PDE_BND_PERIODIC, &
                              GRID_PDE_BND_DIRICHLET,&
                              GRID_PDE_BND_NEUMANN, &
                              Grid_getLocalNumBlks, &
                              Grid_getListOfBlocks, &
                              Grid_getBlkPtr, &
                              Grid_releaseBlkPtr

!use gr_isoInterface, ONLY: gr_isoImageMass, gr_isoImageBdry
use Timers_interface, ONLY : Timers_start, Timers_stop
use Driver_interface, ONLY : Driver_abortFlash


implicit none

integer :: isoln, isrc
integer, dimension(6) :: bc_types
real, dimension(2,6) :: bc_values
real    :: poisfact

integer       :: lb, lnblocks2, MyPE2, MasterPE2
integer       :: i, j, k

integer, save :: bcro,bcri,bcvi,bcpi,bcsi,bczi,bcyi


!===============================================================================

! Get key numbers from the database for the temporary variables we need.

bcro = BIRO_VAR
bcri = BIRI_VAR
bcvi = BIVI_VAR
bcpi = BIPI_VAR
bcsi = BISI_VAR
bczi = BIZI_VAR
bcyi = BIYI_VAR


!===============================================================================

! Call the multigrid Poisson solver for different types of boundary conditions.

#ifdef IMGP_VAR

! This is a hack to make this version of poisson.F90's interface
! consistent with 'blessed' version.


select case (bc_types(1))

!-------------------------------------------------------------------------------
  print*,'isolated'
  case (GRID_PDE_BND_ISOLATED) ! isolated boundary conditions

    call Driver_abortFlash("[poisson]  Isolated BCs are not supported on BPCGStab.")

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_PERIODIC) ! periodic boundary conditions

    call BiPCGStab (isrc,isoln,poisfact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bc_types,bc_values)

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_DIRICHLET) ! Dirichlet boundary conditions

    call BiPCGStab (isrc,isoln,poisfact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bc_types,bc_values)

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_NEUMANN) ! Neumann boundary conditions

    call BiPCGStab (isrc,isoln,poisfact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bc_types,bc_values)

!-------------------------------------------------------------------------------

  case default
    call Driver_abortFlash("[poisson]  invalid boundary condition type!")


!-------------------------------------------------------------------------------

end select


#else

    call Timers_start("BiPCGStab_solve")

    call BiPCGStab (isrc,isoln,poisfact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bc_types,bc_values)

    call Timers_stop("BiPCGStab_solve")
    
#endif


!===============================================================================

return
end subroutine Grid_SolvePoisson
