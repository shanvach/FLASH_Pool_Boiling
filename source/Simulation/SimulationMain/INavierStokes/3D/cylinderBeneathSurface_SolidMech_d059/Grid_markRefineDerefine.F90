!!****if* source/Simulation/SimulationMain/INavierStokes/3D/Snorkel_mcHYPRE_VD_wFS/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var, gr_meshComm

  use tree, ONLY : newchild, refine, derefine, stay,lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells
  use Simulation_data


  implicit none

#include "constants.h"
!#include "MHD.h"
#include "Flash_mpi.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  integer :: ierr


  logical :: gcMask(NUNK_VARS)

  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  real, dimension(7) :: specsA
  real, dimension(7) :: specsB
  real, dimension(7) :: specsC
  real, dimension(7) :: specsD
  real, dimension(3) :: specs2
  real, dimension(6) :: specs4
  real, dimension(4) :: specs3
  integer :: lref,specsSize, specsSize2, specsSize3,specsSize4
  integer :: specsSizeA, specsSizeB, specsSizeC,specsSizeD
  !! End of special refinement treatment ---------

  call Grid_fillGuardCells(CENTER,ALLDIR)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

#define SPECIAL_REFINEMENT 1

!!#ifdef SPECIAL_REFINEMENT
!!  specsSizeC=3
!!  specsC(1) = real(OMGZ_VAR)
!!  specsC(2) =   3.0
!!  specsC(3) =  -3.0
!!
!!  !lref = lrefine_max-2
!!  !lref = lrefine_max-1
!!!  lref = 4!lrefine_max
!!  lref = 5!lrefine_max
!!
!!  call Grid_markRefineSpecialized_KPD (VORTICITY,specsSizeC,specsC,lref)
!!#endif

#ifdef SPECIAL_REFINEMENT
  specsSizeC=3
  specsC(1) = real(OMGZ_VAR)
  specsC(2) =   4.0
  specsC(3) =  -4.0

  !lref = lrefine_max-2
  !lref = lrefine_max-1
!  lref = 4!lrefine_max
  lref = 5!lrefine_max

  call Grid_markRefineSpecialized_KPD (VORTwLS,specsSizeC,specsC,lref)
#endif

     call MPI_BARRIER(gr_meshComm,ierr)

!!#ifdef SPECIAL_REFINEMENT
!!  specsSize2=3
!!  specs2(1) =  real(DFUN_VAR)
!!  specs2(2) =  -0.25
!!  specs2(3) =   0.25
!!  lref = 5!lrefine_max-2
!!  !lref = lrefine_max-1
!!  !lref = lrefine_max
!!  call Grid_markRefineSpecialized_KPD (THRESHOLD,specsSize2,specs2,lref)
!!#endif
!!
!!     call MPI_BARRIER(gr_meshComm,ierr)

#ifdef SPECIAL_REFINEMENT
  specsSizeA=7
  specsA(1) = 0.0
  specsA(2) =-1.09
  specsA(3) = 0.0
  specsA(4) = 0.49
  lref = lrefine_max
  call Grid_markRefineSpecialized (WITHRADIUS,specsSizeA,specsA,lref)
#endif

#ifdef SPECIAL_REFINEMENT
  specsSizeA=7
  specsA(1) = 0.0
  specsA(2) =-1.09
  specsA(3) = 0.0
  specsA(4) = 0.50
  lref = lrefine_max
  call Grid_markRefineSpecialized (WITHRADIUS,specsSizeA,specsA,lref)
#endif

#ifdef SPECIAL_REFINEMENT
  specsSizeA=7
  specsA(1) = 0.0
  specsA(2) =-1.09
  specsA(3) = 0.0
  specsA(4) = 0.51
  lref = lrefine_max
  call Grid_markRefineSpecialized (WITHRADIUS,specsSizeA,specsA,lref)
#endif


#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
