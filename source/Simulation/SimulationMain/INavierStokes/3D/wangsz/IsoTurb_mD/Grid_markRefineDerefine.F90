!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
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
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()


#ifdef FLASH_GRID_PARAMESH

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_eosModeNow
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min
!!$  use physicaldata, ONLY : force_consistency
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells
  use Particles_interface, only: Particles_sinkMarkRefineDerefine
  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask

  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  integer :: lref,specsSize
  !! End of special refinement treatment ---------


  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.
!!$  gcMask(NUNK_VARS+1:maskSize) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

!!$  force_consistency = .FALSE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
     gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! For PARAMESH2, call gr_markRefineDerefine here if it hasn't been called above.
  ! This is necessary to make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

!!#define SPECIAL_REFINEMENT 1

#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -12.5  + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =   12.5  - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -12.5  + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   12.5  - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  -1.0  + 0.005 !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =   1.0  - 0.005 !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !! Bring all qualifying blocks to this level of refinement
  lref = lrefine_min+1

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
#endif



  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  
  call Particles_sinkMarkRefineDerefine()
#endif
  
  return
end subroutine Grid_markRefineDerefine

