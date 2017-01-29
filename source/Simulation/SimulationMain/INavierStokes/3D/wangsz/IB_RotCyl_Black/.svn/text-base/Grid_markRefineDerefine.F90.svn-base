!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Grid_markRefineDerefine
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
                        gr_numRefineVars,gr_refine_var
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min, lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells, &
                             Grid_getListOfBlocks, Grid_getBlkIndexLimits,   &
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_data

  use ins_interface, only :  ins_velomg2center
  use ImBound_data,  only :  ib_BlockMarker_flag

  use Driver_data, ONLY: dr_nstep

  implicit none

#include "constants.h"


  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref

  logical :: gcMask(NUNK_VARS)

  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  integer :: lref,specsSize
  !! End of special refinement treatment ---------


  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,faxezData

  real :: refval,meanRefVar

  integer :: blockID,j,k,lb

  real :: del(MDIM)

  logical, save :: firstcall = .true.

  call Grid_fillGuardCells(CENTER,ALLDIR)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

! do l = 1,gr_numRefineVars
!     iref = gr_refine_var(l)
!     ref_cut = gr_refine_cutoff(l)
!     deref_cut = gr_derefine_cutoff(l)
!     ref_filter = gr_refine_filter(l)
!     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
! end do

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Average Velocities and Vorticity to cell-centers
  call ins_velomg2center(blocklist,blockcount)

  ! Now fill refine - derefine flags:
  do l = 1,gr_numRefineVars

     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)

     do lb = 1,blockCount

         blockID = blockList(lb)

         ! Get Blocks internal limits indexes:
         call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

         ! Get block Center data pointer:
         call Grid_getBlkPtr(blockID,solnData,CENTER)

         refval = 0.
         do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
             do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                refval = max(refval,abs(solnData(iref,i,j,k)))
             enddo
          enddo
         enddo

         ! Refine Criteria, ImBound
         if(ib_BlockMarker_flag(blockID)) refine(blockID) = .TRUE.

         ! Refine criteria, by absolute value of variable:
         if (refval .ge. ref_cut) then
            refine(blockID) = .TRUE.
         endif
         ! Derefine criteria
         if ((refval .le. deref_cut) .and. &
             (refine(blockID) .neqv. .TRUE.)) then
            derefine(blockID) = .TRUE.
         endif

         ! Release block pointer:
         call Grid_releaseBlkPtr(blockID,solnData,CENTER)

      enddo

  enddo

if (dr_nstep .gt. CONSTANT_ONE) firstcall= .false.

if (firstcall) then
#define SPECIAL_REFINEMENT 1

#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -2.5  + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =   2.5  - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -2.5  + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   2.5  - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  -3.0  + 0.005 !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =   3.0  - 0.005 !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !! Bring all qualifying blocks to this level of refinement
  lref = lrefine_max !lrefine_min+1 

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
#endif

endif



#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
