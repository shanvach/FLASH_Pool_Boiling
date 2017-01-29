!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl/Grid_markRefineDerefine
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
  use tree, ONLY : newchild, refine, derefine, stay,lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells, &
                             Grid_getListOfBlocks,Grid_getBlkPtr,Grid_releaseBlkPtr,&
                             Grid_getBlkIndexLimits
  use Simulation_data
  use gr_interface, ONLY : gr_findMean

  use ImBound_data, only : ib_bflags

  implicit none

#include "constants.h"
!#include "MHD.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref


  logical :: gcMask(NUNK_VARS)

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: refval,meanRefVar

  integer :: blockID,j,k,lb

  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  integer :: lref,specsSize
  !! End of special refinement treatment ---------

  real :: del(MDIM)
 
  logical, save :: specref=.true. 

  call Grid_fillGuardCells(CENTER,ALLDIR)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

!!$  do l = 1,gr_numRefineVars
!!$     iref = gr_refine_var(l)
!!$     ref_cut = gr_refine_cutoff(l)
!!$     deref_cut = gr_derefine_cutoff(l)
!!$     ref_filter = gr_refine_filter(l)
!!$     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
!!$  end do

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  if (.not. specref) then ! First refine don't use vort magnitude. Initialize the grid with specified refinement.
  do l = 1,gr_numRefineVars

     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)

     !call gr_findMean(iref,2,.false.,meanRefVar)

     do lb = 1,blockCount

         blockID = blockList(lb)

         ! Get Blocks internal limits indexes:
         call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

         call Grid_getBlkPtr(blockID,solnData,CENTER)
         call Grid_getBlkPtr(blockID,facexData,FACEX)
         call Grid_getBlkPtr(blockID,faceyData,FACEY)

         call Grid_getDeltas(blockID,del)

         !write(*,*) 'OMGM=',OMGM_VAR,'ivar=',iref,ib_bflags(:,blockID)

         solnData(iref,NGUARD:NGUARD+NXB+1,NGUARD:NGUARD+NXB+1,1) = &
         (faceyData(VELC_FACE_VAR,NGUARD+1:NGUARD+NXB+2,NGUARD+1:NGUARD+NYB+2,1)- &
          faceyData(VELC_FACE_VAR,NGUARD:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB+2,1))/del(1) - &
         (facexData(VELC_FACE_VAR,NGUARD+1:NGUARD+NXB+2,NGUARD+1:NGUARD+NYB+2,1)-&
          facexData(VELC_FACE_VAR,NGUARD+1:NGUARD+NXB+2,NGUARD:NGUARD+NYB+1,1))/del(2)


         refval = 0.
         do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
             do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                refval = max(refval,abs(solnData(iref,i,j,k))) !-meanRefVar))

             enddo
          enddo
         enddo

         ! Refine Criteria, ImBound
         if( any(ib_bflags(1:NDIM,blockID) .gt. 0)) refine(blockID) = .TRUE.
     

         ! Refine criteria, by absolute value of variable:
         if (refval .ge. ref_cut) then
            refine(blockID) = .TRUE.
         endif
         ! Derefine criteria
         if ((refval .le. deref_cut) .and. &
             (refine(blockID) .neqv. .TRUE.)) then
            derefine(blockID) = .TRUE.
         endif


         call Grid_releaseBlkPtr(blockID,solnData,CENTER)
         call Grid_releaseBlkPtr(blockID,facexData,FACEX)
         call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


      enddo

  enddo
  endif


!!$#define SPECIAL_REFINEMENT 1

!!$#ifdef SPECIAL_REFINEMENT
!!$  !! Call for the specialized refinement
!!$  specsSize=7
!!$  !! Coordinate information --------------------------------------
!!$  !! define a range of coordinates of the rectangle in x-direction
!!$  specs(1) =  sim_xMin + 1./4.*(sim_xMax - sim_xMin) +.005
!!$  specs(2) =  sim_xMin + 3./4.*(sim_xMax - sim_xMin) -.005
!!$
!!$  !! define a range of coordinates of the rectangle in y-direction
!!$  specs(3) =  sim_yMin + 1./4.*(sim_yMax - sim_yMin) +.005
!!$  specs(4) =  sim_yMin + 3./4.*(sim_yMax - sim_yMin) -.005
!!$
!!$  !! define a range of coordinates of the rectangle in z-direction
!!$  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.005
!!$  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.005
!!$  !! End of coordinate information -------------------------------
!!$
!!$  !! Decide wheather or not we refine only blocks completely 
!!$  !! contained within the rectangle (specs(7) .NE. 0.0)
!!$  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
!!$  specs(7) = 0.0
!!$
!!$  !write(*,*) 'Specs=',specs(1:7)
!!$
!!$  !! Bring all qualifying blocks to this level of refinement
!!$  lref = lrefine_max
!!$
!!$  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif
!!$
!!$
!!$#ifdef SPECIAL_REFINEMENT ! Quadrant II
!!$  !! Call for the specialized refinement
!!$  specsSize=7
!!$  !! Coordinate information --------------------------------------
!!$  !! define a range of coordinates of the rectangle in x-direction
!!$  specs(1) =  sim_xMin +.005
!!$  specs(2) =  sim_xMin + 2./4.*(sim_xMax - sim_xMin) -.005
!!$
!!$  !! define a range of coordinates of the rectangle in y-direction
!!$  specs(3) =  sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
!!$  specs(4) =  sim_yMax -.005
!!$
!!$  !! define a range of coordinates of the rectangle in z-direction
!!$  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
!!$  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
!!$  !! End of coordinate information -------------------------------
!!$
!!$  !! Decide wheather or not we refine only blocks completely 
!!$  !! contained within the rectangle (specs(7) .NE. 0.0)
!!$  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
!!$  specs(7) = 0.0
!!$
!!$  !write(*,*) 'Specs=',specs(1:7)
!!$
!!$  !! Bring all qualifying blocks to this level of refinement
!!$  lref = lrefine_max
!!$
!!$  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif


!!$#ifdef SPECIAL_REFINEMENT
!  if (specref) then
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -3. + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =  10.0 - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -4. + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   4. - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !write(*,*) 'Specs=',specs(1:7)

  !! Bring all qualifying blocks to this level of refinement
  lref = 2 !lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif


!!$#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -2. + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =   5.0 - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -2. + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   2. - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !write(*,*) 'Specs=',specs(1:7)

  !! Bring all qualifying blocks to this level of refinement
  lref = 3 !lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif

!!$#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -1.5 + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =   2.5 - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -1.5 + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   1.5 - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !write(*,*) 'Specs=',specs(1:7)

  !! Bring all qualifying blocks to this level of refinement
  lref = 4 !lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif


!!$#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  -1.   + 0.005 ! sim_xMin + 0./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =   1.25 - 0.005 !sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  -1. + 0.005 !sim_yMin + 2./4.*(sim_yMax - sim_yMin) +.005
  specs(4) =   1. - 0.005 !sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !write(*,*) 'Specs=',specs(1:7)

  !! Bring all qualifying blocks to this level of refinement
  lref = 5 !lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
!!$#endif

  specref = .false.

!  endif







#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
