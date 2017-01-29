!!****if* source/Grid/GridMain/Samrai/Grid_updateRefinement
!!
!! NAME
!!
!!  Grid_updateRefinement
!!
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!  
!!  call Grid_updateRefinement()
!!
!!
!! DESCRIPTION
!!
!!  Apply the user-defined refinment critera to determine which blocks need 
!!  to be refined and derefined.  Once the blocks are marked, call 
!!  amr_refine_derefine to actually carry out the refinements.  During this
!!  stage, the blocks are redistributed across processors (if needed).  
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step uses 
!!  a user-defined prolongation routine (these can be picked in the Modules
!!  file).
!!
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent.
!!
!!
!! PARAMETERS 
!!
!!  nrefs      The number of steps between refinements
!!
!!
!!***

#define DEBUG_POSITIVITY 0

subroutine Grid_updateRefinement(myPe, nstep,time, gridChanged)

  use Grid_data, ONLY : gr_blkList, gr_convertToConsvdForMeshCalls,gr_nrefs
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_markRefineDerefine
  
  implicit none


#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(in) :: nstep,myPe
  real, intent(in) :: time
  logical,intent(out),OPTIONAL :: gridChanged

  integer :: i

  real, pointer, dimension(:,:,:,:) :: solnData

  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: mode = MODE_DENS_EI
  logical :: allBlocks
  integer :: count, level=0 !DEV: level just a dummy var passed in

!=============================================================================

! We only consider refinements every nrefs timesteps.

  if (mod(nstep, gr_nrefs) == 0) then

     call Timers_start("tree")

     blkLimitsGC(LOW,:)=1
     blkLimits(LOW,:)=NGUARD+1
     blkLimits(HIGH,IAXIS)=GRID_IHI
     blkLimits(HIGH,JAXIS)=GRID_JHI
     blkLimits(HIGH,KAXIS)=GRID_KHI
     
     blkLimitsGC(HIGH,IAXIS)=GRID_IHI_GC
     blkLimitsGC(HIGH,JAXIS)=GRID_JHI_GC
     blkLimitsGC(HIGH,KAXIS)=GRID_KHI_GC
     
     call Timers_start("guardcell tree")
     call Grid_fillGuardCells( CENTER, ALLDIR)
     call Timers_stop("guardcell tree")

     call Grid_getListOfBlocks(LEAF, gr_blkList,count)

     do i = 1, count
        call Grid_getBlkPtr(gr_blkList(i),solnData)
        call Eos_inPlacePtr(mode,blkLimitsGC,solnData)
        call Grid_releaseBlkPtr(gr_blkList(i),solnData)
     end do

     call Timers_start("markRefineDerefine")
     call Grid_markRefineDefine()
     call Timers_stop("markRefineDerefine")
     
! Now perform the indicated refinements and derefinements.  Once blocks are
! destroyed and new children are created, the redistribution is done.  During
! this step, the block interiors and a single perimeter of guardcells are moved.
    
     call Timers_start("amr_refine_derefine")
     call samrai_re-grid()
     call Timers_stop("amr_refine_derefine")
     if (present(gridChanged)) gridChanged = .TRUE. ! Maybe re-gridding happened.

! update the grid coordinates for the new mesh
     !!DEV: not sure that Samrai needs to do this.  Not actually sure
     !!why paramesh does either
     !!call Timers_start("updateData")
     !!call gr_updateData()
     !!call Timers_stop("updateData")

! If using conserved varaibles, convert primitive variables to conserved
! in blocks that are not new children. The prolonging will stuff the new
! children.
         
!!$     if (gr_convertToConsvdForMeshCalls) then
!!$        count = 0
!!$        allBlocks = .false.
!!$        do i = 1, lnblocks
!!$           if ( .not. newchild(i)) then
!!$              count = count+1
!!$              gr_blkList(count)=i
!!$           endif
!!$        enddo
!!$        call gr_primitiveToConserve(allBlocks,gr_blkList,count)
!!$     endif
!!$     
     
! Initialize the data in the newly created children by prolonging the data
! from the parent to the children.  Important: the prolongation stencil can
! only include a single coarse zone on either side of the parent of the new
! children, since only one guardcell is guaranteed to be valid.


!     call amr_prolong (MyPE, 1, NGUARD)

     if (gr_convertToConsvdForMeshCalls) then
        allBlocks=.true.
        !call gr_conserveToPrimitive(allBlocks,gr_blkList,lnblocks)
     endif

! Now fill all of the guardcells for the leafs and then parents of children

!!! Dev :: for paramesh 3, and maybe in general, this is not needed.
!!the guardcell fill is done at the beginnging of hydro sweep
!!DEV: guardcell fill may be happening unnecessarily here!
     call Timers_start("guardcell, needed?")
     call Grid_fillGuardCells( CENTER, ALLDIR)
     call Timers_stop("guardcell, needed?")


! Call the EOS to make sure the energy and pressure are consistent in the new 
! blocks in their interiors.
     call Timers_stop("tree")

     call Timers_start("eos")
     call Grid_getListOfBlocks(LEAF, gr_blkList,count)

     do i = 1, count
        call Grid_getBlkPtr(gr_blkList(i),solnData)
        call Eos_inPlacePtr(mode,blkLimitsGC,solnData)
        call Grid_releaseBlkPtr(gr_blkList(i),solnData)
     end do

     call Timers_stop("eos")

  else
     if (present(gridChanged)) gridChanged = .FALSE.
  endif
         
  return
end subroutine Grid_updateRefinement
