!!****if* source/Simulation/SimulationMain/ShockCyl/Grid_markRefineDerefine
!!
!! NAME
!! 
!!  Grid_markRefineDerefine
!!
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine(integer(IN) :: MyPE)
!!
!!
!! DESCRIPTION
!!  This routine is specific to ShockCyl.  It overrides the default
!!  routine located in Grid/common/paramesh
!!
!!  This routine is the main interface to the refining criteria in FLASH.
!!  By default it will check the second derivative condition for the variables
!!  specified by refine_var_[1-4], and set the refine and derefine logical
!!  flags accordingly.  
!! 
!!  After the second derivatives are checked, it loops through all of the 
!!  blocks and makes sure that they will still have a refinement level in
!!  the range lrefine_min < block's refinement level < lrefine_max.
!!
!!  Any user defined refinement should be put in this routine.  For instance,
!!  to refine at circular region of the grid, you can loop over all the 
!!  blocks, compute the distance from the center of the circle, and set
!!  refine(blockID) = .TRUE. and derefine(blockID) = .FALSE. if the block
!!  falls within the circle.  Much of this functionality is provided by the
!!  MarkGridRef module, whose routines should be called from here.
!!
!!***

subroutine Grid_markRefineDerefine(MyPE)

  use Grid_data, ONLY : gr_blkList, gr_numRefineVars
  use Driver_interface, ONLY : Driver_getSimTime
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getLocalNumBlks
  use Grid_data, ONLY : gr_refine_cutoff,gr_derefine_cutoff,gr_refine_filter,&
       gr_numRefineVars, gr_refine_var

  use Tree, ONLY : derefine, refine, lrefine, nodetype, stay, lnblocks
  use Tree, ONLY : newchild, lrefine_min, lrefine_max


  use Simulation_data, ONLY: mach, ref_rect_x, ref_rect_y, xctr, yctr, ymin, ymax
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"


  integer, intent(IN) :: MyPE 

  integer :: idir=ALLDIR
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: mode = MODE_DENS_EI


  integer :: blockCount


  real          ::  rmid, time
!!  real, dimension(2*MDIM) :: xlb
  real, pointer, dimension(:,:,:,:) :: solnData

  integer       :: i, j, k, l, iref, lb
  real          :: d_, dj,  errh, errl, ref_cut, deref_cut, ref_filter
  integer       :: sf6_ind
  real          :: ddens, dpres, dsf6, dref

  ddens = 0.05
  dpres = 0.10
  dsf6  = 0.01

  call Grid_fillGuardCells( CENTER, idir)
  call Grid_getListOfBlocks(LEAF,gr_blkList,blockCount)


  do i = 1, blockCount
     call Grid_getBlkIndexLimits(gr_blkList(i),blkLimits,blkLimitsGC)
     call Eos_wrapped(mode,blkLimitsGC,gr_blkList(i))
  end do

  !!DEV: believe I have to make these calls here ...
  !!DEV: i think we need an accessor function for time

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter)
  end do

  do lb = 1, blockCount

     call Grid_getBlkPtr(gr_blkList(lb),solnData)

     call Grid_getBlkIndexLimits(gr_blkList(lb),blkLimits,blkLimitsGC)

     errh = 0.
     errl = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     refine on density
!!!!!!!!!!!!!!!!!!!!!!!!!

     iref = DENS_VAR
     dref = ddens
     dj   = 0.

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)


              d_ = abs(solnData(iref,i+1,j,k)-solnData(iref,i-1,j,k)) &
                   /min(solnData(iref,i+1,j,k),solnData(iref,i-1,j,k))
              dj = max(d_,dj)
#if N_DIM >= 2

              d_ = abs(solnData(iref,i,j+1,k)-solnData(iref,i,j-1,k)) &
                   /min(solnData(iref,i,j+1,k),solnData(iref,i,j-1,k))
              dj = max(d_,dj)
#endif
#if N_DIM == 3

              d_ = abs(solnData(iref,i,j,k+1)-solnData(iref,i,j,k-1)) &
                   /min(solnData(iref,i,j,k+1),solnData(iref,i,j,k-1))
              dj = max(d_,dj)
#endif

           end do
        end do
     end do
     errh = errh + max(0., dj-dref)
     errl = errl + max(0., dj-0.25*dref)


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     refine on pressure
!!!!!!!!!!!!!!!!!!!!!!!!!

     iref = PRES_VAR
     dref = dpres
     dj   = 0.

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)


              d_ = abs(solnData(iref,i+1,j,k)-solnData(iref,i-1,j,k)) &
                   /min(solnData(iref,i+1,j,k),solnData(iref,i-1,j,k))
              dj = max(d_,dj)
#if N_DIM >= 2

              d_ = abs(solnData(iref,i,j+1,k)-solnData(iref,i,j-1,k)) &
                   /min(solnData(iref,i,j+1,k),solnData(iref,i,j-1,k))
              dj = max(d_,dj)
#endif
#if N_DIM == 3

              d_ = abs(solnData(iref,i,j,k+1)-solnData(iref,i,j,k-1)) &
                   /min(solnData(iref,i,j,k+1),solnData(iref,i,j,k-1))
              dj = max(d_,dj)
#endif

           end do
        end do
     end do
     errh = errh + max(0., dj-dref)
     errl = errl + max(0., dj-0.25*dref)

!!!!!!!!!!!!!!!!!!!!!!!
!!!     check all jumps
!!!!!!!!!!!!!!!!!!!!!!!

     if ( errh.gt.0. ) then
        derefine(gr_blkList(lb)) = .FALSE.
        refine(gr_blkList(lb))   = .TRUE.
     else if ( errl.le.0. ) then
        derefine(gr_blkList(lb)) = .TRUE.
     end if

!!!!!!!!!!!!!!!!!!!!!
!!!     refine on SF6
!!!!!!!!!!!!!!!!!!!!!

     iref = SF6_SPEC
     dref = dsf6
     dj   = 0.

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              dj = max(dj, solnData(iref,i,j,k) )

           end do
        end do
     end do

     if ( dj.gt.dref ) then
        derefine(gr_blkList(lb)) = .FALSE.
        refine(gr_blkList(lb))   = .TRUE.
     end if

     call Grid_releaseBlkPtr(gr_blkList(lb),solnData)
  end do

  !------------------------------------------------------------------------------

  call Driver_getSimTime(time)

  if ( min(ref_rect_x,ref_rect_y) > 0.e0 ) then

     rmid = xctr + 5858.2e0*time + 2.66344e6*time**2

     if ( abs(mach-3.e0) < 0.01 ) rmid = max(xctr, 20.5e0 + 66666.66667e0*time)
     
!!$     xlb=0.0
!!$     xlb(1) = rmid - 0.5e0*ref_rect_x
!!$     xlb(2) = rmid + 0.5e0*ref_rect_x
!!$     xlb(3) = max(ymin, yctr - 0.5e0*ref_rect_y)
!!$     xlb(4) = min(ymax, yctr + 0.5e0*ref_rect_y)
!!$
!!$     call Grid_markRefineSpecialized(RECTANLE, xlb,7, lrefine_max)

  end if

  !------------------------------------------------------------------------------

  ! Make sure that the refinement level of the blocks is between lrefine_min
  ! and lrefine_max

  call Grid_getLocalNumBlks(lnblocks)
  do i = 1, lnblocks

     if ((lrefine(i) == lrefine_min) .and. (nodetype(i) == 1)) & 
          derefine(i) = .FALSE.

     if ((lrefine(i) < lrefine_min) .and. (nodetype(i) == 1)) then
        refine(i) = .TRUE.
        derefine(i) = .FALSE.
     end if

     if ((lrefine(i) == lrefine_max) .and. (nodetype(i) == 1)) & 
          refine(i) = .FALSE.

  end do

  return
end subroutine Grid_markRefineDerefine
