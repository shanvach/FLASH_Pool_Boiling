!*******************************************************************************

!  Routine:     mg_bndry()

!  Description: Update boundary zones on a given level for a given variable.

!  Parameters:  level       Level to update, or 0 for all levels.
!               ivar        Index of variable to update.
!               leaf_only   If nonzero, variable is only defined on leaf-node
!                           blocks; obtain boundary information by restriction
!                           from finer neighbors.
!               iopt        COPY_UNK_TO_WORK: copy unk to work, exchange work's
!                             guardcells
!                           UPDATE_UNK:  as above, then copy guardcells to unk
!                           EXCHANGE_WORK: just exchange work's guardcells
!                             without unk copy
!               call        BEGIN_SERIES: first call in a series
!                           CONTINUE_SERIES: intermediate call in a series
!                           END_SERIES: cleanup call (don't exchange)
!                           STANDALONE: do first call, intermediate, and
!                             cleanup stuff
!               extrap      If nonzero, use boundary extrapolation rather
!                           than homogeneous Dirichlet boundaries for the
!                           exterior of the mesh.

subroutine mg_bndry (level, ivar, nlayers, leaf_only, iopt, call, extrap)

!===============================================================================

use mg_common, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui, &
     ile, jle, kle, iue, jue, kue, &
     BEGIN_SERIES, CONTINUE_SERIES, END_SERIES, STANDALONE, &
     COPY_UNK_TO_WORK, EXCHANGE_WORK, UPDATE_UNK
     

use dBase, ONLY: dBaseGetDataPtrAllBlocks, dBasePropertyInteger,     &
                 dBasePropertyReal, dBaseRefinementLevel,            &
                 dBaseTreePtrNodeType, dBaseNeighborBlockList,       &
                 nguard_work, maxblocks_tr, ndim, k2d, k3d

use perfmon

!temporary until we can think of better way
use dBaseDeclarations, ONLY: work

implicit none
integer :: level, ivar, nlayers, leaf_only, iopt, call, extrap

integer, save :: myPE, numPEs
integer :: i, j, k, lb, lnblocks, lrefine
real    :: time

real,    pointer, dimension(:,:,:,:,:), save :: solnData
integer, pointer, dimension(:), save         :: nodetype

logical, save :: first_call=.true.

!===============================================================================

call timer_start("mg_bndry")

!call timer_start("in_mg_bndry_before_gc")

if (first_call) then 
   first_call = .false.
   solnData => dBaseGetDataPtrAllBlocks()
   myPE     = dBasePropertyInteger("MyProcessor")
   numPEs   = dBasePropertyInteger("NumberOfProcessors")
   nodetype =>dBaseTreePtrNodeType()
end if

time     = dBasePropertyReal("Time")
lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

if ((call == BEGIN_SERIES) .or. (call == CONTINUE_SERIES) .or. &
    (call == STANDALONE)) then
!   call timer_start("work copy")
   
   if ((iopt == COPY_UNK_TO_WORK) .or. (iopt == UPDATE_UNK)) then
      do lb = 1, lnblocks
         lrefine = dBaseRefinementLevel(lb)
         if (((level == 0) .and. (nodetype(lb) == 1)) .or. &
             (lrefine == level) .or. &
             (lrefine == level+1) .or. (lrefine == level-1)) then
            do k = kli, kui
               do j = jli, jui
                  do i = ili, iui
                     work(i,j,k,lb) = solnData(ivar,i,j,k,lb)
                  enddo
               enddo
            enddo
         endif
      enddo
   endif

!   call timer_stop("work copy")

endif

! Temporarily re-mark the blocks' node types.  The PARAMESH guard cell routine
! works by first restricting leaf-node data to parents, then trading boundary
! data at the same level of refinement, then finally interpolating boundary
! data for nodetype 2 blocks to their leaf-node children.  We alter the
! nodetype and neigh_type arrays (with the help of get_tree_nodetypes) in order
! to fool amr_guardcell into doing what we want.

! What we want is:  if the ivar variable is defined everywhere on a level (not
! just on leaf nodes), we do not want to clobber valid data in the restriction
! step.  So we treat all blocks on the chosen level as leaf nodes, all blocks
! above it as temporarily nonexistent, and blocks below it as leaf or parent,
! depending on what they normally are.  This will clobber valid data in
! 'parent' blocks below the chosen level.

if ((level /= 0) .and. ((call == BEGIN_SERIES) .or. (call == STANDALONE))) then
   if (leaf_only == 0) then
      
      do lb = 1, lnblocks
         nodetype(lb) = nodetype_save(lb)
         if (dBaseRefinementLevel(lb) > level)  nodetype(lb) = -1
         if (dBaseRefinementLevel(lb) == level) nodetype(lb) = 1
         if ((dBaseRefinementLevel(lb) == level-1) .and. (nodetype(lb) /= 1)) & 
              nodetype(lb) = 2
      enddo

!      call timer_start("get_tree_nodetypes")
      call get_tree_nodetypes (NumPEs, MyPE)
!      call timer_stop("get_tree_nodetypes")

   ! If the ivar variable is defined only on leaf nodes, then it's OK to clobber
   ! parent data.  We probably don't need to do anything here (since at most a
   ! jump of one level is supposed to be guaranteed), but we mark all blocks at
   ! the next higher level as leaf-node just to make sure they provide boundary
   ! data to our chosen level.
   
   else
   
      do lb = 1, lnblocks
         nodetype(lb) = nodetype_save(lb)
         if (dBaseRefinementLevel(lb) > level+1)  nodetype(lb) = -1
         if (dBaseRefinementLevel(lb) == level+1) nodetype(lb) = 1
         if ((dBaseRefinementLevel(lb) == level) .and. (nodetype(lb) /= 1)) &
              nodetype(lb) = 2
      enddo
!      call timer_start("get_tree_nodetypes")
      call get_tree_nodetypes (NumPEs, MyPE)
!      call timer_stop("get_tree_nodetypes")
      
   endif
end if

!call timer_stop("in_mg_bndry_before_gc")

! Now call the PARAMESH guard cell routine (on the work array).

!call timer_start("mg_guardcell_in_mg_bndry")

if (call == STANDALONE) then
   call mg_guardcell(MyPE, ivar, nlayers, time, 1, 0, extrap)
else if ((call == BEGIN_SERIES).or.(call == CONTINUE_SERIES)) then
   call mg_guardcell(MyPE, ivar, nlayers, time, 0, 0, extrap)
end if

if ((level /= 0) .and. ((call == STANDALONE) .or. (call == END_SERIES))) then
   call mg_restore_nodetypes()
endif

!call timer_stop("mg_guardcell_in_mg_bndry")

!   call timer_start("work copy")
   
if (iopt == UPDATE_UNK) then
  do lb = 1, lnblocks
    lrefine = dBaseRefinementLevel(lb)
    if (((level == 0) .and. (nodetype(lb) == 1)) .or. &
        (lrefine == level) .or. &
        (lrefine == level+1) .or. (lrefine == level-1)) then
      solnData(ivar,ile:ili-1,:,:,lb) = work(ile:ili-1,:,:,lb)
      solnData(ivar,iui+1:iue,:,:,lb) = work(iui+1:iue,:,:,lb)
      if (ndim >= 2) then
        solnData(ivar,:,jle:jli-1,:,lb) = work(:,jle:jli-1,:,lb)
        solnData(ivar,:,jui+1:jue,:,lb) = work(:,jui+1:jue,:,lb)
      endif
      if (ndim == 3) then
        solnData(ivar,:,:,kle:kli-1,lb) = work(:,:,kle:kli-1,lb)
        solnData(ivar,:,:,kui+1:kue,lb) = work(:,:,kui+1:kue,lb)
      endif
    endif
  enddo
endif

!   call timer_stop("work copy")
call timer_stop("mg_bndry")
!===============================================================================

return
end
