!*******************************************************************************

!  Routine:     mg_cycle()

!  Description: Perform one multigrid cycle beginning at a specified mesh
!               level.  This version implements a V-cycle for adaptively
!               refined meshes (Martin, D. and Cartwright, K.  "Solving
!               Poisson's Equation using Adaptive Mesh Refinement," 1996).


  subroutine gr_mgCycleHYPRE (levelmax, levelmin, img_soln, img_src, & 
     &                       img_res, img_corr, img_temp, img_temp2, & 
     &                       mg_solve, mg_residual, mg_residualMG, mg_relax_RBGS, &
                             mg_relax_HYPRE, &
                             img_denx,mg_relax)

!===============================================================================

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs


  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use mg_common, only : solvelevel, nodetype_save

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getListOfBlocks,   &
                             Grid_getLocalNumBlks

  use tree, only : lrefine,nodetype,grid_changed

  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data,      ONLY : dr_nstep 


  implicit none

#include "Flash.h"
#include "constants.h"


  integer :: levelmax, levelmin, img_soln, img_src,      & 
                img_res, img_corr, img_temp, img_temp2, &
                img_denx

  external mg_solve, mg_residual, mg_residualMG, mg_relax_RBGS, &
           mg_relax_HYPRE,mg_relax

  integer :: i,j 

  real, pointer, dimension(:,:,:,:,:) :: unk
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData


  integer, save :: mgrid_npresmooth,mgrid_npossmooth
  logical, save :: first_call = .true.

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer lb,blockcount,NblockX,ii,jj,kk,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)
  real, pointer, dimension(:,:,:,:) :: solnData

  integer lnblocks,pDim

  character*70        :: internalFile
  integer             :: comm, globalTopLevel
  logical             :: requestMap
  logical             :: suppressPfft
  integer             :: gridChanged
  integer,save        :: prevSolveLevel = -1

  real :: xcell, ycell, t_start, t_stop
  real del(MDIM)
  integer :: iii


#include "Flash_mpi.h"

  !==========================================================================

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  !==========================================================================
  !==========================================================================
  !- kpd - Test loop to check block communication
  !==========================================================================
  !do iii=1,gr_meshNumProcs
  do iii=1,1
  !print*,"========================"
  do lb = 1, blockCount                                                        
      blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkCenterCoords(blockId,coord)
     call Grid_getBlkBoundBox(blockId,boundBox)
     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

           !if (blockID .eq. 23 .AND. iii .eq. 1 ) then
           !   print*,"Block3",blockID,i,j,"XY",xcell,ycell, &
           !            facexData(RH1F_FACE_VAR,i,j,1)+facexData(RH2F_FACE_VAR,i,j,1), &
           !            faceyData(RH1F_FACE_VAR,i,j,1)+faceyData(RH2F_FACE_VAR,i,j,1)
           !end if
        end do
     end do
     !print*,"KPD Inside mgCycle",gr_meshNumProcs,blockCount,iii,blockID,xcell,ycell
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo
  !print*,"========================"
  end do

  !==========================================================================
  !==========================================================================

  if (first_call) then
     call RuntimeParameters_get('mgrid_npresmooth',    mgrid_npresmooth)
     call RuntimeParameters_get('mgrid_npossmooth',    mgrid_npossmooth)
     first_call = .false.
  end if

  if (gr_meshMe .eq. 0) then
     print*,"KPD - Entering gr_mgCycleHYPRE V-Cycles with Min/Max no. of levels: ",levelmin,levelmax
  end if

  !------------------------------------------------------------------------
  !- kpd - Initialize (set = 0,zero out) multigrid variable on a particular 
  !           mesh level (here the residual correction, e)
  !------------------------------------------------------------------------
  !- kpd - zero out img_corr on 5,6,7,8,10,11,12,... (includes non-leafs)
  call gr_mgZero (levelmax, img_corr, 0)

  do i = levelmax, levelmin
           call mg_restore_nodetypes (i)
  end do
  !- kpd - Assign mixture density to all leaf blocks
  call gr_mgAssignFden (levelmax, img_denx, 0)

      !-----------------------------------------------------------------|
      !- kpd - Work your way down the V-cycle...Fine -> coarse leg of V.|
      !-----------------------------------------------------------------|
      !- kpd - Only looped for more than 1 difference in levels         |
      !-----------------------------------------------------------------|
      if (levelmax > levelmin) then

        do i = levelmax, levelmin+1,  -1

           !print*,"Restrict at Level",i,"with max/min levels:",levelmax,levelmin
           !print*,"Restrict at Level",i,"with proc:",gr_meshMe

           !--------------------------------------------------------
           !- kpd - Copy img_soln to img_temp on the current level
           !--------------------------------------------------------
           !- kpd - img_temp receives img_soln on 5,6,7,8,10,11,12,... (leaf level 4's)
           !        Works on leaf blocks at "i" level
           call gr_mgCopy (i, img_soln, img_temp, 1)

           !-----------------------------------------------------------------
           !- kpd - Initialize (zero out) the img_corr array on the i-1 level
           !-----------------------------------------------------------------
           !- kpd - zero out img_corr on 3,4,9,14,20,... (includes non-leafs)
           call gr_mgZero (i-1, img_corr, 0)

           !--------------------------------------------------
           !- kpd - Solve the Poisson Eqn on the current level
           !        ==========================================
           !        A(h) * u(h)     = f(h)
           !        A    * img_corr = img_res
           !--------------------------------------------------
           !- kpd - img_corr receives update on 5,6,7,8,10,11,12,... (leaf level 4's)
           !                          RHS      LHS
           !print *,'Before mg_relax_RBGS ',gr_meshMe,i
           call mg_relax_RBGS   (i, img_res, img_corr, mgrid_npresmooth,img_denx,levelmax)
           !print*,"Using mg_relax instead of mg_relax_RBGS in mgCycleHYPRE"
           !call mg_relax (i, img_res, img_corr, mgrid_npresmooth)
           !print *,'After mg_relax_RBGS ',gr_meshMe,i

           !------------------------------------------------------------------
           !- kpd - Correct the current value of the solution on a given level
           !           using the residual on that level.
           !        Soln = Soln + Corr
           !        r(h) = f(h) - A(h)*u(h)
           !------------------------------------------------------------------
           !- kpd - img_soln (DELP) receives updated solution on 5,6,7,8,10,11,12,... (leaf level 4's)
           call gr_mgCorrect (i, img_soln, img_corr, 1)

           !!!if (gr_meshMe == 0) print *,' DONE mg_correct '
           !print *,' DONE mg_correct ',gr_meshMe,i

           !print*,"KPD - Entering mg_residual FIRST in gr_mgCycle.F90 with level # ",i
           !------------------------------------------------------------------------
           !- kpd - Compute the residual of the Poisson eqn to be solved using
           !           multigrid (The residual of the residual correction equation).
           !       img_res=RHS, img_corr=LHS, img_temp2=index of residual variable
           !------------------------------------------------------------------------
           !- kpd - *** img_temp2 receives residual on leafs 
           !                        RHS      LHS       RES
           call mg_restore_nodetypes (i)
           call mg_residualMG (i, img_res, img_corr, img_temp2, 0, 0, img_denx, levelmax)
           !print*,"Using mg_residual instead of mg_residualMG"
           !call mg_residual (i, img_res, img_corr, img_temp2, 0, 0)
           !print*,"KPD - Leaving mg_residual FIRST in gr_mgCycle.F90 with level # ",i

           !!!if (Mype == 0) print *,' DONE mg_residual '

           !----------------------------------------------------------------------------------
           !- kpd - Restrict the residual (RHS of residual correction eqn) from 
           !           one mesh level to the next coarser level. 
           !        Restrict from: img_temp2 @ fine --> img_res @ coarse for coarse NON_LEAF's
           !        img_temp2 is previously only defined for LEAFS and finer blocks ONLY
           !----------------------------------------------------------------------------------
           !if (gr_meshMe .eq. 0) print*,"KPD - Entering gr_mgRestrictFaces",i
           call gr_mgRestrictFaces (i, img_denx, img_denx)
           !if (gr_meshMe .eq. 0) print*,"KPD - Entering gr_mgRestrict",i
           call gr_mgRestrict (i, img_temp2, img_res)
           !if (gr_meshMe .eq. 0) print*,"KPD - Leaving gr_mgRestrict",i

           !------------------------------------------------------------------------------
           !- kpd - Compute the residual of the Poisson eqn to be solved using
           !           multigrid (The residual of the residual correction equation).
           !        Now you have the residual at the next coarsest level.
           !       img_src=DUST=RHS, img_soln=DELP=LHS, img_res=index of residual variable
           !------------------------------------------------------------------------------
           !print*,"KPD - Entering mg_residual SECOND in gr_mgCycle.F90 with level # ",i
           !- kpd - This only updates img_res on 3,25,47,69 (leaf blocks at coarser level)
           call mg_restore_nodetypes (i)
           !                        RHS      LHS       RES
           call mg_residualMG (i-1, img_src, img_soln, img_res, 1, 1, img_denx, levelmax)
           !call mg_residual (i-1, img_src, img_soln, img_res, 1, 1)
           !print*,"KPD - Leaving mg_residual SECOND in gr_mgCycle.F90 with level # ",i

        enddo

      endif

      !if (gr_meshMe == 0) print *,' DONE up sweep '

!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------

      !------------------------------------------------------------------
      !- kpd - Solve at the bottom of the V-cycle (coarsest level)
      !------------------------------------------------------------------

      !---------------------------------------------------------------------
      !- kpd - Solve the Poisson Eqn on the current level, using the
      !           Residual Correction Method. Thus, Ae=r is solved, where
      !           e is the error correction and r is the residual. The final
      !           pressure at the new time step will be P(n+1) = P(n) + e.
      !        img_res is the residual calculated previously.
      !        img_corr is the LHS solution, always 0.0 heading in, and
      !           is as "e" for final pressure P(n+1) = P(n) + e.
      !----------------------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------------------
      !- kpd - img_corr receives final solution on all bottom/coarse level LEAF's (3's, 3,4,9,14,...)
      !----------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------
      !        A(h) * u(h)     = f(h)
      !        A    * img_corr = img_res
      !---------------------------------------------------------------------

      !if (gr_meshMe .eq. 0) print*,"Bottom of Vcycle, proc="
      !print*,"Bottom of Vcycle, proc="

      !                               RHS      LHS
      !print*,"Using RBGS at bottom of multigrid with steps = 5500"
      !call mg_relax_RBGS   (levelmin, img_res, img_corr,5500,img_denx,levelmin)
      !call mg_relax   (levelmin, img_res, img_corr,5500)

      call cpu_time(t_start)
      call mg_relax_HYPRE  (levelmin, img_res, img_corr, levelmin)
      call cpu_time(t_stop)
      if (gr_meshMe .eq. 0) print*,"HYPRE Solver Time         ",t_stop-t_start

      !if (gr_meshMe .eq. 0) print*,"Using PFFT Solver"
      !call gr_mgPfftSolveLevel(img_res, img_corr,levelmin)

      !print*,"Leaving mg_relax_UMF at bottom of V-cycle"

      !------------------------------------------------------------------
      !- kpd - Correct the pressure from the previous time step with the new
      !           error correction calculated from the "Ae=r" Poisson solve 
      !           at the current time step. 
      !        This is based on the Residual Correction Method 
      !        P(n+1) = P(n) + e = img_soln + img_corr 
      !------------------------------------------------------------------
      !- kpd - img_soln (DELP) receives updated solution for leaf blocks at
      !           the bottom (coarsest) level only
      !------------------------------------------------------------------
      call gr_mgCorrect (levelmin, img_soln, img_corr, 1)

      !print*,"Leaving mgCorrect at bottom of V-cycle"
      !!!if (gr_meshMe == 0) print *,' DONE coarse solve '

!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------

      !---------------------------------------------------------------------
      !- kpd - Work your way back up the V-cycle... Coarse -> fine leg of V.
      !---------------------------------------------------------------------
      !- kpd - Only looped for more than 1 level.
      !---------------------------------------------------------------------

      !if (level > SolveLevel) then
      if (levelmax > levelmin) then

        !do i = SolveLevel+1, level
        do i = levelmin+1, levelmax

           !print*,"Prolong at Level",i

           !------------------------------------------------------
           !- kpd - Prolong img_corr (e) from coarse to fine level
           !------------------------------------------------------
           call gr_mgProlong (i-1, img_corr, img_corr, 1)

           !---------------------------------------------------------
           !- kpd - Calculate the residual of the residual correction
           !           equation at the finer level
           !---------------------------------------------------------
           !                        RHS      LHS       RES
           call mg_restore_nodetypes (i)
           call mg_residualMG (i, img_res, img_corr, img_res, 0, 0, img_denx, levelmax)
           !call mg_residual (i, img_res, img_corr, img_res, 0, 0)


           call gr_mgZero (i  , img_temp2, 0)
           call gr_mgZero (i-1, img_temp2, 0)

           !--------------------------------------------------
           !- kpd - Solve the Poisson Eqn on the current level
           !--------------------------------------------------
           !                        RHS      LHS
           call mg_relax_RBGS (i, img_res, img_temp2, mgrid_npossmooth,img_denx,levelmax)
           !print*,"Using mg_relax instead of mg_relax_RBGS in mgCycleHYPRE"
           !call mg_relax (i, img_res, img_temp2, mgrid_npossmooth)


           !- kpd-  img_corr is corrected by img_temp2 on all of level 4
           call gr_mgCorrect (i, img_corr, img_temp2, 0)

           !- kpd - img_soln receives img_temp on 5,6,7,8,10,11,12,... (leaf level 4's)
           call gr_mgCopy (i, img_temp, img_soln, 1)

           !- kpd-  img_soln is corrected by img_corr at level i
           call gr_mgCorrect (i, img_soln, img_corr, 1)

        enddo

        !if (gr_meshMe == 0) print *,' DONE down sweep '

      endif

      !print*,"Done with Multigrid Cycle",gr_meshMe


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

 2    return
      end
