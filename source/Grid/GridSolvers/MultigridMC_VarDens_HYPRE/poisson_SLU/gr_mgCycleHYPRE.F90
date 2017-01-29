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
                             img_denx,mg_relax, icycle)

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
                img_denx, &
                icycle

  external mg_solve, mg_residual, mg_residualMG, mg_relax_RBGS, &
           mg_relax_HYPRE,mg_relax

  integer :: i,j 

  real, pointer, dimension(:,:,:,:,:) :: unk
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData


  integer, save :: mgrid_npresmooth,mgrid_npossmooth, mgrid_solveLevelKPD
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
  integer :: iii, SolveLevelKPD


#include "Flash_mpi.h"

  !==========================================================================

!  if (levelmin .eq. levelmax) then
!     SolveLevelKPD = levelmin
!  else
!     SolveLevelKPD = 3
!  end if

  !==========================================================================

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  call Grid_getLocalNumBlks(lnblocks)

  !==========================================================================
  !==========================================================================
  !==========================================================================
  !==========================================================================
  !- kpd - Test loop to check block communication
  !==========================================================================
  !do iii=1,gr_meshNumProcs
  do iii=1,1
  !print*,"========================"
  !do lb = 1, blockCount                                                        
  !    blockID = blockList(lb)
  do lb = 1, lnblocks                                                     
      blockID = lb

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

           !if (blockID .eq. 18 .AND. iii .eq. 1 ) then
           !   print*,"Block3",blockID,i,j, &
           !            facexData(1,i,j,1),faceyData(1,i,j,1)
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
  !==========================================================================
  !==========================================================================

  if (first_call) then
     call RuntimeParameters_get('mgrid_npresmooth',    mgrid_npresmooth)
     call RuntimeParameters_get('mgrid_npossmooth',    mgrid_npossmooth)
     call RuntimeParameters_get('mgrid_solveLevelKPD', mgrid_solveLevelKPD)
     first_call = .false.
  end if

  !if (gr_meshMe .eq. 0) then
  !   print*,"KPD - Entering gr_mgCycleHYPRE V-Cycles with Min/Max no. of levels: ",levelmin,levelmax
  !end if

  !------------------------------------------------------------------------
  !- kpd - Initialize (set = 0,zero out) multigrid variable on a particular 
  !           mesh level (here the residual correction, e)
  !------------------------------------------------------------------------
  !- kpd - Only top level here, other leaf levels in V-cycle loop 
  call gr_mgZero (levelmax, img_corr, 0)


  !- kpd - Assign mixture density to all leaf blocks
  !do i = levelmax, SolveLevelKPD,-1
  do i = levelmax, mgrid_solveLevelKPD,-1
           call mg_restore_nodetypes (i)
  end do
  call gr_mgAssignFden (levelmax, img_denx, 0)
!  do i = levelmax, mgrid_solveLevelKPD,-1
!           call gr_mgRestrictFaces (i, img_denx, img_denx)
!!           call mg_restore_nodetypes (i)
!  end do

      !-----------------------------------------------------------------|
      !- kpd - Work your way down the V-cycle...Fine -> coarse leg of V.|
      !-----------------------------------------------------------------|
      !- kpd - Only looped for more than 1 difference in levels         |
      !-----------------------------------------------------------------|
      !if (levelmax > levelmin) then
      !if (levelmax > SolveLevelKPD) then
      if (levelmax > mgrid_solveLevelKPD) then

        !do i = levelmax, levelmin+1,  -1
        !do i = levelmax, SolveLevelKPD+1,  -1
        do i = levelmax, mgrid_solveLevelKPD+1,  -1

           !print*,"Restrict at Level",i,"with max/min levels:",levelmax,levelmin
           !print*,"Restrict at Level",i,"with proc:",gr_meshMe

           !------------------------------------------------------------------
           !- kpd - Copy img_soln (LHS) to img_temp on the current level leafs
           !------------------------------------------------------------------
           call gr_mgCopy (i, img_soln, img_temp, 1)
           !------------------------------------------------------------------

           !-----------------------------------------------------------------------------
           !- kpd - Initialize (zero out) img_corr  on the i-1 level (includes non-leafs)
           !-----------------------------------------------------------------------------
           call gr_mgZero (i-1, img_corr, 0)
           !-----------------------------------------------------------------------------

           !------------------------------------------------------------------------------
           !- kpd - Solve Poisson Eqn for img_corr on current level (includes non-leafs)
           !        ====================================================================
           !        A(h) * u(h)     = f(h)
           !        A    * e        = r        --> Pnew=Pold+e
           !        A    * Pcorr    = Res  
           !        A    * img_corr = img_res  --> A = div(1/rho*del)
           !---------------------------------------------------------
           !                          RHS      LHS
           call mg_relax_RBGS   (i, img_res, img_corr, mgrid_npresmooth,img_denx,levelmax)
           !------------------------------------------------------------------------------

           !------------------------------------------------------------------
           !- kpd - Correct the pressure with Pcorr just calculated (LEAFS)
           !        Soln = Soln + Corr
           !        r(h) = f(h) - A(h)*u(h)
           !        DELP = DELP + Pcorr
           !------------------------------------------------------------------
           call gr_mgCorrect (i, img_soln, img_corr, 1)
           !------------------------------------------------------------------

           !-------------------------------------------------------------------------------------
           !- kpd - Compute the residual of the residual correction equation (includes non-leafs)
           !        img_temp2 = A*img_corr - img_res
           !        img_temp2 receives residual of Ae=r on LEAFS 
           !-------------------------------------------------------------------------------------
           call mg_restore_nodetypes (i)
           !                        RHS      LHS       RES
           call mg_residualMG (i, img_res, img_corr, img_temp2, 0, 0, img_denx, levelmax)
           !-------------------------------------------------------------------------------------


           !----------------------------------------------------------------------------------
           !- kpd - Restrict img_temp2 = A*img_corr - img_res
           !        Restrict img_temp2 @ fine LEAF --> img_res @ coarse NON_LEAF's
           !        img_temp2 is previously only defined for LEAFS and finer blocks ONLY
           !     ***Now you have img_res defined on ALL NON LEAF blocks at i-1 level
           !----------------------------------------------------------------------------------
           if (icycle .eq. 1) then
              if (gr_meshMe .eq. 0) print*,"Restricting Density At Level: ",i
              call gr_mgRestrictFaces (i, img_denx, img_denx)
           end if
           call gr_mgRestrict (i, img_temp2, img_res)

           !------------------------------------------------------------------------------
           !- kpd - Compute img_res = A*img_soln - img_src for i-1 LEAFS
           !        img_src=DUST=RHS, img_soln=DELP=LHS, img_res=A*P-RHS
           !     ***Now you have img_res defined on ALL blocks at i-1 level
           !------------------------------------------------------------------------------
           call mg_restore_nodetypes (i)
           !                        RHS      LHS       RES
           call mg_residualMG (i-1, img_src, img_soln, img_res, 1, 1, img_denx, levelmax)

        enddo

      endif

!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !- kpd - Solve at the bottom of the V-cycle (coarsest level)
      !--------------------------------------------------------------------------
      !- kpd - Solve the Poisson Eqn on the current level, using the
      !           Residual Correction Method. Thus, Ae=r is solved, where
      !           e is the error correction and r is the residual. The final
      !           pressure at the new time step will be P(n+1) = P(n) + e.
      !--------------------------------------------------------------------------
      !- kpd - Solve for img_corr on all bottom level Blocks (includes NON LEAFS)
      !--------------------------------------------------------------------------
      !        A(h) * u(h)     = f(h)
      !        A    * img_corr = img_res
      !        A    * e        = r
      !--------------------------------------------------------------------------

      !print*,"Using RBGS at bottom of multigrid with steps = 5500"
      !call mg_relax_RBGS   (SolveLevelKPD, img_res, img_corr,5500,img_denx,SolveLevelKPD)
!      call mg_relax_RBGS   (mgrid_solveLevelKPD, img_res, img_corr,1000,img_denx,mgrid_solveLevelKPD)

      call cpu_time(t_start)
      !                                     RHS       LHS
      !call mg_relax_HYPRE  (SolveLevelKPD, img_res, img_corr, SolveLevelKPD)
      call mg_relax_HYPRE  (mgrid_solveLevelKPD, img_res, img_corr, mgrid_solveLevelKPD)
      call cpu_time(t_stop)
      if (gr_meshMe .eq. 0) print*,"HYPRE All Solver Time         ",t_stop-t_start

      !------------------------------------------------------------------
      !- kpd - Correct the pressure (DELP) on BOTTOM level LEAFS only 
      !        P(n+1) = P(n) + e = img_soln + img_corr 
      !------------------------------------------------------------------
      !call gr_mgCorrect (SolveLevelKPD, img_soln, img_corr, 1)
      call gr_mgCorrect (mgrid_solveLevelKPD, img_soln, img_corr, 1)
      !------------------------------------------------------------------

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
      !if (levelmax > levelmin) then
      !if (levelmax > SolveLevelKPD) then
      if (levelmax > mgrid_solveLevelKPD) then

        !do i = SolveLevel+1, level
        !do i = levelmin+1, levelmax
        !do i = SolveLevelKPD+1, levelmax
        do i = mgrid_solveLevelKPD+1, levelmax

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

      endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

 2    return
      end
