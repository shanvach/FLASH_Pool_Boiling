!!****if* source/Grid/GridSolvers/MultigridMC_VarDens_HYPRE/poisson/poisson_mg_relax_HYPRE
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Driver routine for poisson solvers in the Grid
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!  
!!  Most implementations only support a limited subset of boundary condition types.
!!  Some implementations require that all significant elements of bcTypes are the same.
!!  (That is always the case for GRID_PDE_BND_ISOLATED.)
!!  Even if an implementation supports combinations of different boundary conditions
!!  on different sides of the domain, the types at left and right sides for the same
!!  axis direction will usually have to be the samme.
!!
!!  Suport in some implementations provided with FLASH3.3:
!!
!!   GridSolvers/Multipole:                  GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/Pfft:                       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.), 
!!                                            GRID_PDE_BND_NEUMANN(hom.)
!!                                            in various combinations,
!!                                            depending on GridSolvers/Pfft subdirectory 
!!                                            (i.e. implementation configured in)
!!                                           (requires UG in pencil shape or Paramesh as Grid)
!!   GridSolvers/Multigrid hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations,
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid)
!!   
!!***

!!!                             LHS    RHS     BCs      BCs       BCs
!!subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)
!!subroutine poisson_mg_relax_HYPRE (level, irhs, ilhs, nsmooth, idenvar, levelmax)

!                                        rhs   lhs
subroutine poisson_mg_relax_HYPRE (level,iSrc, iSoln, levelmax)

  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Simulation_data,  ONLY : sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax 
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  !use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreCreateMatrix_KPD, &
  use gr_interface,     ONLY : gr_hypreCreateMatrix, &
                               gr_hypreComputeB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel, Grid_getLocalNumBlks, &
                               Grid_getBlkBoundBox,                       &
                               GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET, Grid_getBlkBC

  use Driver_data,    ONLY : dr_nstep

  
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hypreSetup
 
  use tree, only : maxblocks_tr, lrefine, grid_changed, nodetype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get 

  Use paramesh_mpi_interfaces, only : amr_morton_process

  use tree, only : newchild

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)     :: iSoln, iSrc, level, levelmax
  !integer, intent(in)    :: bcTypes(6)
  !real, intent(in)       :: bcValues(2,6)
  integer                 :: bcTypes(6)
  real                    :: bcValues(2,6)
  !real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW  
  
  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec,facexData,faceyData,facezData
  integer :: blockID, lb, lb_AMR
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  logical :: mask(NUNK_VARS)

  integer :: iFactorB
  integer :: iFactorA
  integer :: iFactorC  
  real    :: dt, theta
  real    :: t_start,t_stop


  integer :: blockCount, lnblocks, blockCount_AMR, iii
  integer :: blockList(MAXBLOCKS), blockList_AMR(MAXBLOCKS)

  !kpd - for density matching
  integer :: nodetype_perm(maxblocks_tr)
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer :: localbcTypes(6) 
  character(len=32) :: matfile

  real, pointer, dimension(:,:,:,:), save :: unk

  real, dimension(2,MDIM) :: boundBox
  real :: rXmin,rYmin,rZmin

  !- kpd - Fix for INTEL compiler
  integer, dimension(NDIM) :: temp_gr_hypreLower,temp_gr_hypreUpper 
  real, dimension(NXB*NYB*NZB) :: temp_RHSVal

  logical :: newchild_old(maxblocks_tr) !Shizhao
  
  call Timers_start("Grid_solvePoisson")    

  call RuntimeParameters_get('xmin', rXmin)
  call RuntimeParameters_get('ymin', rYmin)
  call RuntimeParameters_get('zmin', rZmin)

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call Grid_getLocalNumBlks(lnblocks)

  !!-----------------------------------------------------------------------
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------

  !nodetype_perm(:) = nodetype(:)
  nodetype_perm(1:lnblocks) = nodetype(1:lnblocks)

#ifdef PREVIOUS_NODETYPE
  blockCount_AMR = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then

        blockCount_AMR = blockCount_AMR + 1
        blockList_AMR(blockCount_AMR) = lb_AMR

        nodetype(lb_AMR) = 1
     else
        nodetype(lb_AMR) = 2
     end if
  end do
#endif

  !---------------------------------
  ! Shizhao
  newchild_old(1:lnblocks) = newchild(1:lnblocks)

  blockCount_AMR = 0
  newchild = .FALSE.
  do lb_AMR = 1, lnblocks
    nodetype(lb_AMR) = nodetype_perm(lb_AMR)
    if (lrefine(lb_AMR) > level) nodetype(lb_AMR) = -1
    if (lrefine(lb_AMR) == level) then
      nodetype(lb_AMR) = 1
      blockCount_AMR = blockCount_AMR + 1
      blockList_AMR(blockCount_AMR) = lb_AMR
      newchild(lb_AMR) = .TRUE.
    end if
    if (lrefine(lb_AMR) == level-1 .and. nodetype(lb_AMR) /= 1) &
      nodetype(lb_AMR) = 2
  end do
  call amr_morton_process()
  !call mg_restore_nodetypes (level)  ! Shizhao.
  !---------------------------------

  call gr_hypreGridStatus (blockCount_AMR, blockList_AMR)  

  !KPD - Print which cores are being solved on and how many blocks they have...
  !if (dr_nstep .eq. 1 .AND. blockCount_AMR .gt. 0) print*,"Bottom Level Dist:",gr_meshMe,blockCount_AMR

#ifdef NOTETYPR_DEBUG
  call mg_restore_nodetypes (level)  ! Shizhao.
  ! Shizhao, debug
  blockCount_AMR = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then
        blockCount_AMR = blockCount_AMR + 1
     end if
  end do
  if (dr_nstep .eq. 1 .AND. blockCount_AMR .gt. 0) print*,"Bottom Level Dist 2:",gr_meshMe,blockCount_AMR
  ! Shizhao, debug
  blockCount_AMR = 0
  do lb_AMR = 1, lnblocks
     if (nodetype(lb_AMR) == 1) then

        blockCount_AMR = blockCount_AMR + 1
     end if
  end do
  if (dr_nstep .eq. 1 .AND. blockCount_AMR .gt. 0) print*,"Bottom Level Dist 3:",gr_meshMe,blockCount_AMR

  ! Shizhao, debug
  blockCount_AMR = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then

        blockCount_AMR = blockCount_AMR + 1
        blockList_AMR(blockCount_AMR) = lb_AMR

        nodetype(lb_AMR) = 1
     elseif(lrefine(lb_AMR) == level-1) then
        nodetype(lb_AMR) = 2
     else
        nodetype(lb_AMR) = -1
     end if
  end do
  if (dr_nstep .eq. 1 .AND. blockCount_AMR .gt. 0) print*,"Bottom Level Dist 4:",gr_meshMe,blockCount_AMR

#endif
  !!-----------------------------------------------------------------------

  !- kpd - This CAN NOT have the same integer value as a Flash.h Variable,
  !           otherwise this will overwrite it !!!  
  !        iFactorB is commented out in my implementation, bc it is only 
  !           needed for using HYPRE with AMR > 1 level. I am only using 
  !           HYPRE at the bottom of the V-Cycle where the grid is uniform.
  !iFactorB = 1
  !iFactorB = 1

  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.

  !!---------------------------------------------------------------------------------
  !- kpd - I edited these out, these are artifacts fromwhen gr_hypreComputeB is used
  !!---------------------------------------------------------------------------------
  !call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)   
  !call HYPRE_SStructVectorGather(gr_hypreVecB, ierr)     
  !!-----------------------------------------------------------------------

  !************************************************************************
  !!-----------------------------------------------------------------------
  !!-----------------------------------------------------------------------
  !do lb = 1, blockCount 
  !   blockID = blockList(lb)
  lb = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then
   
     blockID = lb_AMR
     lb = lb+1

     !print*,"BOTTOM LEVEL BLOCK/NODETYPE  I",blockID,nodetype(blockID),nodetype_perm(blockID)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getBlkPtr(blockID, unk, CENTER)
     call Grid_getBlkRefineLevel(blockID,mylevel)

     call Grid_getBlkBoundBox(blockId,boundBox)
     
     !print*,"BOUND BOX",gr_meshMe,blockID,boundBox(1,1),boundBox(1,2),boundBox(1,3)
     !print*,"SIMULATIO",gr_meshMe,rXmin, rYmin, rZmin

     !- kpd - mypart should = 0 for uniform grid 
     mypart = mylevel - gr_hypreRefineMIN
     !print*,"KPD MY PART",mypart
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
     
     !if (gr_meshMe .eq. 0) print*,"RHS Size: ",lb, product(dataSize(1:NDIM))

     allocate(RHSVal(product(dataSize(1:NDIM))))     
     
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  

              !- kpd - solnVec-->unk
              !RHSVal(ii) = -solnVec(iSrc,i,j,k)
              !RHSVal(ii) = -unk(iSrc,i,j,k)
              temp_RHSVal(ii) = -unk(iSrc,i,j,k)

              !*********************************
              !- kpd - For FLASH4 implementation
              !*********************************
              if (lb .eq. 1 .AND.                                    &
                                  boundBox(1,1).eq.rXmin .AND.       &
                                  boundBox(1,2).eq.rYmin .AND.       &  
                                  boundBox(1,3).eq.rZmin .AND.       &  
                                  i .eq. blkLimits(LOW, IAXIS) .AND. &
                                  j .eq. blkLimits(LOW, JAXIS) .AND. &
                                  k .eq. blkLimits(LOW, KAXIS)) then 
                 !RHSVal(ii) = 0.0
                 temp_RHSVal(ii) = 0.0
                 !print*,"REFERENCE PRESSURE:",gr_meshMe,blockID,i,j,k
              end if

              !print*,"Bvec",lb,blockID,i,j,ii,solnVec(iSrc,i,j,k)
              
           end do
        end do
     end do
     

     !- kpd - Initial Conditions, iFactorB is not used for my implementation,
     !           and the initial solution vector is zeroed out
     !solnVec(iFactorB,:,:,:) = 1.0

     !- kpd - solnVec-->unk
     !solnVec(iSoln,:,:,:)    = 0.0
     unk(iSoln,:,:,:)    = 0.0

     temp_gr_hypreLower = 0.0
     temp_gr_hypreUpper = 0.0
!     temp_RHSVal = 0.0

     do iii=1,NDIM
        temp_gr_hypreLower(iii) = gr_hypreLower(lb, iii)
        temp_gr_hypreUpper(iii) = gr_hypreUpper(lb, iii)
     end do

!     do iii=1,(NXB*NYB*NZB)
!        temp_RHSVal(iii) = RHSVal(iii)
!     end do
          
     !print*,"SIZES",size(temp_gr_hypreLower),size(temp_gr_hypreUpper),size(temp_RHSVal),size(RHSVal)

     !!-----------------------------------------------------------------------
     !call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
     !     gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)
     !call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
     !     gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)  
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, temp_RHSVal, ierr)  
     !!-----------------------------------------------------------------------
     
     call Grid_releaseBlkPtr(blockID, solnVec)
     call Grid_releaseBlkPtr(blockID, unk, CENTER)
     
     deallocate (RHSVal)
     
     end if  !- kpd - End lref=level if
  end do     !- kpd - End blocklist loop


     !---------------------------------------------------
     !- kpd - This is where the b Vector is ASSEMBLED!
     !---------------------------------------------------
     call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)
     !---------------------------------------------------
  !!-----------------------------------------------------------------------
  !!-----------------------------------------------------------------------


  !***********************************************************************************************
  !***********************************************************************************************
  !***********************************************************************************************
  !KPD - DENSITY MATCHING...
  !***********************************************************************************************
  !***********************************************************************************************
  !***********************************************************************************************
  gcMask = .FALSE.
  gcMask(NUNK_VARS+MGW8_FACE_VAR) = .TRUE.                 ! X-dens
  gcMask(NUNK_VARS+1*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! Y-dens
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! Z-dens
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=LEAF)
  !***********************************************************************************************
  !nodetype(:) = nodetype_perm(:)
  nodetype(1:lnblocks) = nodetype_perm(1:lnblocks) ! Shizhao add 1:lnblocks
  newchild = newchild_old ! Shizhao
  call amr_morton_process() ! Shizhao
 
  !***********************************************************************************************
  lb = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then

     blockID = lb_AMR
     lb = lb+1

     ! Get blocks BCs:
     call Grid_getBlkBC(blockID,faces,onBoundary)

     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     !KPD - X-Dir Sweep for Density Mismatch
     do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
           !do i=GRID_ILO_GC,GRID_IHI_GC

           !KPD - Low-Side
           if ( (ABS(facexData(MGW8_FACE_VAR,1,j,k) - facexData(MGW8_FACE_VAR,NGUARD+2,j,k)) .GT. 1e-8 ) .AND. &
                (onBoundary(LOW,IAXIS) .ne. PERIODIC) ) then
              facexData(RH1F_FACE_VAR,NGUARD+1,j,k) = 0.0
              facexData(RH2F_FACE_VAR,NGUARD+1,j,k) = (facexData(MGW8_FACE_VAR,1,j,k) + &
                                                       facexData(MGW8_FACE_VAR,NGUARD+2,j,k))/2.
              !print*,"X MisMatch Low ",blockID,NGUARD+1,j,k,facexData(MGW8_FACE_VAR,1,j,k),facexData(MGW8_FACE_VAR,NGUARD+2,j,k)
           end if

           !KPD - High-Side
           if ((ABS(facexData(MGW8_FACE_VAR,GRID_IHI_GC,j,k) - facexData(MGW8_FACE_VAR,GRID_IHI-2,j,k)) .GT. 1e-8 ) .AND. &
              (onBoundary(HIGH,IAXIS) .ne. PERIODIC) ) then
              facexData(RH1F_FACE_VAR,GRID_IHI+1,j,k) = 0.0
              facexData(RH2F_FACE_VAR,GRID_IHI+1,j,k) = (facexData(MGW8_FACE_VAR,GRID_IHI_GC,j,k) + &
                                                facexData(MGW8_FACE_VAR,GRID_IHI-2 ,j,k))/2.
              !print*,"X MisMatch High",blockID,GRID_IHI+1,j,k,facexData(MGW8_FACE_VAR,GRID_IHI_GC,j,k),facexData(MGW8_FACE_VAR,GRID_IHI-2,j,k)
           end if
           !end do
        end do
     end do

     !KPD - Y-Dir Sweep for Density Mismatch
     do k=GRID_KLO_GC,GRID_KHI_GC
        !do j=GRID_JLO_GC,GRID_JHI_GC
           do i=GRID_ILO_GC,GRID_IHI_GC

           !KPD - Low-Side
           if ((ABS(faceyData(MGW8_FACE_VAR,i,1,k) - faceyData(MGW8_FACE_VAR,i,NGUARD+2,k)) .GT. 1e-8) .AND. &
               (onBoundary(LOW,JAXIS) .ne. PERIODIC) ) then
              faceyData(RH1F_FACE_VAR,i,NGUARD+1,k) = 0.0
              faceyData(RH2F_FACE_VAR,i,NGUARD+1,k) = (faceyData(MGW8_FACE_VAR,i,1,k) + &
                                                       faceyData(MGW8_FACE_VAR,i,NGUARD+2,k))/2.
              !print*,"Y MisMatch Low ",blockID,i,NGUARD+1,k,faceyData(MGW8_FACE_VAR,i,1,k),faceyData(MGW8_FACE_VAR,i,NGUARD+2,k)
           end if

           !KPD - High-Side
           if ((ABS(faceyData(MGW8_FACE_VAR,i,GRID_JHI_GC,k) - faceyData(MGW8_FACE_VAR,i,GRID_JHI-2,k)) .GT. 1e-8 ) .AND. &
               (onBoundary(HIGH,JAXIS) .ne. PERIODIC) ) then
              faceyData(RH1F_FACE_VAR,i,GRID_JHI+1,k) = 0.0
              faceyData(RH2F_FACE_VAR,i,GRID_JHI+1,k) = (faceyData(MGW8_FACE_VAR,i,GRID_JHI_GC,k) + &
                                                         faceyData(MGW8_FACE_VAR,i,GRID_JHI-2 ,k))/2.
              !print*,"Y MisMatch High",blockID,i,GRID_JHI+1,k,faceyData(MGW8_FACE_VAR,i,GRID_JHI_GC,k),faceyData(MGW8_FACE_VAR,i,GRID_JHI-2,k)
           end if
           end do
        !end do
     end do

#if NDIM == 3
     !KPD - Z-Dir Sweep for Density Mismatch
     !do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
           do i=GRID_ILO_GC,GRID_IHI_GC

           !KPD - Low-Side
           if ((ABS(facezData(MGW8_FACE_VAR,i,j,1) - facezData(MGW8_FACE_VAR,i,j,NGUARD+2)) .GT. 1e-8) .AND. &
               (onBoundary(LOW,KAXIS) .ne. PERIODIC) ) then
              facezData(RH1F_FACE_VAR,i,j,NGUARD+1) = 0.0
              facezData(RH2F_FACE_VAR,i,j,NGUARD+1) = (facezData(MGW8_FACE_VAR,i,j,1) + &
                                                       facezData(MGW8_FACE_VAR,i,j,NGUARD+2))/2.
              !print*,"Z MisMatch Low ",blockID,i,j,NGUARD+1,facezData(MGW8_FACE_VAR,i,j,1),facezData(MGW8_FACE_VAR,i,j,NGUARD+2)
           end if

           !KPD - High-Side
           if ((ABS(facezData(MGW8_FACE_VAR,i,j,GRID_KHI_GC) - facezData(MGW8_FACE_VAR,i,j,GRID_KHI-2)) .GT. 1e-8 ) .AND. &
               (onBoundary(HIGH,KAXIS) .ne. PERIODIC) ) then
              facezData(RH1F_FACE_VAR,i,j,GRID_KHI+1) = 0.0
              facezData(RH2F_FACE_VAR,i,j,GRID_KHI+1) = (facezData(MGW8_FACE_VAR,i,j,GRID_KHI_GC) + &
                                                         facezData(MGW8_FACE_VAR,i,j,GRID_KHI-2 ))/2.
              !print*,"Z MisMatch High",blockID,i,j,GRID_KHI+1,facezData(MGW8_FACE_VAR,i,j,GRID_KHI_GC),facezData(MGW8_FACE_VAR,i,j,GRID_KHI-2)
           end if
           end do
        end do
     !end do
#endif

     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     end if
  end do
  !***********************************************************************************************
  !***********************************************************************************************
  !***********************************************************************************************
  !KPD - End Density Matching ^^^
  !***********************************************************************************************
  !***********************************************************************************************
  !***********************************************************************************************
  
  !!-----------------------------------------------------------------------
  !!     2. Exchange iFactorB across processors. Needs to be done only in
  !!        PARAMESH / AMR, if UG the function will return without any action.
  !!-----------------------------------------------------------------------
  !!call gr_hypreExchangeFacB (iFactorB, blockCount, blockList, level)  
  
  !!-----------------------------------------------------------------------
  !!     3. Set initial guess for solver typically 
  !!        iVar at previous time step is used.
  !- kpd -  The x Vector is Assembled in here
  !!-----------------------------------------------------------------------
  !call gr_hypreSetIniGuess (iSoln, blockCount, blockList, level)      
  !call gr_hypreSetIniGuess (iSoln, blockCount_AMR, blockList_AMR, level)      
  call gr_hypreSetIniGuess (iSoln, blockCount_AMR, blockList_AMR)      

  !!-----------------------------------------------------------------------
  !!     4. Compute the actual A matrix in AX=B
  !!-----------------------------------------------------------------------  

  !! iFactorA, iFactorC not used in routine.
  iFactorA = 1
  iFactorC = 1
  dt       = 1.0
  theta    = 1.0 
  
  localbcTypes = GRID_PDE_BND_NEUMANN
  !localbcTypes = GRID_PDE_BND_PERIODIC    !- kpd - Does this even change the BC from Neumann ???
  bcValues(1,:) = 0.0                      !        The faces array is set from the flash.par global BCs!!!
  bcValues(2,:) = 0.0
  

  !- kpd - This Only works on the A matrix (Interior and Boundary)
  !             --------------------------
  !        The A matrix is also assembled in here!
  !            -----------------------------------
  !print*,"USING gr_hypreCreateMatrix NOT gr_hypreCreateMatrix_KPD"
  !call gr_hypreCreateMatrix(iSoln, iFactorB, iFactorA, localbcTypes, bcValues, dt, theta,  &
  !     blockCount, blockList, .TRUE., iFactorC)       
  !call gr_hypreCreateMatrix_KPD(iSoln, localbcTypes, bcValues, dt, theta,  &
  !     blockCount, blockList, .TRUE., level)       
  call gr_hypreCreateMatrix_KPD(iSoln, localbcTypes, bcValues, dt, theta,  &
       blockCount_AMR, blockList_AMR, .TRUE., level)       
 
  !!$  call gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, &
  !!$       bcTypes, bcValues, iFactorD)  
  
  
  !if (dr_nstep .gt. 1) then
  !           call Driver_abortFlash("HYPRE stop")
  !end if
  !!-----------------------------------------------------------------------
  !!     7. Solve AX = B
  !!-----------------------------------------------------------------------
  !- kpd - The b Vector is assembled in here
  call cpu_time(t_start)
  call gr_hypreSolve ()
  call cpu_time(t_stop)
  if (gr_meshMe .eq. 0) print*,"gr_hypreSolve Solver Time         ",t_stop-t_start
  
  !!-----------------------------------------------------------------------
  !!     8. Update unk variable using X.
  !!-----------------------------------------------------------------------
  !                        ilhs
  call gr_hypreUpdateSoln (iSoln, blockCount_AMR, blockList_AMR)  
  
  
  call Timers_stop("Grid_solvePoisson") 
  
  return
end subroutine poisson_mg_relax_HYPRE

