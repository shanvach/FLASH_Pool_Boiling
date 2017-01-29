!!****if* source/Grid/GridSolvers/MultigridMC_VarDens_HYPRE/poisson/2D_MUMPS_poisson_mg_relax_HYPRE
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
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreCreateMatrix_KPD, &
                               gr_hypreComputeB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel, Grid_getLocalNumBlks

  use Driver_data,    ONLY : dr_nstep

  
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hypreSetup, gr_hypreSolver
 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET

  use tree, only : lrefine

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)     :: iSoln, iSrc, level, levelmax
  !integer, intent(in)    :: bcTypes(6)
  !real, intent(in)       :: bcValues(2,6)
  integer                 :: bcTypes(6)
  real                    :: bcValues(2,6)
  !real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW  
  
  integer :: mylevel, mypart, var, ii,jj,i,j,k,m,n, ierr
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID, lb, lb_AMR
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  logical :: mask(NUNK_VARS)

  integer :: iFactorB
  integer :: iFactorA
  integer :: iFactorC  
  real    :: dt, theta

  real :: RHSVal2(1024)
  real :: AA,RR
  real, dimension(MDIM)     :: del, coord, bsize
  integer :: nentries
  integer, dimension(2,MDIM):: faces
  real, allocatable :: kpdA(:)

  integer :: blockCount, lnblocks, blockCount_AMR, iii
  integer :: blockList(MAXBLOCKS), blockList_AMR(MAXBLOCKS)

  integer :: localbcTypes(6) 
  character(len=32) :: matfile

  real, pointer, dimension(:,:,:,:), save :: unk

  integer, allocatable :: blockOrder(:),zPlane(:)
  real, allocatable  :: blockX(:),blockY(:),blockZ(:)
  real, allocatable  :: yzPlane(:,:)
  real, dimension(2,MDIM) :: boundBox
  integer :: nBlocksInXDir,nBlocksInYDir,nBlocksInZDir,iCount,jCount,iFlag, NN, NZ
  real :: rx, ry, rz

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  call Timers_start("Grid_solvePoisson")    

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call Grid_getLocalNumBlks(lnblocks)

  !!-----------------------------------------------------------------------
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------
  blockCount_AMR = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then

        blockCount_AMR = blockCount_AMR + 1
        blockList_AMR(blockCount_AMR) = lb_AMR

     end if
  end do
  call gr_hypreGridStatus (blockCount_AMR, blockList_AMR)  

  blockCount=blockCount_AMR
  allocate(blockOrder(lnblocks))
  allocate(zPlane(lnblocks))
  allocate(blockX(lnblocks))
  allocate(blockY(lnblocks))
  allocate(blockZ(lnblocks))
  allocate(yzPlane(blockCount,lnblocks))

  !!-----------------------------------------------------------------------

  !- kpd - This CAN NOT have the same integer value as a Flash.h Variable,
  !           otherwise this will overwrite it !!!  
  !        iFactorB is commented out in my implementation, bc it is only 
  !           needed for using HYPRE with AMR > 1 level. I am only using 
  !           HYPRE at the bottom of the V-Cycle where the grid is uniform.
  !iFactorB = 1
  iFactorB = 123

  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.

  !!---------------------------------------------------------------------------------
  !- kpd - I edited these out, these are artifacts from when gr_hypreComputeB is used
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

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getBlkPtr(blockID, unk, CENTER)
     call Grid_getBlkRefineLevel(blockID,mylevel)
     

     !- kpd - mypart should = 0 for uniform grid 
     mypart = mylevel - gr_hypreRefineMIN
     !print*,"KPD MY PART",mypart
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
     
     allocate(RHSVal(product(dataSize(1:NDIM))))     
     
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  

              !- kpd - solnVec-->unk
              !RHSVal(ii) = -solnVec(iSrc,i,j,k)
              RHSVal(ii) = -unk(iSrc,i,j,k)           !- kpd - iSrc is the input RHS

              !*********************************
              !- kpd - For FLASH4 implementation
              !*********************************
              if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                  i .eq. blkLimits(LOW, IAXIS) .AND. &
                                  j .eq. blkLimits(LOW, JAXIS) .AND. &
                                  k .eq. blkLimits(LOW, KAXIS)) then 
                 RHSVal(ii) = 0.0
              end if

              !print*,"BLOCKID",blockID,i,j,k,ii,gr_meshMe
              !if (gr_meshMe .eq. 1 .AND. blockID .eq. 2 .AND. i .eq. 6 .AND. j.eq.10) then
              !   RHSVal(ii) = 999999.999999
              !end if


              !print*,"RHS",lb,blockID,i,j,ii,RHSVal(ii)
              
           end do
        end do
     end do
     
     !- kpd - Initial Conditions, iFactorB is not used for my implementation,
     !           and the initial solution vector is zeroed out
     !solnVec(iFactorB,:,:,:) = 1.0
     unk(iFactorB,:,:,:) = 1.0

     !- kpd - solnVec-->unk
     !solnVec(iSoln,:,:,:)    = 0.0
     unk(iSoln,:,:,:)    = 0.0
          
     !!-----------------------------------------------------------------------
     !call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
     !     gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)  
     !!-----------------------------------------------------------------------
     
     call Grid_releaseBlkPtr(blockID, solnVec)
     call Grid_releaseBlkPtr(blockID, unk, CENTER)
     
     deallocate (RHSVal)
     
     end if  !- kpd - End lref=level if
  end do     !- kpd - End blocklist loop

!**************************************************************************************************
!********************** MUMPS CODE ****************************************************************
!**************************************************************************************************


        !- kpd - Count Blocks in X,Y, and Z directions 
        !---------------------------------------------------------------
        nBlocksInXDir=1
        nBlocksInYDir=1
        nBlocksInZDir=1
        iCount = 0
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
           call Grid_getBlkCenterCoords(lb,coord)
           iCount = iCount+1
           if (iCount .eq. 1) then
              rx = coord(IAXIS)
              ry = coord(JAXIS)
#if NDIM==3
              rz = coord(KAXIS)
#endif
           end if
           if (iCount .gt. 1 .AND. coord(IAXIS) .eq. rx .AND. coord(JAXIS) .eq. ry) nBlocksInXDir = nBlocksInXDir + 1
           if (iCount .gt. 1 .AND. coord(JAXIS) .eq. ry .AND. coord(IAXIS) .eq. rx) nBlocksInYDir = nBlocksInYDir + 1
#if NDIM==3
           if (iCount .gt. 1 .AND. coord(KAXIS) .eq. rz .AND. coord(JAXIS) .eq. ry ) nBlocksInZDir = nBlocksInZDir + 1
#endif
           end if
        end do

        !- kpd - Calculate the lower left corner of each block
        !-----------------------------------------------------
        do lb = 1, lnblocks

           blockX(lb) = 10000.0
           blockY(lb) = 10000.0
           blockZ(lb) = 10000.0

           if (lrefine(lb) == level) then

              !- kpd - Grab block pointers 
              call Grid_getBlkCenterCoords(lb,coord)

              blockX(lb) = coord(IAXIS)
              blockY(lb) = coord(JAXIS)
#if NDIM==3
              blockZ(lb) = coord(KAXIS)
              print*,"Block",lb,coord(IAXIS),coord(JAXIS),coord(KAXIS)
#endif

              !print*,"Block",lb,coord(IAXIS),coord(JAXIS),blockX(lb),blockY(lb)

           end if
        end do

print*,"Blocks in X,Y,Z Dir",nBlocksInXDir,nBlocksInYDir,nBlocksInZDir,"BlockCount",blockCount

        !- kpd - Find out what z-plane each block lies in
        !---------------------------------------------------------------
        iCount = 0
        do jj=1,nBlocksInZDir
           iCount = iCount+1
           rz = MINVAL(blockZ)
           do lb = 1, lnblocks
              if (lrefine(lb) == level) then
              !- kpd - Grab block pointers 
              call Grid_getBlkCenterCoords(lb,coord)
              if (blockZ(lb) .eq. rz ) then
                 zPlane(lb) = iCount
                 blockZ(lb) = 100000.0
                 !print*,"Zplane",lb,coord(KAXIS),zPlane(lb)
              end if
              end if
           end do
        end do

        do jj=1,nBlocksInZDir
           do lb = 1, lnblocks
              if (lrefine(lb) == level) then
              !- kpd - Grab block pointers 
              call Grid_getBlkCenterCoords(lb,coord)

                 yzPlane(jj) = coord(JAXIS)

              end if
           end do
        end do

        !- kpd - Put the XY plane blocks in left-to-right --> top-to-bottom order
        !------------------------------------------------------------------------
        iCount = 0
        jCount = 0
        do jj=1,1!nBlocksInZDir
           jCount = jCount+1


        !iCount = 0
        do j=1,blockCount 
           iFlag = 0


           do lb = 1, lnblocks
              if (lrefine(lb) == level) then

print*,lb,jCount,zPlane(lb),"Y",blockY(lb),MINVAL(blockY)

                 if (blockY(lb) .eq. MINVAL(blockY) .AND. iFlag .eq. 0 .AND. zPlane(lb) .eq. jCount) then
                 !if (blockY(lb) .eq. MINVAL(blockY) .AND. iFlag .eq. 0) then
                    iCount = iCount+1
                    iFlag  = 1
                    blockOrder(iCount) = lb
                    blockY(lb) = 10000.0
                 end if

              end if
           end do
        end do

        end do

        !!- kpd - Count Blocks in X and Y directions 
        !!---------------------------------------------------------------
        !nBlocksInXDir=1
        !nBlocksInYDir=1
        !do lb=1,blockCount
        !   blockID = blockOrder(lb)
        !   call Grid_getBlkCenterCoords(blockID,coord)
        !   if (lb .eq. 1) then
        !      rx = coord(IAXIS)
        !      ry = coord(JAXIS)
        !   end if
        !   if (lb .gt. 1 .AND. coord(IAXIS) .eq. rx ) nBlocksInXDir = nBlocksInXDir + 1 
        !   if (lb .gt. 1 .AND. coord(JAXIS) .eq. ry ) nBlocksInYDir = nBlocksInYDir + 1 
        !end do

        !- kpd - Test output
        do i = 1, blockCount
           print*,"ORDER",i,nBlocksInXDir,nBlocksInYDir,blockOrder(i)
        end do
        call Driver_abortFlash("HYPRE stop")

!**************************************************************************************************
!kpd - MUMPS Test Loop 

  open(unit=3,file="HEADkpdTest.dat")
  open(unit=5,file="BkpdTest.dat")
  open(unit=7,file="AkpdTest.dat")

  nentries = 2*NDIM + 1

  allocate(kpdA(nentries))

  ii = 0
  lb = 0
  NZ = 0

  do m = 1, blockCount

     blockID = blockOrder(m)
     lb = m

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID, unk, CENTER)
     call Grid_getBlkBC (blockID, faces)
     call Grid_getDeltas(blockID, del)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

              ii = ii + 1

              !========================= A MATRIX ==========================================
              do iii=1,nentries
                 kpdA(iii)=0.0
              end do

              ! LHS Boundary (i=LOW)
              if ((i /= blkLimits(LOW , IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then

                    kpdA(2) = -1.d0 / (del(IAXIS)**2.0)

                    !- kpd - Normalize the pressure field
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                    kpdA(2) = 0.0
                    end if

              else
                    kpdA(2) = 0.0
                    !print*,"BD1",blockID,i,j,k
              end if

              ! RHS Boundary (i=HIGH)
              if ((i /= blkLimits(HIGH, IAXIS)) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then

                    kpdA(3) = -1.d0 / (del(IAXIS)**2.0)

                    !- kpd - Normalize the pressure field
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                    kpdA(3) = 0.0
                    end if
              else
                    kpdA(3) = 0.0
                    !print*,"BD2",blockID,i,j,k
              end if

              ! BOTTOM Boundary (j=LOW)
              if ((j /= blkLimits(LOW , JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then

                    kpdA(4) = -1.d0 / (del(IAXIS)**2.0)

                    !- kpd - Normalize the pressure field
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                    kpdA(4) = 0.0
                    end if
              else
                    kpdA(4) = 0.0
                    !print*,"BD3",blockID,i,j,k
              end if

              ! BOTTOM Boundary (j=LOW)
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then

                    kpdA(5) = -1.d0 / (del(IAXIS)**2.0)

                    !- kpd - Normalize the pressure field
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                    kpdA(5) = 0.0
                    end if
              else
                    kpdA(5) = 0.0
                    !print*,"BD4",blockID,i,j,k
              end if

              kpdA(1) = -1.0*(kpdA(2)+kpdA(3)+kpdA(4)+kpdA(5))
              NZ = NZ+1

              !- kpd - Normalize the pressure field
              if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                  i .eq. blkLimits(LOW, IAXIS) .AND. &
                                  j .eq. blkLimits(LOW, JAXIS) .AND. &
                                  k .eq. blkLimits(LOW, KAXIS)) then
              kpdA(1) = 1.0
              end if

                    !-----------------------------------------------------------------------------------------
                    !- kpd - Write A matrix in MUMPS Format...
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBOR",ii,kpdA(1)
                       write(7,*) ii,ii,kpdA(1)

                    if (i .eq. NGUARD+1 .AND. abs(kpdA(2)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBXL",ii-(NXB*NXB-(NXB-1)),kpdA(2)
                       write(7,*) ii,ii-(NXB*NXB-(NXB-1)),kpdA(2)
                       NZ = NZ+1
                    else if (abs(kpdA(2)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBLL",ii-1,kpdA(2)
                       write(7,*) ii,ii-1,kpdA(2)
                       NZ = NZ+1
                    end if

                    if (i .eq. NXB+NGUARD .AND. abs(kpdA(3)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBXR",ii+(NXB*NXB-(NXB-1)),kpdA(3)
                       write(7,*) ii,ii+(NXB*NXB-(NXB-1)),kpdA(3)
                       NZ = NZ+1
                    else if (abs(kpdA(3)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBRR",ii+1,kpdA(3)
                       write(7,*) ii,ii+1,kpdA(3)
                       NZ = NZ+1
                    end if

                    if ( m .eq. blockCount .AND. i .eq. blkLimits(HIGH, IAXIS) .AND. &
                                                 j .eq. blkLimits(HIGH, JAXIS) .AND. &
                                                 k .eq. blkLimits(HIGH, KAXIS)) then
                       if (j .eq. NGUARD+1 .AND. abs(kpdA(4)) .gt. 0.00001) then
                          !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBYL",ii-((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(4)
                          write(7,*) ii,ii-((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(4),"        :values"
                          NZ = NZ+1
                       else if (abs(kpdA(4)) .gt. 0.00001) then
                          !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBBB",ii-NXB,kpdA(4)
                          write(7,*) ii,ii-NXB,kpdA(4),"        :values"
                          NZ = NZ+1
                       end if

                    else

                       if (j .eq. NGUARD+1 .AND. abs(kpdA(4)) .gt. 0.00001) then
                          !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBYL",ii-((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(4)
                          write(7,*) ii,ii-((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(4)
                          NZ = NZ+1
                       else if (abs(kpdA(4)) .gt. 0.00001) then
                          !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBBB",ii-NXB,kpdA(4)
                          write(7,*) ii,ii-NXB,kpdA(4)
                          NZ = NZ+1
                       end if
             
                    end if

                    if (j .eq. NXB+NGUARD .AND. abs(kpdA(5)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBYR",ii+((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(5)
                       write(7,*) ii,ii+((nBlocksInXDir-1)*NXB*NXB+NXB),kpdA(5)
                       NZ = NZ+1
                    else if (abs(kpdA(5)) .gt. 0.00001) then
                       !write(7,*) blockID,i,j,"OWNER",ii,"NEIGHBTT",ii+NXB,kpdA(5)
                       write(7,*) ii,ii+NXB,kpdA(5)
                       NZ = NZ+1
                    end if
                    !-----------------------------------------------------------------------------------------

           end do  !i
        end do     !j
     end do        !k

  end do           !m = 1, blockCount

  print*,"Number of Cells",ii,"NonZeros",NZ
  NN = ii
     write(3,*) NN,"              :N"
     write(3,*) NZ,"              :NZ"

  ii = 0
  lb = 0

  do m = 1, blockCount

     blockID = blockOrder(m)
     lb = m

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID, unk, CENTER)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

              ii = ii + 1

              !========================= b Vector ==========================================

              if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                  i .eq. blkLimits(LOW, IAXIS) .AND. &
                                  j .eq. blkLimits(LOW, JAXIS) .AND. &
                                  k .eq. blkLimits(LOW, KAXIS)) then
                 RR = 0.0
              else
                 RR = -unk(iSrc,i,j,k)
              end if

              !write(5,*) blockID,i,j,k,ii,RR
              if (ii .eq. NN) then
                 write(5,*) RR,"        :RHS"
              else
                 write(5,*) RR
              end if

           end do  !i
        end do     !j
     end do        !k

  end do           !m = 1, blockCount


  deallocate(kpdA)

  close(5)
  close(7)

  deallocate(blockOrder)
  deallocate(zPlane)
  deallocate(blockX)
  deallocate(blockY)
  deallocate(blockZ)
  deallocate(yzPlane)


  !call Driver_abortFlash("HYPRE stop")

!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************

     !---------------------------------------------------
     !- kpd - This is where the b Vector is ASSEMBLED!
     !---------------------------------------------------
     call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)
     !---------------------------------------------------
  !!-----------------------------------------------------------------------
  !!-----------------------------------------------------------------------

  !!mask = .false.  
  !!mask(iSoln) = .true.  
  !!mask(iSrc)  = .true.
  !!mask(iFactorB) = .false.  

  
  !?????????????????????????????????????????????????????????????????????????
  !- kpd - ??????  What is this for  ???????
  !call Timers_start("Grid_fillGuardCells")   
  !call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, mask=mask,selectBlockType=LEAF)  
  !call Timers_stop("Grid_fillGuardCells")   
  !?????????????????????????????????????????????????????????????????????????

  
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
  call gr_hypreSetIniGuess (iSoln, blockCount, blockList, level)      
  
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
  !call gr_hypreCreateMatrix(iSoln, iFactorB, iFactorA, localbcTypes, bcValues, dt, theta,  &
  !     blockCount, blockList, .TRUE., iFactorC, level, idenvar)       
  !call gr_hypreCreateMatrix(iSoln, iFactorB, iFactorA, localbcTypes, bcValues, dt, theta,  &
  !     blockCount, blockList, .TRUE., iFactorC, level )       
  call gr_hypreCreateMatrix_KPD(iSoln, localbcTypes, bcValues, dt, theta,  &
       blockCount, blockList, .TRUE., level)       
 
  !!$  call gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, &
  !!$       bcTypes, bcValues, iFactorD)  
  
  
  !if (dr_nstep .gt. 1) then
  !           call Driver_abortFlash("HYPRE stop")
  !end if
  !!-----------------------------------------------------------------------
  !!     7. Solve AX = B
  !!-----------------------------------------------------------------------
  !- kpd - The b Vector is assembled in here
  call gr_hypreSolve ()
  
  !!-----------------------------------------------------------------------
  !!     8. Update unk variable using X.
  !!-----------------------------------------------------------------------
  !                        ilhs
  call gr_hypreUpdateSoln (iSoln, blockCount, blockList, level)  
  
  
  call Timers_stop("Grid_solvePoisson") 
  
  return
end subroutine poisson_mg_relax_HYPRE

