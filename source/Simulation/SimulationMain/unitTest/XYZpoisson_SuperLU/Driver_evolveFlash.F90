!!****if* source/Simulation/SimulationMain/unitTest/XYZpoisson_SuperLU/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is the main global driver for simulations that are:
!!      Spatially refined, State form, strang split
!!
!!  DOC: Driver_evolveFlash needs more explanation 
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_myPE or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

#include "Flash.h"

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_globalComm, dr_nbegin,       &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime

  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime, Driver_abortFlash
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump

  use Grid_interface,    ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                                GRID_PDE_BND_DIRICHLET, Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement, &
                                Grid_getBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_solvePoisson, &
                                Grid_releaseBlkPtr, &
                                Grid_getBlkCornerID, &
                                Grid_getMaxCommonRefinement, &
                                Grid_getBlkCenterCoords

  use Grid_data, only : gr_meshMe,gr_meshNumProcs,gr_meshComm

  use gr_interface, ONLY : gr_findMean

  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt


  use ut_qsortInterface, ONLY : ut_qsort
  use gr_interfaceTypeDecl
  
#ifdef FLASH_GRID_PARAMESH
  use tree, only : lrefine,lrefine_max,surr_blks
  !use gr_mgPfftData, ONLY : gr_mgbcTypes
#endif
 
  
  implicit none

#include "constants.h"
 include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: endRun  

  integer, dimension(6) :: bc_types, gr_PoissonBcTypes 
  real, dimension(2,6)  :: bc_values
  real poisfact

  integer blkpoints, blkpointsaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux
  
  integer blockID,lb,i,ierr

  integer TA(2),count_rate
  real   ET
  
  real meanPHI,meanPHIaux,meanANL,meanANLaux 

  integer :: eachBoundary
  character(len=MAX_STRING_LENGTH), dimension(2*MDIM) :: bcTypeStr


  !!! slu data module
  integer :: gr_sluLower(MAXBLOCKS,MDIM),gr_sluUpper(MAXBLOCKS,MDIM)
  integer :: slu_max_FullyRefinedLevel


  !!! slu init data
  integer :: blockCount_set
  integer :: blockList_set(MAXBLOCKS)  
  integer :: stride (MDIM), lb_stride(MAXBLOCKS, MDIM)
  integer :: datasize(MDIM)
  integer :: ilower(MDIM), iupper(MDIM)
  

  integer :: aux, blockCount_start_idx, neighProcsCount
  integer, allocatable, dimension(:,:) :: neighProc_blkcnt
  integer, allocatable, dimension(:,:) :: neighProc_blkList
  integer, allocatable, dimension(:,:) :: neighproc_blklist_idx 
  integer :: neighProc_idx,neighProc_blocks


  !!! Find Processors data
  type (AllBlockRegions_t) :: surrBlksSummary
  integer,dimension(BLKNO:TYPENO) :: negh_prop
  integer :: vec_ij(2),vec_jk(2),vec_ik(2)
  integer :: ii,j,k, ENDK, ENDJ, ENDI, numNegh, procCount, procSize, procID
  integer, allocatable, dimension(:) :: procs
  integer, allocatable, dimension(:) :: neighProcsList
  integer :: p,b,u,n,allCenters
  integer :: recvs(gr_meshNumProcs),sts(MPI_STATUS_SIZE,gr_meshNumProcs),vecBuff(2),TAG


  real :: coord(MDIM)

  meanPHI = 0.0
  meanPHIaux =0.0 
  meanANL = 0.0 
  meanANLaux = 0.0

  call Logfile_stamp(dr_globalMe, 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  !Initialise to -1 to help us spot errors.
  bc_types = -1

  ! Retrieve Poisson solution Boundary Conditions for each face: 
  call RuntimeParameters_get("xl_boundary_type", bcTypeStr(1))
  call RuntimeParameters_get("xr_boundary_type", bcTypeStr(2))
  if (NDIM >= 2) then
     call RuntimeParameters_get("yl_boundary_type", bcTypeStr(3))
     call RuntimeParameters_get("yr_boundary_type", bcTypeStr(4))
  endif
  if (NDIM == 3) then
     call RuntimeParameters_get("zl_boundary_type", bcTypeStr(5))
     call RuntimeParameters_get("zr_boundary_type", bcTypeStr(6))
  endif

  ! Boundary Conditions:
!!$#ifdef FLASH_GRID_UG
!!$
!!$
!!$  do eachBoundary = 1, 2*NDIM
!!$
!!$     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), bc_types(eachBoundary))
!!$
!!$  enddo
!!$
!!$
!!$  bc_values = 0.
!!$
!!$#else


  do eachBoundary = 1, 2*NDIM

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), gr_PoissonBcTypes(eachBoundary))

     select case(gr_PoissonBcTypes(eachBoundary))
     case (PERIODIC)
        bc_types(eachBoundary)  = GRID_PDE_BND_PERIODIC
        write(*,*) eachBoundary,'GRID_PDE_BND_PERIODIC'
     case (OUTFLOW)
        bc_types(eachBoundary)  = GRID_PDE_BND_NEUMANN
        write(*,*) eachBoundary,'GRID_PDE_BND_NEUMANN'
     case (DIRICHLET)
        bc_types(eachBoundary)  = GRID_PDE_BND_DIRICHLET
        write(*,*) eachBoundary,'GRID_PDE_BND_DIRICHLET'
     case default
           write(*,*) 'In DriverEvolveFlash: Poisson Problem Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')
     end select
  enddo


  bc_values(1,:) = 0.
  bc_values(2,:) = -1.

!!$#endif





  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (dr_globalMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  


!!$!!! 
!!$  ! Initialization:  
!!$
!!$  !! Global indexing:
!!$  !! ------ --------
!!$
!!$  ! Define a refinement level to do the solve: 
!!$  call Grid_getMaxCommonRefinement(gr_MeshComm, slu_max_FullyRefinedLevel)
!!$  !slu_max_FullyRefinedLevel = 1
!!$  
!!$  ! Make a list of blocks that belong to the level:
!!$  call Grid_getListOfBlocks(REFINEMENT,blockList_set,blockCount_set, &
!!$                            refinementLevel=slu_max_FullyRefinedLevel)
!!$
!!$  ! Set Global numeration of cells in each block, this is dimension by dimension.
!!$  do lb=1, blockCount_set
!!$     blockID = blockList_set(lb)
!!$     ! Get Corner ID and Stride, blkLimits:
!!$     call Grid_getBlkCornerID(blockId, ilower(1:MDIM), stride)
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!$
!!$     write(*,*) ilower(1:MDIM)
!!$
!!$     ilower(1:MDIM) =  ceiling(real(ilower(1:MDIM)) / real(stride(1:MDIM)))
!!$     datasize  (1:MDIM)= max(blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM), datasize(1:MDIM))
!!$     iupper(1:MDIM) = ilower(1:MDIM) + datasize (1:MDIM)
!!$
!!$
!!$     do i=1, NDIM
!!$        gr_sluLower(lb,i) =  ilower(NDIM-i+1)
!!$        gr_sluUpper(lb,i) =  iupper(NDIM-i+1)
!!$     end do
!!$
!!$     lb_stride(lb,1:MDIM) = stride(1:MDIM)
!!$
!!$  end do
!!$  do i=NDIM+1, MDIM
!!$     gr_sluLower(:,i) = 1
!!$     gr_sluUpper(:,i) = 1
!!$  end do
!!$  do lb=1, blockCount_set
!!$     blockID = blockList_set(lb)
!!$     
!!$     write(*,*) "Mype=",dr_globalMe,blockID
!!$     write(*,*) gr_sluLower(lb,1),gr_sluLower(lb,2),gr_sluLower(lb,3)
!!$     write(*,*) gr_sluUpper(lb,1),gr_sluUpper(lb,2),gr_sluUpper(lb,3)
!!$
!!$  enddo
!!$  write(*,*) 'blockCount_set=',blockCount_set,slu_max_FullyRefinedLevel  
!!$  call Driver_abortFlash('Done Global Numeration');
!!$
!!$
!!$  ! Global numeration, by processor and block in the level:
!!$  ! Exchange information on how many blocks have unknowns in
!!$  ! global numeration before mine:
!!$  CALL mpi_scan(blockCount_set,aux,CONSTANT_ONE,FLASH_INTEGER,MPI_SUM,dr_globalComm,ierr)
!!$
!!$  blockCount_start_idx = aux - blockCount_set
!!$  !write(*,*) 'Mype=',dr_globalMe,blockCount_set, blockCount_start_idx+1, &
!!$  !                               blockCount_start_idx+ blockCount_set
!!$  !call Driver_abortFlash('Done Global Numeration');
!!$  ! Loop of blocks goes from  blockCount_start_idx+1 to  blockCount_start_idx+blockCount_set
!!$  
!!$
!!$  ! Find which are the processors that share block boundaries with me:
!!$  ! gives me neighProcsCount  
!!$  ! call gr_pdsluFindNeighs(blockCount_set,blockList_set,...)
!!$  ENDI = 3
!!$  ENDJ = max(1,K2D*3)
!!$  ENDK = max(1,K3D*3)
!!$  allCenters = 2**NDIM
!!$
!!$  ! Up to "2**(NDIM-1)" neighbors for each guard cell region
!!$  ! Exactly "((ENDI*ENDJ*ENDK)-1)" guard cell regions in each block
!!$  ! Exactly "blockCount_set" blocks in this MPI rank.
!!$  ! Add "1" to ensure my MPI rank is in the procs array.
!!$  procSize = (2**(NDIM-1) * ((ENDI*ENDJ*ENDK)-1) * blockCount_set) + 1
!!$  allocate(procs(procSize))
!!$
!!$  procCount = 0
!!$  do b = 1, blockCount_set
!!$     blockID = blockList_set(b)
!!$  
!!$     call Grid_getBlkCenterCoords(blockID,coord)
!!$
!!$     !write(*,*) 'Mype=',gr_meshMe,coord(1:MDIM)
!!$
!!$     ! The call to this routine assumes the block is a leaf block,
!!$     ! It should work for lrefinements that cover the whole computaitonal domain
!!$     ! but will not work for parents.
!!$     ! call gr_findAllNeghID(blockID, surrBlksSummary)
!!$     do k = 1, ENDK
!!$        do j = 1, ENDJ
!!$           do i = 1, ENDI
!!$
!!$           vec_ij = (/ i , j /)
!!$           vec_jk = (/ j , k /)
!!$           vec_ik = (/ i , k /)
!!$            
!!$
!!$           !Initialize all neighbor details in region (i,j,k) to NONEXISTENT.
!!$           surrBlksSummary % regionInfo(i,j,k) % numNegh = NONEXISTENT
!!$           surrBlksSummary % regionInfo(i,j,k) &
!!$                % details(BLKNO:TYPENO, 1:2**(NDIM-1)) = NONEXISTENT
!!$
!!$           if ( ((i*j*k) /= allCenters) .and.  &
!!$                    (all(vec_ij .eq. 2)  .or.  &
!!$                     all(vec_jk .eq. 2)  .or.  &
!!$                     all(vec_ik .eq. 2)) ) then ! No Edge or corner Blocks needed 1 level compute.
!!$
!!$              !May not be a block at this position.  We may be at 
!!$              !the edge of the domain, and so surr_blks at i,j,k contains 
!!$              !the external boundary conditions.
!!$              negh_prop(:) = surr_blks(:,i,j,k,blockID)
!!$
!!$              if(gr_meshMe .eq. MASTER_PE) then
!!$                write(*,*) 'M_PE 0=',i,j,k,negh_prop(PROCNO) 
!!$              endif
!!$
!!$              !First check for an external boundary.
!!$              if (negh_prop(BLKNO) <= PARAMESH_PHYSICAL_BOUNDARY) then
!!$
!!$
!!$                 !write(*,*) 'Mype=',gr_meshMe,'PhysicalBoundary'
!!$                 !This is the case for e.g. an external reflecting boundary.
!!$                 !...Just copy the external boundary conditions.
!!$                 surrBlksSummary % regionInfo(i,j,k) % numNegh = 0
!!$                 surrBlksSummary % regionInfo(i,j,k) &
!!$                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)
!!$
!!$              !Now check for wacky values in surr_blks.
!!$              else if ( &
!!$                   (negh_prop(BLKNO) < -1) .or. &
!!$                   (negh_prop(BLKNO) == 0) .or. &
!!$                   (negh_prop(BLKNO) > MAXBLOCKS) .or. &
!!$                   (negh_prop(PROCNO) < -1) .or. &
!!$                   (negh_prop(PROCNO) >= gr_meshNumProcs) &
!!$                   ) then
!!$
!!$                 print *, gr_meshMe, negh_prop(:)
!!$                 call Driver_abortFlash("Unexpected surr_blks values")
!!$
!!$
!!$              else
!!$
!!$                 !! This is the situation when the neighbor is at the
!!$                 !! same level of resulution. There is only one neighbor
!!$                 !! very simply found.
!!$                 surrBlksSummary % regionInfo(i,j,k) % numNegh = 1
!!$                 surrBlksSummary % regionInfo(i,j,k) &
!!$                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)
!!$
!!$              endif
!!$
!!$
!!$              numNegh = surrBlksSummary % regionInfo(i,j,k) % numNegh
!!$              do n = 1, numNegh
!!$                 procID = surrBlksSummary % regionInfo(i,j,k) % &
!!$                          details(PROCNO,n)
!!$                 if (procID /= gr_meshMe) then
!!$                    procCount = procCount + 1
!!$                    procs(procCount) = procID
!!$                 end if
!!$              end do
!!$           end if
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  if(procCount > 0) then
!!$     call ut_qsort(procs, procCount)
!!$     allocate(neighProcsList(procCount+1))
!!$     neighProcsList = -1
!!$     neighProcsList(1) = gr_meshMe
!!$     u = 1
!!$     do p = 1, procCount
!!$        if (.not.(any(procs(p)==neighProcsList))) then
!!$           u = u + 1
!!$           neighProcsList(u) = procs(p)
!!$        end if
!!$     end do
!!$     deallocate(procs)
!!$  end if
!!$  neighProcsCount = u
!!$
!!$  
!!$  !do p=1,neighProcsCount
!!$     write(*,*) 'Mype=',gr_meshMe,neighProcsList(1:neighProcsCount)
!!$!                        neighProc_blkList_idx(neighProc_blkList(1:neighProc_blocks,p),p)
!!$  !enddo
!!$!  call Driver_abortFlash("DONE 111")
!!$
!!$  ! Allocate data arrays:
!!$  allocate(neighProc_blkcnt(CONSTANT_TWO,neighProcsCount)) 
!!$  allocate(neighProc_blkList(MAXBLOCKS,neighProcsCount)); neighProc_blkList = 0
!!$  
!!$
!!$  ! Send and receive blockCount_set among list of neighProcs => fill neighProc_blkcnt
!!$  recvs  = MPI_REQUEST_NULL
!!$  TAG    = CONSTANT_ONE
!!$  if (neighProcsCount > 0) then 
!!$     do p=1,neighProcsCount
!!$        call MPI_Irecv(neighProc_blkcnt(1:2,p),2,FLASH_INTEGER,neighProcsList(p),TAG,&
!!$                       gr_meshComm,recvs(p),ierr)
!!$     enddo
!!$     vecBuff(1) = blockCount_start_idx
!!$     vecBuff(2) = blockCount_set
!!$     do p=1,neighProcsCount
!!$        call MPI_Ssend(vecBuff(1:2),2,FLASH_INTEGER,neighProcsList(p),TAG,&
!!$                       gr_meshComm,ierr) !sends(p)
!!$     enddo
!!$  end if
!!$
!!$  call MPI_WaitAll(neighProcsCount, recvs, sts, ierr)
!!$  if (ierr /= MPI_SUCCESS) then
!!$     call Driver_abortFlash("Send MPI_Waitall error for neighProc_blkcnt")
!!$  endif
!!$
!!$  ! Send and receive blockList_set among list of neighProcs  => fill neighProc_blkidx
!!$  ! With these we find in lo
!!$  recvs  = MPI_REQUEST_NULL
!!$  TAG    = CONSTANT_TWO
!!$  if (neighProcsCount > 0) then 
!!$     do p=1,neighProcsCount
!!$        neighProc_blocks = neighProc_blkcnt(2,p)
!!$        call MPI_Irecv(neighProc_blkList(1:neighProc_blocks,p),neighProc_blocks,FLASH_INTEGER,&
!!$                       neighProcsList(p),TAG,gr_meshComm,recvs(p),ierr)
!!$     enddo
!!$     do p=1,neighProcsCount
!!$        call MPI_Ssend(blockList_set(1:blockCount_set),blockCount_set,&
!!$                       FLASH_INTEGER,neighProcsList(p),TAG,gr_meshComm,ierr) !sends(p)
!!$     enddo
!!$  end if
!!$
!!$  call MPI_WaitAll(neighProcsCount, recvs, sts, ierr)
!!$  if (ierr /= MPI_SUCCESS) then
!!$     call Driver_abortFlash("Send MPI_Waitall error for neighProc_blkidx")
!!$  endif
!!$
!!$  do p=1,neighProcsCount
!!$     neighProc_blocks = neighProc_blkcnt(2,p)
!!$     write(*,*) 'Mype=',gr_meshMe,neighProcsList(p), &
!!$                        neighProc_blkList(1:neighProc_blocks,p)
!!$  enddo
!!$  call Driver_abortFlash('Done Finding Neighbor Procs')
!!$
!!$  ! Now build Local to global block indexes of gr_meshMe and neighProcs.
!!$  ! Make it simple, stupid!!
!!$  ! gr_meshMe is number 1 on the list
!!$  allocate(neighProc_blkList_idx(MAXBLOCKS,neighProcsCount)); neighProc_blkList_idx = 0
!!$
!!$  do p=1,neighProcsCount
!!$
!!$     neighProc_idx    = neighProc_blkcnt(1,p) 
!!$     neighProc_blocks = neighProc_blkcnt(2,p)
!!$
!!$     do b = 1,neighProc_blocks
!!$
!!$        blockID = neighProc_blkList(b,p)
!!$        neighProc_blkList_idx(blockID,p) = neighProc_idx + b
!!$
!!$     enddo
!!$  enddo
!!$
!!$ 
!!$  do i = 1,gr_meshNumProcs
!!$     if (gr_meshMe .eq. (i-1)) then
!!$        do p=1,neighProcsCount
!!$           neighProc_blocks = neighProc_blkcnt(2,p)
!!$           write(*,*) 'Mype=',gr_meshMe,neighProcsList(p), &
!!$                neighProc_blkList_idx(neighProc_blkList(1:neighProc_blocks,p),p)
!!$        enddo
!!$     endif
!!$     call mpi_barrier(gr_meshComm,ierr)
!!$  enddo
!!$  call Driver_abortFlash('Done Finding Neighbor Procs')
!!$
!!$  ! Build Supermatrix Matrix A in Parallel:



  ! Call Poisson Solver
  if (dr_globalMe .eq. 0) write(*,*) 'Into Grid Solve Poisson ..'

  poisfact = 1.0 
  call Grid_solvePoisson(VPHI_VAR, VSRC_VAR, bc_types, bc_values, poisfact)


  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (dr_globalMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2),count_rate)
      ET=REAL(TA(2)-TA(1),8)/count_rate
      write(*,*) 'Poisson Solver time =',ET
   endif


  

  call Timers_stop("evolution")
  call Logfile_stamp(dr_globalMe, 'Exiting evolution loop' , '[Driver_evolveFlash]')

  ! Mean Phi
  call gr_findMean(VPHI_VAR,2,.false.,meanPHI)
  call gr_findMean(VANL_VAR,2,.false.,meanANL)

  write(*,*) 'Mean PHI=',meanPHI
  write(*,*) 'Mean ANL=',meanANL


  !! Errors.
  call Grid_getLocalNumBlks(localNumBlocks)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Check error in the solution:
  L2_err = 0.
  blkpoints = 0
  Linf_err = 0.
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   =       &
     abs(  solnData(VPHI_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   - &
           meanPHI                                                         - &
           solnData(VANL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   )


#ifdef FLASH_GRID_PARAMESH
     if (lrefine(blockID) .eq. lrefine_max) then
#endif
 
     blkpoints = blkpoints + &
                 (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
                 (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
                 (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     ! L2 norm of error:
     L2_err = L2_err + sum( solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) ))
     
#ifdef FLASH_GRID_PARAMESH
     endif
#endif

     ! Release Pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  enddo  


  ! Sum processors meanPHI
  meanPHIaux = meanPHI
  call MPI_Allreduce(meanPHIaux, meanPHI, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Sum processors meanANL
  meanANLaux = meanANL
  call MPI_Allreduce(meanANLaux, meanANL, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Sum processors L2 norm of error squared
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Compute L2 norm for whole domain
  L2_err = sqrt(1./real(blkpoints) * L2_err)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( dr_globalMe .eq. 0) then
     write(*,*) 'Mean Anl, Num=',meanANL,meanPHI     
     write(*,*) 'L2 error = ',L2_err
     write(*,*) 'Linf error = ',Linf_err
     write(*,*) '############################################'

  endif

  ! Export to Tecplot:
  call outtotecplot(dr_globalMe,0.0,1.,1,0,0.0,blockList,blockCount,0)



!  if(.NOT.endRun) call IO_outputFinal(dr_myPE, dr_numProcs)
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



