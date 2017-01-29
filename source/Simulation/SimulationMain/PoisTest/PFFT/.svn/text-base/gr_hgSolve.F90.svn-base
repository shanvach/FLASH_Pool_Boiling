!!****if* source/Simulation/SimulationMain/PoisTest/PFFT/gr_hgSolve
!!
!! NAME
!!  gr_hgSolve
!!
!! SYNOPSIS
!!  gr_hgSolve(integer, intent(in) :: gr_iSource,
!!             integer, intent(in) :: gr_iSoln,
!!             integer, intent(in) :: gr_iSls,
!!             integer, intent(in) :: gr_iCorr,
!!             external            :: SolveBlock,
!!             integer, intent(in) :: bndTypes(6),
!!             real, intent(in)    :: src_fact)
!!
!! DESCRIPTION
!!  This is the main Poisson solve routine for the Huang & Greengard
!!  (2000, SIAM J. Sci. Comput., 21, 1551) algorithm.  This routine 
!!  defines the multigrid cycle as expressed in the article.
!!
!! coarse ^
!!      o---->o                o---->o
!!      ^    s \               ^    s \ c
!!     r|     o o              |r    o o o
!!     e|      l \             |e     l \ r
!!     s|       v o            |s      v o r
!!     t|        e \           |t       e \ e
!!     r|           o          |r          o c
!!     i|            \         |i           \ t
!!     c|             o        |c            o
!!     t|              \       |t             \
!!  src |               o----->o gr_iSls       o------>...until |r| < rtol
!! fine |______________take residual____________take new residual__>
!!
!! ARGUMENTS
!!  gr_iSource - the source variable
!!  gr_iSoln - the solution variable
!!  gr_iSls       - the residual variable
!!  gr_iCorr      - the correction variable
!!  SolveBlock - the name of the function used for individual block solves
!!  bndTypes   - the BC types as defined in Multigrid.h
!!  src_fact   - The factor by which the solution is multiplied
!!  
!! EXAMPLE
!!
!!  gr_hgSolve(idens, gpot, hgwk1, hgwk2, gr_hgPoissonSolveBlock,
!!             (/(MG_BND_DIRICHLET,i=1,6)/), newton)
!!
!!  This call to this equation takes the density variable and puts
!!  a 0-boundary-valued solution to poisson's equation in gpot, using
!!  hgwk1 and hgwk2 as work variables (for the residual and correction)
!!  
!!  USED BY
!!
!!   Primarily Grid_solvePoisson.
!!
!!  SEE ALSO
!!   
!!   Grid_solvePoisson
!!   gr_hgSolveLevel
!!    
!!***

#include "constants.h"
#include "Flash.h"


!Instead of doing the multigrid solve, we just move data to pfft grid and back again.
!This is only tests the data movement to and from pfft grid when starting 
!with a mixed refinement PARAMESH grid. 

!We test data movement at all possible refinement levels.


subroutine gr_hgSolve(gr_iSource, gr_iSoln, gr_iSls, gr_iCorr, SolveBlock, bndTypes, src_fact, dt, chi)

  !=================================================================

  use Grid_data, ONLY: gr_myPE, gr_nblockX, gr_nblockY, gr_nblockZ
  use Logfile_interface, ONLY: Logfile_open, Logfile_close
  use tree, ONLY : lrefine
  use Grid_interface, ONLY : Grid_getMaxCommonRefinement, Grid_getListOfBlocks, &
       Grid_pfftMapToInput, Grid_pfftMapFromOutput
  use gr_hgInterface, ONLY : gr_hgPfftInit, gr_hgPfftFinalize
  use gr_hgPfftData, ONLY : gr_hgPfftInArray
  implicit none

  include "Flash_mpi.h"

  interface
     subroutine LabelBlocksAtRefineLevel(blkList, blkCount, iSource, logUnit)
       implicit none
       integer, dimension(MAXBLOCKS), intent(IN) :: blkList
       integer, intent(IN) :: blkCount, iSource, logUnit
     end subroutine LabelBlocksAtRefineLevel

     subroutine CheckSourceAndSoln(iSource, iSoln, blkList, blkCount, logUnit)
       implicit none
       integer, intent(IN) :: iSource, iSoln
       integer, dimension(MAXBLOCKS), intent(IN) :: blkList
       integer, intent(IN) :: blkCount, logUnit
     end subroutine CheckSourceAndSoln

     subroutine PlaceGarbageInGrid(iSource, iSoln)
       implicit none
       integer, intent(IN) :: iSource, iSoln
     end subroutine PlaceGarbageInGrid
  end interface

  integer, intent(in) :: gr_iSource, gr_iSoln, gr_iSls, gr_iCorr
  integer, intent(in) :: bndTypes(6)
  real, intent(IN)    :: src_fact

  real, intent(IN),OPTIONAL :: dt, chi
  
  external               SolveBlock

  logical, parameter :: logUnitLocal = .true.
  integer, dimension(MAXBLOCKS) :: blkList
  integer             :: comm, globalTopLevel, blkCount, ierr, refineLevel
  integer :: logUnit, lowerBound, upperBound
  integer, dimension(MDIM) :: initialBlocks

  initialBlocks(1:MDIM) = 1
  initialBlocks(IAXIS) = gr_nblockX
  if (NDIM >= 2) initialBlocks(JAXIS) = gr_nblockY
  if (NDIM == 3) initialBlocks(KAXIS) = gr_nblockZ


  !Determine the maximum common refinement level in the grid.
  comm = MPI_COMM_WORLD
  call Grid_getMaxCommonRefinement(comm, globalTopLevel)


  !Just to see the range of leaf block refinement levels.
  !call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  !print *, "Processor", gr_myPE, "LEAF min refinement:", &
  !     minval(lrefine(blkList(1:blkCount))), "LEAF max refinement:", &
  !     maxval(lrefine(blkList(1:blkCount)))


  if (gr_myPE == MASTER_PE) then
     print *, "Complete block coverage up to refinement level:", &
          globalTopLevel, " NDIM:", NDIM, &
          " Number of initial blocks(1:NDIM):", initialBlocks(1:NDIM)
  end if

  call Logfile_open(logUnit,logUnitLocal)


  !These are the only levels at which there are blocks at 
  !this refinement that completely cover the domain.
  do refineLevel = globalTopLevel, 1, -1

     if (gr_myPE == MASTER_PE) then
        print *, ""
        print *, "Testing data movement at refinement level:", refineLevel, " Contains:", &
             product(initialBlocks(1:NDIM)) * ((2**(refineLevel-1))**NDIM), "block(s)."
     end if

     !Before we start ensure that there is garbage in 
     !the specified grid elements, and distinct garbage too!
     call PlaceGarbageInGrid(gr_iSource, gr_iSoln)


     call Grid_getListOfBlocks(REFINEMENT, blkList, blkCount, refineLevel)
     !print *, "Processor", gr_myPE, "has", blkCount, "blocks at level", refineLevel
     write(logUnit,*) ""
     write(logUnit,*) "----------------------------------------------------"
     write(logUnit,*) "Refinement level:", refineLevel, "at which I own", &
          blkCount, "blocks."


     !Initialise the gr_iSource grid variable at refinement 
     !level "refineLevel" with a globally unique block ID.
     call LabelBlocksAtRefineLevel(blkList, blkCount, gr_iSource, logUnit)


     !Create arrays in PFFT space that are sized appropriately for 
     !each chosen refinement level.
     call gr_hgPfftInit(refineLevel)
     gr_hgPfftInArray(:) = -7.5 !fill with garbage.


     !Map from gr_iSource grid variable to PFFT array.
     !gr_iSource contains the global block identifier of each block.
     call Grid_pfftMapToInput(gr_iSource, gr_hgPfftInArray)


     !Output PFFT array to logfile.
     if (blkCount > 0) then
        write(logUnit,*) "*** Printing PFFT array:"

        !Last 2 elements are FFT specific and can be ignored.
        lowerBound = lbound(gr_hgPfftInArray,1)
        upperBound = ubound(gr_hgPfftInArray,1)-2
        write(logUnit,*) int(gr_hgPfftInArray(lowerBound:upperBound))
     else
        write(logUnit,*) "No data placed in PFFT array!"
     end if


     !As we do not operate on the array in PFFT space, we
     !just copy the PFFT array to the gr_iSoln grid variable.
     call Grid_pfftMapFromOutput(gr_iSoln, gr_hgPfftInArray)
     

     !If the map back and forth was successful, we will have the
     !same value in gr_iSource and gr_iSoln grid variables.
     call CheckSourceAndSoln(gr_iSource, gr_iSoln, blkList, blkCount, logUnit)


     !Clean up PFFT allocations at this refinement level.
     call gr_hgPfftFinalize()
      
  end do
  call Logfile_close(logUnitLocal)


  !We have tested the data movement now so end the simulation. 
  !Ensure all processes have successfully reached this point before 
  !assuming test was a success.
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (gr_myPE == MASTER_PE) then
     print *, ""
     print *, "TEST WAS SUCCESSFUL AT ALL REFINEMENT LEVELS!!!!"
     print *, "(see process specific log files for details of data movement)"
  end if
  call MPI_Finalize(ierr)
  stop

end subroutine gr_hgSolve



subroutine LabelBlocksAtRefineLevel(blkList, blkCount, iSource, logUnit)

  use Grid_data, ONLY: gr_myPE, gr_numProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_open, Logfile_close
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none

  include "Flash_mpi.h"

  integer, dimension(MAXBLOCKS), intent(IN) :: blkList
  integer, intent(IN) :: blkCount, iSource, logUnit

  real, dimension(:,:,:,:), pointer :: solnData
  integer, dimension(:), allocatable :: globalCountAtRefine, globalCountAtRefineTmp
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: blockID, expectedValue, startBlockID
  integer :: ierr, i, error

  !Define a global block identifier for each block at the requested refinement level.  
  !This is different from the normal block ID which is local to each processor.
  !To define a global block identifier, we create an array on each processor which 
  !has size = no.processors (one element per processor).  Each processor places the 
  !number of blocks it owns at the requested refinement level in its designated 
  !element in the array.  We then perform a global allreduce so that each processor 
  !knows the number of blocks on each processor.  We can use this information to 
  !work out a unique block identifier across all processors.
  !---------------------------------------------------------------------
  allocate(globalCountAtRefine(gr_numProcs), globalCountAtRefineTmp(gr_numProcs), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if

  globalCountAtRefine = 0; globalCountAtRefineTmp = 0;
  globalCountAtRefineTmp(gr_myPE+1) = blkCount

  !MPI_IN_PLACE is an MPI-2 construct, hence we use 2 arrays.
  call MPI_Allreduce(globalCountAtRefineTmp,globalCountAtRefine,gr_numProcs,&
       FLASH_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)  

  if (gr_myPE == 0) then
     startBlockID = 0
  else
     startBlockID = sum(globalCountAtRefine(1:gr_myPE))
  end if


  !Loop over each block at this refinement level, and label the internal grid 
  !points with the global block ID.
  !---------------------------------------------------------------------
  if (blkCount > 0) then
     write(logUnit,*) "*** Printing iSource grid variable:"

     do i = 1, blkCount
        blockID = blkList(i)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        
        expectedValue = startBlockID + i - 1
        solnData(iSource,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = real(expectedValue)

        write(logUnit,*) "local blockID:", blockID, ", iSource grid data:", &
             int(solnData(iSource,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     end do

  else
     write(logUnit,*) "No blocks to print (iSource)"
  end if


  deallocate(globalCountAtRefine, globalCountAtRefineTmp, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be deallocated!")
  end if

end subroutine LabelBlocksAtRefineLevel




subroutine CheckSourceAndSoln(iSource, iSoln, blkList, blkCount, logUnit)

  use Grid_data, ONLY: gr_myPE, gr_numProcs
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: iSource, iSoln
  integer, dimension(MAXBLOCKS), intent(IN) :: blkList
  integer, intent(IN) :: blkCount, logUnit

  real, dimension(:,:,:,:), pointer :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID, expectedValue, i


  if (blkCount > 0) then
     write(logUnit,*) "*** Printing iSoln grid variable:"

     do i = 1, blkCount
        blockID = blkList(i)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        expectedValue = solnData(iSource,blkLimits(LOW,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS))

        write(logUnit,*) "local blockID:", blockID, ", iSoln grid data:", &
             int(solnData(iSoln,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

        if(all(int(solnData(iSoln,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) == expectedValue)) then
        else
           print *, "Failure!!! - Processor", gr_myPE, "block", blockID, &
                "messed up during the map."
           call Driver_abortFlash("FAILURE DURING MAP")
        end if

        call Grid_releaseBlkPtr(blockID,solnData,CENTER)     
     end do

  else
     write(logUnit,*) "No blocks to print (iSoln)"
  end if

end subroutine CheckSourceAndSoln




subroutine PlaceGarbageInGrid(iSource, iSoln)

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none
  integer, intent(IN) :: iSource, iSoln

  real, dimension(:,:,:,:), pointer :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkCount, blockID, i

  call Grid_getListOfBlocks(ALL_BLKS, blkList, blkCount)

  do i = 1, blkCount
     blockID = blkList(i)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     solnData(iSource,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = -2.5
     solnData(iSoln,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = -4.5

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  end do

end subroutine PlaceGarbageInGrid
