!*******************************************************************************

!  Routine:     poisson()

!  Description: Driver routine for the SUPERLU Poisson solver.  

!  Parameters:  isrc            Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               isoln           Index for solution array.  The solution is
!                               written directly into this variable.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: The following have been added to the interface only to get this
! version of the paramesh3.x poisson solver to work.  These have no effect. 
!               bc_types(6)     Boundary condition types array:
!                                 GRID_PDE_BND_PERIODIC
!                                 GRID_PDE_BND_DIRICHLET
!                                 GRID_PDE_BND_NEUMANN
!
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z  6 = +z
!               bc_values(2,6)  Values for dirichlet and neumann boundary
!                               conditions. All 0 for now (homogeneous). 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine Grid_SolvePoisson (isoln, isrc, bc_types, bc_values, poisfact)

  !===============================================================================
#include "Flash.h"
#include "Superlu.h"

  use superlu_common

  use Grid_data, only : gr_meshMe,gr_meshComm,gr_meshNumProcs

  use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED, &
                                GRID_PDE_BND_PERIODIC, &
                                GRID_PDE_BND_DIRICHLET,&
                                GRID_PDE_BND_NEUMANN, &
                                Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr, &
                                Grid_releaseBlkPtr

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash


  use slu_interface, only : gr_slugetSubsetBlks,       &
                            gr_slugetSubsetStartidx,   &
                            gr_slugetSubsetNeighProcs, &
                            gr_slugetSubsetNeighLists, &
                            gr_sluBuildAdist,          &
                            gr_sluBuildrhsAdist,       &
                            gr_sluDumpSolndist

  use superlu_mod

  implicit none
#include "Flash_mpi.h"
  integer :: isoln, isrc
  integer, dimension(6)   :: bc_types
  real,    dimension(2,6) :: bc_values
  real    :: poisfact

  ! Local Variables:
  integer :: blockList_set(MAXBLOCKS),blockCount_set
  integer :: slu_flag , diagblk_flag

  ! Super LU Variables:
  integer(superlu_ptr) :: grid
  integer(superlu_ptr) :: options
  integer(superlu_ptr) :: ScalePermstruct
  integer(superlu_ptr) :: LUstruct
  integer(superlu_ptr) :: SOLVEstruct
  integer(superlu_ptr) :: A
  integer(superlu_ptr) :: stat
  integer, save :: nprow, npcol
  integer :: nprowK
  integer :: nrhs, init, info,ierr
  integer :: n
  real    :: berr
  real, save, allocatable, dimension(:) :: rhsA


  logical, save :: firstcall = .true.
  integer :: i,p, neighProc_blocks

  integer :: ldb
  integer*4 :: ldb4


  if (firstcall .or. gr_sluGridChanged) then !!! If first time or grid changed

     if (gr_meshMe .eq. 0) print*,"Setting up SuperLU data structures."

     !KPD  - Initialize solving core count
     nprowK = 0

     ! Call routine to define blocks that belong to specified subset
     slu_flag = MAXUNIF_LEVEL
     call gr_slugetSubsetBlks(slu_flag,blockList_set,blockCount_set)

     ! Obtain the blocks Global Numeration index
     call gr_slugetSubsetStartidx(blockCount_set) 

     ! Find neighbor Blocks: Fills neighProcsCount and neighProcsList of superlu_common
     diagblk_flag = DIAGBLKS_NO
     call gr_slugetSubsetNeighProcs(blockList_set,blockCount_set,diagblk_flag)

     ! Get Neighbor Procs-Block indexes, Build Global block lists:
     ! Sets arrays neighProc_blkcnt,neighProc_blkList,neighProc_blkList_idx
     ! defined in superlu_common
     call gr_slugetSubsetNeighLists(blockList_set,blockCount_set)

     ! Initialize the SuperLu Solver vars:
     ! Create Fortran handles for the C structures used in SuperLU_DIST
     call f_create_gridinfo_handle(grid)
     call f_create_options_handle(options)
     call f_create_ScalePerm_handle(ScalePermstruct)
     call f_create_LUstruct_handle(LUstruct)
     call f_create_SOLVEstruct_handle(SOLVEstruct)
     call f_create_SuperMatrix_handle(A)
     call f_create_SuperLUStat_handle(stat)

     ! Initialize the SuperLU_DIST process grid
     !KPD - nprow*npcol MUST equal the number of solving cores!!!
     !nprow = gr_meshNumProcs 
     nprow = 0
     npcol = 1

     if (nloc_dofs .gt. 0) then 
        nprowK = nprowK + 1
     end if

     call MPI_Allreduce(nprowK, nprow, 1, FLASH_INTEGER,&
                        MPI_SUM, MPI_COMM_WORLD, ierr)

     print*,"Final nprow on",gr_meshMe,"is:",nprow

     call f_superlu_gridinit(gr_meshComm, nprow, npcol, grid)

     if (nloc_dofs .eq. 0) then
        print*,"SuperLU skipping core",gr_meshMe
        go to 200
     end if

     ! Here m_glob n_blob should be set.
     ! Build Supermatrix Matrix A in Parallel:
     call gr_sluBuildAdist(A, m_glob, n_glob)

     ! Apply LU decomposition in parallel:
     ! Set the default input options
     call f_set_default_options(options)

     ! Modify one or more options
     !call set_superlu_options(options,ColPerm=NATURAL)
     !call set_superlu_options(options,RowPerm=NATURAL)
     !call set_superlu_options(options,Fact=FACTORED)     

     ! Initialize ScalePermstruct and LUstruct
     call get_SuperMatrix(A,nrow=m_glob,ncol=n_glob)
     call f_ScalePermstructInit(m_glob, n_glob, ScalePermstruct)
     call f_LUstructInit(m_glob, n_glob, LUstruct)
     
     ! Initialize the statistics variables
     call f_PStatInit(stat)
     
     firstcall = .false.
     gr_sluGridChanged = .false.

  end if  !!! End first time or grid changed


  ! Build Right hand side isrc in Parallel:
  nrhs = 1
  allocate(rhsA(nloc_dofs))
  call gr_sluBuildrhsAdist(isrc,poisfact,nloc_dofs,rhsA)


  if (gr_meshMe == 0) write(*,*) 'rhsA(1)=',rhsA(1)

  !KPD - nloc_dofs is the local number of rows
  ! Call the solver:
  call f_pdgssvx(options, A, ScalePermstruct,rhsA, nloc_dofs, nrhs, &
                 grid, LUstruct, SOLVEstruct, berr, stat, info)

  if (gr_meshMe == 0) then
     write (*,*) 'Backward error: ', berr
     write(*,*) 'INFO from f_pdgssvx = ', info
  endif


  if (gr_meshMe == 0) write(*,*) 'rhsA(1) soln=',rhsA(1)

  ! Dump Solution in isoln:  
  call gr_sluDumpSolndist(isoln,nloc_dofs,rhsA)

  ! Deallocate rhsA:
  deallocate(rhsA)  

200 return
!  return

end subroutine Grid_SolvePoisson
