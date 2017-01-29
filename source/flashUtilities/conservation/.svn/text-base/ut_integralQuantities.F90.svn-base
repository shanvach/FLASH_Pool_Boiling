#include "constants.h"
#include "Flash.h"

!!REORDER(4):solnData

#define CHECK_CONSERVATION

module ut_integralQuantities

  implicit none
  real, save :: energy_t1, energy_t2, mass_t1, mass_t2

#ifdef CHECK_CONSERVATION
  logical, parameter :: isStubFunction = .false.
#else
  logical, parameter :: isStubFunction = .true.
#endif

  private
  public :: save_integral_quantities, verify_conservation


contains

  subroutine save_integral_quantities(firstOrSecond, isSolnInPrimitiveForm)
    use Grid_interface, ONLY : Grid_getListOfBlocks, &
         Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
         Grid_releaseBlkPtr, Grid_restrictAllLevels
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none
    include "Flash_mpi.h"

    integer, intent(in) :: firstOrSecond
    logical, intent(in) :: isSolnInPrimitiveForm

    integer :: lb, count
    integer :: error
    integer :: blockList(MAXBLOCKS)
    integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

    integer, parameter ::  nGlobalSum = 2 ! Number of globally-summed quantities
    real :: gsum(nGlobalSum) !Global summed quantities
    real :: lsum(nGlobalSum)

    integer :: i, j, k
    real :: dvol
    real, dimension(:,:,:,:), POINTER :: solnData

    integer :: point(MDIM)

    if (isStubFunction) then
       return
    end if


    if (firstOrSecond /= 1 .and. firstOrSecond /= 2) then
       call Driver_abortFlash("[save_integral_quantities]: Invalid argument")
    end if


    ! Sum quantities over all locally held leaf-node blocks.
    gsum = 0.
    lsum = 0.

    call Grid_restrictAllLevels() 
    call Grid_getListOfBlocks(REFINEMENT, blockList, count, refinementLevel=1)

    do lb = 1, count
       !get the index limits of the block
       call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

       ! get a pointer to the current block of data
       call Grid_getBlkPtr(blockList(lb), solnData)

       ! Sum contributions from the indicated blkLimits of cells.
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                point(IAXIS) = i
                point(JAXIS) = j
                point(KAXIS) = k

                !! Get the cell volume for a single cell
                call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)

                ! mass   
#ifdef DENS_VAR
                lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k) * dvol
#endif


                ! total energy
#ifdef ENER_VAR
                if (isSolnInPrimitiveForm) then
                   lsum(2) = lsum(2) + solnData(ENER_VAR,i,j,k) * & 
                        solnData(DENS_VAR,i,j,k)*dvol
                else
                   lsum(2) = lsum(2) + solnData(ENER_VAR,i,j,k) * dvol
                end if
#endif


             enddo
          enddo
       enddo
       call Grid_releaseBlkPtr(blockList(lb), solnData)

    enddo



    ! Now the MASTER_PE sums the local contributions from all of
    ! the processors and writes the total to a file.
    call MPI_AllReduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
         &               MPI_COMM_WORLD, error)


    call MPI_Barrier (MPI_COMM_WORLD, error)


    if (firstOrSecond == 1) then
       mass_t1 = gsum(1)
       energy_t1 = gsum(2)
    else
       mass_t2 = gsum(1)
       energy_t2 = gsum(2)
    end if

    !=============================================================================

  end subroutine save_integral_quantities


  subroutine verify_conservation(myPE, callingLocation)
    implicit none
    integer, intent(in) :: myPE
    character(len=*), intent(IN) :: callingLocation
    character(len=2000) :: userMsg
    character (len=*), parameter  :: logStr = &
         "(a, / a,3es25.18, / a,3es25.18)"
    real, parameter :: tol = 1E-12
    real :: err
    logical :: inErrorState

    if (isStubFunction) then
       return
    end if
    inErrorState = .false.

    !Check conservation of mass.
    err = relative_error(mass_t1, mass_t2)
    if (err > tol) then
       inErrorState = .true.
    end if

    !Check conservation of total energy.
    err = relative_error(energy_t1, energy_t2)
    if (err > tol) then
       inErrorState = .true.
    end if

    if (inErrorState) then
       userMsg = "Non conservation in section "//trim(callingLocation)
       if (myPE == 0) then
          write(6,logStr) trim(userMsg), &
          "Mass:", mass_t1, mass_t2, mass_t2-mass_t1, &
          "Energy:", energy_t1, energy_t2, energy_t2-energy_t1
       end if
       !call Driver_abortFlash(trim(userMsg))
    end if

  end subroutine verify_conservation


  function relative_error(term1, term2) result(err)
    implicit none
    real :: err
    real, intent(IN) :: term1, term2
    err = abs((term1 - term2) / term1)
  end function relative_error

end module ut_integralQuantities
