!!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/Driver_evolveFlash
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
!!  This is the global driver specifically for tesing the transposing of
!!  data within the pfft unit.  It opens a file for each of the processors 
!!  that it is running and calls Grid_unitTest with a file unit. Grid_unitTest
!!  carries out the testing of the unit, and reports success or
!!  failure through a logical argument "perfect". Upon return from 
!!  Grid_unitTest, the Driver_evolveFlash function makes appropriate notification 
!!  in the file which the test suite can parse to determine if the test was
!!  successful.
!!
!!***

#include "constants.h"
#include "Pfft.h"  
#include "Flash.h"

subroutine Driver_evolveFlash()
  use Driver_data, ONLY: dr_globalMe
  use Grid_interface, ONLY : Grid_pfftMapToInput, &
       Grid_pfftMapFromOutput, Grid_unitTest
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftData, ONLY : pfft_inLen, pfft_midLen, pfft_outLen, &
       pfft_work1, pfft_work2, pfft_globalLen, &
       pfft_comm, pfft_procGrid
  use gr_pfftInterface, ONLY : gr_pfftTranspose
  use Simulation_data, ONLY : sim_outputGridData
  use sim_module, ONLY : sim_printLocalPencilData

#ifdef FLASH_GRID_UG
  use Grid_data, ONLY : gr_axisNumProcs
#endif

  implicit none

  real, dimension(:), allocatable :: inArray, outArray
  integer, dimension(MDIM) :: origDimOrder
  integer, dimension(4) :: prNum
  integer, parameter :: fileUnit = 2
  integer :: i, temp, error
  character (len=20) :: fileName
  logical :: perfect

  perfect = .true.
  temp=dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
       char(48+prNum(2))//char(48+prNum(1))
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe


#ifdef FLASH_GRID_UG
  !Limitation with the mapping routines.  Something for me to investigate
  !when I get time.  This condition stops the cases that break, but I need
  !to think more about this to be sure that this is the reason for the failures.
  do i = 1, NDIM
     if (mod(pfft_globalLen(i),gr_axisNumProcs(IAXIS)) /= 0) then
        call Driver_abortFlash("[Driver_evolveFlash]: Grid cannot be handled.")
     end if
  end do
#endif


  !---------------------------------------------------------------------
  ! START SUMMARY INFO PRINTS.
  ! We are not concerned with any FFT related stuff so only need
  ! three shape arrays: pfft_inLen, pfft_midLen, pfft_outLen.
  !---------------------------------------------------------------------
  if(dr_globalMe == 0) then
     write(6,*) ""
     if (NDIM == 1) then

        write(6,'(a,i4,a,i5,a)') " The global grid shape is", &
             pfft_globalLen(1:NDIM), &
             "  (product is", product(pfft_globalLen(1:NDIM)), &
             " and grid orientation is [IAXIS])"
        write(6,'(a,i4)') " The pencil processor grid is", &
             pfft_procGrid(1:NDIM)
        write(6,'(a,i4,a,i5,a)') " The initial local grid shape is", &
             pfft_inLen(1:NDIM), &
             "  (product is", product(pfft_inLen(1:NDIM)), &
             " and grid orientation is [IAXIS])"

     else if (NDIM == 2) then

        write(6,'(a,2i4,a,i5,a)') " The global grid shape is", &
             pfft_globalLen(1:NDIM), &
             "  (product is", product(pfft_globalLen(1:NDIM)), &
             " and grid orientation is [IAXIS,JAXIS])"
        write(6,'(a,2i4)') " The pencil processor grid is", &
             pfft_procGrid(1:NDIM)
        write(6,'(a,2i4,a,i5,a)') " The initial local grid shape is", &
             pfft_inLen(1:NDIM), &
             "  (product is", product(pfft_inLen(1:NDIM)), &
             " and grid orientation is [IAXIS,JAXIS])"
        write(6,'(a,2i4,a,i5,a)') " After transpose 1 the local grid shape is", &
             pfft_outLen(1:NDIM), &
             "  (product is", product(pfft_outLen(1:NDIM)), &
             " and grid orientation is [JAXIS,IAXIS])"

     else if (NDIM == 3) then

        write(6,'(a,3i4,a,i5,a)') " The global grid shape is", &
             pfft_globalLen(1:NDIM), &
             "  (product is", product(pfft_globalLen(1:NDIM)), &
             " and grid orientation is [IAXIS,JAXIS,KAXIS])"
        write(6,'(a,3i4)') " The pencil processor grid is", &
             pfft_procGrid(1:NDIM)
        write(6,'(a,3i4,a,i5,a)') " The initial local grid shape is", &
             pfft_inLen(1:NDIM), &
             "  (product is", product(pfft_inLen(1:NDIM)), &
             " and grid orientation is [IAXIS,JAXIS,KAXIS])"
        write(6,'(a,3i4,a,i5,a)') " After transpose 1 the local grid shape is", &
             pfft_midLen(1:NDIM), &
             "  (product is", product(pfft_midLen(1:NDIM)), &
             " and grid orientation is [JAXIS,KAXIS,IAXIS])"
        write(6,'(a,3i4,a,i5,a)') " After transpose 2 the local grid shape is", &
             pfft_outLen(1:NDIM), &
             "  (product is", product(pfft_outLen(1:NDIM)), &
             " and grid orientation is [KAXIS,IAXIS,JAXIS])"

     end if
     write(6,*) ""
  end if
  !---------------------------------------------------------------------
  ! END SUMMARY INFO PRINTS.
  !---------------------------------------------------------------------

  !Use the oversized pfft_work1 and pfft_work2 in all transpose operations so
  !transposes work for a wider range of grid configurations (see check-in 8269).
  pfft_work1(:)= -1.0; pfft_work2(:)= -1.0

  allocate(inArray(product(pfft_inLen(1:NDIM))), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Driver_evolveFlash]: inArray cannot be allocated!")
  end if
  allocate(outArray(product(pfft_outLen(1:NDIM))), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Driver_evolveFlash]: outArray cannot be allocated!")
  end if


  !origDimOrder is just [IAXIS,JAXIS,KAXIS] (in 3d).
  do i = 1, MDIM
     if (i <= NDIM) then
        origDimOrder(i) = i
     else
        origDimOrder(i) = 0
     end if
  end do


  !Map from FLASH grid DENS variable to a pencil grid.
  call Grid_pfftMapToInput(DENS_VAR, inArray)

  !Copy the original input array into our oversized work array.
  pfft_work1(1:size(inArray)) = inArray(:) 

  !Print the original pencil grid.
  call sim_printLocalPencilData(pfft_work1, origDimOrder, &
       pfft_inLen,  "output_initial_")

  !---------------------------------------------------------------------
  ! START TRANSPOSES.
  ! Perform the forward and back transposes, and create an output file
  ! showing the current local grid and portion of global grid belonging
  ! to each processor after each forward transpose.
  !---------------------------------------------------------------------
  if (NDIM == 1) then
     !NO-OP
  else if (NDIM == 2) then
     !--------------------------- Forward ------------------------------
     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,pfft_work1,&
          pfft_work2,pfft_inLen,pfft_outLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     pfft_work1 = -1.0

     call sim_printLocalPencilData(pfft_work2, (/JAXIS,IAXIS,0/), &
          pfft_outLen,  "output_transpose_1_")
     !------------------------------------------------------------------

     !--------------------------- Inverse ------------------------------
     call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work2,&
          pfft_work1,pfft_outLen,pfft_inLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     pfft_work2 = -1.0
     !------------------------------------------------------------------
  else if (NDIM == 3) then
     !--------------------------- Forward ------------------------------
     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,pfft_work1,&
          pfft_work2,pfft_inLen,pfft_midLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     pfft_work1 = -1.0

     call sim_printLocalPencilData(pfft_work2, (/JAXIS,KAXIS,IAXIS/), &
          pfft_midLen, "output_transpose_1_")

     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,pfft_work2,&
          pfft_work1,pfft_midLen,pfft_outLen,&
          pfft_procGrid(KAXIS),pfft_comm(KAXIS))
     pfft_work2 = -1.0
     !------------------------------------------------------------------

     call sim_printLocalPencilData(pfft_work1, (/KAXIS,IAXIS,JAXIS/), &
          pfft_outLen,  "output_transpose_2_")

     !--------------------------- Inverse ------------------------------
     call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work1,&
          pfft_work2,pfft_outLen,pfft_midLen,&
          pfft_procGrid(KAXIS),pfft_comm(KAXIS))
     pfft_work1 = -1.0
     call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work2,&
          pfft_work1,pfft_midLen,pfft_inLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     pfft_work2 = -1.0
     !------------------------------------------------------------------
  else
     call Driver_abortFlash("[Driver_evolveFlash]: Dimensionality not handled")
  end if
  !---------------------------------------------------------------------
  ! END TRANSPOSES.
  !---------------------------------------------------------------------


  !Print the final pencil grid (should be the same as the orginal pencil grid).
  call sim_printLocalPencilData(pfft_work1, origDimOrder, &
       pfft_inLen,  "output_final_")

  !Copy from oversized array into our output array.
  outArray(:) = pfft_work1(1:size(outArray))

  !Map from pencil grid to FLASH grid PDEN_VAR variable.
  call Grid_pfftMapFromOutput(PDEN_VAR, outArray)


  deallocate(inArray, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Driver_evolveFlash]: inArray cannot be deallocated!")
  end if
  deallocate(outArray, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Driver_evolveFlash]: outArray cannot be deallocated!")
  end if


  !Grid_unitTest checks that DENS_VAR contains the same as PDEN_VAR.
  call Grid_unitTest(fileUnit,perfect) 
  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
    if(dr_globalMe == 0) then
       write(6,*) "Success on processor 0 - see output files"
    end if
  else
    write(fileUnit,'("test failed.")')
    if(dr_globalMe == 0) then
       write(6,*) "Failure on processor 0 - see output files"
    end if
  endif
  close(fileUnit)
end subroutine Driver_evolveFlash
