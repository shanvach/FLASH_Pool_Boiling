!!****if* source/Simulation/SimulationMain/unitTest/RayPath/Driver_evolveFlash
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
!!  This is the global driver specifically for tesing pfft
!!  It opens a file for each of the processor that it is running
!!  on and call the Grid_unitTest with the file unit. The Grid_unitTest
!!  functions carries out the testing of the unit, and reports success or
!!  failure through a logical argument "perfect". Upon return from 
!!  Grid_unitTest, the Driver_evolveFlash function makes appropriate notification 
!!  in the file which the test suite can parse to determine if the test was
!!  successful.
!!
!!
!!***

subroutine Driver_evolveFlash()

#include "constants.h"
#include "Pfft.h"  
#include "Flash.h"


  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use GridRayTrace_interface, ONLY :Grid_unitTest 

  implicit none

  logical ::  perfect = .true., needMap = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer :: blockID=1
  integer :: i,j,k,nbegin,nend,mbegin,mend,fullSize,temp

  real :: err
  real,dimension(:,:,:,:),pointer :: solnData

  integer :: iTimes, error, ii, jj, kk

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

  call sim_analyticEndpoint()
  call Grid_unitTest( fileUnit, perfect)
  
  
  print*,'the error is',err
  if(err < 1e-05) then
     write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)

  return

end subroutine Driver_evolveFlash


!realSpacing = 1 for real data.
!realSpacing = 2 for complex data.
subroutine myIO_globalArray(localArray, realSpacing, globalSize)

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real, pointer, dimension(:) :: localArray
  integer, intent(IN) :: realSpacing
  integer, dimension(MDIM), intent(IN) :: globalSize
  integer :: dataStart, dataEnd, ierr, i, originalProc, error
  real, allocatable, dimension(:) :: globalArrayTmp, globalArray

  if(NDIM > 2) then
     if(dr_globalMe == 0) then
        print *, "[Driver_evolveFlash] Can't print information in 3D"
     end if
     return
  end if


  allocate(globalArray(size(localArray) * dr_globalNumProcs), &
       globalArrayTmp(size(localArray) * dr_globalNumProcs), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if

  dataStart = (dr_globalMe * size(localArray,1)) + 1
  dataEnd = dataStart + size(localArray,1) - 1

  globalArrayTmp(:) = 0.0
  globalArrayTmp(dataStart:dataEnd) = localArray(:)

  call MPI_ALLREDUCE(globalArrayTmp, globalArray, size(localArray,1)*dr_globalNumProcs, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)


  if(dr_globalMe == 0) then
     do i = 1, globalSize(JAXIS)
        originalProc = (i-1) / (globalSize(JAXIS) / dr_globalNumProcs)

        dataStart = ((i-1) * globalSize(IAXIS) * realSpacing) + 1
        dataEnd= i * globalSize(IAXIS) * realSpacing
        print *, "[Driver_evolveFlash] Processor", originalProc, ", global array (", &
             dataStart, ":", dataEnd, ":", realSpacing, ") = ", & 
             int(globalArray(dataStart:dataEnd:realSpacing))
     end do
  end if

  deallocate(globalArray, globalArrayTmp, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be deallocated!")
  end if

end subroutine myIO_globalArray



subroutine testTranspose(inputArray, globalSize)

#include "constants.h"
#include "Pfft.h"  
#include "Flash.h"

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use gr_pfftData, ONLY : pfft_midLen,pfft_t1Len, pfft_transformType, &
       pfft_comm, pfft_me, pfft_procGrid,&
       pfft_work1
  use gr_pfftInterface, ONLY : gr_pfftTranspose


  implicit none
  real, dimension(:), intent(IN) :: inputArray
  integer, dimension(MDIM), intent(IN) :: globalSize
  integer, dimension(MDIM) :: temp_t1Len, temp_midLen
  real, allocatable, dimension(:) ::  tempRealArray, tempComplexArray
  integer :: i, j

  interface
     subroutine myIO_globalArray(localArray, realSpacing, globalSize)       
       implicit none
       real, dimension(:) :: localArray
       integer, intent(IN) :: realSpacing
       integer, dimension(MDIM), intent(IN) :: globalSize
     end subroutine myIO_globalArray
  end interface


  !------------------------------------------------------
  ! TEST WITH A REAL TO REAL TRANSPOSE
  if(dr_globalMe == 0) then
     print *, "Testing real to real transpose"
     print *, ""
  end if
  allocate(tempRealArray(size(inputArray)))
  tempRealArray(:) = -1
  pfft_work1(:)= -1

  call myIO_globalArray(inputArray, 1, globalSize)

  tempRealArray(:) = inputArray(:)

  call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,tempRealArray,&
       pfft_work1,pfft_t1Len,pfft_midLen,&
       pfft_procGrid(JAXIS),pfft_comm(JAXIS))

  tempRealArray(:) = -1
  tempRealArray(:) = pfft_work1(1:size(inputArray))

  if(dr_globalMe == 0) then
     print *, ""
  end if
  call myIO_globalArray(tempRealArray, 1, globalSize)

  deallocate(tempRealArray)


  !------------------------------------------------------
  ! TEST WITH A COMPLEX TO COMPLEX TRANSPOSE
  if(dr_globalMe == 0) then
     print *, "Testing real to complex, complex to complex transpose"
     print *, ""
  end if

  allocate(tempComplexArray(size(inputArray) * 2))
  tempComplexArray(:) = -1
  pfft_work1(:)= -1

  do i = 1, size(inputArray)
     j = (i*2) - 1
     tempComplexArray(j) = inputArray(i)  !Fill in the real component only.
  end do

  call myIO_globalArray(tempComplexArray, 2, globalSize)

  call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_COMPLEX,tempComplexArray,&
       pfft_work1,pfft_t1Len,pfft_midLen,&
       pfft_procGrid(JAXIS),pfft_comm(JAXIS))

  tempComplexArray(:) = -1
  tempComplexArray(:) = pfft_work1(:)

  if(dr_globalMe == 0) then
     print *, ""
  end if
  call myIO_globalArray(tempComplexArray, 2, globalSize)

  deallocate(tempComplexArray)

  !print *, "Processor", dr_globalMe, "finished in testTranspose"

end subroutine testTranspose



!!$  nbegin=1
!!$  mbegin= blkLimits(LOW,IAXIS)
!!$  mend=blkLimits(HIGH,IAXIS)
!!$  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!$     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!!$        nend=nbegin+localSize(IAXIS)-1
!!$        inArray(nbegin:nend)=solnData(1,mbegin:mend,j,k)
!!$        nbegin=nend+1
!!$     end do
!!$  end do
