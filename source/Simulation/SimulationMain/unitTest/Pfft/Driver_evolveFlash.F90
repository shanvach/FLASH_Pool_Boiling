!!****if* source/Simulation/SimulationMain/unitTest/Pfft/Driver_evolveFlash
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
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits,Grid_getBlkIndexLimits,&
       Grid_getBlkPtr,&
       Grid_releaseBlkPtr, Grid_pfft, Grid_pfftInit,&
       Grid_pfftFinalize, Grid_pfftMapToInput
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftData, ONLY : pfft_inLen,  pfft_work1,  pfft_work2

  implicit none

  logical ::  perfect = .true., needMap = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer :: blockID=1
  integer :: i,j,k,nbegin,nend,mbegin,mend,fullSize,temp
  integer :: sim_jprocs, sim_kprocs
  ! stays true if no errors are found

  integer, dimension(MDIM) :: comm,procGrid,me
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MDIM) :: globalSize, localSize,transformType,blkLen
  integer, dimension(4) :: prNum
  real :: err
  real, dimension(:), allocatable :: inArray, outArray, transArray
  real,dimension(:,:,:,:),pointer :: solnData
  integer,dimension(LOW:HIGH,MDIM) :: configLimits,phaseLimits
  integer :: iTimes, error, ii, jj, kk


  interface
     subroutine myIO_globalArray(localArray, realSpacing, globalSize)       
       implicit none
       real, dimension(:) :: localArray
       integer, intent(IN) :: realSpacing
       integer, dimension(MDIM), intent(IN) :: globalSize
     end subroutine myIO_globalArray
  end interface

  interface 
     subroutine testTranspose(inputArray, globalSize)
       implicit none
       real, dimension(:) :: inputArray
       integer, dimension(MDIM), intent(IN) :: globalSize
     end subroutine testTranspose
  end interface


  call RuntimeParameters_get('sim_jprocs',sim_jprocs)
  call RuntimeParameters_get('sim_kprocs',sim_kprocs)

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
  call Grid_getGlobalIndexLimits(globalSize)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  transformType(IAXIS)=PFFT_REAL2C
  transformType(JAXIS:KAXIS)=PFFT_COMPLEX


  !Chris temp
  !I label all internal grid points with the processor ID.  
  !After all the data transfers, I can scan the grid points in 
  !the constructed pfft array to check if communication 
  !implemented correctly.
  !--------------------------------------------------------
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           solnData(DENS_VAR,ii,jj,kk) = real(dr_globalMe)
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  !--------------------------------------------------------

  fullSize = pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)

  !inArray is used to hold real data from the mesh.
  allocate(inArray(fullSize), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if


  !transArray and outArray participate in complex-to-complex transforms 
  !so are allocated twice the size as the real array, inArray.
  allocate(transArray(2*fullSize), outArray(2*fullSize), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if


  call Grid_pfftMapToInput(DENS_VAR, inArray)
  
  !if(NDIM == 2) then
  !   call testTranspose(inArray, globalSize)
  !end if

  call Grid_pfft(PFFT_FORWARD,inArray,transArray)

  call Grid_pfft(PFFT_INVERSE,transArray,outArray)


  err = maxval(abs(outArray(1:fullSize)-inArray(1:fullSize)))

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
