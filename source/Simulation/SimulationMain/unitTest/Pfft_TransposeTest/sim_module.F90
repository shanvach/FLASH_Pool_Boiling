!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/sim_module
!!
!! NAME
!!
!!  sim_module
!!
!!***

#include "constants.h"
#include "Pfft.h"
#include "Flash.h"

module sim_module
  implicit none

contains
  subroutine sim_printLocalPencilData(arr, dimOrder, localShape, &
       fileNamePrefix)
    use Driver_interface, ONLY : Driver_abortFlash
    use Simulation_data, ONLY : sim_outputGridData, sim_meshMe
    implicit none
    real, dimension(:), intent(IN) :: arr
    integer, dimension(MDIM), intent(IN) :: dimOrder, localShape
    character (len=*), intent(IN) :: fileNamePrefix  

    integer, parameter :: STR_LEN = 200
    integer, dimension(LOW:HIGH,MDIM) :: localLimits
    integer :: dataStart, dataEnd, i, numColumns
    character (len=STR_LEN) :: fileName, uniqueName, formatStr

    if (sim_outputGridData .eqv. .true.) then
       if ( len(fileNamePrefix) > (STR_LEN-15) ) then
          call Driver_abortFlash("[sim_printLocalPencilData]: "//&
               "Possible buffer overflow.")
       end if

       !Create a Fortran format string something like: '(016i3)'
       !This works with gfortran, but other compilers will complain?
       write (formatStr, fmt='(''('',i3.3,''i3)'')') localShape(IAXIS)

       !Create a full file name.
       write(uniqueName, fmt='(i3.3,''.dat'')') sim_meshMe
       fileName = fileNamePrefix // uniqueName

       !Get the pencil limits owned by this processor.
       !Assert that the limit range is equal to the shape range.
       localLimits = -1 !For simulations with NDIM < MDIM.
       do i = 1,NDIM
          call sim_pfftGetLocalLimits(i,dimOrder(i),localShape,localLimits)
       end do

       print*,sim_meshMe,' LL:',localLimits
       print*,sim_meshMe,' LS:',localShape

       if ( product(localLimits(HIGH,1:NDIM)-localLimits(LOW,1:NDIM)+1) /= &
            product(localShape(1:NDIM)) ) then
          call Driver_abortFlash("[sim_printLocalPencilData]: "//&
               "Inconsistent limits and shape.")
       end if


       open(unit=12, file=fileName)
       write(12,'(a,2i4,a,2i4,a,2i4)') "  X:", localLimits(LOW:HIGH,IAXIS), &
            ", Y:", localLimits(LOW:HIGH,JAXIS), &
            ", Z:", localLimits(LOW:HIGH,KAXIS)
       !The IAXIS is fully enclosed, so we print the remaining dimensions
       !together.  As such, when we view all the files together with 
       !>>> more output_*
       !we will see how the global grid is distributed.

       !Works for standard distributions.
       do i = 1, localShape(JAXIS) * localShape(KAXIS)
          dataStart = ((i-1) * localShape(IAXIS)) + 1
          dataEnd = i * localShape(IAXIS)
          write(12,formatStr) int(arr(dataStart:dataEnd))
       end do
       close(12)
    end if
  end subroutine sim_printLocalPencilData


  !We really need to create a unified version of gr_pfftGetLocalLimits.
  !I have just thrown this one together and it seems to give the right
  !limits at each step.  This version is completely decoupled from all
  !the FFT stuff that doubles certain dimensions e.t.c.
  subroutine sim_pfftGetLocalLimits(axis1,axis2,localPortion,localLimits)
    use gr_pfftData, ONLY : pfft_me,pfft_globalLen,pfft_localLimits,pfft_outLen
    use gr_pfftInterface, ONLY : gr_pfftGetLocalLimitsAnytime
    implicit none
    integer,intent(IN) :: axis1, axis2
    integer,dimension(MDIM),intent(IN) :: localPortion
    integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: localLimits
    integer :: len,s

    call gr_pfftGetLocalLimitsAnytime(axis1,axis2,axis1,localPortion,PFFT_PCLDATA_REAL,localLimits)
    localLimits(:,axis1) = localLimits(:,axis1) + 1
    return                      !!!!!!
    len=pfft_globalLen(axis2) 
    s=pfft_me(axis1)*localPortion(axis1)
    if(s>len) then
       localLimits(LOW:HIGH,axis1)=1
    else
       localLimits(LOW,axis1)=s+1
       localLimits(HIGH,axis1)=s+localPortion(axis1)
    end if
  end subroutine sim_pfftGetLocalLimits

end module sim_module
