!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/save/op_initLowTemp
!!
!! NAME
!!
!!  op_initLowTemp
!!
!! SYNOPSIS
!!
!!  call op_initLowTemp ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the section of the low temperature opacity unit. Several arrays are
!!  allocated here and initialized with all the data. The number of atomic elements
!!  and their corresponding atomic numbers for the current opacity run must be known
!!  at this stage, otherwise the program will stop with a message.
!!
!!  The goal of this routine is to get and hold only the data corresponding to the
!!  atomic elements needed in memeory. The big arrays containing the data for all the
!!  elements will only be alive during the time spent in this routine.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initLowTemp ()

  use Opacity_data,     ONLY : op_totalElements,        &
                               op_element2AtomicNumber

  use op_lowTempData,   ONLY : op_Aij4,                 &
                               op_Jmax,                 &
                               op_PEenergyRange,        &
                               op_elementAij4,          &
                               op_elementJmax,          &
                               op_elementPEenergyRange
                                  
  use op_interface,     ONLY : op_setPEcoeffsAij4,      &
                               op_setPEarrayJmax,       &
                               op_setPEenergyRange

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

# include "Opacity.h"

  integer :: n
  integer :: status
  integer :: Z
!
!
!   ...Allocate the big arrays:
!
!           1) photoelectric cross section A(i,j,4) coefficient array
!           2) j-index delimiter array
!           3) photoelectronic energy range array
!
!
  allocate (op_Aij4 (1:4,1:12,1:100), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_Aij4 () allocate failed')
  end if

  allocate (op_Jmax (1:100), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_Jmax () allocate failed')
  end if

  allocate (op_PEenergyRange (LOW:HIGH,1:12,1:100), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_PEenergyRange () allocate failed')
  end if
!
!
!   ...Load the data into the big arrays.
!
!
  call op_setPEcoeffsAij4  ()
  call op_setPEarrayJmax   ()
  call op_setPEenergyRange ()
!
!
!   ...Optional printout of the big arrays. (For checking purposes only)
!
!
  call op_writePEdata ()
!
!
!   ...Allocate the small arrays:
!
!           1) photoelectric cross section A(i,j,4) coefficient array
!           2) j-index delimiter array
!           3) photoelectronic energy range array
!
!
  if ((op_totalElements < 1) .or. (op_totalElements > 100)) then
       call Driver_abortFlash ('[op_initLowTemp] ERROR: # of atomic elements < 1 or > 100')
  end if

  allocate (op_elementsAij4 (1:4,1:12,1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementsAij4 () allocate failed')
  end if

  allocate (op_elementsJmax (1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementsJmax () allocate failed')
  end if

  allocate (op_elementsPEenergyRange (LOW:HIGH,1:12,1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementsPEenergyRange () allocate failed')
  end if
!
!
!   ...Load the data into the small arrays.
!
!
  do n = 1,op_totalElements

     Z = op_element2AtomicNumber (n)

     if (Z < 1 .or. Z > 100) then
         call Driver_abortFlash ('[op_initLowTemp] ERROR: Atomic number Z < 1 or > 100')
     end if

     op_elementsJmax                        (n) = op_Jmax                        (Z)
     op_elementsAij4               (1:4,1:12,n) = op_Aij4               (1:4,1:12,Z)
     op_elementsPEenergyRange (LOW:HIGH,1:12,n) = op_PEenergyRange (LOW:HIGH,1:12,Z)

  end do
!
!
!   ...Deallocate the big arrays.
!
!
  deallocate (op_Aij4)
  deallocate (op_Jmax)
  deallocate (op_PEenergyRange)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initLowTemp
