!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/save/op_finalizeLowTemp
!!
!! NAME
!!
!!  op_finalizeLowTemp
!!
!! SYNOPSIS
!!
!!  call op_finalizeLowTemp ()
!!
!! DESCRIPTION
!!
!!  Finalizes the low temperature section of the opacity unit.
!!  All temporary arrays are deallocated to save memory.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeLowTemp ()

  use Opacity_dataLowTemp, ONLY : op_atomName,         &
                                  op_Aij4,             &
                                  op_Jmax,             &
                                  op_PEenergyRange

  implicit none
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_atomName)     ) deallocate (op_atomName)
  if (allocated (op_Aij4)         ) deallocate (op_Aij4)
  if (allocated (op_Jmax)         ) deallocate (op_Jmax)
  if (allocated (op_PEenergyRange)) deallocate (op_PEenergyRange)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeLowTemp
