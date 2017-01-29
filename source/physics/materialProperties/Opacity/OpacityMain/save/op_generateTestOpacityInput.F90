!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/op_generateTestOpacityInput
!!
!! NAME
!!
!!  op_generateTestOpacityInput
!!
!! SYNOPSIS
!!
!!  call op_generateTestOpacityInput ()
!!
!! DESCRIPTION
!!
!!  Generates the opacity input file for testing purposes. This file contains
!!  info about what kind of opacities are wanted for each species.
!!
!!  Warning:
!!
!!  The structure of this routine is closely tied to how the input is read in
!!  the Opacity_init routine. Close inspection of that routine is necessary
!!  to avoid surprises. The input must be read in in exactly the same way.
!!
!! ARGUMENTS
!!
!!***
subroutine op_generateTestOpacityInput ()

  use Opacity_dataUnitTest, ONLY : op_totalSpecies,      &
                                   op_absorptionKind,    &
                                   op_emissionKind,      &
                                   op_transportKind,     &
                                   zero

  implicit none

  character (len=40) :: dummyLine  = '----------------------------------------'
  character (len= 1) :: tabledot   = '.'
  character (len=10) :: tableBase  = 'OPACITIES_'
  character (len= 3) :: tableLabel = 'imx'
  character (len= 6) :: tableKind  = 'IONMIX'
  character (len=17) :: tableName

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: species
!
!
!    ...Open the 'opacity_input.txt' file and proceed from there.
!
!
  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit, file = "opacity_input.txt")
!
!
!   ...Loop over all species. Since we only want to test the correctness of reading
!      the tables, we set each kind (absorption, emission and transport) equal among
!      themselves. For the same reason we do not include the option of constant
!      opacities in the test.
!
!
  do species = 1,op_totalSpecies

     write (tableName,'(A10,I3.3,A1,A3)') tableBase,species,tabledot,tableLabel

     write (fileUnit,'(A20)'  ) op_absorptionKind (species)
     write (fileUnit,'(A20)'  ) op_emissionKind   (species)
     write (fileUnit,'(A20)'  ) op_transportKind  (species)
     write (fileUnit,'(E12.6)') zero
     write (fileUnit,'(E12.6)') zero
     write (fileUnit,'(E12.6)') zero
     write (fileUnit,'(A6)'   ) tableKind
     write (fileUnit,'(A17)'  ) tableName
     write (fileUnit,'(A40)'  ) dummyLine

  end do
!
!
!    ...Close the 'opacity_input.txt' file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!
  return
end subroutine op_generateTestOpacityInput
