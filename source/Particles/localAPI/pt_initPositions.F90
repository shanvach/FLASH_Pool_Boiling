!!****if* source/Particles/localAPI/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    Initializes particle locations for one block in the grid.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!!***


subroutine pt_initPositions (blockID,success)


  implicit none

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  success = .true. ! DEV: returns true because this stub creates no particles,
                   ! therefore all of those zero particles were created successfully
  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositions


