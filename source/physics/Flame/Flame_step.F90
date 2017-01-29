!!****f* source/physics/Flame/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  call Flame_step ( integer, intent(in) :: num_blocks,
!!                    integer(:), intent(in) :: blockList,
!!                    real, intent(in) :: dt )
!!
!! DESCRIPTION
!!
!!  Performs a flame step.
!!
!! ARGUMENTS
!!
!!   num_blocks - number of blocks.
!!    blockList - the block list.
!!           dt - the time step.
!!
!!***

subroutine Flame_step( num_blocks, blockList, dt  )    

  implicit none
    
! args
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList
  real,    INTENT(in)                        :: dt

  return
  
end subroutine Flame_step

