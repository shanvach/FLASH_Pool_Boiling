!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhErfc
!!
!! NAME
!!
!!  gr_bhErfc
!!
!!
!! SYNOPSIS
!!
!!   real var = gr_bhErfc(real:: x)
!!
!! DESCRIPTION
!!
!!   A dummy implementation to use in place of the complementary error
!!   function erfc().  Always returns 1.0.
!!
!! ARGUMENTS
!!
!!   x : argument of erfc(x)
!!
!! NOTES
!!
!!  Use this to be able to build FLASH on a system / with a compiler
!!  that does not provide an erfc implementation.  Define
!!  USER_ERFC for compiling FLASH. Two simple ways to do this are:
!!  Variant (1): Add
!!                    -defines=USER_ERFC=gr_bhErfc
!!               to the setup command line.
!!  Variant (2): Add
!!                    PPDEFINE USER_ERFC gr_bhErfc
!!               to your simulation's Config file.
!!
!!
!!  If you use this fake implementation and also turn on gr_bhUseEwaldDecomp,
!!  you get what you deserve!
!!
!! SEE ALSO
!!
!!  README.erfc
!!***

function gr_bhErfc(x)
! Totally fake Chebyshev polynomial approximation
  implicit none
  real, intent(in) :: x
  real             :: gr_bhErfc

  gr_bhErfc = 1.0

  return
end function gr_bhErfc
