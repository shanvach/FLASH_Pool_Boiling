!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/save/op_PlanckMeanOpacity
!!
!! NAME
!!
!!  op_PlanckMeanOpacity
!!
!! SYNOPSIS
!!
!!  call op_PlanckMeanOpacity ()
!!
!! DESCRIPTION
!!
!!  This routine constructs the group Planck mean opacity for a particular species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_PlanckMeanOpacity ()

  use Opacity_dataLowTemp, ONLY : op_atomName,             &
                                  op_Aij4,                 &
                                  op_Jmax,                 &
                                  op_PEenergyRange,        &
                                  op_eBaseLowestExponent,  &
                                  op_eBaseLargestExponent, &
                                  op_eBaseRescaleExponent

  use Driver_interface,    ONLY : Driver_abortFlash

  use Timers_interface,    ONLY : Timers_start,            &
                                  Timers_stop

  implicit none

  integer :: n,nSec

  integer, parameter :: maxSec = 10

  real    :: A1Sec     (1:maxSec)
  real    :: A2Sec     (1:maxSec)
  real    :: A3Sec     (1:maxSec)
  real    :: A4Sec     (1:maxSec)
  real    :: intLimits (1:maxSec+1)
  real    :: BiggsPlanckIntegral
  real    :: PlanckIntegral
!
!
!   ...Proceed.
!
!
!  write (*,*) ' # of sections (< 10) ? '
!  read  (*,*) nSec
!
!  do n = 1,nSec
!     write (*,*) ' A1 (',n,') ? '
!     read  (*,*) A1Sec (n)
!     write (*,*) ' A2 (',n,') ? '
!     read  (*,*) A2Sec (n)
!     write (*,*) ' A3 (',n,') ? '
!     read  (*,*) A3Sec (n)
!     write (*,*) ' A4 (',n,') ? '
!     read  (*,*) A4Sec (n)
!     if (n==1) then
!         write (*,*) ' Lower integration limit (1) ? '
!         read  (*,*) intLimits (1)
!     end if
!     write (*,*) ' Upper integration limit (',n,') ? '
!     read  (*,*) intLimits (n+1)
!  end do

  nSec = 4

  A1Sec (1) = 1.E8
  A2Sec (1) = -1.E7
  A3Sec (1) = 1.E6
  A4Sec (1) = -1.E9

  A1Sec (2) = 2.E5
  A2Sec (2) = 4.E6
  A3Sec (2) = 1.E8
  A4Sec (2) = -2.E10

  A1Sec (3) = -1.E-3
  A2Sec (3) = 9.E-1
  A3Sec (3) = -3.E2
  A4Sec (3) = 7.E5

  A1Sec (4) = 7.E0
  A2Sec (4) = 1.E-5
  A3Sec (4) = -8.E4
  A4Sec (4) = 2.E-3

  intLimits (1) = 1
  intLimits (2) = 5
  intLimits (3) = 20
  intLimits (4) = 100
  intLimits (5) = 2000

  call Timers_start  ("Integrator")

  do n = 1,1

  call op_BiggsPlanckGroupIntegrate (nSec,                                          &
                                     A1Sec,A2Sec,A3Sec,A4Sec,                       &
                                     intLimits,                                     &
                                     op_eBaseLowestExponent,                        &
                                     op_eBaseLargestExponent,                       &
                                     op_eBaseRescaleExponent,                       &
                                                              BiggsPlanckIntegral,  &
                                                              PlanckIntegral        )

  end do

  call Timers_stop  ("Integrator")

  write (*,*) ' Biggs-Planck Integral = ',BiggsPlanckIntegral
  write (*,*) '       Planck Integral = ',PlanckIntegral
!
!
!   ...Ready! 
!
!
  return
end subroutine op_PlanckMeanOpacity
