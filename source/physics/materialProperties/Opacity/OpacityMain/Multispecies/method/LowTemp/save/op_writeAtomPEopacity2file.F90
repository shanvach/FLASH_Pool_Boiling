!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/save/op_writeAtomPEopacity2file
!!
!! NAME
!!
!!  op_writeAtomPEopacity2file
!!
!! SYNOPSIS
!!
!!  call op_writeAtomPEopacity2file (integer (in) :: Z,
!!                                   real    (in) :: Elower,
!!                                   real    (in) :: Eupper,
!!                                   integer (in) :: nPoints)
!!
!! DESCRIPTION
!!
!!  This routine writes the photoelectron opacities for element Z between the
!!  energy range Elower - Eupper in nPoints points to a specific file.
!!  The routine is meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeAtomPEopacity2file (Z,Elower,Eupper,nPoints)

  use Opacity_dataLowTemp, ONLY : op_Aij4,             &
                                  op_Jmax,             &
                                  op_PEenergyRange

  implicit none

# include "Opacity.h"

  integer, intent (in) :: nPoints
  integer, intent (in) :: Z

  real,    intent (in) :: Elower
  real,    intent (in) :: Eupper

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: i,j,n
  integer :: jmax
  integer :: nSteps

  real    :: Elow,Ehigh
  real    :: Energy
  real    :: Estep
  real    :: Opacity
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "Atom_opacities_plot.txt", &
        form = 'formatted')

  jmax = op_Jmax (Z)

  nSteps = nPoints - 1
  Estep  = (Eupper - Elower) / real (nSteps)

  do n = 1,nPoints

     Energy = Elower + real (n - 1) * Estep

     write (*,*) ' point #, energy = ',n,Energy


     do j = 1,jmax

        Elow  = op_PEenergyRange (LOW ,j,Z)
        Ehigh = op_PEenergyRange (HIGH,j,Z)

        write (*,*) ' Elow, Ehigh = ',Elow, Ehigh

        if ((Elow <= Energy) .and. (Energy < Ehigh)) then

             Opacity = 0.0
             do i = 1,4
                Opacity = Opacity + op_Aij4 (i,j,Z) / (Energy ** i)
             end do

             write (fileUnit,'(2ES10.3)') Energy,Opacity

             exit

        end if

     end do
  end do

  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeAtomPEopacity2file
