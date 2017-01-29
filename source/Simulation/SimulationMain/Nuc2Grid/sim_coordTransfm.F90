!!****if* source/Simulation/SimulationMain/Nuc2Grid/sim_coordTransfm
!!
!! NAME
!!
!!  sim_coordTransfm
!!
!! SYNOPSIS
!!
!!  call sim_coordTransfm(real(IN)    :: x,
!!                        real(IN)    :: y,
!!                        real(IN)    :: z,
!!                        real(OUT)   :: xout,
!!                        real(OUT)   :: yout,
!!                        real(OUT)   :: zout,
!!                        integer(IN) :: ndim,
!!                        real(IN)    :: velI,
!!                        real(IN)    :: velJ,
!!                        real(IN)    :: velK,
!!                        real(OUT)   :: velIOut,
!!                        real(OUT)   :: velJOut,
!!                        real(OUT)   :: velKOut)
!!
!!
!! DESCRIPTION
!!
!!  Convert from Cartesian to other coordinates.
!!
!! ARGUMENTS
!!
!!  x - First Cartesian coordinate.
!!  y - Second Cartesian coordinate.
!!  z - Third Cartesian coordinate.
!!  r - First output coordinate.
!!  theta - Second output coordinate.
!!  phi - Third output coordinate.
!!  ndim - dimensionality
!!
!!
!!***

#include "constants.h"


subroutine sim_coordTransfm(x,y,z, xout,yout,zout, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
  use Driver_interface, ONLY: Driver_abortFlash
  use Simulation_data, ONLY: sim_ptInGeometry, sim_geometry, sim_meshMe
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: xout,yout,zout
  integer,intent(IN) :: ndim
  real,intent(IN) :: velI,velJ,velK
  real,intent(OUT) :: velIOut,velJOut,velKOut
  logical,save :: first_call = .TRUE.

  
  if (first_call) then
     print*,'sim_ptInGeometry=',sim_ptInGeometry
     print*,'sim_geometry=',sim_geometry
  end if

  velIout = velI
  velJOut = velJ
  velKOut = velK
  if (sim_geometry==sim_ptInGeometry) then
     if (sim_meshMe==MASTER_PE) then
        if (first_call) print*,'Input geometry ok.'
     end if
     xout = x
     yout = y
     zout = z
  else if (sim_ptInGeometry==CARTESIAN) then
     if (sim_geometry==SPHERICAL) then
        call sim_coordTransfmSph(x,y,z, xout,yout,zout, ndim)
     else if (sim_geometry==CYLINDRICAL) then
        call sim_coordTransfmCyl(x,y,z, xout,yout,zout, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
     else
        call Driver_abortFlash("invalid particlesInputGeometry")
     end if
  else
     call Driver_abortFlash("invalid geometry")
  end if

  first_call = .FALSE.


end subroutine sim_coordTransfm

!!****if* source/Simulation/SimulationMain/PhoenixInputKeepNames/sim_coordTransfm
!!
!! NAME
!!
!!  sim_coordTransfmSph
!!
!! SYNOPSIS
!!
!!  call sim_coordTransfmSph(real(IN)    :: x,
!!                        real(IN)    :: y,
!!                        real(IN)    :: z,
!!                        real(OUT)   :: r,
!!                        real(OUT)   :: theta,
!!                        real(OUT)   :: phi,
!!                        integer(IN) :: ndim)
!!
!!
!! DESCRIPTION
!!
!!  Convert from Cartesian to spherical coordinates.
!!
!! ARGUMENTS
!!
!!  x - First Cartesian coordinate.
!!  y - Second Cartesian coordinate.
!!  z - Third Cartesian coordinate.
!!  r - First spherical coordinate.
!!  theta - Second spherical coordinate.
!!  phi - Third spherical coordinate.
!!  ndim - dimensionality
!!
!!
!!***
subroutine sim_coordTransfmSph(x,y,z, r,theta,phi, ndim)
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: r,theta,phi
  integer,intent(IN) :: ndim

  real :: rsq

  Theta = PI*0.5
  phi = 0.0

  rsq = x**2
  if (ndim > 1) rsq = rsq + y**2
  if (ndim > 2) rsq = rsq + z**2
  r = sqrt(rsq)

  if (r == 0.0) return

  if (ndim > 1) then
     theta = acos(z/r)
  end if

  if (ndim > 2) then
     phi = atan2(y,x)
  end if

end subroutine sim_coordTransfmSph

subroutine sim_coordTransfmCyl(x,y,z, r,zout,phi, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: r,zout,phi
  integer,intent(IN) :: ndim
  real,intent(IN) :: velI,velJ,velK
  real,intent(INOUT) :: velIOut,velJOut,velKOut

  real :: rsq, velSq, velIJ

  zout = 0.0
  phi = 0.0

  rsq = x**2
  if (ndim > 1) rsq = rsq + y**2
  r = sqrt(rsq)

  if (ndim > 1) then
     velSq = velI**2 + velJ**2
     zout = z
     if (velSq.NE.0.0) then
!!$        velIJ = sqrt(velSq)
        if (r .NE. 0.0) then
           velIOut = (velI*x + velJ*y) / r
           velJOut = (velJ*x - velI*y) / r
        end if
     end if
  end if

  if (r == 0.0) return

  if (ndim > 2) then
     phi = atan2(y,x)
  end if

end subroutine sim_coordTransfmCyl

#ifdef DEBUG
program test
  implicit none

  real :: x,y,z,r,theta,phi
  do while(1)
     print*,'Enter x,y,z:'
     read*,x,y,z
     call sim_coordTransfmSph(x,y,z,r,theta,phi,3)
     print*,'R    =',r
     print*,'theta=',theta,' (',theta*180/PI,' degrees)'
     print*,'phi  =',phi  ,' (',phi  *180/PI,' degrees)'
  end do

end program test
#endif
