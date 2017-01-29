!!****if* source/Particles/ParticlesInitialization/FromFile/pt_initNextNParticles
!!
!! NAME
!!    pt_initNextNParticles
!!
!! SYNOPSIS
!!
!!    pt_initNextNParticles(integer(IN) :: pe,
!!                          integer(OUT) :: numReturned,
!!                          integer(OUT) :: pb)
!!
!! DESCRIPTION
!!    
!!   Fill in the next N particles, most often by reading the values from some file
!!
!! ARGUMENTS
!!
!!  pe : pointer to the last particle
!!  numReturned : number of particles read in this invocation
!!  pb : pointer to the first particle to be read in
!!  
!!
!!***
subroutine pt_initNextNParticles(pe,numReturned,pb)
  use Particles_data, ONLY : particles,pt_maxPerProc,pt_numAtOnce
  use Grid_data, ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax

  implicit none
#include "Flash.h"
#include "constants.h"
  
  integer, intent(IN) :: pe
  integer, intent(OUT) :: numReturned,pb

  logical, save:: firstcall=.true.
  real,save :: posx,posy,posz
  real      :: del=0.01
  integer   :: i
  
  if(firstcall) then
     firstcall = .false.
     posx=gr_imin
     posy=gr_jmin
     posz=gr_kmin
  end if

  numReturned=pt_numAtOnce
  pb=pe-numReturned+1

  do i = pb,pe
     particles(POSX_PART_PROP,i)=posx+del
     particles(POSY_PART_PROP,i)=posy+K2D*del
     particles(POSZ_PART_PROP,i)=posz+K3D*del
     posx=posx+del
     posy=posy+K2D*del
     posz=posz+K3D*del
  end do
end subroutine pt_initNextNParticles
