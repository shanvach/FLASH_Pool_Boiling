!!****if* source/Simulation/SimulationMain/unitTest/particleTests/VParticles/pt_advanceCustom
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)    :: dtold,
!!                        real(in)    :: dtnew,
!!                        real(inout) :: particlesunused(NPART_PROPS,p_countUnused),
!!                        integer(in) :: p_countunused)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particlesunused -- particles on which to operate
!!   p_countunused - the number of particles in the list to advance
!!
!!
!!
!!***

subroutine pt_advanceCustom(dtOld,dtNew, particles,p_count, ind)
    
  use Simulation_data, ONLY : sim_deltaMove
  implicit none

#include "Flash.h"
#include "constants.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

  integer :: i,j
  integer,dimension(MDIM):: pos

  pos(IAXIS)=POSX_PART_PROP
  pos(JAXIS)=POSY_PART_PROP
  pos(KAXIS)=POSZ_PART_PROP

  !write(*,*) 'Particle Count=',p_count

  do i = 1,p_count
     do j= 1,NDIM
	!write(*,*)'Old pos= ',particles(pos(j),i)
        particles(pos(j),i)=particles(pos(j),i)+sim_deltaMove(j)
	!write(*,*) 'New pos= ',particles(pos(j),i)
     end do
  end do
     
end subroutine pt_advanceCustom
