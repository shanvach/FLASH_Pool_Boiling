!!****if* source/Simulation/SimulationMain/unitTest/Particles/sim_addParticles
!!
!! NAME
!!
!!  sim_addParticles
!!
!! SYNOPSIS
!!
!!  sim_addParticles
!!
!! DESCRIPTION
!!
!!
!! NOTES
!!
!!
!!
!!***



subroutine sim_addParticles()

  use Simulation_data, ONLY : sim_meshMe, sim_addPartCount, sim_globalBndBox, sim_addPartDisp, &
       sim_grPtRemove
  use Particles_interface, ONLY : Particles_addNew
  use Particles_data, ONLY : particles, pt_numLocal, pt_maxPerProc, pt_indexList, pt_indexCount
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"

  logical,parameter :: extraRemove=.FALSE. !Do an extra attempt to remove particles here if gr_ptRemove is on?

  real, allocatable, dimension(:,:) :: pos
  real,dimension(MDIM):: disp
  integer :: i,j
  logical :: success

  allocate(pos(MDIM,sim_addPartCount))

  pos=0.0
  do i = 1,NDIM
     disp(i)=(sim_globalBndBox(HIGH,i)-sim_globalBndBox(LOW,i))*sim_addPartDisp
     pos(i,1)=sim_globalBndBox(LOW,i)+disp(i)
  end do

  do j = 2,sim_addPartCount
     do i = 1,NDIM
        pos(i,j)=pos(i,j-1)+disp(i)
     end do
  end do

  99 format('pt_numLocal on',I5,' before adding:',I15)
  print 99,sim_meshMe,pt_numLocal
  call Particles_addNew(sim_addPartCount,pos,success)
  if(.not.success) then
#ifdef DEBUG_ALL
     print*,'Tags',particles(TAG_PART_PROP,1:pt_numLocal), pt_numLocal
#endif
     if (sim_grPtRemove) then
        if (extraRemove) then
           call gr_ptSetIndices(pt_indexList,pt_indexCount)
           call gr_ptHandleExcess(particles, NPART_PROPS, pt_numLocal, pt_maxPerProc)
           call gr_ptResetIndices(pt_indexList,pt_indexCount)
        else
           print*,'Did not add particles for lack of space on', sim_meshMe
           call Logfile_stamp("Did not add particles for lack of space.","[sim_addParticles]")
        end if
     else
        call Driver_abortFlash("no space for adding particles")
     end if
  end if
  deallocate(pos)
  return
  
end subroutine Sim_addParticles



