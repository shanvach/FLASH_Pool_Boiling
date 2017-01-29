!!****ih* source/Particles/ParticlesInitialization/FromFile/pt_initFromFileInterface
!!
!! This is the header file for the initialization subunit of the Particles
!! unit. It defines the interfaces with subunit scope.
!! 
!!***
Module pt_initFromFileInterface

#include "constants.h"

  interface
     subroutine pt_initNumToGet(numToGet)
       integer,intent(OUT) :: numToGet
     end subroutine pt_initNumToGet
  end interface

  interface
     subroutine pt_initNextNParticles(pe,numReturned,pb)
       integer,intent(IN) :: pe
       integer,intent(OUT) :: numReturned,pb
     end subroutine pt_initNextNParticles
  end interface

  interface 
     subroutine pt_initIfInBlock(boundBox,numParticles,putStart,getEnd)
       real, dimension(LOW:HIGH,MDIM), intent(IN) :: boundBox
       integer, intent(IN) :: putStart
       integer, intent(INOUT) :: getEnd,numParticles
     end subroutine pt_initIfInBlock
  end interface
end Module pt_initFromFileInterface
