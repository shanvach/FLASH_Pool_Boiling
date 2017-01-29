!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/pt_initNumToGet
!!
!! NAME
!!    pt_numToGet
!!
!! SYNOPSIS
!!
!!    pt_numToGet(integer(OUT) :: numToGet)
!!
!! DESCRIPTION
!!    
!!   Return the number of particles to be got (read in from the file most likely)
!!
!! ARGUMENTS
!!
!!  numToGet :  the number of particles
!!
!!
!!***
subroutine pt_initNumToGet(numToGet)

 use pt_initData, ONLY : pt_numX, pt_numY, pt_numZ
 use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#include "constants.h"
#include "Flash.h"

implicit none
  
  integer, intent(OUT) :: numToGet

  call pt_init()

  if (NDIM == 1) then
     numToGet = pt_numX
  else if (NDIM == 2) then
     numToGet = pt_numX * pt_numY
  else if (NDIM == 3) then
     numToGet = pt_numX * pt_numY * pt_numZ
  end if

  !If we are using more than 1 processor:
  !numToGet = numToGet / pt_meshNumProcs
  !(However, we must normalise to ensure that the global sum
  !of all numToGet's adds up to total particles).

  !Above is commented out because NBLOCKX = NBLOCKY = NBLOCKZ = 1, 
  !and so the number to get is just the global number of particles.  Choosing 
  !a different value for NBLOCK[X,Y,Z] is not possible because 
  !multigrid requires the number of blocks at the coarsest level to be 1.

end subroutine pt_initNumToGet
