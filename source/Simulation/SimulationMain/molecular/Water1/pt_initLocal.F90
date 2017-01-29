!!****if* source/Simulation/SimulationMain/molecular/Water1/pt_initLocal
!!
!! NAME
!!    pt_initLocal
!!
!! SYNOPSIS
!!
!!    pt_initLocal()
!!
!! DESCRIPTION
!!    Local initialization of  particle locations.  Specific initializations that are
!!      needed only with gas density.  Calculates the total volume and the average density
!!      across all blocks.
!!      Initializes random fields.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pt_pRand:                integer  Something to do with initial random distribution
!!
!! NOTES
!!
!!   This version:
!!    Simplified copy of ParticlesInitialization/WithDensity/ptInitLocal.F90,
!!    density assumed constant. - KW
!!
!!   There is a nice description of the Fortran90 random number routines at
!!         http://www.nsc.liu.se/~boein/f77to90/a5.html#section21c
!!
!!***

!===============================================================================

subroutine pt_initLocal ()

  use Particles_data, ONLY:  pt_geometry,pt_myPE, pt_numProcs,pt_pRand, &
       pt_totalMass, pt_totalVolume, pt_averageDensity, pt_numParticlesWanted

  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
    Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits, &
    Grid_getBlkCenterCoords, Grid_releaseBlkPtr, Grid_getSingleCellVol

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"



  ! random number generator
  integer :: seed_size
  real    :: urp

  real          :: localMass, localVolume, localDensity
  real          :: dvol
  integer       :: thisBlock, bb, i, j, k, ierr, blockCount
  integer       :: numPEs
  integer       :: blockList(MAXBLOCKS)
  integer, dimension(MDIM) :: point
  integer, dimension(2,MDIM):: blkLimits, blkLimitsGC



  !-------------------------------------------------------------------------------

  ! Runtime Parameters


  ! Currently we only support Cartesian and 2D axisymmetric geometries.

  if ( (pt_geometry /= CARTESIAN) .and. &
       (.not. ((pt_geometry == CYLINDRICAL) .and. (NDIM == 2))) ) &
       call Driver_abortFlash ("pt_initLocal:  unsupported geometry for with density particle initialization!")

  ! In this routine, we determine the total volume and average density on the
  ! grid and save it.  Note that this will only work correctly if pt_initPositions
  ! has been called after the mesh has been set up OR the DENS_VAR variable contains
  ! accurate zone-average densities.

  ! In axisymmetric geometry we don't worry about factors of 2*pi or quadrants, since
  ! we are only concerned about the relative amount of mass in each block.

  localMass = 0.
  localVolume = 0.
  localDensity = 0.

  ! loop over all local leaf blocks
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do bb = 1, blockCount

     thisBlock = blockList(bb)

     ! get dimension limits etc.
     call Grid_getBlkIndexLimits(thisBlock,blkLimits,blkLimitsGC)


     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        point(3) = k
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           point(2) = j
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(1) = i
              call Grid_getSingleCellVol(thisBlock,EXTERIOR,point,dvol)
              localMass = localMass + dvol
              localVolume = localVolume + dvol
           enddo
        enddo
     enddo

     !  release the pointer

  enddo  !! of looping over all local leaf blocks
  

  !! get the mass across all processors and store it in data variable pt_totalMass
  call mpi_allreduce(localMass, pt_totalMass, 1, FLASH_REAL, MPI_SUM, &
       MPI_COMM_WORLD, ierr)
  call mpi_allreduce(localVolume, pt_totalVolume, 1, FLASH_REAL, MPI_SUM, &
       MPI_COMM_WORLD, ierr)

  ! now calculate the average density
  pt_averageDensity = pt_totalMass / pt_totalVolume

  !-------------------------------------------------------------------------------

  ! randomize the initial particle positions

  !  returned value seed_size gives the number of integers the processor uses for the
  !    starting value
  call random_seed(SIZE=seed_size)
  
  !  generates a large (from pt_pRand) and unique integer for each processor

  i = int(pt_pRand * pt_myPE) + pt_numProcs

  !  initializes the random number vector with a fixed seed (from i)
  call random_seed(PUT=(/(i, j = 1, seed_size)/))
  
  !  now call the random_number generator the same number of times as has been seeded.
  !  WHY we call the random_number generator so many times this is unknown -- perhaps to 
  !    really mix things up?  If so, then why a fixed seed?
  do j = 1,2*i
     !  returns a single uniform random number urp between zero and 
     call random_number (harvest=urp)
  end do



  !-------------------------------------------------------------------------------

  return

end subroutine pt_initLocal


