!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/InitParticlesUsingDensityDist
!!
!! NAME
!!    InitParticlesUsingDensityDist
!!
!! SYNOPSIS
!!
!!    InitParticlesUsingDensityDist( integer, INTENT(in) :: blockID,
!!                             integer, INTENT(out) :: numParticlesThisBlock, 
!!                             real, INTENT(in) :: maxDensity
!!                             real, INTENT(in) :: minDensity)
!!
!! DESCRIPTION
!!
!! Apply a rejection method to initialize the particle positions.  
!! Here, we choose positions based upon the density distibution in a block.
!! The number of particles to initialize is given by the input argument, 
!! numParticlesThisBlock.
!!
!! A modified version of the passive particle rejection method:
!! (Particles/ParticlesInitialization/WithDensity/RejectionMethod/pt_initPositions.F90)
!! is used.  In this version, however, the maximum and minimum density 
!! in the block is considered, so that the particle clustering follows the density
!! distribution in DENS_VAR.  This is done because the maximum and minimum density 
!! in this distribution are relatively similar, and both are non-zero.
!!
!! ARGUMENTS
!!
!! blockID: ID of block in which we will intialise particles.
!! numParticlesThisBlock:  The number of particles to create for this block.
!! maxDensity: The maximum density stored in a DENS_VAR cell in this block.
!! minDensity: The minimum density stored in a DENS_VAR cell in this block.
!!
!!***

subroutine InitParticlesUsingDensityDist(blockID, numParticlesThisBlock, maxDensity, minDensity)

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getSingleCellVol, Grid_getBlkBoundBox, &
       Grid_getBlkPhysicalSize, Grid_getDeltas

  use Driver_interface, ONLY : Driver_abortFlash
  use Particles_data, ONLY:  pt_numLocal, particles

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: blockID, numParticlesThisBlock
  real, intent(IN) :: maxDensity, minDensity

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(MDIM) :: blockSize, delta
  real, dimension(2,MDIM) :: boundBox
  integer, dimension(2,MDIM):: blkLimits, blkLimitsGC
  integer, parameter :: randomValue = 1
  integer       :: i, j, k, ierr, ii, jj, partID, p
  real          :: rho, xpos, ypos, zpos
  logical       :: accept
  logical, save :: firstCall = .true.

  ! random number generator
  integer :: seed_size
  real    :: urp
  real :: densityRange

  print *, "In InitParticlesUsingDensityDist.  pt_numLocal=", pt_numLocal


  !CD: Taken from Particles/ParticlesInitialization/WithDensity/pt_initLocal.F90
  !----------------------------------------------------------------------------
  if (firstCall) then

     !  randomize the initial particle positions     

     !  returned value seed_size gives the number of integers the processor uses for the
     !  starting value
     call random_seed(SIZE=seed_size)

     !  generates a large (from randomValue) and unique integer for each processor
     ii = int(randomValue * gr_meshMe) + gr_meshNumProcs
     !(Remove gr_meshNumProcs from the above expression to get the same particle 
     !distribution for any no. of procs).

     !  initializes the random number vector with a fixed seed (from i)
     call random_seed(PUT=(/(ii, jj = 1, seed_size)/))


     !  now call the random_number generator the same number of times as has been seeded.
     !  WHY we call the random_number generator so many times this is unknown -- perhaps to 
     !  really mix things up?  If so, then why a fixed seed?
     do jj = 1,2*ii

        !  returns a single uniform random number urp between zero and 
        call random_number (harvest=urp)

     end do

     firstCall = .false.
  end if
  !----------------------------------------------------------------------------



  !CD: Modified from Particles/ParticlesInitialization/WithDensity/RejectionMethod/pt_initPositions.F90
  !----------------------------------------------------------------------------
  call Grid_getBlkPhysicalSize(blockID,blockSize)  ! physical size of the block
  call Grid_getBlkBoundBox(blockID,boundBox)       ! physical corner of the block
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID,delta)

  densityRange = maxDensity - minDensity
  p = 0
  ! Initialize the particles for this block.
  do while (p < numParticlesThisBlock)

     ! Initialize particle positions using the rejection method:  choose uniform
     ! random x, y, and z coordinates and density.  If the density lies below the
     ! mesh density for the zone in which the coordinates lie, keep the particle;
     ! otherwise, choose another set of random values.  In determining the random
     ! choice box dimension for the density, use the max density times 1.01 to make
     ! sure we don't have any problems with roundoff in the highest-density zones
     ! in the block.
     ! NOTE that some whole complicated (and LBR doesn't understand it) procedure
     ! has been called in Particles_initWithDensity to initialize the random seed.

     call random_number (harvest=urp)

     !This expression is not used because it does not generate the required 
     !particle clustering.  The reason is because minDensity is similar to 
     !maxDensity, and both are larger than zero.
     !rho  = urp * maxDensity * 1.01     ! generate a random density less than maximum over this block

     rho = minDensity + (urp * densityRange)

     ! Calculate X position and index, needed in all geometries
     call random_number (harvest=urp)

     xpos = urp * blockSize(IAXIS)             
     i    = int(xpos/delta(IAXIS)) + blkLimits(LOW,IAXIS)
     xpos = xpos + boundBox(LOW,IAXIS)

     if (NDIM >= 2) then
        call random_number (harvest=urp)
        ypos = urp * blockSize(JAXIS)
        j    = int(ypos/delta(JAXIS)) + blkLimits(LOW,JAXIS)
        ypos = ypos + boundBox(LOW,JAXIS)
     else
        ypos = 0.
        j    = blkLimits(LOW,JAXIS)
     endif

     if (NDIM == 3) then
        call random_number (harvest=urp)
        zpos = urp * blockSize(KAXIS)
        k    = int(zpos/delta(KAXIS)) + blkLimits(LOW,KAXIS)
        zpos = zpos + boundBox(LOW,KAXIS)
     else
        zpos = 0.
        k    = blkLimits(LOW,KAXIS)
     endif


     if( (i<blkLimits(LOW,IAXIS)) .or. (i>blkLimits(HIGH,IAXIS)) ) then
        print *, "i not in internal range"
        call Driver_abortFlash ("out of range!")
     end if
     if( (j<blkLimits(LOW,JAXIS)) .or. (j>blkLimits(HIGH,JAXIS)) ) then
        print *, "j not in internal range"
        call Driver_abortFlash ("out of range!")
     end if
     if( (k<blkLimits(LOW,KAXIS)) .or. (k>blkLimits(HIGH,KAXIS)) ) then
        print *, "k not in internal range"
        call Driver_abortFlash ("out of range!")
     end if


     accept = (rho <= solnData(DENS_VAR,i,j,k))


     !If we choose to accept the particle then initilise using 
     !the guided random particle positions.
     !Then add one to the PDEN grid variable so that we 
     !can perform particle based refinement.

     if (accept) then

        p = p + 1
        partID = pt_numLocal + p

        particles(BLK_PART_PROP,partID) = real(blockID)
        particles(PROC_PART_PROP,partID) = real(gr_meshMe)
        particles(POSX_PART_PROP,partID) = xpos
        particles(POSY_PART_PROP,partID) = ypos
        particles(POSZ_PART_PROP,partID) = zpos
        particles(MASS_PART_PROP,partID) = 1.0

        particles(VELX_PART_PROP,partID) = 0.0
        particles(VELY_PART_PROP,partID) = 0.0
        particles(VELZ_PART_PROP,partID) = 0.0

        !****************************************************!
        ! This is how we perform particle based refinement   !
        !****************************************************!

        solnData(PDEN_VAR,i,j,k) = solnData(PDEN_VAR,i,j,k) + 1
        !print *, "Adding a PDEN_VAR hit on block:", blockID, "cell:", i, j, k
        !****************************************************!
        ! This is how we perform particle based refinement   !
        !****************************************************!

     endif

  end do  !loop over number of particles on each block.

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  !----------------------------------------------------------------------------  


end subroutine InitParticlesUsingDensityDist
