!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/InitializeParticles
!!
!! NAME
!!    InitializeParticles
!!
!! SYNOPSIS
!!
!!    InitializeParticles()
!!
!! DESCRIPTION
!!
!!    Determine the total number of particles to initialize, and 
!!    then estimate the number of particles that should
!!    be intialised on each block.  Our estimate is obtained from the 
!!    relative density on the block, which is obtained from the 
!!    Gaussian distribution which is stored temporarily 
!!    in the grid element DENS_VAR.
!!
!!***

subroutine InitializeParticles()
  
  use Particles_data, ONLY:  pt_numLocal
  use Simulation_interface, ONLY : Simulation_initBlock
  use Grid_data, ONLY : gr_meshMe
  use pt_initFromFileInterface, ONLY : pt_initNumToGet, pt_initNextNParticles
  use Particles_data, ONLY:  particles, pt_maxPerProc, pt_numAtOnce
  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkPtr,Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getPointData

 implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer       :: i, j, k, ipos,jpos,kpos
  integer       :: p,pb,pe
  integer       :: blockID
  logical       :: notDone
  integer       :: blkCount, numToGet,numReturned, iBlk, ierr, numParticlesThisBlock
  integer,dimension(MAXBLOCKS) :: blkList
  real,dimension(MAXBLOCKS) :: densPerBlock, maxDensPerBlock , minDensPerBlock
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(MDIM):: del,delInv
  integer, dimension(MDIM) :: point
  real :: localDensity, cellDensity, globalDensity
  character (len=*), PARAMETER :: baseName="particle_dist"
  character (len=200) :: fileName, uniqueName!** NOTE 200 characters!!!
  integer, PARAMETER :: FILE_NUMB = 16
  integer :: eachP

  interface
     subroutine InitParticlesUsingDensityDist(blockID, numParticlesThisBlock, maxDensity, minDensity)
       integer, intent(IN) :: blockID, numParticlesThisBlock
       real, intent(IN) :: maxDensity, minDensity
     end subroutine InitParticlesUsingDensityDist
  end interface


  pt_numLocal=0


  !First determine the number of particles we will initialize in the simulation.
  call pt_initNumToGet(numToGet)


  !We must estimate how many particles should be placed in each top 
  !level block.
  localDensity = 0.0
  densPerBlock(:) = 0.0
  maxDensPerBlock(:) = 0.0
  minDensPerBlock(:) = 0.0

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  print *, "Processor", gr_meshMe, "has", blkCount, "blocks"


  !Sum the local density for each block on a processor.
  do iBlk = 1, blkCount

     blockID = blkList(iBlk)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        point(3) = k
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           point(2) = j
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(1) = i

              !Extract the density that has been placed in each cell.
              call Grid_getPointData(blockID, CENTER, DENS_VAR, EXTERIOR, point, cellDensity)

              densPerBlock(iBlk) = densPerBlock(iBlk) + cellDensity

              !Find the minimum density per block.
              if((i == blkLimits(LOW,IAXIS)).and.(j == blkLimits(LOW,JAXIS)).and.(k == blkLimits(LOW,KAXIS))) then
                 minDensPerBlock(iBlk) = cellDensity
              else
                 if (minDensPerBlock(iBlk) > cellDensity) then
                    minDensPerBlock(iBlk) = cellDensity
                 end if
              end if

              !Find the maximum density per block.
              if(maxDensPerBlock(iBlk) < cellDensity) then
                 maxDensPerBlock(iBlk) = cellDensity
              end if

           end do
        end do
     end do
     
     print *, "Density for block", blockID, "is:", densPerBlock(iBlk), "max density is:", maxDensPerBlock(iBlk), &
          "min density is:", minDensPerBlock(iBlk)
     localDensity = localDensity + densPerBlock(iBlk)

  end do  !! loop over all local leaf blocks


  call mpi_allreduce(localDensity, globalDensity, 1, FLASH_REAL, MPI_SUM, &
       MPI_COMM_WORLD, ierr)


  !Determine the number of particles that need to be initialized on 
  !each block based upon the localDensity.
  do iBlk = 1, blkCount

     blockID = blkList(iBlk)
     numParticlesThisBlock = (densPerBlock(iBlk) / globalDensity) * numToGet
     print *, "number of particles on block", iBlk, "is", numParticlesThisBlock

     !numParticlesThisBlock needs to be normalised so that 
     !the global sum of particles adds to the total number 
     !of particles.

     !Stops buffer overflow.  (Will not happen though because there is 
     !only one initial block in the computational domain.)
     !numParticlesThisBlock = numParticlesThisBlock * 0.7

     call InitParticlesUsingDensityDist(blockID, numParticlesThisBlock, maxDensPerBlock(iBlk), minDensPerBlock(iBlk))

     !--------------------------------------------------------------
     !Print the x and y position of each particle to file.
     write(uniqueName, fmt='(i3.3,''_'',i3.3,''.dat'')') gr_meshMe, iBlk
     fileName = baseName // uniqueName
     open(unit=FILE_NUMB, file=filename)

     do eachP = max(1,pt_numLocal), pt_numLocal + numParticlesThisBlock

        write(FILE_NUMB, *) particles(POSX_PART_PROP,eachP), particles(POSY_PART_PROP,eachP)

     end do
     close(FILE_NUMB)
     !--------------------------------------------------------------

     pt_numLocal = pt_numLocal + numParticlesThisBlock

  end do


  !Lets see how clustered the initial distribution is:
  if(blkCount == 1) then

     call Grid_getBlkPtr(blockID,solnData,CENTER)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              !print *, "i:", i, ",j:", j, ", PDEN_VAR:", solnData(PDEN_VAR,i,j,k)

           end do
        end do
     end do

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end if

  return

end subroutine InitializeParticles
