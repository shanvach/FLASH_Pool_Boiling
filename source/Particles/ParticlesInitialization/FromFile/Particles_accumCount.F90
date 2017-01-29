!!****if* source/Particles/ParticlesInitialization/FromFile/Particles_accumCount
!!
!! NAME
!!    Particles_accumCount
!!
!! SYNOPSIS
!!
!!    Particles_accumCount(integer, intent(IN) :: var)
!!
!! DESCRIPTION
!!    
!!   This routine is to be used in the refinement based upon the number of
!!   particles in a block. It adds a weight to the cell of the grid in the
!!   specified grid variable if a particle is found to be contained in the 
!!   cell 
!!
!! ARGUMENTS
!!
!!  var:        the grid variable to add the weights if particle found in the cell
!!
!!
!!***


subroutine Particles_accumCount(var)

  use Particles_data, ONLY:  particles, pt_maxPerProc, pt_numAtOnce
  
  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getDeltas,Grid_getBlkIndexLimits,&
                             Grid_getBlkBoundBox,Grid_getBlkPtr,Grid_releaseBlkPtr
  use pt_interface, ONLY : pt_initNumToGet, pt_initNextNParticles

  use Particles_data, ONLY:  pt_numLocal


  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: var
  integer       :: i, j, k, ipos,jpos,kpos
  integer       :: p,pb,pe
  integer       :: blockID
  logical       :: notDone
  integer       :: blkCount, numInFile,numReturned
  integer,dimension(MAXBLOCKS) :: blkList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(MDIM):: del,delInv

!----------------------------------------------------------------------

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  pt_numLocal=0

  !call pt_initNumInFile(numInFile)

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  

  

  jpos=1; kpos=1
  do i = 1,numInFile,pt_numAtOnce
     pe = pt_maxPerProc

     call pt_initNextNParticles(pe,numReturned,pb)

     notDone=.true.
     
     j=0
     do while(notDone.and.j<blkCount)
        blockID=blkList(j+1)
        call Grid_getBlkBoundBox(blockID,boundBox)
        k=1
        call pt_initIfInBlock(boundBox,numReturned,k,pe)
        if(numReturned>0) then
           call Grid_getBlkPtr(blockID,solnData,CENTER)
           call Grid_getDeltas(blockID,del)
           call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
           delInv(1:NDIM)=1.0/del(1:NDIM)
           do k = 1,numReturned
              ipos=(particles(POSX_PART_PROP,k)-boundBox(LOW,IAXIS))*delInv(IAXIS)&
                   +blkLimits(LOW,IAXIS)
              if(NDIM>1)jpos=(particles(POSY_PART_PROP,k)-boundBox(LOW,JAXIS))*delInv(JAXIS)&
                   +blkLimits(LOW,JAXIS)
              if(NDIM>2)kpos=(particles(POSZ_PART_PROP,k)-boundBox(LOW,KAXIS))*delInv(KAXIS)&
                   +blkLimits(LOW,KAXIS)
              solnData(var,ipos,jpos,kpos)=solnData(var,ipos,jpos,kpos)+1
           end do
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        end if
        notDone= (pe>=pb)
        j=j+1
     end do
  end do
  
  return
  
  !----------------------------------------------------------------------
  
end subroutine Particles_accumCount


