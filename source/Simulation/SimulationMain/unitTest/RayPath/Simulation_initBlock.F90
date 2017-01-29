!!****if* source/Simulation/SimulationMain/unitTest/RayPath/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a slab somewhere in the domain, which can 
!!  refract the light rays incident on it. The refractive index itself can
!!  be chosen to be either constant throughout the slab, or it can chosen 
!!  such that its gradient is constant throughout the slab.
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
       Grid_getBlkCoordinates, Grid_getBlkBoundBox
  use Simulation_data, ONLY : sim_slabBoundBox, sim_refract, sim_refractType

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer, dimension(MDIM) :: axis, globalSize

  integer :: i, j, k, i1, var
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer, dimension(MDIM) :: blkSize
  real,pointer,dimension(:,:,:,:)::solnData
  real, boundBox(LOW:HIGH,IAXIS:JAXIS)
  real, dimension(MDIM) :: pos
  logical :: overlap
  real,allocatable,dimension(:) :: iCoord,jCoord,kCoord

  integer, dimension(MDIM) :: lowEnd,highEnd


  call Grid_getBlkPtr(blockID,solnData,CENTER)
  
  solnData(:,:,:,:)=1.0
  
  call Grid_getBlkBoundBox(blockId, boundBox)

  call sim_overlap(boundBox,sim_slabBndBox,overlap)
  
  if(overlap) then
     call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
     lengths(1:NDIM)=blkLimitsGC(HIGH,1:NDIM)-blkLimits(LOW,1:NDIM)+1
     allocate(iCoord(lengths(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,CENTER, gcell, &
          iCoord, lengths(IAXIS))
     if(NDIM >1) then
        allocate(jCoord(lengths(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,CENTER, gcell, &
             jCoord, lengths(JAXIS))
     end if
     if(NDIM >2) then
        allocate(kCoord(lengths(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,CENTER, gcell, &
             kCoord, lengths(KAXIS))
     end if
     
     lowEnd(1:MDIM)=1
     highEnd(1:MDIM)=1
     do i = 1, NDIM
        if(sim_slabBndBox(LOW,i)<boundBox(LOW,i)) then
           lowEnd(i)=blkLimits(LOW,i)
        else
           lowEnd(i)=(sim_slabBndBox(LOW,i)-boundBox(LOW,i))/delta(i)
        end if
        if(sim_slabBndBox(HIGH,i)>boundBox(HIGH,i)) then
           highEnd(i)=blkLimits(HIGH,i)
        else
           highEnd(i)=(sim_slabBndBox(HIGH,i)-boundBox(Low,i))/delta(i)
        end if
     end do
     do k=lowEnd(KAXIS),highEnd(KAXIS)
        do j=lowEnd(JAXIS),highEnd(JAXIS)
           do i=lowEnd(IAXIS),highEnd(IAXIS)
              solnData(i,j,k,REFR_VAR)=sim_refract
           end do
        end do
     end do
  end if
  
  deallocate(iCoord)
  if(NDIM>1)deallocate(jCoord)
  if(NDIM>2)deallocate(kCoord)

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  
  return
end subroutine Simulation_initBlock
