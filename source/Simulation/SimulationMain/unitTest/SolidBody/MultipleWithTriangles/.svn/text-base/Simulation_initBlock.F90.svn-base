!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(integer(in) :: blockID)
!!                         
!!
!! DESCRIPTION   
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!! PARAMETERS
!!
!!***


 subroutine Simulation_initBlock(blockId)

!============================================================================
  use Simulation_data, ONLY : sim_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  implicit none 

#include "constants.h"
#include "Flash.h"


  integer, intent(IN) :: blockID
  real, pointer, dimension(:,:,:,:) :: fcx, fcy, fcz
  real, dimension(:), allocatable :: coord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: facecent,  faces
  integer :: axis, coordSize, i

  call Grid_getBlkPtr(blockID, fcx, FACEX)
  call Grid_getBlkPtr(blockID, fcy, FACEY)
  call Grid_getBlkPtr(blockID, fcz, FACEZ)

  do axis = IAXIS, KAXIS
     if (axis .eq. IAXIS) then
        facecent = (/ FACES, CENTER, CENTER/)
        faces = (/FACEX,CENTER,CENTER/)
     elseif (axis .eq. JAXIS) then
        facecent = (/ CENTER, FACES, CENTER/)
        faces = (/CENTER,FACEY,CENTER/)
     elseif (axis .eq. JAXIS) then
        facecent = (/ CENTER, CENTER, FACES/)
        faces = (/CENTER,CENTER,FACEZ/)
     endif
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, faces(axis))
     coordSize = blkLimitsGC(HIGH,axis)-blkLimitsGC(LOW,axis)+1
     allocate(coord(coordSize))
     
     call Grid_getCellCoords(axis, blockID, facecent(axis), .true., coord, coordSize)
     i = 0
     do i=blkLimitsGC(LOW,axis),blkLimitsGC(HIGH,axis)-1
        if (axis .eq. IAXIS) then
           fcx(VELC_FACE_VAR,i,:,:) = coord(i)
        elseif (axis .eq. JAXIS) then
           fcy(VELC_FACE_VAR,:,i,:) = coord(i)
        elseif (axis .eq. KAXIS) then
           fcz(VELC_FACE_VAR,:,:,i) = coord(i)
        endif
     end do

     deallocate(coord)
  enddo

  call Grid_releaseBlkPtr(blockID,fcx)
  call Grid_releaseBlkPtr(blockID,fcy)
  call Grid_releaseBlkPtr(blockID,fcz)

end subroutine Simulation_initBlock
