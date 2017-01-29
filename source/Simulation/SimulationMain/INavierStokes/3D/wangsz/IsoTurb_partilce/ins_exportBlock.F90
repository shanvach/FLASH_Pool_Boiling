! Modified from subroutine ins_exportslices_z
! Explorts flow information for  turbulence simulations:
! One file will be written per bloci.
! Shizhao Wang
! Dec 05, 2014
! --------------------------------------------------------------------------

  subroutine ins_exportBlock(istep,iwrite,time,mvisc,blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkCenterCoords,&
                             Grid_getBlkPhysicalSize,&
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkCornerID,    &
                             Grid_getGlobalIndexLimits

  
!  use ins_interface, only  :  ins_velomg2center
  
  use Grid_data, only : gr_meshMe, gr_meshNumProcs, gr_meshComm, gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "IncompNS.h"

  integer :: istep, iwrite
  real    :: time, mvisc
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real :: del(MDIM),coord(MDIM)

  integer i,j,k,lb,blockID

  character(80) :: filename
  integer :: ixs, ixe, iys, iye, izs, ize

  integer :: iFile

#ifdef FLASH_GRID_PARAMESH
  write(*,*) 'The routine ins_exportZone has not been tested for the ParaMesh.'
  write(*,*) 'The routine ins_exportZone just allows one block per procs..'
  stop
#endif

  ixs = NGUARD + 1
  ixe = NGUARD + NXB
  iys = NGUARD + 1
  iye = NGUARD + NYB
  izs = NGUARD + 1
  ize = NGUARD + NZB

  do lb = 1,blockCount
      blockID = blockList(lb)
 
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
 
      ! Point to blocks center and face vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER)
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
      call Grid_getBlkPtr(blockID,facezData,FACEZ)

      write(filename,'("./IOData/stat_block.",i4.4,".",i6.6,".bin")') &
                     iwrite, gr_meshMe  
      iFile = 51000000+gr_meshMe
      open(iFile,file=trim(filename),status='replace',form='unformatted')

        write(iFile) time,mvisc
        write(iFile) istep, NXB, NYB, NZB
        write(iFile) del
        write(iFile) gr_imin,gr_imax,gr_jmin,gr_jmax, gr_kmin, gr_kmax

        write(iFile)  (((solnData(PRES_VAR,i,j,k), i=ixs,ixe),j=iys,iye),k=izs,ize) 
        write(iFile)  (((facexData(VELC_FACE_VAR,i,j,k), i=ixs,ixe+1),j=iys,iye),k=izs,ize) 
        write(iFile)  (((faceyData(VELC_FACE_VAR,i,j,k), i=ixs,ixe),j=iys,iye+1),k=izs,ize) 
        write(iFile)  (((facezData(VELC_FACE_VAR,i,j,k), i=ixs,ixe),j=iys,iye),k=izs,ize+1) 

      close(iFile)

      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
 
   enddo

  return
  end subroutine ins_exportBlock



