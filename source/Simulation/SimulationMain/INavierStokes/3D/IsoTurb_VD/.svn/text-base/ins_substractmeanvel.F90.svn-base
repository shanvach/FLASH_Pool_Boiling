! --------------------------------------------------------------------------

  subroutine ins_substractmeanvel(mype,istep,time,mvisc, &
                                  blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits

  use Grid_data, only : gr_meshMe
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  include 'mpif.h'

  integer :: mype, istep
  real :: time, mvisc
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  integer, parameter ::  ng = NGUARD
  integer, parameter ::  nxi= NGUARD + NXB
  integer, parameter ::  nyj= NGUARD + NYB
  integer, parameter ::  nzk= NGUARD + NZB
  integer, parameter ::  nxc= NGUARD + NXB + 1
  integer, parameter ::  nyc= NGUARD + NYB + 1
  integer, parameter ::  nzc= NGUARD + NZB + 1

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real del(3),dx,dy,dz
  real turbkinaux,volcell,tvol,tvolaux
  real varl,varlaux,vart,vartaux,dudx_v,dudy_v
  real lambdaf,lambdag
  real uavg,vavg,wavg,uavgaux,vavgaux,wavgaux,fconstn
  real dissavg,dissavgaux

  integer lb,blockID,ierr

  ! -------------------------------------------------------------------
  ! Average velocities, Check there is no drift:
  uavg = 0.
  vavg = 0.
  wavg = 0.
  tvol = 0.

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Cell volume:
     volcell = del(DIR_X)*del(DIR_Y)*del(DIR_Z) ! dx*dy*dz

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! U average:
     uavg = uavg + &
      sum(facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk))*volcell

     ! V average:
     vavg = vavg + &
      sum(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk))*volcell

     ! W average:
     wavg = wavg + &
      sum(facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc))*volcell

     ! Sum to processors total volume:
     tvol = tvol + real(NXB*NYB*NZB)*volcell

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo
  

  ! Sum processors average velocities:
  uavgaux = uavg
  call MPI_Allreduce(uavgaux, uavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 

  vavgaux = vavg
  call MPI_Allreduce(vavgaux, vavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 

  wavgaux = wavg
  call MPI_Allreduce(wavgaux, wavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 


  ! Sum processors volumes
  tvolaux = tvol  
  call MPI_Allreduce(tvolaux, tvol, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)


  ! Average velocities:
  uavg = uavg/tvol
  vavg = vavg/tvol
  wavg = wavg/tvol


  ! Substract Mean Velocities: 
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     
     facexData(VELC_FACE_VAR,:,:,:) = facexData(VELC_FACE_VAR,:,:,:) - uavg
     faceyData(VELC_FACE_VAR,:,:,:) = faceyData(VELC_FACE_VAR,:,:,:) - vavg
     facezData(VELC_FACE_VAR,:,:,:) = facezData(VELC_FACE_VAR,:,:,:) - wavg

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo

  !if (gr_meshMe .eq. MASTER_PE) write(*,*) 'Average Velocities=',tvol,uavg,vavg,wavg

  return
  end subroutine


