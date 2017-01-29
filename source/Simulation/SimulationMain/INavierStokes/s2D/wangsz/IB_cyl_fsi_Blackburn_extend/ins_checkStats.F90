




subroutine ins_checkStats( blockCount, blockList, var_flg, note)

#include "Flash.h"
  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits

  use Grid_data, only : gr_domainBC,gr_meshComm

  use IncompNS_data, only : ins_meshMe, ins_Qin, ins_meanDivUstar, ins_deltamass

  use gr_interface, ONLY : gr_findMean

  use ins_interface, only : ins_divergence 

  implicit none

#include "constants.h"
  include "Flash_mpi.h"
  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
  integer, INTENT(IN) :: var_flg
  character(len=*), INTENT(IN) :: note
  !! -----------------------------------------------------


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,solnData

  integer :: lb,blockID,ierr

  real :: del(MDIM), dx,dy,dz,dxdydz, delta_massi
  
  real :: max_var, min_var, mean_var, max_var_tmp, min_var_tmp


  max_var = 0.0
  min_var = 0.0
  mean_var = 0.0

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif


#if NDIM == 2
      max_var_tmp = maxval(solnData(var_flg,GRID_ILO:GRID_IHI, & 
                  & GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI))

      min_var_tmp = minval(solnData(var_flg,GRID_ILO:GRID_IHI, & 
                  & GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI))

#endif
 
    ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif            
  enddo

   call MPI_Allreduce(max_var_tmp, max_var, 1, FLASH_REAL,&
                      MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(min_var_tmp, min_var, 1, FLASH_REAL,&
                      MPI_MIN, MPI_COMM_WORLD, ierr)


  ! Gather total delta mass:
!  call MPI_Allreduce(delta_massi,ins_deltamass,1,FLASH_REAL,    &
!                     FLASH_SUM, gr_meshComm, ierr)

  ! Mean Div of Ustar:
  call gr_findMean(var_flg, 2,.false.,mean_var)


  if ((ins_meshMe .eq. MASTER_PE))  then
  write(*,*) '====================='
  write(*,*) note
  write(*,*) 'var id:', var_flg
  write(*,*) 'max var', max_var
  write(*,*) 'min var', min_var
  write(*,*) 'mean var', mean_var
  write(*,*) '====================='
  endif

  return

end subroutine ins_checkStats
