subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_invRe,ins_meshMe
   use Multiphase_data, only: mph_thco2,mph_vis2,mph_cp2,mph_rho2,&
                              mph_thco1,mph_vis1,mph_cp1,mph_rho1

   use Grid_interface, ONLY : Grid_getDeltas,         &
                                 Grid_getBlkIndexLimits, &
                                 Grid_getCellCoords,     &
                                 Grid_getBlkPtr,         &
                                 Grid_releaseBlkPtr,     &
                                 Grid_getBlkBoundBox,    &
                                 Grid_getBlkCenterCoords, &
                                 Grid_fillGuardCells

   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 
   use Heat_AD_interface, ONLY: Heat_getQmicro

   implicit none

#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb,i,j,k
   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   integer :: ierr,iter
   real :: maxdfun_local, maxdfun_global
   real :: beta, chi, soln, a_I, b_I, x1, x2, f1, f2, h
   real :: del(MDIM)
   real :: thermalBL_dt
   integer :: thermalBL_nt
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

   real :: dxmin

   call RuntimeParameters_get("Pr",ht_Pr)
   call RuntimeParameters_get("St",ht_St)
   call RuntimeParameters_get("hfit",ht_hfit)
   call RuntimeParameters_get("Ab",ht_Ab)
   call RuntimeParameters_get("Cb",ht_Cb)
   call RuntimeParameters_get("Bb",ht_Bb)
   call RuntimeParameters_get("Ra",ht_Ra)

   if (ins_meshMe .eq. MASTER_PE) then

     write(*,*) 'ht_Pr   =',ht_Pr
     write(*,*) 'ht_St   =',ht_St
     write(*,*) 'ht_hfit =',ht_hfit
     write(*,*) 'ht_Ab   =',ht_Ab
     write(*,*) 'ht_Bb   =',ht_Bb
     write(*,*) 'ht_Cb   =',ht_Cb
     write(*,*) 'ht_Ra   =',ht_Ra

   end if

   ht_Twall_low    =  1.0
   ht_Twall_high   =  0.0
   ht_Tsat         =  0.0
   ht_AMR_specs(:) =  0.0

   dxmin    = 1e10

   do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getDeltas(blockID,del)

     dxmin = min(dxmin,del(JAXIS))

   end do

   call MPI_ALLREDUCE(dxmin,ht_dxmin,1,FLASH_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)

   if (ins_meshMe .eq. MASTER_PE) call Heat_getQmicro(ht_qmic,ht_fmic,ht_dxmin)

   call MPI_BCAST(ht_qmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(ht_fmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

   print *,"qmic,fmic: ",ht_qmic,ht_fmic

   !thermalBL_dt = 0.0001/0.0101

   !thermalBL_nt = 0

   !do while(thermalBL_dt*thermalBL_nt .lt. (0.08/0.0101))

   !do lb = 1,blockCount

   !  blockID = blockList(lb)

   !  call Grid_getDeltas(blockID,del)

   !  call Grid_getBlkPtr(blockID,solnData,CENTER)

   !  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

   !  call Grid_getBlkPtr(blockID,solnData,CENTER)
   !  call Grid_getBlkPtr(blockID,facexData,FACEX)
   !  call Grid_getBlkPtr(blockID,faceyData,FACEY)

   !  solnData(TOLD_VAR,:,:,:) = solnData(TEMP_VAR,:,:,:)

   !  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
   !   do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
   !     do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

   !       solnData(TEMP_VAR,i,j,k) = solnData(TOLD_VAR,i,j,k) + &
   !                                  thermalBL_dt*ins_invRe/ht_Pr*&
   !                                  (solnData(TOLD_VAR,i-1,j,k)+solnData(TOLD_VAR,i+1,j,k)-2*solnData(TOLD_VAR,i,j,k))/(del(DIR_X)**2) + &
   !                                   thermalBL_dt*ins_invRe/ht_Pr*&
   !                                  (solnData(TOLD_VAR,i,j-1,k)+solnData(TOLD_VAR,i,j+1,k)-2*solnData(TOLD_VAR,i,j,k))/(del(DIR_Y)**2) 
                                     


   !     end do
   !   end do
   !  end do

   !  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
   !  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
   !  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

   ! gcMask = .FALSE.
   ! gcMask(TEMP_VAR)=.TRUE.

   ! call Grid_fillGuardCells(CENTER,ALLDIR,&
   !      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

   !end do

  !thermalBL_nt = thermalBL_nt + 1

  !end do

end subroutine Heat_AD_init
