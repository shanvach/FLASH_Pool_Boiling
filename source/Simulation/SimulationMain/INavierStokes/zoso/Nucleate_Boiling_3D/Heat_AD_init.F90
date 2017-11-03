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
                                 Grid_getBlkCenterCoords

   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Heat_AD_interface, only: Heat_getQmicro
   use Driver_data, only: dr_restart

   implicit none

#include "constants.h"
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
   real :: dxmin
   real :: del(MDIM)

   call RuntimeParameters_get("Pr",ht_Pr)
   call RuntimeParameters_get("St",ht_St)
   call RuntimeParameters_get("hfit",ht_hfit)
   call RuntimeParameters_get("Ab",ht_Ab)
   call RuntimeParameters_get("Cb",ht_Cb)
   call RuntimeParameters_get("Bb",ht_Bb)

   if (ins_meshMe .eq. MASTER_PE) then

     write(*,*) 'ht_Pr   =',ht_Pr
     write(*,*) 'ht_St   =',ht_St
     write(*,*) 'ht_hfit =',ht_hfit
     write(*,*) 'ht_Ab   =',ht_Ab
     write(*,*) 'ht_Bb   =',ht_Bb
     write(*,*) 'ht_Cb   =',ht_Cb

   end if

   ht_Twall_low  = 1.0
   ht_Twall_high = 0.0
   ht_Tsat       = 0.0
   ht_AMR_specs  = 0.0

   if(dr_restart .eqv. .FALSE.) then

        !dxmin    = 1e10

        !do lb = 1,blockCount

        !  blockID = blockList(lb)
        !  call Grid_getDeltas(blockID,del)
        !  dxmin = min(dxmin,del(JAXIS))

        !end do

        !call MPI_ALLREDUCE(dxmin,ht_dxmin,1,FLASH_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)

        ht_dxmin = 0.01

        if (ins_meshMe .eq. MASTER_PE) call Heat_getQmicro(ht_qmic,ht_fmic,ht_dxmin)

        call MPI_BCAST(ht_qmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ht_fmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

   end if

   print *,"qmic,fmic: ",ht_qmic,ht_fmic

end subroutine Heat_AD_init
