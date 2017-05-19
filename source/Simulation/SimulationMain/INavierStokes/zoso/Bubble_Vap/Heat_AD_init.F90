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
   real :: beta, chi, soln, a_I, b_I, x1, x2, f1, f2, h,cps

   call RuntimeParameters_get("Pr",ht_Pr)
   call RuntimeParameters_get("St",ht_St)
   call RuntimeParameters_get("hfit",ht_hfit)
   call RuntimeParameters_get("beta",beta)

   if (ins_meshMe .eq. MASTER_PE) then

     write(*,*) 'ht_Pr   =',ht_Pr
     write(*,*) 'ht_St   =',ht_St
     write(*,*) 'ht_hfit =',ht_hfit
     write(*,*) 'beta    =',beta

   end if

   ht_Twall_low  = 1.0
   ht_Twall_high = 1.0
   ht_Tsat       = 0.0

   !beta = sqrt(3/acos(-1.0))*(ht_St)*(mph_rho2/mph_rho1)

   !if((mph_rho2 .gt. 10.0)  .and. (mph_rho2 .lt. 30.0))  !beta =
   !if((mph_rho2 .gt. 50.0)  .and. (mph_rho2 .lt. 70.0))  !beta =
   !if((mph_rho2 .gt. 100.0) .and. (mph_rho2 .lt. 150.0)) beta = 6.8074241158925353901506643989232
   !if((mph_rho2 .gt. 1600.0)) beta = 5.0880841483201690539025583264782

   chi  = 1.0 - (mph_rho1/mph_rho2)
   cps  = 1.0 - ((mph_cp1/mph_rho1)/(mph_cp2/mph_rho2))
   soln = 0.0

   do lb = 1,blockCount
     blockID = blockList(lb)

    ! call Grid_getDeltas(blockID,del)

    ! call Grid_getBlkBoundBox(blockId,boundBox)
    ! bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

    ! call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)


     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     !maxdfun_local = maxval(abs(solnData(DFUN_VAR,:,:,:)))

     !call MPI_Allreduce(maxdfun_local, maxdfun_global, 1, FLASH_REAL,&
     !                MPI_MAX, MPI_COMM_WORLD, ierr)

     do k=1,blkLimitsGC(HIGH,KAXIS)
      do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS) 

        if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) then
     
        solnData(TEMP_VAR,i,j,k) = ht_Tsat
        !solnData(TEMP_VAR,i,j,k) = ht_Twall_low

        else 


        a_I   = 1 - 0.005/(0.005+abs(solnData(DFUN_VAR,i,j,k)))
        b_I   = 0.99
        soln   = 0.
        h     = (b_I-a_I)/200

        do iter=1,200

             x1 = a_I + h*(iter-1)
             x2 = a_I + h*(iter)

             !f1 = ((1-x1)**(1-2*chi*beta*beta))*exp((-beta*beta)/((1-x1)**2))*(1./((1-x1)**2));
             !f2 = ((1-x2)**(1-2*chi*beta*beta))*exp((-beta*beta)/((1-x2)**2))*(1./((1-x2)**2));

             f1 = exp(2*chi*beta*beta*x1)*exp((-beta*beta)/((1-x1)**2));
             f2 = exp(2*chi*beta*beta*x2)*exp((-beta*beta)/((1-x2)**2));

            soln = soln + (h/2)*(f1+f2)

        end do

        solnData(TEMP_VAR,i,j,k) = ht_Twall_low - 2.0*beta*beta*(mph_rho1/mph_rho2)*((1/ht_St)+cps)*exp(beta*beta)*soln;
        !solnData(TEMP_VAR,i,j,k) = ht_Tsat

        end if 


        if(solnData(TEMP_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = 0.

         end do
       end do
     end do

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

   end do

end subroutine Heat_AD_init
