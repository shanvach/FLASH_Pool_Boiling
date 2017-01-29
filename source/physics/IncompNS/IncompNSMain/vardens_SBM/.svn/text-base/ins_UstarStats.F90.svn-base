




subroutine ins_UstarStats( blockCount, blockList, print_flg, qin_flg)

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
  logical, INTENT(IN) :: print_flg, qin_flg
  !! -----------------------------------------------------


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,solnData

  integer :: lb,blockID,ierr

  real :: del(MDIM), dx,dy,dz,dxdydz, delta_massi



  ! Initialize values of div Ustar - delta_mass:
  ins_meanDivUstar = 0.
  ins_deltamass    = 0.

  ! Compute mean divergence of Ustar + Delta_mass.
  delta_massi = 0.
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS); dy = del(JAXIS);  
#if NDIM == MDIM
     dz = del(KAXIS)
#else
     dz = 1.
#endif
     dxdydz = dx*dy*dz

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! compute divergence of intermediate velocities
     call ins_divergence(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       del(IAXIS),del(JAXIS),del(KAXIS),&
                       solnData(DUST_VAR,:,:,:) )

     ! Delta_mass = -int_Omg(div(Ustar))dOmg (- sign because is the 
     ! mass injected in the domain, i,e, Qin) Remember : density =1.
     delta_massi = delta_massi - dxdydz * &
      sum(solnData(DUST_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif            
  enddo

  ! Gather total delta mass:
  call MPI_Allreduce(delta_massi,ins_deltamass,1,FLASH_REAL,    &
                     FLASH_SUM, gr_meshComm, ierr)

  ! Mean Div of Ustar:
  call gr_findMean(DUST_VAR,2,.false.,ins_meanDivUstar)

  ! Add residual delta_mass to inflow volume. 
  if (qin_flg) ins_Qin = ins_Qin + ins_deltamass

  if ((ins_meshMe .eq. MASTER_PE) .and. print_flg) &
  write(*,*) 'Mean DivUstar, DeltaMass=', ins_meanDivUstar, ins_deltamass

  return

end subroutine ins_UstarStats
