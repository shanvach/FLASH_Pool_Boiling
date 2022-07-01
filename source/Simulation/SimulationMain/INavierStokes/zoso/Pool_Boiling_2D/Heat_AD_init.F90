subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_invRe,ins_meshMe
   use Multiphase_data, only: mph_thco2,mph_vis2,mph_cp2,mph_rho2,&
                              mph_thco1,mph_vis1,mph_cp1,mph_rho1,&
                              mph_timeStampAll

   use Grid_interface, ONLY : Grid_getDeltas,         &
                                 Grid_getBlkIndexLimits, &
                                 Grid_getCellCoords,     &
                                 Grid_getBlkPtr,         &
                                 Grid_releaseBlkPtr,     &
                                 Grid_getBlkBoundBox,    &
                                 Grid_getBlkCenterCoords

   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Heat_AD_interface, only: Heat_getQmicro

   use Simulation_data, only: sim_nucSiteDens, sim_nuc_radii, &
                              sim_nuc_site_x, sim_nuc_site_y, &
                              sim_nuc_site_z, sim_sinkB

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
   integer :: nuc_index,tSI
   real, dimension(MDIM)  :: coord,bsize
   real ::  boundBox(2,MDIM),ycell

   call RuntimeParameters_get("Pr",ht_Pr)
   call RuntimeParameters_get("St",ht_St)
   call RuntimeParameters_get("hfit",ht_hfit)
   call RuntimeParameters_get("Ab",ht_Ab)
   call RuntimeParameters_get("Cb",ht_Cb)
   call RuntimeParameters_get("Bb",ht_Bb)
   call RuntimeParameters_get("tsat",ht_Tsat)
   call RuntimeParameters_get("Ra",ht_Ra)
!---------------------------------------------------------------------------------------------Shantanu
   call RuntimeParameters_get("twait_hydrophobic",ht_tWait_hydrophobic)
   call RuntimeParameters_get("twait_hydrophilic",ht_tWait_hydrophilic)
!-------------------------------------------------------
   call RuntimeParameters_get("tnuc",ht_Tnuc)
   call RuntimeParameters_get("qmic",ht_qmic)
   call RuntimeParameters_get("fmic",ht_fmic)

   if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) 'ht_Pr     =',ht_Pr
     write(*,*) 'ht_St     =',ht_St
     write(*,*) 'ht_hfit   =',ht_hfit
     write(*,*) 'ht_Ab     =',ht_Ab
     write(*,*) 'ht_Bb     =',ht_Bb
     write(*,*) 'ht_Cb     =',ht_Cb
     write(*,*) 'ht_Tsat   =',ht_Tsat
     write(*,*) 'ht_Ra     =',ht_Ra
!--------------------------------------------------------------------------------Shantanu
     write(*,*) 'ht_tWait_hydrophilic  =',ht_tWait_hydrophilic
     write(*,*) 'ht_tWait_hydrophobic  =',ht_tWait_hydrophobic
!---------------------------------------------------------------------------------
     write(*,*) 'sim_sinkB =',sim_sinkB
     write(*,*) 'ht_Tnuc   =',ht_Tnuc
   end if

   ht_Twall_low  = 1.0
   ht_Twall_high = 0.0
   ht_AMR_specs  = 0.0

   ht_microFlg   = .TRUE.

   if(dr_restart .eqv. .TRUE.) then

     if (ins_meshMe .eq. MASTER_PE) print *,"Entering heat restart 1"

     ht_psi          = (45.0/180.0)*acos(-1.0)

   end if

end subroutine Heat_AD_init
