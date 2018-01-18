subroutine Plasma_init(blockCount,blockList,restart)

! Read runtime parameters and initialize constants declared in Plasma_data
! script

   use Plasma_data
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                                Driver_getComm, Driver_getNstep




   implicit none

  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"
#include "Plasma.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   logical, INTENT(IN)    :: restart


   call Driver_getMype(MESH_COMM, pls_meshMe)
   call Driver_getNumProcs(MESH_COMM, pls_meshNumProcs)
   call Driver_getComm(MESH_COMM, pls_meshComm)

   call RuntimeParameters_get("cflflg", pls_cflflg)
   call RuntimeParameters_get("cfl", pls_cfl)
   call RuntimeParameters_get("sigma",pls_sigma)
   call RuntimeParameters_get("dtspec",pls_dtspec)
   call RuntimeParameters_get("vel_prolong_method",pls_prol_method)

   call Driver_getNstep(pls_nstep)
   pls_restart=restart

   pls_dcoeff = 1e3 ! test electron diffusion coefficient
   
   pls_Ckb = 1.38064852e-23
   pls_Cme = 9.10938356e-31
   pls_Ce  = 1.60217662e-19
   pls_gam = EXP(0.577) 
   pls_Cpi = 3.14159265359
   
   ! collision diameters for all heavy species
   pls_RSCD(1) = 1.
   pls_RSCD(2) = 1.
   pls_RSCD(3) = 2*1.55
   pls_RSCD(4) = 2*1.55
   pls_RSCD(5) = 1.55
   pls_RSCD(6) = 2*1.52
   pls_RSCD(7) = 1.52
   pls_RSCD(8) = 2*1.52
   pls_RSCD(9) = 1.52
   pls_RSCD(10) = 1.52+1.55

   ! molar binary masses
   pls_MHSP(1) = 4.0
   pls_MHSP(2) = 4.0
   pls_MHSP(3) = 28.0
   pls_MHSP(4) = 28.0
   pls_MHSP(5) = 14.0
   pls_MHSP(6) = 32.0
   pls_MHSP(7) = 16.0
   pls_MHSP(8) = 32.0
   pls_MHSP(9) = 16.0 
   pls_MHSP(10) = 14.0 + 16.0
   
  if (pls_meshMe .eq. MASTER_PE) then
     write(*,*) 'pls_cfl   =',pls_cfl
     write(*,*) 'pls_sigma =',pls_sigma
     write(*,*) 'pls_dtspec=',pls_dtspec
  endif

end subroutine Plasma_init
