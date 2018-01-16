subroutine Plasma_init(blockCount,blockList,restart)

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

   pls_dcoeff = 1.0

  if (pls_meshMe .eq. MASTER_PE) then
     write(*,*) 'pls_cfl   =',pls_cfl
     write(*,*) 'pls_sigma =',pls_sigma
     write(*,*) 'pls_dtspec=',pls_dtspec
  endif

end subroutine Plasma_init
