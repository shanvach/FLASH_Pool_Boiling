subroutine Plasma_computeDtLocal(blockID,   &
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              solnData,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc )

  use Plasma_data, ONLY : pls_cflflg, pls_cfl, pls_sigma, pls_dtspec, pls_dcoeff

  use Grid_interface, ONLY : Grid_getBlkCenterCoords

  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash.h"
#include "constants.h"
  integer, intent(IN) :: blockID
  integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  integer, intent(IN) :: isize,jsize,ksize
  real, intent(IN) :: dx, dy, dz
  real, pointer,dimension(:,:,:,:)  :: solnData,facexData,faceyData,facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  ! Local variables:
  real, parameter :: eps = 1.e-12
  real :: dtl,pi
  real :: dcoeff

  pi = acos(-1.0)

  if (pls_cflflg .eq. 0) then
     dtLocal    = pls_dtspec
     lminloc(:) = 0
     return
  endif

  dcoeff = MAX(MAXVAL(ABS(solnData(DFEL_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH0_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH1_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH2_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH3_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH4_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFH5_VAR,:,:,:))),&
               MAXVAL(ABS(solnData(DFEL_VAR,:,:,:))))
               !MAXVAL(ABS(solnData(DFH6_VAR,:,:,:))),&
               !MAXVAL(ABS(solnData(DFH7_VAR,:,:,:))),&
               !MAXVAL(ABS(solnData(DFH8_VAR,:,:,:))),&
               !MAXVAL(ABS(solnData(DFH9_VAR,:,:,:))))
       
  dtl = (pls_cfl/dcoeff)*min(dx**2,dy**2) !adaptive time step

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
     lminloc(1) = NGUARD + NXB/2
     lminloc(2) = NGUARD + NYB/2
     lminloc(3) = NGUARD + NZB/2
     lminloc(4) = blockID
     lminloc(5) = gr_meshMe

  endif

end subroutine Plasma_computeDtLocal
