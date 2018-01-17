subroutine Plasma_computeDtLocal(blockID,   &
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
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
  real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  ! Local variables:
  real, parameter :: eps = 1.e-12
  real :: dtl,pi


  pi = acos(-1.0)

  if (pls_cflflg .eq. 0) then
     dtLocal    = pls_dtspec
     lminloc(:) = 0
     return
  endif

  dtl = (pls_cfl/pls_dcoeff)*min(dx**2,dy**2) !adaptive time step

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
     lminloc(1) = NGUARD + NXB/2
     lminloc(2) = NGUARD + NYB/2
     lminloc(3) = NGUARD + NZB/2
     lminloc(4) = blockID
     lminloc(5) = gr_meshMe

  endif

end subroutine Plasma_computeDtLocal
