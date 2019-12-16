!!****if* source/physics/IncompNS/IncompNSMain/boussinesq_stretched/ins_computeDtLocal
!!
!! NAME
!!
!!  ins_computeDtLocal
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***

subroutine ins_computeDtLocal(blockID,   & 
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc )

  use IncompNS_data, ONLY : ins_cflflg, ins_cfl, ins_sigma, ins_invsqrtRa_Pr, ins_dtspec

  use Grid_interface, ONLY : Grid_getBlkCenterCoords
  
  use Grid_data, ONLY : gr_meshMe

  use Heat_AD_data, only: ht_invsqrtRaPr

  implicit none

#include "Flash.h"
#include "constants.h"
  integer, intent(IN) :: blockID
  integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  integer, intent(IN) :: isize,jsize,ksize
  real, dimension(:,:),      intent(IN) :: dx, dy, dz
  real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  ! Local variables:
  real, parameter :: eps = 1.e-12
  real :: dtc,dtv,dtl,velcoeff
  integer :: j,k

  if (ins_cflflg .eq. 0) then
     dtlocal    = ins_dtspec
     lminloc(:) = 0
     return
  endif



#if NDIM == MDIM
  velcoeff = eps

  do k = GRID_KLO, GRID_KHI
    do j = GRID_JLO, GRID_JHI 

      velcoeff = MAX(velcoeff,                                                                                    &
                 MAXVAL(ABS(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,j,k))*dx(GRID_ILO:GRID_IHI+1,LEFT_EDGE)), &
                 MAXVAL(ABS(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI  ,j,k))*dy(j,LEFT_EDGE)),                   &
                 MAXVAL(ABS(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI  ,j,k))*dz(k,LEFT_EDGE)))
    end do
  end do

  velcoeff = MAX(velcoeff,                                                                                                  &
             MAXVAL(ABS(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI+1,GRID_KLO:GRID_KHI))*dy(GRID_JHI+1,LEFT_EDGE)), &
             MAXVAL(ABS(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KHI+1))*dz(GRID_KHI+1,LEFT_EDGE)))

  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps
  endif
  
  dtv = ins_sigma / (MAX(ins_invsqrtRa_Pr, ht_invsqrtRaPr)*MAX( 2.*MAXVAL(dx(:,CENTER))**2.,  &
                                                                2.*MAXVAL(dy(:,CENTER))**2.,  &
                                                                2.*MAXVAL(dz(:,CENTER))**2. ))

#else
  velcoeff = eps

  do j = GRID_JLO, GRID_JHI
    velcoeff = MAX(velcoeff,                                                                                    &
               MAXVAL(ABS(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,j,1))*dx(GRID_ILO:GRID_IHI+1,LEFT_EDGE)), &
               MAXVAL(ABS(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,j,1))*dy(j,LEFT_EDGE)))
  end do

  velcoeff = MAX(velcoeff,MAXVAL(ABS(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI+1,1))*dy(GRID_JHI+1,LEFT_EDGE)))
                  
  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps  
  endif
  
  dtv = ins_sigma / (MAX(ins_invsqrtRa_Pr, ht_invsqrtRaPr)*MAX( 2.*MAXVAL(dx(:,CENTER))**2.,  &
                                                                2.*MAXVAL(dy(:,CENTER))**2. ))
#endif

  dtl = MIN(dtc,dtv)

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
     ! Cell located at center of Block - Used to define block location. 
     lminloc(IAXIS) = NGUARD + NXB/2
     lminloc(JAXIS) = NGUARD + NYB/2
#if NDIM == MDIM
     lminloc(KAXIS) = NGUARD + NZB/2
#else
     lminloc(KAXIS) = CONSTANT_ONE
#endif
     lminloc(4) = blockID
     lminloc(5) = gr_meshMe
  endif

  return

end subroutine ins_computeDtLocal
