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
  use Heat_AD_data, ONLY  : ht_invsqrtRaPr 
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID
  integer,dimension(2,MDIM), intent(IN) :: blkLimits, blkLimitsGC
  integer, intent(IN) :: isize, jsize, ksize
  real, dimension(:,:), intent(IN) :: dx, dy, dz
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData, facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  real, parameter :: eps = 1.e-12
  real :: dtl, velcoeff, nu
  integer :: imin, imax, jmin, jmax, kmin, kmax
  integer :: i, j, k

  if (ins_cflflg .eq. 0) then
     dtlocal    = ins_dtspec
     lminloc(:) = 0
     return
  endif

  imin = blkLimits(LOW, IAXIS)
  jmin = blkLimits(LOW, JAXIS)
  kmin = blkLimits(LOW, KAXIS)
  imax = blkLimits(HIGH, IAXIS)
  jmax = blkLimits(HIGH, JAXIS)
  kmax = blkLimits(HIGH, KAXIS)

  nu = max(ins_invsqrtRa_Pr, ht_invsqrtRaPr) / ins_sigma

#if NDIM == MDIM
  velcoeff = eps

  do k = kmin, kmax 
    do j = jmin, jmax
      do i = imin, imax

        velcoeff = max(velcoeff,                                       &
            abs(facexData(VELC_FACE_VAR,i+1,j,k)) * dx(i,RIGHT_EDGE) + &
            abs(faceyData(VELC_FACE_VAR,i,j+1,k)) * dy(j,RIGHT_EDGE) + &
            abs(facezData(VELC_FACE_VAR,i,j,k+1)) * dz(k,RIGHT_EDGE) + &
            nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 + dz(k,CENTER)**2 ) )

      end do
    end do
  end do

  do k = kmin, kmax
    do j = jmin, jmax
      i = imin
      velcoeff = max(velcoeff,                                       &
          abs(facexData(VELC_FACE_VAR,i  ,j,k)) * dx(i,LEFT_EDGE ) + &
          abs(faceyData(VELC_FACE_VAR,i,j+1,k)) * dy(j,RIGHT_EDGE) + &
          abs(facezData(VELC_FACE_VAR,i,j,k+1)) * dz(k,RIGHT_EDGE) + &
          nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 + dz(k,CENTER)**2 ) )
    end do
  end do

  do k = kmin, kmax
    do i = imin, imax
      j = jmin
      velcoeff = max(velcoeff,                                       &
          abs(facexData(VELC_FACE_VAR,i+1,j,k)) * dx(i,RIGHT_EDGE) + &
          abs(faceyData(VELC_FACE_VAR,i,j  ,k)) * dy(j,LEFT_EDGE ) + &
          abs(facezData(VELC_FACE_VAR,i,j,k+1)) * dz(k,RIGHT_EDGE) + &
          nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 + dz(k,CENTER)**2 ) )
    end do
  end do

  do j = jmin, jmax
    do i = imin, imax
      k = kmin
      velcoeff = max(velcoeff,                                       &
          abs(facexData(VELC_FACE_VAR,i+1,j,k)) * dx(i,RIGHT_EDGE) + &
          abs(faceyData(VELC_FACE_VAR,i,j+1,k)) * dy(j,RIGHT_EDGE) + &
          abs(facezData(VELC_FACE_VAR,i,j,k  )) * dz(k,LEFT_EDGE ) + &
          nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 + dz(k,CENTER)**2 ) )
    end do
  end do

  dtl = ins_cfl / velcoeff
  
#else
  velcoeff = eps

  do j = jmin, jmax
    do i = imin, imax

      velcoeff = max(velcoeff,                                       &
          abs(facexData(VELC_FACE_VAR,i+1,j,1)) * dx(i,RIGHT_EDGE) + &
          abs(faceyData(VELC_FACE_VAR,i,j+1,1)) * dy(j,RIGHT_EDGE) + &  
          nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 ) )

    end do
  end do

  do j = jmin, jmax
    i = imin
    velcoeff = max(velcoeff,                                       &
        abs(facexData(VELC_FACE_VAR,i  ,j,1)) * dx(i,LEFT_EDGE ) + &
        abs(faceyData(VELC_FACE_VAR,i,j+1,1)) * dy(j,RIGHT_EDGE) + &
        nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 ) )
  end do

  do i = imin, imax
    j = jmin
    velcoeff = max(velcoeff,                                       &
        abs(facexData(VELC_FACE_VAR,i+1,j,1)) * dx(i,RIGHT_EDGE) + &
        abs(faceyData(VELC_FACE_VAR,i,j  ,1)) * dy(j,LEFT_EDGE ) + &
        nu * ( dx(i,CENTER)**2 + dy(j,CENTER)**2 ) )
  end do

  dtl = ins_cfl / velcoeff
  
#endif

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
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
