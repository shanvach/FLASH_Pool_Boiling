!    
! File:  
! 
! Modified from the file sm_assemble_IntForce_rigid.F90.
! To calculate the repulsive forces among rigid bodies.
! The model in the work of Wan and Turek (2007) is used.
!  Wan and Turek, Journal of Computational and Applied Mathematics, 2007
!  (203):561-580
! Shizhao Wang
! Nov 05, 2014

subroutine sm_assemble_repulForce_rigid(ibd)

  use SolidMechanics_Data, only : sm_BodyInfo, sm_pos, sm_structure, sm_numBodies, sm_cm_rhou, sm_cm_eps1, sm_cm_eps2
!  use sm_Misc_interface, only : sm_crossProd
  use Driver_interface, only : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas

  implicit none
#include "SolidMechanics.h"
#include "Flash.h"
#include "constants.h"    
  !IO Variables
  integer, intent(in)    :: ibd ! body number
  
  ! Local variables:
  type(sm_structure), pointer :: body
  
  real :: rhou, eps1, eps2
  real :: rCir, rpos(NDIM), rr, rc1, rc2, fRpl(NDIM)

  real :: del(MDIM), dh
  integer :: ind, idx, i, imaster

  ! Get the body
  body => sm_BodyInfo(ibd)

  imaster = BODYFRAME_NODE

  ! Set the paramters in the model
  rhou = sm_cm_rhou
  eps1 = sm_cm_eps1
  eps2 = sm_cm_eps2

  rCir = 0.5

  rc1 = 2.*rCir
  rc2 = rc1 + rhou


  fRpl = 0.
  do ind = 1, ibd-1
    rpos(:) = sm_pos(:,ibd) - sm_pos(:,ind)
    rr = sqrt(dot_product(rpos,rpos))
    if(rr > rc2) then
      cycle
    elseif(rr>rc1) then
      fRpl(:) = fRpl(:) + rpos(:)*(rc2-rr)*(rc2-rr)/eps1
    else
      fRpl(:) = fRpl(:) + rpos(:)*(rc1-rr)/eps2
    endif
  enddo

  do ind = ibd+1, sm_numBodies
    rpos(:) = sm_pos(:,ibd) - sm_pos(:,ind)
    rr = sqrt(dot_product(rpos,rpos))
    if(rr > rc2) then
      cycle
    elseif(rr>rc1) then
      fRpl(:) = fRpl(:) + rpos(:)*(rc2-rr)*(rc2-rr)/eps1
    else
      fRpl(:) = fRpl(:) + rpos(:)*(rc1-rr)/eps2
    endif
  enddo

  ! Add to whatever is in Rpl:
  ! Forces:
  do i = body%ix,body%ex
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Rpl(idx)  = fRpl(i)
     endif
  end do

  ! Orientation:
  do i = body%ia,body%ea
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Rpl(idx)  = 0.
     endif
  end do 

  ! Moments: 
  do i = body%iw,body%ew
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Rpl(idx)  = 0.
     endif
  end do


!  write(*,*) 'rc1:', rc1
!  write(*,*) 'rc2:', rc2
!  write(*,*) 'sm_pos:', sm_pos
!  write(*,*) 'fRpl:', fRpl, 'ibd', ibd
!  write(*,*) 'Rpl:', body%Rpl, 'ibd', ibd

  return
    
end subroutine sm_assemble_repulForce_rigid
