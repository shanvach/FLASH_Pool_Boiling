!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_EstDT_rigid
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!! 
!! Compute DT for body
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

subroutine sm_EstDT_rigid(ibd, dt)

  use IncompNS_data, ONLY : ins_invRe
  use SolidMechanics_data, only: sm_bodyInfo, sm_structure

  implicit none
  ! IO Variables
  integer, intent(in)  :: ibd
  real,    intent(out) :: dt
 
  ! Internal variables:
   type(sm_structure),  pointer :: body
   integer :: imaster,i, idx, i_dim
   integer :: ix,ex,ia,ea,iw,ew
   real :: x_bod,y_bod,z_bod, xd_bod,yd_bod,zd_bod,xdd_bod,ydd_bod,zdd_bod
   real :: umag

   ! Get the body
   body => sm_BodyInfo(ibd)
 
   ! Extract info from Master:
   imaster = body % borigin_node
 
   ix = body%ix; ex = body%ex;
 
   xd_bod  = body%qdn(body%ID(ix,imaster))
   yd_bod  = body%qdn(body%ID(ix+1,imaster))
   zd_bod  = body%qdn(body%ID(ix+2,imaster))

   umag = sqrt(xd_bod**2+yd_bod**2+zd_bod**2)
 
   dt = 1.e12

   if(umag > 1.0d-12 ) then
     dt = 1.0d-4/ins_invRe*umag
   endif
   if(dt < 1.0d-4) dt = 1.0d-4

  return

end subroutine sm_EstDT_rigid

