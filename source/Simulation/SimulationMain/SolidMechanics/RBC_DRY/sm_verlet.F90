subroutine sm_verlet_partone(ibd,dt)
  
#include "Flash.h"
#include "SolidMechanics.h"
  
  
  use SolidMechanics_rbc_data, ONLY: sm_rbc_Lambda,sm_rbc_mi,Fext_M
  use SolidMechanics_data, Only : sm_bodyInfo, sm_structure
  use Simulation_data, Only : Nplus, Nminus
  implicit none
  
  ! Argument list
  integer,INTENT(IN) :: ibd
  real, intent(IN) :: dt
  type(sm_structure),  pointer :: body
  real :: dtt
  integer :: i,nnp
  integer,allocatable :: inds(:)

  body => sm_bodyinfo(ibd)
  
 ! write(*,*) 'Fluid forces'
 ! write(*,*) 'body%Hs= '
 ! do i=1,3*body%nnp
 ! write(*,*) body%Hs(i)
 ! end do 
      nnp=body%nnp 
      allocate(inds(nnp))
      write(*,*) 'nnp',nnp
      inds(:)=(/(i,i=1,nnp*3,3)/)
      !write(*,*) inds(1:10)
      !stop
      if (ibd==1) then 
         body%Hs(inds)=0.
      elseif (ibd==2) then
         body%Hs(:)=0.
      end if

      
      !write(*,*) 'body%Hs(1:3:nnp-2)',body%Hs(1:9),' dt=',dt
      dtt=dt;
      body%qn(:)= body%qms(:,1)  + body%qdms(:,1) * dtt + 0.5 * (body%qddms(:,1) + body%Hs(:)+ &
           body%Cddms(:,1) ) * dtt *dtt / 1.;
      body%qdi(:)= body%qdms(:,1) + (sm_rbc_Lambda * dtt * (body%qddms(:,1) + body%Hs(:))+ body%Cddms(:,1) ) /1.;
      
      
!  Pos2 = Pos1 + (vel1 * dt) + (0.5 * Fn * dt * dt) / 1.; 
!  vel2 = vel1 + (sm_rbc_Lambda * dt * Fn) /1.;
    
  if (maxval(body%qn(:))>100) then
     print*,'There is a problem with the Pos2'
     stop
  end if
end subroutine sm_verlet_partone

subroutine sm_verlet_parttwo(ibd,dt)
  
#include "Flash.h"
#include "SolidMechanics.h"
  
  
  use SolidMechanics_rbc_data, ONLY: sm_rbc_mi
  use SolidMechanics_data, Only : sm_bodyInfo, sm_structure
  
  implicit none
  
  ! Argument list
  integer,INTENT(IN) :: ibd
  real, intent(IN) :: dt
  type(sm_structure),  pointer :: body
  integer :: i
  real ::ddt
  ddt=dt;
  body => sm_bodyInfo(ibd)
  
  !body%qdms(:,1),body%qdn(:),body%qddms(:,1),body%qddn(:),body%nnp

  body%qdn(:)  = body%qdms(:,1) + (0.5 * (body%qddms(:,1)+body%qddn(:)) + body%Hs(:) + &
        body%Cddms(:,1))*ddt /1.;
  !vel2 = vel1 + (0.5 * dt *(Fn+Fnp1))/1.

  body%qms(:,1)   = body%qn(:);      ! Positions
  body%qdi(:)     = body%qdms(:,1);  ! Old velocities to intermediate velocities
  body%qdms(:,1)  = body%qdn(:);     ! Velocities
  body%qddms(:,1) = body%qddn(:);    ! Forces


end subroutine sm_verlet_parttwo
