!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/2D/sm_el01_mapParticles
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

subroutine sm_el01_mapParticles_inout(body, e, ptelem,  &
                                xpos,ypos,        &
                                xvel,yvel,        &
                                xacc,yacc,        &
                                xnrm,ynrm,        &
                                areai, loc_num,   &
                                outloc )
  use SolidMechanics_data, only: sm_structure
  use Driver_interface, only: Driver_abortFlash
  use Driver_data, only: dr_simTime

  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  integer, intent(in)  :: ptelem
  real, dimension(ptelem) :: xpos, ypos, xvel, yvel, xacc, yacc, xnrm, ynrm, areai, outloc
  integer, dimension(ptelem) :: loc_num

  ! Internal variables
  integer, parameter :: nen_e = TWO_NODES
  integer :: a, i, idx1, idx2, nEta,ind
  integer :: idx_in, idx_out
  real, dimension(nen_e,NDIM) :: u, ud, udd
  real, dimension(NDIM) :: norm,v12
  real :: del,area,area_sub
  real :: ei,N1,N2

  !
  ! Get absolute position of each node in the element
  !

  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     u(a,1) = body%x(idx1)
     u(a,2) = body%y(idx1)
    
     do i = 1,NDIM

        if(i==1) then

        idx2     = body%ID(i,idx1)

        body%qn(idx2)  = 0.0*dr_simTime
        body%qdn(idx2) = 0.0

        u(a,i)   = u(a,i) + body%qn(idx2)
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)
 
        else
        idx2     = body%ID(i,idx1)
        u(a,i)   = u(a,i) + body%qn(idx2)
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)

        end if

     end do

     !Inflow-Outflow condition
     idx_in = body%ws_IEN_inflow(a,e)
     if(idx_in == 1) then
       ind = e+(a-1)
       ud(a,IAXIS) =  0.0
       ud(a,JAXIS) = -1.5*(body%x(ind)*(1.0-body%x(ind)))*4.0      
       !print*,"this is ",a,body%xB(ind),body%x(ind),ud(a,2)
     endif

     !idx_out = body%ws_IEN_outflow(a,e)
     !if(idx_out == 1) then
     !  ind = e+(a-1)
       !ud(a,JAXIS) =  0.0
       !ud(a,IAXIS) = 1.5*(body%y(ind)*(1.0-body%y(ind)))*4.0  
       !print*,"this is ",ind,body%y(ind),ud(a,IAXIS)
     !endif
 
  enddo
 
  ! nXi -> nEta
  nEta  = body%ws_nXi(e)
  
  !Get the normal:
  v12(IAXIS)=u(2,IAXIS)-u(1,IAXIS);
  v12(JAXIS)=u(2,JAXIS)-u(1,JAXIS);

  ! Segment Length
  area = sqrt( v12(IAXIS)**2. + v12(JAXIS)**2. )

  ! Normal outside: Counter clockwise node numbering in 2D 
  ! Therefore n = v12 x k 
  norm(IAXIS) =  v12(JAXIS)/area
  norm(JAXIS) = -v12(IAXIS)/area

  area_sub = area/real(ptelem);
  del=1./real(nEta);

  idx_out = body%ws_IEN_outflow(1,e)

  ! Map nEta particles:
  do i = 1, nEta

     ei  = (0.5+real(i-1))*del

     ! Shape functions:
     N1 = 1. - ei
     N2 = ei  

     if(idx_out == 1) then
       outloc(i) = 1.0       
     else
       outloc(i) = 0.0
     endif

     ! x,y positions of internal particles:
     xpos(i) = N1*u(1,IAXIS) + N2*u(2,IAXIS);
     ypos(i) = N1*u(1,JAXIS) + N2*u(2,JAXIS);

     ! Vel
     xvel(i) = N1*ud(1,IAXIS) + N2*ud(2,IAXIS);
     yvel(i) = N1*ud(1,JAXIS) + N2*ud(2,JAXIS);

     ! Acc
     xacc(i) = N1*udd(1,IAXIS) + N2*udd(2,IAXIS);
     yacc(i) = N1*udd(1,JAXIS) + N2*udd(2,JAXIS);

     xnrm(i) = norm(IAXIS) ;
     ynrm(i) = norm(JAXIS) ;

     loc_num(i)=i;
     areai(i)  = area_sub;

     !if(idx_out == 1) print*,xpos(i),ypos(i)

  enddo

  return

end subroutine sm_el01_mapParticles_inout
