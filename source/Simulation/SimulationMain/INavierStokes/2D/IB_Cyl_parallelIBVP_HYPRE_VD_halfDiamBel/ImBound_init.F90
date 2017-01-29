!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/ImBound_init
!!
!! NAME
!!
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init(restart)
!!  
!! ARGUMENTS
!!
!! restart - restart flag, logical.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)

  use Driver_data, ONLY: dr_simTime, dr_meshMe

  use Driver_interface, ONLY : Driver_abortFlash, &
       Driver_getMype, Driver_getComm

  use ImBound_data, only : ib_meshMe, ib_meshComm, ib_globalMe
  use ImBound_data
  use gr_sbData, only : gr_sbNumBodies, gr_sbBodyInfo, NodesPerElem

  use Grid_interface, only :  Grid_sbSelectMaster, &
    Grid_sbBroadcastParticles, Grid_getBoundboxCentroids

  use gr_sbInterface, ONLY: gr_sbCreateGroups, gr_sbCreateParticles

  use gr_ptVPData, only : gr_ptVPBufferFactor

  implicit none
#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

  logical, INTENT(IN) :: restart

  ! Local variables:

  real :: xpt, ypt, dtheta, angle, tita, dsb
  integer :: i, b, ibd, nodelocpos, NumVertices, NumAelem

  real :: L1,L2,n1(MDIM),n2(MDIM)


  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)

  ! Set the buffer factor for Virtual particles on guardcell region to 1.5:
  gr_ptVPBufferFactor=1.5


  tita  = omg*dr_simTime 

!#define DATAFLG 1

#ifdef DATAFLG 
!DATAFLG == OURS

  !========= KPD Does Not Go In Here ! ====================

  ib_nbd = 1

  ib_nela  = 628   !200 !380
  ib_nnoda = ib_nela + 1

  allocate(ib_ABODY(ib_nbd))

  !print*,"DATAFLG ALLOCATING XO,YO,ZO !!!!!!!!!!"
  allocate(xo(ib_nbd),yo(ib_nbd))  


  ibd = 1 
  ib_ABODY(ibd) % lb = 0
  ib_ABODY(ibd) % mb = ib_nnoda
  ib_ABODY(ibd) % elb = 0
  ib_ABODY(ibd) % emb = ib_nela
  ib_ABODY(ibd) % memflag = 1
  ib_ABODY(ibd) % ron = 0.

  !========= KPD Does Not Go In Here ! ====================
  xo(ibd) = 0. !0.5
  yo(ibd) = -0.75 !0.5

  ! Define positions and normals
  dtheta = 2.*PI/ib_nela

  angle = 0.
  ib_sb(1) = 0.
  do i=1,ib_nnoda


    ! Positions of markers, normals and arclength var
    ib_xbus(i) = Ro*cos(angle)
    ib_ybus(i) = Ro*sin(angle)

    ib_xb(i) = xo + ib_xbus(i)*cos(tita) - ib_ybus(i)*sin(tita)
    ib_yb(i) = yo + ib_xbus(i)*sin(tita) + ib_ybus(i)*cos(tita)

    ib_nxL(i) = ib_xb(i)/Ro
    ib_nyL(i) = ib_yb(i)/Ro

    if (i .gt. 1) then
       !dsb = sqrt( (ib_xb(i)-ib_xb(i-1))**2 + (ib_yb(i)-ib_yb(i-1))**2 )
       dsb = dtheta*Ro
       ib_sb(i) = ib_sb(i-1) + dsb

       ib_AELEM(1,i-1) = 2
       ib_AELEM(2,i-1) = i-1
       ib_AELEM(3,i-1) = i
    endif

    ! Set velocities and accelerations
    ib_ubd(i) = -omg*Ro*sin(angle) ! 0.
    ib_vbd(i) =  omg*Ro*cos(angle) ! 0.
    ib_ubdd(i)= 0.
    ib_vbdd(i)= 0.
 
    ! Update angle
    angle = angle + dtheta

  end do

!#define DATAFLG

#else 
!if defined DATAFLG 

  !KPD... Goes Here for Parallel Implementation

  do b=1,gr_sbNumBodies

     !---------------------------------------------------
     !-------- KPD added for ALL procs ------------------
     !---------------------------------------------------
     !print*,"NO DATAFLG ALLOCATING XO,YO,ZO !!!!!!!!!!"
     allocate(xo(gr_sbNumBodies),yo(gr_sbNumBodies))
     xo(b) =  0.0
     yo(b) = -0.75
     gr_sbBodyInfo(b)%xo = xo(b)
     gr_sbBodyInfo(b)%yo = yo(b)
     !---------------------------------------------------

     if (ib_meshMe .eq. gr_sbBodyInfo(b)%bodyMaster) then 

     NumVertices = gr_sbBodyInfo(b)%NumVertices
     NumAelem = gr_sbBodyInfo(b)%NumAelem

     gr_sbBodyInfo(b) % AELTYPE(1:NumAelem) =  TWO_NODE_SEGMENT_VERT !TWO_NODE_SEGMENT_VERT

!     !KPD
!     allocate(xo(gr_sbNumBodies),yo(gr_sbNumBodies))
!     xo(b) =  0.0
!     yo(b) = -0.75
!     gr_sbBodyInfo(b)%xo = xo(b)
!     gr_sbBodyInfo(b)%yo = yo(b)


     gr_sbBodyInfo(b) % sb(1) = 0.
     do i=1,NumVertices

        gr_sbBodyInfo(b)%xb(i) = xo(b) + gr_sbBodyInfo(b)%xbus(i)*cos(tita) - gr_sbBodyInfo(b)%ybus(i)*sin(tita)
        gr_sbBodyInfo(b)%yb(i) = yo(b) + gr_sbBodyInfo(b)%xbus(i)*sin(tita) + gr_sbBodyInfo(b)%ybus(i)*cos(tita) + ao*sin(2.*PI*freq_t*dr_simTime) 

        ! Set velocities and accelerations
        gr_sbBodyInfo(b)%ubd(i) = -omg*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*cos(tita) 
        gr_sbBodyInfo(b)%vbd(i) =  omg*gr_sbBodyInfo(b)%xbus(i)*cos(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*sin(tita) + (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime)
        gr_sbBodyInfo(b)%ubdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*cos(tita) + omg**2.*gr_sbBodyInfo(b)%ybus(i)*sin(tita)
        gr_sbBodyInfo(b)%vbdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg**2.*gr_sbBodyInfo(b)%ybus(i)*cos(tita) - (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)

       if (i .gt. 1) then
          dsb = sqrt( (gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2 + (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2 )
          gr_sbBodyInfo(b)%sb(i) = gr_sbBodyInfo(b)%sb(i-1) + dsb
       endif

     enddo

     ! Normals:
     i = 1
     L1    = sqrt((gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))**2. + &
                  (gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1))**2.)
     L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
                  (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2.) 

     n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1))
     n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))

     n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
     n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))

     gr_sbBodyInfo(b) % nxl(NumVertices) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(NumVertices) = 0.5*(n1(2)+n2(2))

     gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))

     do i =2,gr_sbBodyInfo(b)%NumVertices-1

       L1    = sqrt((gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2. + &
                    (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2.)
       L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
                    (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2.) 

       n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))
       n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))

       n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
       n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))

       gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
       gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))
 
     enddo
     
     ! Write .dat file:
     open(unit=113,file='./IOData/CYL2.dat',form='formatted')
     write(113,'(I8)') gr_sbBodyInfo(b) % NumVertices
     write(113,'(2I8)') gr_sbBodyInfo(b) % NumAelem,NodesPerElem

     do i = 1,gr_sbBodyInfo(b) % NumVertices
        write(113,'(2E22.12)') gr_sbBodyInfo(b)%xbus(i),gr_sbBodyInfo(b)%ybus(i)
     enddo

     do i = 1,gr_sbBodyInfo(b) % NumAelem
        write(113,'(2I8)') gr_sbBodyInfo(b)%AELEM(2,i),gr_sbBodyInfo(b)%AELEM(3,i)
     enddo
     close(113)

     endif


  enddo


  call Grid_getBoundboxCentroids()
  call Grid_sbSelectMaster()
  call gr_sbCreateParticles()
  call Grid_sbBroadcastParticles()


#endif






  if(ib_globalMe==MASTER_PE)print*,'Immersed_Boundaries initialized'
  return

end subroutine ImBound_init
