!!****if* source/Simulation/SimulationMain/INavierStokes/3D/Snorkel_mcHYPRE_VD_wFS/ImBound_init
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
!! VARIABLES
!!
!! restart = restart flag, logical.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)

  use Driver_data, ONLY: dr_simTime

  use Driver_interface, ONLY : Driver_abortFlash, &
       Driver_getMype, Driver_getComm

  use ImBound_data, only : ib_meshMe, ib_meshComm, ib_globalMe
  use ImBound_data
  use gr_sbData, only : gr_sbNumBodies, gr_sbBodyInfo, NodesPerElem

  use Grid_interface, only :  Grid_sbSelectMaster, &
    Grid_sbBroadcastParticles, Grid_getBoundboxCentroids

  use gr_sbInterface, ONLY: gr_sbCreateGroups, gr_sbCreateParticles

  use gr_ptVPData, only : gr_ptVPBufferFactor

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

  logical, INTENT(IN) :: restart

  ! Local variables:

  real :: dtheta, angle, tita, dsb
  integer :: i, b, ibd, nodelocpos, NumVertices, NumAelem

  real :: L1,L2,n1(MDIM),n2(MDIM)
  real :: nrm

  real, save :: cylDepth

  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)

  ! Set the buffer factor for Virtual particles on guardcell region to 1.5:
  gr_ptVPBufferFactor=1.5


  tita  = omg*dr_simTime 

!#define DATAFLG 1

   call RuntimeParameters_get('cylDepth',    cylDepth)
print*,"CYLINDER dEPTH",cylDepth

#ifdef DATAFLG 

   !_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_
   !_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_

#else 

  !KPD... Goes Here for Parallel Implementation

  allocate(xo(gr_sbNumBodies))
  allocate(yo(gr_sbNumBodies))
  allocate(zo(gr_sbNumBodies))

  do b=1,gr_sbNumBodies

     !KPD
     xo(b) = 0.0
     yo(b) = 0.0
     zo(b) = cylDepth 

     gr_sbBodyInfo(b)%xo = xo(b)
     gr_sbBodyInfo(b)%yo = yo(b)
!#if NDIM == 3
     gr_sbBodyInfo(b)%zo = zo(b)
!#endif


     if (ib_meshMe .eq. gr_sbBodyInfo(b)%bodyMaster) then 

     NumVertices = gr_sbBodyInfo(b)%NumVertices
     NumAelem = gr_sbBodyInfo(b)%NumAelem

     !gr_sbBodyInfo(b) % AELTYPE(1:NumAelem) =  TWO_NODE_SEGMENT_VERT 
     gr_sbBodyInfo(b) % AELTYPE(1:NumAelem) = THREE_NODE_TRIANG_VERT 

#if NDIM == 2
     gr_sbBodyInfo(b)%sb(1) = 0.
#endif

     do i=1,NumVertices

!        gr_sbBodyInfo(b)%xb(i) = xo(b) + gr_sbBodyInfo(b)%xbus(i)*cos(tita) - gr_sbBodyInfo(b)%ybus(i)*sin(tita)
!        gr_sbBodyInfo(b)%yb(i) = yo(b) + gr_sbBodyInfo(b)%xbus(i)*sin(tita) + gr_sbBodyInfo(b)%ybus(i)*cos(tita) + ao*sin(2.*PI*freq_t*dr_simTime) 
!!#if NDIM == 3
!        gr_sbBodyInfo(b)%zb(i) = zo(b)
        gr_sbBodyInfo(b)%xb(i) =  gr_sbBodyInfo(b)%xo + gr_sbBodyInfo(b)%xbus(i)
        gr_sbBodyInfo(b)%yb(i) =  gr_sbBodyInfo(b)%yo + gr_sbBodyInfo(b)%ybus(i)
        gr_sbBodyInfo(b)%zb(i) =  gr_sbBodyInfo(b)%zo + gr_sbBodyInfo(b)%zbus(i)
!#endif

        ! Set velocities and accelerations
        gr_sbBodyInfo(b)%ubd(i) = -omg*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*cos(tita) 
        gr_sbBodyInfo(b)%vbd(i) =  omg*gr_sbBodyInfo(b)%xbus(i)*cos(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*sin(tita) + (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime)
        gr_sbBodyInfo(b)%ubdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*cos(tita) + omg**2.*gr_sbBodyInfo(b)%ybus(i)*sin(tita)
        gr_sbBodyInfo(b)%vbdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg**2.*gr_sbBodyInfo(b)%ybus(i)*cos(tita) - (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)

       if (i .gt. 1) then
          dsb = sqrt( (gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2 + (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2 )

#if NDIM == 2
          gr_sbBodyInfo(b)%sb(i) = gr_sbBodyInfo(b)%sb(i-1) + dsb
#endif
       endif

     enddo

     ! Normals:
     i = 1
     L1    = sqrt((gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))**2. + &
                  (gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1)) + &
                  (gr_sbBodyInfo(b)%zb(NumVertices)-gr_sbBodyInfo(b)%zb(NumVertices-1))**2.)
     L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
                  (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2. + &
                  (gr_sbBodyInfo(b)%zb(i+1)-gr_sbBodyInfo(b)%zb(i))**2.) 

     n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1))
     n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))

     n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
     n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))

     gr_sbBodyInfo(b) % nxl(NumVertices) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(NumVertices) = 0.5*(n1(2)+n2(2))

     gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))

     do i =2,gr_sbBodyInfo(b)%NumVertices-1

!       !L1    = sqrt((gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2. + &
!       !             (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2.)
!       !L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
!       !             (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2.) 
!       L1    = sqrt((gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2. + &
!                    (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2. + &
!                    (gr_sbBodyInfo(b)%zb(i)-gr_sbBodyInfo(b)%zb(i-1))**2.)
!       L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
!                    (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2. + &
!                    (gr_sbBodyInfo(b)%zb(i+1)-gr_sbBodyInfo(b)%zb(i))**2.) 
!
!       n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))
!       n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))
!
!       n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
!       n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))
!
!       !KPD - x,y,z Normals
!       gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
!       gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))

       !KPD - x,y,z Normals
       gr_sbBodyInfo(b)%nxL(i) =  (gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xo)!/(Ro)
       gr_sbBodyInfo(b)%nyL(i) =  (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yo)!/(Ro)
       gr_sbBodyInfo(b)%nzL(i) =  (gr_sbBodyInfo(b)%zb(i)-gr_sbBodyInfo(b)%zo)!/(Ro)
       nrm = sqrt(gr_sbBodyInfo(b)%nxL(i)**2 + gr_sbBodyInfo(b)%nyL(i)**2 + gr_sbBodyInfo(b)%nzL(i)**2)
       gr_sbBodyInfo(b)%nxL(i) = gr_sbBodyInfo(b)%nxL(i)/nrm
       gr_sbBodyInfo(b)%nyL(i) = gr_sbBodyInfo(b)%nyL(i)/nrm
       gr_sbBodyInfo(b)%nzL(i) = gr_sbBodyInfo(b)%nzL(i)/nrm


 
     enddo

     ! Write .dat file:
     open(unit=113,file='./IOData/splashCyl_LxThree.surf',form='formatted')
     write(113,'(I8)') gr_sbBodyInfo(b) % NumVertices
     write(113,'(2I8)') gr_sbBodyInfo(b) % NumAelem,NodesPerElem

     do i = 1,gr_sbBodyInfo(b) % NumVertices
       !write(113,'(2E22.12)') gr_sbBodyInfo(b)%xbus(i),gr_sbBodyInfo(b)%ybus(i)
        write(113,'(3E22.12)') gr_sbBodyInfo(b)%xbus(i),gr_sbBodyInfo(b)%ybus(i),gr_sbBodyInfo(b)%zbus(i)
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
