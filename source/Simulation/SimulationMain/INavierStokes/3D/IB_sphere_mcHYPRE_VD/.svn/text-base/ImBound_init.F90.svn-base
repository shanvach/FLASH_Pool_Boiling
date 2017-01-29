!!****if* source/Simulation/SimulationMain/INavierStokes/3D/IB_sphere_mcHYPRE_VD/ImBound_init
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
    Grid_sbBroadcastParticles

  use gr_sbInterface, ONLY: gr_sbCreateParticles, gr_sbFinalize

  use gr_ptVPData, only : gr_ptVPBufferFactor

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


  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)

  ! Set the buffer factor for Virtual particles on guardcell region to 1.5:
  gr_ptVPBufferFactor=1.5


!#define DATAFLG 1

#ifdef DATAFLG 

   !_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_
   !_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_

#else 

  !KPD... Goes Here for Parallel Implementation

  do b=1,gr_sbNumBodies

     if (ib_meshMe .eq. gr_sbBodyInfo(b)%bodyMaster) then 

     allocate(xo(gr_sbNumBodies))
     allocate(yo(gr_sbNumBodies))
     allocate(zo(gr_sbNumBodies))

     NumVertices = gr_sbBodyInfo(b)%NumVertices
     NumAelem = gr_sbBodyInfo(b)%NumAelem

     gr_sbBodyInfo(b) % AELTYPE(1:NumAelem) = THREE_NODE_TRIANG_VERT 

     !KPD
     xo(b) =  0.00
     yo(b) =  0.00
     zo(b) =  0.00

     gr_sbBodyInfo(b)%xo = xo(b)
     gr_sbBodyInfo(b)%yo = yo(b)
     gr_sbBodyInfo(b)%zo = zo(b)


     do i=1,gr_sbBodyInfo(b)%NumVertices

        write(*,*) i,gr_sbBodyInfo(b)%xbus(i),gr_sbBodyInfo(b)%ybus(i),gr_sbBodyInfo(b)%zbus(i)

        ! Positions of markers, normals and arclength var
        gr_sbBodyInfo(b)%xb(i) =  gr_sbBodyInfo(b)%xo + gr_sbBodyInfo(b)%xbus(i)   
        gr_sbBodyInfo(b)%yb(i) =  gr_sbBodyInfo(b)%yo + gr_sbBodyInfo(b)%ybus(i)
        gr_sbBodyInfo(b)%zb(i) =  gr_sbBodyInfo(b)%zo + gr_sbBodyInfo(b)%zbus(i)
    
        gr_sbBodyInfo(b)%nxL(i) =  (gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xo)/(Ro)
        gr_sbBodyInfo(b)%nyL(i) =  (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yo)/(Ro)
        gr_sbBodyInfo(b)%nzL(i) =  (gr_sbBodyInfo(b)%zb(i)-gr_sbBodyInfo(b)%zo)/(Ro)

        nrm = sqrt(gr_sbBodyInfo(b)%nxL(i)**2 + gr_sbBodyInfo(b)%nyL(i)**2 + gr_sbBodyInfo(b)%nzL(i)**2)
 
        gr_sbBodyInfo(b)%nxL(i) = gr_sbBodyInfo(b)%nxL(i)/nrm
        gr_sbBodyInfo(b)%nyL(i) = gr_sbBodyInfo(b)%nyL(i)/nrm
        gr_sbBodyInfo(b)%nzL(i) = gr_sbBodyInfo(b)%nzL(i)/nrm


        ! Set velocities and accelerations
        gr_sbBodyInfo(b)%ubd(i) = 0.
        gr_sbBodyInfo(b)%vbd(i) = 0.
        gr_sbBodyInfo(b)%wbd(i) = 0.
        gr_sbBodyInfo(b)%ubdd(i)= 0.
        gr_sbBodyInfo(b)%vbdd(i)= 0.
        gr_sbBodyInfo(b)%wbdd(i)= 0.

     enddo



     ! Write .dat file:
     open(unit=113,file="sphere_coarse_altered_write.surf",form='formatted')
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


!  call Grid_getBoundboxCentroids()
!  call Grid_sbSelectMaster()
!  call gr_sbCreateParticles()
!  call Grid_sbBroadcastParticles()
!  call gr_sbFinalize()

  call Grid_getBoundboxCentroids()
  call Grid_sbSelectMaster()
  call gr_sbCreateParticles()
  call Grid_sbBroadcastParticles()




!  call gr_sbFinalize()

#endif

  if(ib_globalMe==MASTER_PE)print*,'Immersed_Boundaries initialized'
  return

end subroutine ImBound_init
